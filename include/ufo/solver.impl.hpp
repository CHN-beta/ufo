# pragma once
# include <ufo/solver.hpp>

namespace Eigen
{
  constexpr inline auto serialize(auto & archive, Eigen::Matrix3d& matrix)
    { return archive(std::span(matrix.data(), matrix.size())); }
  constexpr inline auto serialize(auto & archive, const Eigen::Matrix3d& matrix)
    { return archive(std::span(matrix.data(), matrix.size())); }
  constexpr inline auto serialize(auto & archive, Eigen::Vector3d& vector)
    { return archive(std::span(vector.data(), vector.size())); }
  constexpr inline auto serialize(auto & archive, const Eigen::Vector3d& vector)
    { return archive(std::span(vector.data(), vector.size())); }
}

namespace ufo
{
  inline HighFive::CompoundType create_compound_complex()
    { return {{ "r", HighFive::AtomicType<double>{}}, {"i", HighFive::AtomicType<double>{}}}; }

  inline concurrencpp::generator<std::pair<Eigen::Vector<unsigned, 3>, unsigned>> Solver::triplet_sequence
    (Eigen::Vector<unsigned, 3> range)
  {
    for (unsigned x = 0; x < range[0]; x++)
      for (unsigned y = 0; y < range[1]; y++)
        for (unsigned z = 0; z < range[2]; z++)
          co_yield
          {
            Eigen::Vector<unsigned, 3>{{x}, {y}, {z}},
            x * range[1] * range[2] + y * range[2] + z
          };
  }

  template <ZppSerializable T> inline void Solver::zpp_write(const T& object, std::string filename)
  {
    auto [data, out] = zpp::bits::data_out();
    out(object).or_throw();
    static_assert(sizeof(char) == sizeof(std::byte));
    std::ofstream file(filename, std::ios::binary | std::ios::out);
    file.exceptions(std::ios::badbit | std::ios::failbit);
    file.write(reinterpret_cast<const char*>(data.data()), data.size());
  }
  template <ZppSerializable T> inline T Solver::zpp_read(std::string filename)
  {
    auto input = std::ifstream(filename, std::ios::binary | std::ios::in);
    input.exceptions(std::ios::badbit | std::ios::failbit);
    static_assert(sizeof(std::byte) == sizeof(char));
    std::vector<std::byte> data;
    {
      std::vector<char> string(std::istreambuf_iterator<char>(input), {});
      data.assign
      (
        reinterpret_cast<std::byte*>(string.data()),
        reinterpret_cast<std::byte*>(string.data() + string.size())
      );
    }
    auto in = zpp::bits::in(data);
    T output;
    in(output).or_throw();
    return output;
  }

  inline Solver::Hdf5file& Solver::Hdf5file::open_for_read(std::string filename)
  {
    File_ = HighFive::File(filename, HighFive::File::ReadOnly);
    return *this;
  }
  inline Solver::Hdf5file& Solver::Hdf5file::open_for_write(std::string filename)
  {
    File_ = HighFive::File(filename, HighFive::File::ReadWrite | HighFive::File::Create
      | HighFive::File::Truncate);
    return *this;
  }
  template <typename T> inline Solver::Hdf5file& Solver::Hdf5file::read(T& object, std::string name)
  {
    object = File_->getDataSet(name).read<std::remove_cvref_t<decltype(object)>>();
    return *this;
  }
  template <typename T> inline Solver::Hdf5file& Solver::Hdf5file::write(const T& object, std::string name)
  {
    File_->createDataSet(name, object);
    return *this;
  }

  inline Solver::DataFile::DataFile
    (YAML::Node node, std::set<std::string> supported_format, std::string config_file, bool allow_same_as_config_file)
  {
    if (auto _ = node["SameAsConfigFile"])
    {
      auto __ = _.as<bool>();
      if (__ && !allow_same_as_config_file)
        throw std::runtime_error("\"SameAsConfigFile: true\" is not allowed here.");
      ExtraParameters["SameAsConfigFile"] = __;
      if (__)
      {
        Filename = config_file;
        Format = "yaml";
        return;
      }
    }
    Filename = node["Filename"].as<std::string>();
    Format = node["Format"].as<std::string>();
    if (!supported_format.contains(Format))
      throw std::runtime_error(fmt::format("Unsupported format: \"{}\"", Format));
    if (auto _ = node["RelativeToConfigFile"])
    {
      auto __ = _.as<bool>();
      ExtraParameters["RelativeToConfigFile"] = __;
      if (__)
        Filename = std::filesystem::path(config_file).parent_path() / Filename;
    }
  };

  template <typename T> inline Solver::PipizedFunction<T>::PipizedFunction(T&& function)
    : Function(std::forward<T>(function)) {}
  template <typename T> inline decltype(auto) operator|(auto&& input, Solver::PipizedFunction<T>&& function)
    { return function.Function(std::forward<decltype(input)>(input)); }
  template <typename T> inline Solver::PipizedFunction<T> Solver::pipize(T&& function)
    { return PipizedFunction(std::forward<T>(function)); }

  inline constexpr UnmoveHealerType unmove;
  template <typename T> inline decltype(auto) operator|(T&& input, const UnmoveHealerType&)
    { return static_cast<T&>(input); }

  template <std::size_t N = std::dynamic_extent> inline constexpr ToSpanHelperType<N> to_span;
  template <int N> inline decltype(auto) operator|(auto&& input, const ToSpanHelperType<N>&)
  {
    if constexpr (N == std::dynamic_extent)
      return std::span(input.data(), input.size());
    else
      return std::span(input.data(), N);
  }

  inline constexpr ToEigenHelper to_eigen;
  template <typename T, std::size_t N> inline decltype(auto) operator|
    (const std::span<T, N>& input, const ToEigenHelper&)
  {
    if constexpr (N == std::dynamic_extent)
      return Eigen::VectorX<T>(input.data(), input.size());
    else
      return Eigen::Vector<T, N>(input.data());
  }
  template <typename T, std::size_t M, std::size_t N> inline decltype(auto) operator|
    (const std::span<std::span<T, N>, M>& input, const ToEigenHelper&)
  {
    if constexpr (M == std::dynamic_extent)
      if constexpr (N == std::dynamic_extent)
      {
        Eigen::Matrix
        return Eigen::MatrixX<T>(input.data()->data(), input.size(), input.data()->size());
      }
      else
        return Eigen::Matrix<T, Eigen::Dynamic, N>(input.data()->data(), input.size());
  }
}
