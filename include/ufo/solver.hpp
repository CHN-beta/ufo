# pragma once
# include <iostream>
# include <array>
# include <numbers>
# include <numeric>
# include <fstream>
# include <optional>
# include <array>
# include <utility>
# include <execution>
# include <syncstream>
# include <any>
# include <map>
# include <vector>
# include <span>
# include <ranges>
# include <yaml-cpp/yaml.h>
# include <Eigen/Dense>
# include <concurrencpp/concurrencpp.h>
# include <fmt/format.h>
# include <fmt/std.h>
# include <fmt/ranges.h>
# include <highfive/H5File.hpp>
# include <zpp_bits.h>
# include <matplot/matplot.h>
# include <matplot/backend/opengl.h>
# include <range/v3/all.hpp>
# include <biu.hpp>

// 在相位中, 约定为使用 $\exp (2 \pi i \vec{q} \cdot \vec{r})$ 来表示原子的运动状态
//  (而不是 $\exp (-2 \pi i \vec{q} \cdot \vec{r})$)
// 一些书定义的倒格矢中包含了 $2 \pi$ 的部分, 我们这里约定不包含这部分.
//  也就是说, 正格子与倒格子的转置相乘, 得到单位矩阵.

namespace Eigen
{
  constexpr auto serialize(auto & archive, Eigen::Matrix3d& matrix);
  constexpr auto serialize(auto & archive, const Eigen::Matrix3d& matrix);
  constexpr auto serialize(auto & archive, Eigen::Vector3d& vector);
  constexpr auto serialize(auto & archive, const Eigen::Vector3d& vector);
}

namespace ufo
{
  using namespace std::literals;
  struct PhonopyComplex { double r, i; };
  HighFive::CompoundType create_compound_complex();

	namespace detail_
	{
		template <typename T> struct SpecializationOfBitsMembersHelper : std::false_type {};
		template <std::size_t N> struct SpecializationOfBitsMembersHelper<zpp::bits::members<N>> : std::true_type {};
	}
	template <typename T> concept ZppSerializable
    = requires() { detail_::SpecializationOfBitsMembersHelper<T>::value == true; };

  struct Solver
  {
    Solver() = delete;

    static concurrencpp::generator<std::pair<Eigen::Vector<unsigned, 3>, unsigned>>
      triplet_sequence(Eigen::Vector<unsigned, 3> range);

    template <ZppSerializable T> static void zpp_write(const T& object, std::string filename);
    template <ZppSerializable T> static T zpp_read(std::string filename);

    class Hdf5file
    {
      public:
        Hdf5file& open_for_read(std::string filename);
        Hdf5file& open_for_write(std::string filename);
        template <typename T> Hdf5file& read(T& object, std::string name);
        template <typename T> Hdf5file& write(const T& object, std::string name);
      protected:
        std::optional<HighFive::File> File_;
    };

    struct DataFile
    {
      std::string Filename;
      std::string Format;
      std::map<std::string, std::any> ExtraParameters;
      inline DataFile() = default;
      DataFile
      (
        YAML::Node node, std::set<std::string> supported_format,
        std::string config_file, bool allow_same_as_config_file = false
      );
    };

    template <typename T> struct PipizedFunction
    {
      T Function;
      PipizedFunction(T&& function);
    };
    template <typename T> static PipizedFunction<T> pipize(T&& function);
  };
  template <typename T> static decltype(auto) operator|(auto&& input, Solver::PipizedFunction<T>&& function);

  struct UnmoveHealerType {};
  template <typename T> decltype(auto) operator|(T&& input, const UnmoveHealerType&);

  template <std::size_t N> struct ToSpanHelperType {};
  template <std::size_t N> decltype(auto) operator|(auto&& input, const ToSpanHelperType<N>&);
}

HIGHFIVE_REGISTER_TYPE(ufo::PhonopyComplex, ufo::create_compound_complex)
