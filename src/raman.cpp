# include <ufo.hpp>

namespace ufo
{
  struct RamanData
  {
    std::map<std::pair<std::size_t, std::size_t>, double> ModeRatio;
    double MaxDisplacement;
    enum { Primative, Super } Cell;
    using serialize = zpp::bits::members<3>;
  };

  void raman_create_displacement(std::string config_file)
  {
    struct Config
    {
      // 要计算的是原胞还是超胞
      decltype(RamanData::Cell) Cell;
      // 要计算哪个 q 点的哪些模式，总是假定是 gamma 点的模式
      std::map<std::size_t, std::set<std::size_t>> SelectedModes;
      // 原子最大位移大小，单位为埃
      double MaxDisplacement;
      // 输出的POSCAR所在的目录
      std::string OutputPoscarDirectory;
      // 输入的数据文件名
      std::string InputDataFile;
      std::string OutputDataFile;
    };

    auto generate_poscar = []
    (
      Eigen::Matrix3d SuperCell, Eigen::MatrixX3d AtomPosition,
      std::vector<std::pair<std::string, std::size_t>> AtomType
    )
    {
      std::stringstream ss;
      ss << "some random comment to make VASP happy\n1.0\n";
      for (std::size_t i = 0; i < 3; i++)
      {
        for (std::size_t j = 0; j < 3; j++) ss << SuperCell(i, j) << " ";
        ss << std::endl;
      }
      ss << "{}\n"_f(ranges::accumulate
      (
        AtomType | ranges::views::transform([](auto&& atom) { return atom.first; }),
        ""s, [](auto&& a, auto&& b) { return a + " " + b; }
      ));
      ss << "{}\n"_f(ranges::accumulate
      (
        AtomType | ranges::views::transform([](auto&& atom) { return atom.second; }),
        ""s, [](auto&& a, auto&& b) { return a + " " + std::to_string(b); }
      ));
      ss << "Direct\n";
      for (const auto& position : AtomPosition.rowwise())
      {
        for (std::size_t i = 0; i < 3; i++) ss << position(i) << " ";
        ss << std::endl;
      }
      return ss.str();
    };

    biu::Logger::Guard log(config_file);
    auto config = YAML::LoadFile(config_file).as<Config>();
    auto input = biu::deserialize<CommonData>
      (biu::read<std::byte>(config.InputDataFile));
    RamanData output;
    output.MaxDisplacement = config.MaxDisplacement;
    output.Cell = config.Cell;

    auto process = [&](auto& cell)
    {
      auto mass = cell.AtomType
        | ranges::views::transform([&](auto&& atom)
          { return ranges::views::repeat_n(input.AtomMass[atom.first], atom.second); })
        | ranges::views::join | ranges::to_vector | biu::toEigen<>;
      for (auto i_of_qpoint : config.SelectedModes | ranges::views::keys)
        for (auto i_of_mode : config.SelectedModes[i_of_qpoint])
        {
          // 假定虚部总是为零
          auto move_eigenvector = cell.Qpoint[i_of_qpoint].Mode[i_of_mode].EigenVector.real()
            .cwiseProduct(mass.cwiseSqrt().cwiseInverse().rowwise().replicate(3)).eval();
          // 归一化
          auto ratio = config.MaxDisplacement / move_eigenvector.rowwise().norm().maxCoeff();
          auto atom_movement = (move_eigenvector * ratio).eval();
          // 输出
          for (auto direction : {1, -1})
          {
            auto path = "{}/{}/{}/{}"_f
              (config.OutputPoscarDirectory, i_of_qpoint, i_of_mode, direction > 0 ? "+" : "-");
            std::filesystem::create_directories(path);
            std::ofstream("{}/POSCAR"_f(path)) << generate_poscar
              (cell.Cell, cell.AtomPosition + atom_movement * cell.Cell.inverse() * direction, cell.AtomType);
          }
          output.ModeRatio[{i_of_qpoint, i_of_mode}] = ratio;
          log.debug("Write mode {} {} {}"_f(i_of_qpoint, i_of_mode, atom_movement.rowwise().norm().eval()));
        }
    };
    if (config.Cell == RamanData::Primative) process(input.Primative);
    else process(input.Super);

    std::filesystem::create_directories("{}/{}"_f(config.OutputPoscarDirectory, "origin"));
    std::ofstream("{}/origin/POSCAR"_f(config.OutputPoscarDirectory)) << generate_poscar
      (input.Super.Cell, input.Super.AtomPosition, input.Super.AtomType);
    std::ofstream(config.OutputDataFile, std::ios::binary) << biu::serialize<char>(output);
  }

  void raman_extract(std::string path)
  {
    biu::Logger::Guard log(path);
    auto get_dielectric_tensor = [](std::filesystem::path file)
    {
      biu::Logger::Guard log(file);
      auto in_stream = std::ifstream(file);
      std::string line;
      // search line containing: MACROSCOPIC STATIC DIELECTRIC TENSOR
      while (std::getline(in_stream, line))
      {
        if (line.find("MACROSCOPIC STATIC DIELECTRIC TENSOR") != std::string::npos)
        {
          log.debug("find in {}"_f(file));
          // skip 1 line
          std::getline(in_stream, line);
          Eigen::Matrix3d dielectric_tensor;
          for (std::size_t i = 0; i < 3; i++) for (std::size_t j = 0; j < 3; j++)
            in_stream >> dielectric_tensor(i, j);
          log.debug("read tensor {}"_f(dielectric_tensor));
          return dielectric_tensor;
        }
      }
      // test if the file is read correctly
      throw std::runtime_error("Error reading file: {}"_f(file));
    };
    auto search_path = [&](this auto self, std::size_t indent, std::string path) -> void
    {
      if (std::filesystem::exists(path + "/OUTCAR"))
      {
        auto e = get_dielectric_tensor(path + "/OUTCAR");
        for (std::size_t i = 0; i < 3; i++) std::cout << "{}- [ {}, {}, {} ]\n"_f
          (std::string(indent * 2, ' '), e(i, 0), e(i, 1), e(i, 2));
      }
      else for (auto entry : std::filesystem::directory_iterator(path)) if (entry.is_directory())
      {
        std::cout << "{}{}:\n"_f
        (
          std::string(indent * 2, ' '),
          (std::set{"+"s, "-"s}.contains(entry.path().filename()))
            ? "\"{}\""_f(entry.path().filename()) : entry.path().filename().string()
        );
        self(indent + 1, entry.path());
      }
    };
    search_path(0, path);
  }

  void raman_apply_contribution(std::string config_file)
  {
    struct Config
    {
      Eigen::Matrix3d OriginalSusceptibility;
      std::map<std::size_t, std::map<std::size_t, std::map<std::string, Eigen::Matrix3d>>> Susceptibility;
      std::array<Eigen::Vector3d, 2> Polarization;
      std::string InputDataFile;
      std::string RamanInputDataFile;
      std::string OutputDataFile;
    };

    biu::Logger::Guard log(config_file);
    auto config = YAML::LoadFile(config_file).as<Config>();
    auto input = biu::deserialize<CommonData>
      (biu::read<std::byte>(config.InputDataFile));
    auto raman_input = biu::deserialize<RamanData>
      (biu::read<std::byte>(config.RamanInputDataFile));

    input.RamanPolarization = config.Polarization;
    auto process = [&](auto& cell)
    {
      for (auto&& [i, ratio] : raman_input.ModeRatio)
      {
        auto&& [i_of_qpoint, i_of_mode] = i;
        auto&& _ = cell.Qpoint[i_of_qpoint].Mode[i_of_mode];
        auto&& susceptibility = config.Susceptibility[i_of_qpoint][i_of_mode];
        Eigen::Matrix3d raman_tensor = (susceptibility["+"] - susceptibility["-"]) / 2 / ratio;
        _.RamanTensor = raman_tensor | biu::fromEigen;
        _.WeightOnRaman = config.Polarization[0].transpose() * raman_tensor * config.Polarization[1];
        log.info("{}:{:.2f}:{}:{:.2f} {}"_f
        (
          i_of_qpoint, fmt::join(cell.Qpoint[i_of_qpoint].Qpoint, ", "),
          i_of_mode, _.Frequency, *_.WeightOnRaman
        ));
      }
    };
    if (raman_input.Cell == RamanData::Primative) process(input.Primative);
    else process(input.Super);

    std::ofstream(config.OutputDataFile, std::ios::binary) << biu::serialize<char>(input);
  }
}
