
# include <ufo.hpp>

namespace ufo
{
  struct RamanData
  {
    struct ModeType
    {
      std::size_t QpointIndex;
      std::size_t ModeIndex;
      double Ratio;
      using serialize = zpp::bits::members<3>;
    };
    std::vector<ModeType> Mode;
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
      // 要计算的模式，总是假定是 gamma 点的模式
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
      std::size_t i_of_poscar = 0;
      auto mass = cell.AtomType
        | ranges::views::transform([&](auto&& atom)
          { return ranges::views::repeat_n(input.AtomMass[atom.first], atom.second); })
        | ranges::views::join | ranges::to_vector | biu::toEigen<>;
      for (auto i_of_qpoint : config.SelectedModes | ranges::views::keys)
        for (auto i_of_mode : config.SelectedModes[i_of_qpoint])
        {
          // 假定虚部总是为零
          auto atom_movement = cell.Qpoint[i_of_qpoint].Mode[i_of_mode].EigenVector.real()
            .cwiseProduct(mass.cwiseSqrt().cwiseInverse().rowwise().replicate(3)).eval();
          // 归一化
          auto ratio = config.MaxDisplacement / atom_movement.rowwise().norm().maxCoeff();
          atom_movement *= ratio;
          // 输出
          auto path = "{}/{}"_f(config.OutputPoscarDirectory, i_of_poscar);
          std::filesystem::create_directories(path);
          std::ofstream("{}/POSCAR"_f(path)) << generate_poscar
          (
            cell.Cell,
            cell.AtomPosition + atom_movement * cell.Cell.inverse(),
            cell.AtomType
          );
          output.Mode.push_back({i_of_qpoint, i_of_mode, ratio});
          log.debug("Write mode {} {} {}"_f(i_of_qpoint, i_of_mode, atom_movement.rowwise().norm().eval()));
          i_of_poscar++;
        }
    };
    if (config.Cell == RamanData::Primative) process(input.Primative);
    else process(input.Super);

    std::filesystem::create_directories("{}/{}"_f(config.OutputPoscarDirectory, "origin"));
    std::ofstream("{}/origin/POSCAR"_f(config.OutputDataFile)) << generate_poscar
      (input.Super.Cell, input.Super.AtomPosition, input.Super.AtomType);
    std::ofstream(config.OutputDataFile, std::ios::binary) << biu::serialize<char>(output);
  }

  void raman_extract(std::vector<std::string> files)
  {
    biu::Logger::Guard log(files);
    std::vector<Eigen::Matrix3d> electricity_tensors;

    for (const auto& file : files)
    {
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
          electricity_tensors.emplace_back();
          for (std::size_t i = 0; i < 3; i++) for (std::size_t j = 0; j < 3; j++)
            in_stream >> electricity_tensors.back()(i, j);
          log.debug("read tensor {}"_f(electricity_tensors.back()));
          break;
        }
      }
      // test if the file is read correctly
      if (!in_stream) throw std::runtime_error("Error reading file: {}"_f(file));
    }

    // output the result
    for (auto e: electricity_tensors) std::cout <<
R"(- - [ {}, {}, {} ]
  - [ {}, {}, {} ]
  - [ {}, {}, {} ]
)"_f
      (
        e(0, 0), e(0, 1), e(0, 2),
        e(1, 0), e(1, 1), e(1, 2),
        e(2, 0), e(2, 1), e(2, 2)
      );
  }

  void raman_apply_contribution(std::string config_file)
  {
    struct Config
    {
      Eigen::Matrix3d OriginalSusceptibility;
      std::vector<Eigen::Matrix3d> Susceptibilities;
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
      for (auto&& [i_of_mode, mode] : ranges::views::enumerate(raman_input.Mode))
      {
        auto&& _ = cell.Qpoint[mode.QpointIndex].Mode[mode.ModeIndex];
        Eigen::Matrix3d raman_tensor = (config.Susceptibilities[i_of_mode] - config.OriginalSusceptibility)
          / mode.Ratio / raman_input.MaxDisplacement;
        _.RamanTensor = raman_tensor | biu::fromEigen;
        _.WeightOnRaman = config.Polarization[0].transpose() * raman_tensor * config.Polarization[1];
        log.info("{}:{:.2f}:{}:{:.2f} {}"_f
        (
          mode.QpointIndex, fmt::join(cell.Qpoint[mode.QpointIndex].Qpoint, ", "),
          mode.ModeIndex, _.Frequency, *_.WeightOnRaman
        ));
      }
    };
    if (raman_input.Cell == RamanData::Primative) process(input.Primative);
    else process(input.Super);

    std::ofstream(config.OutputDataFile, std::ios::binary) << biu::serialize<char>(input);
  }
}
