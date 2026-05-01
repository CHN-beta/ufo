# include <ufo.hpp>

void ufo::mode_effective_charge(std::string config_file)
{
  struct Config { std::string InputFile, OutputFile, VaspOutputFile; };
  struct Output { std::vector<std::vector<std::array<double, 3>>> Super; };

  biu::Logger::Guard log(config_file);
  auto config = YAML::LoadFile(config_file).as<Config>();
  auto input = biu::deserialize<CommonData>
    (biu::read<std::byte>(config.InputFile));
  auto mass = input.Super.AtomType
    | std::views::transform([&](auto&& atom)
      { return std::views::repeat(input.AtomMass[atom.first], atom.second); })
    | std::views::join
    | std::ranges::to<std::vector<double>>();
  // 下标：[原子序号](位移方向, 电场方向)
  auto bec = biu::Hdf5file(config.VaspOutputFile)
    .read<std::vector<Eigen::Matrix3d>>("/results/linear_response/born_charges");
  assert(bec.size() == mass.size());
  auto calc_mec = [&](const Eigen::MatrixX3cd& mode, const std::vector<Eigen::Matrix3d>& bec)
  {
    Eigen::Vector3d mec = Eigen::Vector3d::Zero();
    for (std::size_t i = 0; i < mode.rows(); i++)
      for (std::size_t j = 0; j < 3; j++) // 电场方向
        for (std::size_t k = 0; k < 3; k++) // 位移方向
          mec[j] += mode(i, k).real() / std::sqrt(mass[i]) * bec[i](k, j);
    return mec;
  };
  Output output { .Super = input.Super.Qpoint
    | std::views::transform([&](auto&& qpoint)
      {
        return qpoint.Mode
          | std::views::transform([&](auto&& mode) { return calc_mec(mode.EigenVector, bec) | biu::fromEigen; })
          | std::ranges::to<std::vector<std::array<double, 3>>>();
      })
    | std::ranges::to<std::vector<std::vector<std::array<double, 3>>>>()
  };
  std::ofstream(config.OutputFile) << YAML::Node(output);
}
