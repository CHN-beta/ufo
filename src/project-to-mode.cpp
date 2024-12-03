# include <ufo.hpp>

void ufo::project_to_mode(std::string config_file)
{
  // 将选定的超胞中的某一个 q 点上的每一个模式到反折叠后的某个 q 点的投影，投射到单胞中的对应 q 点的模式上

  struct Config
  {
    // 在单胞内取几个平面波的基矢
    Eigen::Vector<std::size_t, 3> PrimativeCellBasisNumber;
    // 要选定的超胞中的 q 点，sub qpoint，以及单胞中的 q 点
    std::size_t SuperQpointIndex, SubQpointIndex, PrimativeQpointIndex;
    // 输入输出文件的路径，直接输出 yaml，因为这些数据不再需要后处理
    std::string InputDataFile, OutputFile;
  };
  struct Output
  {
    Eigen::Vector3d SuperQpoint, SubQpoint, PrimativeQpoint;
    std::vector<std::vector<double>> Coefficient;
  };

  biu::Logger::Guard log(config_file);
  auto config = YAML::LoadFile(config_file).as<Config>();
  auto input = biu::deserialize<CommonData>
    (biu::read<std::byte>(config.InputDataFile));
  Output output
  {
    .SuperQpoint = input.Super.Qpoint[config.SuperQpointIndex].Qpoint,
    .SubQpoint = input.Super.Qpoint[config.SuperQpointIndex].SubQpoint[config.SubQpointIndex],
    .PrimativeQpoint = input.Primative.Qpoint[config.PrimativeQpointIndex].Qpoint,
    .Coefficient = std::vector<std::vector<double>>(input.Super.Qpoint[config.SuperQpointIndex].Mode.size(),
      std::vector<double>(input.Primative.Qpoint[config.PrimativeQpointIndex].Mode.size()))
  };
  if ((output.SubQpoint - output.PrimativeQpoint).norm() > 0.01)
    log.error("sub qpoint {} != primative qpoint {}"_f(output.SubQpoint, output.PrimativeQpoint));

  // 基在原胞中的表示
  auto primative_basis = [&]
  {
    biu::Logger::Guard log;
    std::vector<Eigen::VectorXcd> basis(config.PrimativeCellBasisNumber.prod());
    for (auto [xyz_of_basis, i_of_basis]
      : biu::sequence(config.PrimativeCellBasisNumber))
    {
      // 波矢就是 xyz_of_basis, 单位为单胞的倒格矢
      // 将它的单位转换成埃^-1
      auto wavevector = input.Primative.Cell.inverse() * xyz_of_basis.cast<double>();
      log.debug("xyz_of_basis: {}"_f(xyz_of_basis));
      log.debug("wavevector: {:.3f}"_f(fmt::join(wavevector, ", ")));
      // 将原子坐标单位也转换为埃
      auto atom_position = input.Primative.AtomPosition * input.Primative.Cell;
      // 计算基矢
      basis[i_of_basis] = (2i * std::numbers::pi_v<double> * (atom_position * wavevector)).array().exp();
    }
    return basis;
  }();
  // 基在超胞中的表示
  auto super_basis = [&]
  {
    biu::Logger::Guard log;
    std::vector<Eigen::VectorXcd> basis(config.PrimativeCellBasisNumber.prod());
    auto diff_of_sub_qpoint_by_reciprocal_modified_super_cell = [&]
    {
      for
      (
        auto [diff_of_sub_qpoint_by_reciprocal_modified_super_cell, i_of_sub_qpoint]
          : biu::sequence(input.Super.CellMultiplier)
      )
        if (i_of_sub_qpoint == config.SubQpointIndex) return diff_of_sub_qpoint_by_reciprocal_modified_super_cell;
      std::unreachable();
    }();
    log.debug("diff_of_sub_qpoint_by_reciprocal_modified_super_cell: {}"_f
      (diff_of_sub_qpoint_by_reciprocal_modified_super_cell));
    for (auto [xyz_of_basis, i_of_basis]
      : biu::sequence(config.PrimativeCellBasisNumber))
    {
      // 计算波矢, 单位为单胞的倒格矢
      auto wavevector_by_reciprocal_primative_cell = xyz_of_basis.cast<double>()
        + input.Super.CellMultiplier.cast<double>().cwiseInverse().asDiagonal()
        * diff_of_sub_qpoint_by_reciprocal_modified_super_cell.cast<double>();
      // 将单位转换为埃^-1
      auto wavevector = input.Primative.Cell.inverse() * wavevector_by_reciprocal_primative_cell;
      // 将原子坐标单位也转换为埃
      auto atom_position = input.Super.AtomPosition
        * (input.Super.CellDeformation * input.Super.CellMultiplier.cast<double>().asDiagonal() * input.Primative.Cell);
      // 计算基矢
      basis[i_of_basis] = (2i * std::numbers::pi_v<double> * (atom_position * wavevector)).array().exp();
    }
    return basis;
  }();
  // 将所有模式投影到基矢上
  auto primative_projection_coefficient =
    input.Primative.Qpoint[config.PrimativeQpointIndex].Mode
      | ranges::views::transform([&](const auto& mode)
        {
          return primative_basis
            | ranges::views::transform([&](const auto& basis)
              { return (basis.transpose().conjugate() * mode.EigenVector).eval() | biu::fromEigen; })
            | ranges::to_vector
            | biu::toEigen<>;
        })
      | ranges::to_vector;
  for (auto i_of_primative_mode : std::views::iota(0u, primative_projection_coefficient.size()))
    for (auto i_of_direction : std::array{0, 1, 2})
    {
      if (i_of_primative_mode != 0 && i_of_primative_mode != 11) continue;
      log.debug("primitive mode: {} {} {:.3f}"_f
        (i_of_primative_mode, i_of_direction, fmt::join(primative_projection_coefficient[i_of_primative_mode].col(i_of_direction).array().abs(), ", ")));
      log.debug("primitive mode: {} {} {:.3f}"_f
        (i_of_primative_mode, i_of_direction, fmt::join(primative_projection_coefficient[i_of_primative_mode].col(i_of_direction).array().arg(), ", ")));
    }
  log.debug("primitive 0 {}"_f(primative_projection_coefficient[0].col(2).transpose().conjugate() * primative_projection_coefficient[0].col(2)));
  log.debug("primitive 11 {}"_f(primative_projection_coefficient[11].col(2).transpose().conjugate() * primative_projection_coefficient[11].col(2)));
  log.debug("0 x 11 {}"_f(primative_projection_coefficient[11].col(2).transpose().conjugate() * primative_projection_coefficient[0].col(2)));
  auto super_projection_coefficient =
    input.Super.Qpoint[config.SuperQpointIndex].Mode
      | ranges::views::transform([&](const auto& mode)
        {
          return super_basis
            | ranges::views::transform([&](const auto& basis)
              { return (basis.transpose().conjugate() * mode.EigenVector).eval() | biu::fromEigen; })
            | ranges::to_vector
            | biu::toEigen<>;
        })
      | ranges::to_vector;
  log.debug("super 2 x primitive 0 {}"_f(super_projection_coefficient[2].col(2).transpose().conjugate() * primative_projection_coefficient[0].col(2)));
  log.debug("super 2 x primitive 11 {}"_f(super_projection_coefficient[2].col(2).transpose().conjugate() * primative_projection_coefficient[11].col(2)));
  log.debug("super 2 x super 2 {}"_f(super_projection_coefficient[2].col(2).transpose().conjugate() * super_projection_coefficient[2].col(2)));
  for (auto i_of_direction : std::array{0, 1, 2})
  {
    log.debug("super mode: 2 {} {:.3f}"_f
      (i_of_direction, fmt::join(super_projection_coefficient[2].col(i_of_direction).array().abs(), ", ")));
    log.debug("super mode: 2 {} {:.3f}"_f
      (i_of_direction, fmt::join(super_projection_coefficient[2].col(i_of_direction).array().arg(), ", ")));
  }
  // 将它们互相投影，得到结果
  for (auto [i_of_super_mode, super_mode] : ranges::views::enumerate(super_projection_coefficient))
  {
    for (auto [i_of_primative_mode, primative_mode]
      : ranges::views::enumerate(primative_projection_coefficient))
      for (auto i_of_direction : std::array{0, 1, 2})
        output.Coefficient[i_of_super_mode][i_of_primative_mode] +=
          (super_mode.col(i_of_direction).transpose().conjugate() * primative_mode.col(i_of_direction)).norm();
    auto sum = ranges::accumulate(output.Coefficient[i_of_super_mode], 0.);
    for (auto& coefficient : output.Coefficient[i_of_super_mode])
      coefficient *=
        input.Super.Qpoint[config.SuperQpointIndex].Mode[i_of_super_mode].WeightOnUnfold[config.SubQpointIndex] / sum;
  }
  // 输出
  std::ofstream(config.OutputFile) << YAML::Node(output);
}
