# include <ufo.hpp>

void ufo::project_to_mode(std::string config_file)
{
  // 将选定的超胞中的某一个 q 点上的每一个模式到反折叠后的某个 q 点的投影，投射到单胞中的对应 q 点的模式上
  // 与 unfold 的方法类似，只有一个重要的区别：
  // 组成基矢的平面波的波矢的两个部分中，第二个部分不再需要。因为这个部分的结果已经算出来了，只需要将结果乘上 WeightOnUnfold 即可。
  struct Config
  {
    // 在单胞内取几个平面波的基矢
    Eigen::Vector3i PrimativeCellBasisNumber;
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
    std::vector<Eigen::VectorXcd> basis((config.PrimativeCellBasisNumber.array() * 2 - 1).prod());
    for (auto [xyz_of_basis, i_of_basis] : biu::sequence
      ((-config.PrimativeCellBasisNumber + Eigen::Vector3i::Ones()).eval(), config.PrimativeCellBasisNumber))
    {
      // 波矢就是 xyz_of_basis, 单位为单胞的倒格矢
      // 将它的单位转换成埃^-1
      auto wavevector =
        (xyz_of_basis.transpose().cast<double>() * input.Primative.Cell.inverse().transpose()).transpose();
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
    std::vector<Eigen::VectorXcd> basis((config.PrimativeCellBasisNumber.array() * 2 - 1).prod());
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
    for (auto [xyz_of_basis, i_of_basis] : biu::sequence
      ((-config.PrimativeCellBasisNumber + Eigen::Vector3i::Ones()).eval(), config.PrimativeCellBasisNumber))
    {
      // 计算波矢, 单位为单胞的倒格矢
      auto wavevector_by_reciprocal_primative_cell = xyz_of_basis.cast<double>()
        + input.Super.CellMultiplier.cast<double>().cwiseInverse().asDiagonal()
        * diff_of_sub_qpoint_by_reciprocal_modified_super_cell.cast<double>();
      // 将单位转换为埃^-1
      auto wavevector =
        (wavevector_by_reciprocal_primative_cell.transpose() * input.Primative.Cell.inverse().transpose()).transpose();
      // 将原子坐标单位也转换为埃
      auto atom_translation =
        (input.Super.AtomTranslation.value_or(std::array{0., 0., 0.}) | biu::toEigen<>)
        .transpose().replicate(input.Super.AtomPosition.rows(), 1).eval();
      auto atom_position = (input.Super.AtomPosition - atom_translation)
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
  // 将它们互相投影，得到结果
  for (auto [i_of_super_mode, super_mode] : ranges::views::enumerate(super_projection_coefficient))
  {
    for (auto [i_of_primative_mode, primative_mode]
      : ranges::views::enumerate(primative_projection_coefficient))
      output.Coefficient[i_of_super_mode][i_of_primative_mode] =
        std::norm((super_mode.transpose().conjugate() * primative_mode).trace());
    auto sum = ranges::accumulate(output.Coefficient[i_of_super_mode], 0.);
    for (auto& coefficient : output.Coefficient[i_of_super_mode])
      coefficient *=
        input.Super.Qpoint[config.SuperQpointIndex].Mode[i_of_super_mode].WeightOnUnfold[config.SubQpointIndex] / sum;
  }
  // 输出
  std::ofstream(config.OutputFile) << YAML::Node(output);
}
