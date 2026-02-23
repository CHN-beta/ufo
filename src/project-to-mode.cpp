# include <ufo.hpp>
# include <execution>

void ufo::project_to_mode(std::string config_file)
{
  // 将选定的超胞中的某一个 q 点上的每一个模式到反折叠后的某个 q 点的投影，投射到单胞中的对应 q 点的模式上
  // 不使用平面波方法，而是匹配对应的原子然后做内积
  struct Config
  {
    // 容许的原子位置偏差，单位为埃
    double AtomPositionTolerance;
    // q 点位置的容许偏差，单位为埃^-1
    double QpointTolerance;
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
  if (((output.SubQpoint - output.PrimativeQpoint).transpose() * input.Primative.Cell.inverse().transpose()).norm()
    > config.QpointTolerance)
    log.error("sub qpoint {} != primative qpoint {}"_f(output.SubQpoint, output.PrimativeQpoint));

  auto atom_is_near = [&](Eigen::Vector3d atom_in_super, Eigen::Vector3d atom_in_primative) -> bool
  {
    // 如果两个原子经过变换后正好重合，那么应该有这个关系：
    // (atom_in_super - AtomTranslation) * CellDeformation * CellMultiplier * PrimativeCell
    //    = atom_in_primative * PrimativeCell
    // 即:
    // (atom_in_super - AtomTranslation) * CellDeformation * CellMultiplier - atom_in_primative = 0
    // 但这里要考虑两个因素：
    // 首先，两个原子可能会差整数个 PrimativeCell 的格矢，所以需要对上面的结果加或者减一个整数，使得结果落在正负 0.5 范围内
    // 其次，我们允许一定误差，因此上面的数值乘以 PrimativeCell 的模应该小于 AtomPositionTolerance
    biu::Logger::Guard log(atom_in_super.transpose(), atom_in_primative.transpose());
    auto translation = input.Super.AtomTranslation.value_or(std::array{0., 0., 0.})
      | biu::toEigen<>;
    auto diff = ((atom_in_super - translation).transpose() * input.Super.CellDeformation
      * input.Super.CellMultiplier.cast<double>().asDiagonal() - atom_in_primative.transpose()).eval();
    for (auto i = 0; i < 3; i++) diff[i] -= std::round(diff[i]);
    log.debug("Diff: {}"_f(diff.transpose()));
    return (diff * input.Primative.Cell).norm() < config.AtomPositionTolerance;
  };

  // 匹配原子
  int number_of_matched_atoms = 0, number_of_unmatched_atoms = 0;
  auto type_in_super = input.Super.AtomType
    | std::views::transform([](const auto& pair) { return std::views::repeat(pair.first, pair.second); })
    | std::views::join
    | std::ranges::to<std::vector<std::string>>();
  auto type_in_primative = input.Primative.AtomType
    | std::views::transform([](const auto& pair) { return std::views::repeat(pair.first, pair.second); })
    | std::views::join
    | std::ranges::to<std::vector<std::string>>();
  std::vector<std::optional<int>> matched_atoms_in_super(input.Super.AtomPosition.rows());
  for (auto i_of_super_atom : std::views::iota(0, input.Super.AtomPosition.rows()))
  {
    auto result = std::views::iota(0, input.Primative.AtomPosition.rows())
      | std::views::filter([&](auto i_of_primative_atom)
        {
          return type_in_super[i_of_super_atom] == type_in_primative[i_of_primative_atom]
            && atom_is_near
              (input.Super.AtomPosition.row(i_of_super_atom), input.Primative.AtomPosition.row(i_of_primative_atom));
        })
      | std::ranges::to<std::vector<int>>();
    if (result.empty())
    {
      log.debug("Unmatched atom in super cell: {} {} {}"_f
      (
        i_of_super_atom, type_in_super[i_of_super_atom],
        input.Super.AtomPosition.row(i_of_super_atom).transpose()
      ));
      number_of_unmatched_atoms++;
    }
    else if (result.size() > 1)
    {
      log.error("Multiple matched atoms in primative cell: {} {} {}"_f
      (
        i_of_super_atom, type_in_super[i_of_super_atom],
        input.Super.AtomPosition.row(i_of_super_atom).transpose()
      ));
      number_of_unmatched_atoms++;
    }
    else
    {
      matched_atoms_in_super[i_of_super_atom] = result[0];
      number_of_matched_atoms++;
    }
  }
  log.info("Matched {} atoms, unmatched {} atoms"_f(number_of_matched_atoms, number_of_unmatched_atoms));

  // 对每个模式进行投影
  std::for_each_n
  (std::execution::par_unseq,
    std::views::iota(0).begin(),
    input.Super.Qpoint[config.SuperQpointIndex].Mode.size(),
    [&](auto i_of_mode)
    {
      auto& super_mode = input.Super.Qpoint[config.SuperQpointIndex].Mode[i_of_mode];
      for (auto i_of_primative_mode
        : std::views::iota(0u, input.Primative.Qpoint[config.PrimativeQpointIndex].Mode.size()))
      {
        auto& primative_mode = input.Primative.Qpoint[config.PrimativeQpointIndex].Mode[i_of_primative_mode];
        output.Coefficient[i_of_mode][i_of_primative_mode] = std::norm(*std::ranges::fold_left_first
        (
          std::views::iota(0, input.Super.AtomPosition.rows())
            | std::views::filter([&](auto i_of_super_atom)
              { return matched_atoms_in_super[i_of_super_atom].has_value(); })
            | std::views::transform([&](auto i_of_super_atom) -> std::complex<double>
              {
                return super_mode.EigenVector.row(i_of_super_atom)
                  * primative_mode.EigenVector.row(*matched_atoms_in_super[i_of_super_atom]).transpose().conjugate();
              }),
          std::plus{}
        ));
      }
      // 将结果进行放缩，使得系数之和为 input.Super.Qpoint[config.SuperQpointIndex].Mode[i_of_mode]
      //    .WeightOnUnfold[config.SubQpointIndex]
      double ratio;
      if (auto sum = *std::ranges::fold_left_first(output.Coefficient[i_of_mode], std::plus{}); sum == 0)
        ratio = 0;
      else
        ratio = input.Super.Qpoint[config.SuperQpointIndex].Mode[i_of_mode].WeightOnUnfold[config.SubQpointIndex] / sum;
      for (auto& coefficient : output.Coefficient[i_of_mode])
        coefficient *= ratio;
    }
  );
  // 输出
  std::ofstream(config.OutputFile) << YAML::Node(output);
}
