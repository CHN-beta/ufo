# include <ufo/ufo.impl.hpp>

int main(int argc, const char** argv)
{
  if (argc != 2)
    throw std::runtime_error(fmt::format("Usage: {} config.yaml", argv[0]));

  std::cerr << "Reading input file..." << std::flush;
  Input input(argv[1]);
  std::cerr << "Done." << std::endl;

  // 反折叠的原理: 将超胞中的原子运动状态, 投影到一组平面波构成的基矢中.
  // 每一个平面波的波矢由两部分相加得到: 一部分是单胞倒格子的整数倍, 所取的个数有一定任意性, 论文中建议取大约单胞中原子个数那么多个;
  //  对于没有缺陷的情况, 取一个应该就足够了.
  // 另一部分是超胞倒格子的整数倍, 取 n 个, n 为超胞对应的单胞的倍数, 其实也就是倒空间中单胞对应倒格子中超胞的格点.
  // 只要第一部分取得足够多, 那么单胞中原子的状态就可以完全被这些平面波描述.
  // 将超胞中原子的运动状态投影到这些基矢上, 计算出投影的系数, 就可以将超胞的原子运动状态分解到单胞中的多个 q 点上.

  // 构建基
  // 每个 q 点对应的一组 sub qpoint。不同的 q 点所对应的 sub qpoint 是不一样的，但 sub qpoint 与 q 点的相对位置一致。
  // 这里 xyz_of_diff_of_sub_qpoint 即表示这个相对位置。
  // 由于基只与这个相对位置有关（也就是说，不同 q 点的基是一样的），因此可以先计算出所有的基，这样降低计算量。
  // 外层下标对应超胞倒格子的整数倍那部分(第二部分), 也就是不同的 sub qpoint
  // 内层下标对应单胞倒格子的整数倍那部分(第一部分), 也就是 sub qpoint 上的不同平面波（取的数量越多，结果越精确）
  std::cerr << "Calculating basis..." << std::flush;
  std::vector<std::vector<Eigen::VectorXcd>> basis(input.SuperCellMultiplier.prod());
  // 每个 q 点对应的一组 sub qpoint。不同的 q 点所对应的 sub qpoint 是不一样的，但 sub qpoint 与 q 点的相对位置一致。
  // 这里 xyz_of_diff_of_sub_qpoint 即表示这个相对位置，单位为超胞的倒格矢
  for (auto [xyz_of_diff_of_sub_qpoint_by_reciprocal_modified_super_cell, i_of_sub_qpoint]
    : triplet_sequence(input.SuperCellMultiplier))
  {
    basis[i_of_sub_qpoint].resize(input.PrimativeCellBasisNumber.prod());
    for (auto [xyz_of_basis, i_of_basis] : triplet_sequence(input.PrimativeCellBasisNumber))
    {
      // 计算 q 点的坐标, 单位为单胞的倒格矢
      auto diff_of_sub_qpoint_by_reciprocal_primative_cell = xyz_of_basis.cast<double>()
        + input.SuperCellMultiplier.cast<double>().cwiseInverse().asDiagonal()
        * xyz_of_diff_of_sub_qpoint_by_reciprocal_modified_super_cell.cast<double>();
      // 将 q 点坐标转换为埃^-1
      auto qpoint = (diff_of_sub_qpoint_by_reciprocal_primative_cell.transpose()
        * (input.PrimativeCell.transpose().inverse())).transpose();
      // 计算基矢
      basis[i_of_sub_qpoint][i_of_basis]
        = (2i * std::numbers::pi_v<double> * (input.AtomPosition * qpoint)).array().exp();
    }
  }
  std::cerr << "Done." << std::endl;

  // 计算投影的结果
  // 最外层下标对应反折叠前的 q 点, 第二层下标对应不同模式, 第三层下标对应这个模式在反折叠后的 q 点(sub qpoint)
  std::vector<std::vector<std::vector<double>>> projection_coefficient(input.QPointData.size());
  std::atomic<unsigned> finished_qpoint(0);
  // 对每个 q 点并行
  std::transform
  (
    std::execution::par, input.QPointData.begin(), input.QPointData.end(),
    projection_coefficient.begin(), [&](const auto& qpoint_data)
  {
    std::osyncstream(std::cerr) << fmt::format("\rCalculating projection coefficient...({}/{})",
      finished_qpoint, input.QPointData.size()) << std::flush;
    std::vector<std::vector<double>> projection_coefficient(qpoint_data.ModeData.size());
    // 这里, qpoint_data 和 projection_coefficient 均指对应于一个 q 点的数据
    for (unsigned i_of_mode = 0; i_of_mode < qpoint_data.ModeData.size(); i_of_mode++)
    {
      auto& _ = projection_coefficient[i_of_mode];
      _.resize(input.SuperCellMultiplier.prod());
      for (unsigned i_of_sub_qpoint = 0; i_of_sub_qpoint < input.SuperCellMultiplier.prod(); i_of_sub_qpoint++)
        // 对于 basis 中, 对应于单胞倒格子的部分, 以及对应于不同方向的部分, 分别求内积, 然后求模方和
        for (unsigned i_of_basis = 0; i_of_basis < input.PrimativeCellBasisNumber.prod(); i_of_basis++)
          _[i_of_sub_qpoint] +=
            (basis[i_of_sub_qpoint][i_of_basis].transpose().conjugate() * qpoint_data.ModeData[i_of_mode].AtomMovement)
              .array().abs2().sum();
      // 如果是严格地将向量分解到一组完备的基矢上, 那么不需要对计算得到的权重再做归一化处理
      // 但这里并不是这样一个严格的概念. 因此对分解到各个 sub qpoint 上的权重做归一化处理
      auto sum = std::accumulate(_.begin(), _.end(), 0.);
      for (auto& __ : _)
        __ /= sum;
    }
    finished_qpoint++;
    std::osyncstream(std::cerr) << fmt::format("\rCalculating projection coefficient...({}/{})",
      finished_qpoint, input.QPointData.size()) << std::flush;
    return projection_coefficient;
  });
  std::cerr << "Done." << std::endl;

  // 填充输出对象
  std::cerr << "Filling output object..." << std::flush;
  Output output;
  for (unsigned i_of_qpoint = 0; i_of_qpoint < input.QPointData.size(); i_of_qpoint++)
  {
    // 当 SuperCellDeformation 不是单位矩阵时, input.QPointData[i_of_qpoint].QPoint 不一定在 reciprocal_primative_cell 中
    // 需要首先将 q 点平移数个周期, 进入不包含 SuperCellDeformation 的超胞的倒格子中
    auto qpoint_by_reciprocal_modified_super_cell_in_modified_reciprocal_super_cell
      = !input.SuperCellDeformation ? input.QPointData[i_of_qpoint].QPoint : [&]
      {
        auto current_qpoint = input.QPointData[i_of_qpoint].QPoint;
        // 给一个 q 点打分
        // 计算这个 q 点以 modified_reciprocal_supre_cell 为单位的坐标, 依次考虑每个维度, 总分为每个维度之和.
        // 如果这个坐标大于 0 小于 1, 则打 0 分.
        // 如果这个坐标小于 0, 则打这个坐标的相反数分.
        // 如果这个坐标大于 1, 则打这个坐标减去 1 的分.
        auto score = [&](Eigen::Vector3d qpoint_by_reciprocal_super_cell)
        {
          // SuperCell = SuperCellDeformation * SuperCellMultiplier.asDiagonal() * PrimativeCell
          // ModifiedSuperCell = SuperCellMultiplier.asDiagonal() * PrimativeCell
          // ReciprocalSuperCell = SuperCell.inverse().transpose()
          // ReciprocalModifiedSuperCell = ModifiedSuperCell.inverse().transpose()
          // qpoint.transpose() = qpoint_by_reciprocal_super_cell.transpose() * ReciprocalSuperCell
          // qpoint.transpose() = qpoint_by_reciprocal_modified_super_cell.transpose() * ReciprocalModifiedSuperCell
          auto qpoint_by_reciprocal_modified_super_cell =
            (input.SuperCellDeformation->inverse() * qpoint_by_reciprocal_super_cell).eval();
          double score = 0;
          for (unsigned i = 0; i < 3; i++)
          {
            auto coordinate = qpoint_by_reciprocal_modified_super_cell[i];
            if (coordinate < 0)
              score -= coordinate;
            else if (coordinate > 1)
              score += coordinate - 1;
          }
          return score;
        };
        while (score(current_qpoint) > 0)
        {
          double min_score = std::numeric_limits<double>::max();
          Eigen::Vector3d min_score_qpoint;
          for (int x = -1; x <= 1; x++)
            for (int y = -1; y <= 1; y++)
              for (int z = -1; z <= 1; z++)
              {
                auto this_qpoint = (current_qpoint
                  + Eigen::Matrix<int, 3, 1>{{x}, {y}, {z}}.cast<double>()).eval();
                auto this_score = score(this_qpoint);
                if (this_score < min_score)
                {
                  min_score = this_score;
                  min_score_qpoint = this_qpoint;
                }
              }
          current_qpoint = min_score_qpoint;
        }
        return input.SuperCellDeformation->inverse() * current_qpoint;
      }();
    for (auto [xyz_of_diff_of_sub_qpoint_by_reciprocal_modified_super_cell, i_of_sub_qpoint]
      : triplet_sequence(input.SuperCellMultiplier))
    {
      auto& _ = output.QPointData.emplace_back();
      auto reciprocal_modified_super_cell =
        (input.SuperCellMultiplier.cast<double>().asDiagonal() * input.PrimativeCell).inverse().transpose();
      // sub qpoint 的坐标，单位为埃^-1
      auto sub_qpoint = ((xyz_of_diff_of_sub_qpoint_by_reciprocal_modified_super_cell.cast<double>()
        + qpoint_by_reciprocal_modified_super_cell_in_modified_reciprocal_super_cell)
        .transpose() * reciprocal_modified_super_cell).transpose();
      // 将坐标转换为相对于单胞的倒格矢的坐标并写入
      // 由 sub_qpoint.transpose() = sub_qpoint_by_reciprocal_primative_cell.transpose()
      //  * PrimativeCell.transpose().inverse()
      // 得到 sub_qpoint_by_reciprocal_primative_cell = PrimativeCell * sub_qpoint
      _.QPoint = input.PrimativeCell * sub_qpoint;
      _.Source = input.QPointData[i_of_qpoint].QPoint;
      _.SourceIndex_ = i_of_qpoint;
      for (unsigned i_of_mode = 0; i_of_mode < input.QPointData[i_of_qpoint].ModeData.size(); i_of_mode++)
      {
        auto& __ = _.ModeData.emplace_back();
        __.Frequency = input.QPointData[i_of_qpoint].ModeData[i_of_mode].Frequency;
        __.Weight = projection_coefficient[i_of_qpoint][i_of_mode][i_of_sub_qpoint];
      }
    }
  }
  std::cerr << "Done." << std::endl;

  std::cerr << "Writing output file..." << std::flush;
  for (auto& output_file : input.QPointDataOutputFile)
    output.write(output_file.FileName, output_file.Format);
  std::cerr << "Done." << std::endl;
}
