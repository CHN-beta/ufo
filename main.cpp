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
# include <yaml-cpp/yaml.h>
# include <eigen3/Eigen/Dense>
# include <concurrencpp/concurrencpp.h>
# include <fmt/format.h>
# include <highfive/H5File.hpp>

using namespace std::literals;

struct PhonopyComplex { double r, i; };
HighFive::CompoundType create_compound_complex()
  { return {{"r", HighFive::AtomicType<double>{}}, {"i", HighFive::AtomicType<double>{}}}; }
HIGHFIVE_REGISTER_TYPE(PhonopyComplex, create_compound_complex)

// 在相位中, 约定为使用 $\exp (2 \pi i \vec{q} \cdot \vec{r})$ 来表示原子的运动状态
//  (而不是 $\exp (-2 \pi i \vec{q} \cdot \vec{r})$)
// 一些书定义的倒格矢中包含了 $2 \pi$ 的部分, 我们这里约定不包含这部分.
//  也就是说, 正格子与倒格子的转置相乘, 得到单位矩阵.

struct Input
{
  // 单胞的三个格矢，每行表示一个格矢的坐标，单位为埃
  Eigen::Matrix3d PrimativeCell;
  // 单胞到超胞的格矢转换时用到的矩阵
  // SuperCellMultiplier 是一个三维列向量且各个元素都是整数，表示单胞在各个方向扩大到多少倍之后，可以得到和超胞一样的体积
  // SuperCellDeformation 是一个行列式为 1 的矩阵，它表示经过 SuperCellMultiplier 扩大后，还需要怎样的变换才能得到超胞
  // SuperCell = (SuperCellDeformation * SuperCellMultiplier.asDiagonal()) * PrimativeCell
  // ReciprocalPrimativeCell = (SuperCellDeformation * SuperCellMultiplier.asDiagonal()).transpose()
  //  * ReciprocalSuperCell
  // Position = PositionToCell(line vector) * Cell
  // InversePosition = InversePositionToCell(line vector) * ReciprocalCell
  // PositionToSuperCell(line vector) * SuperCell = PositionToPrimativeCell(line vector) * PrimativeCell
  // ReciprocalPositionToSuperCell(line vector) * ReciprocalSuperCell
  //  = ReciprocalPositionToPrimativeCell(line vector) * ReciprocalPrimativeCell
  Eigen::Vector<unsigned, 3> SuperCellMultiplier;
  Eigen::Matrix<double, 3, 3> SuperCellDeformation;
  // 在单胞内取几个平面波的基矢
  Eigen::Vector<unsigned, 3> PrimativeCellBasisNumber;
  // 超胞中原子的坐标，每行表示一个原子的坐标，单位为埃
  Eigen::MatrixX3d AtomPosition;

  // 是否调整输出结果, 使得结果中的模式适合人类阅读. 默认为 true.
  // 这包括合并相近的模式, 去除权重过小的模式, 限制输出的小数位数.
  // 如果想用结果来进一步画图, 则建议关闭.
  std::optional<bool> Filter;

  // 关于各个 Q 点的数据
  struct QPointDataType_
  {
    // Q 点的坐标，单位为超胞的倒格矢
    Eigen::Vector3d QPoint;

    // 关于这个 Q 点上各个模式的数据
    struct ModeDataType_
    {
      // 模式的频率，单位为 THz
      double Frequency;
      // 模式中各个原子的运动状态
      // 这个数据是这样得到的: phonopy 输出的动态矩阵的 eigenvector 乘以 $\exp(-2 \pi i \vec q \cdot \vec r)$
      // 这个数据可以认为是原子位移中, 关于超胞有周期性的那一部分, 再乘以原子质量的开方.
      // 这个数据在读入后不会被立即归一化.
      Eigen::MatrixX3cd AtomMovement;
    };
    std::vector<ModeDataType_> ModeData;
  };
  std::vector<QPointDataType_> QPointData;

  Input(std::string yaml_file, std::optional<std::string> hdf5_file);
};

struct Output
{
  // 关于各个 Q 点的数据
  struct QPointDataType_
  {
    // Q 点的坐标，单位为单胞的倒格矢
    Eigen::Vector3d QPoint;

    // 来源于哪个 Q 点, 单位为超胞的倒格矢
    Eigen::Vector3d Source;

    // 关于这个 Q 点上各个模式的数据
    struct ModeDataType_
    {
      // 模式的频率，单位为 THz
      double Frequency;
      // 模式的权重
      double Weight;
    };
    std::vector<ModeDataType_> ModeData;
  };
  std::vector<QPointDataType_> QPointData;
};

concurrencpp::generator<std::pair<Eigen::Vector<unsigned, 3>, unsigned>>
  triplet_sequence(Eigen::Vector<unsigned, 3> range)
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

int main(int argc, const char** argv)
{
  if (argc < 3 || argc > 5)
    throw std::runtime_error("Usage: " + std::string(argv[0]) + " input.yaml [input.hdf5] output.yaml");

  std::cerr << "Reading input file..." << std::flush;
  Input input(argv[1], argc > 3 ? std::make_optional(argv[2]) : std::nullopt);
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
    auto qpoint_by_reciprocal_super_cell_in_modified_reciprocal_super_cell = [&]
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
          (input.SuperCellDeformation.inverse() * qpoint_by_reciprocal_super_cell).eval();
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
      return current_qpoint;
    }();
    for (auto [xyz_of_diff_of_sub_qpoint_by_reciprocal_modified_super_cell, i_of_sub_qpoint]
      : triplet_sequence(input.SuperCellMultiplier))
    {
      auto& _ = output.QPointData.emplace_back();
      // 这一步推导过程在计算 score 的函数中
      auto qpoint_by_reciprocal_modified_super_cell_in_modified_reciprocal_super_cell =
        input.SuperCellDeformation.inverse() * qpoint_by_reciprocal_super_cell_in_modified_reciprocal_super_cell;
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
      if (input.Filter.value_or(true))
      {
        // 从小到大枚举所有的模式，并将相近的模式（相差小于 0.1 THz）合并
        std::map<double, double> frequency_to_weight;
        for (unsigned i_of_mode = 0; i_of_mode < input.QPointData[i_of_qpoint].ModeData.size(); i_of_mode++)
        {
          auto frequency = input.QPointData[i_of_qpoint].ModeData[i_of_mode].Frequency;
          auto weight = projection_coefficient[i_of_qpoint][i_of_mode][i_of_sub_qpoint];
          auto it_lower = frequency_to_weight.lower_bound(frequency - 0.1);
          auto it_upper = frequency_to_weight.upper_bound(frequency + 0.1);
          if (it_lower == it_upper)
            frequency_to_weight[frequency] = weight;
          else
          {
            auto frequency_sum = std::accumulate(it_lower, it_upper, 0.,
              [](const auto& a, const auto& b) { return a + b.first * b.second; });
            auto weight_sum = std::accumulate(it_lower, it_upper, 0.,
              [](const auto& a, const auto& b) { return a + b.second; });
            frequency_sum += frequency * weight;
            weight_sum += weight;
            frequency_to_weight.erase(it_lower, it_upper);
            frequency_to_weight[frequency_sum / weight_sum] = weight_sum;
          }
        }
        // 仅保留权重大于 0.1 的模式
        for (auto& mode : frequency_to_weight)
          if (mode.second > 0.1)
          {
            auto& __ = _.ModeData.emplace_back();
            __.Frequency = mode.first;
            __.Weight = mode.second;
          }
      }
      else
        for (unsigned i_of_mode = 0; i_of_mode < input.QPointData[i_of_qpoint].ModeData.size(); i_of_mode++)
        {
          auto& __ = _.ModeData.emplace_back();
          __.Frequency = input.QPointData[i_of_qpoint].ModeData[i_of_mode].Frequency;
          __.Weight = projection_coefficient[i_of_qpoint][i_of_mode][i_of_sub_qpoint];
        }
    }
  }
  std::cerr << "Done." << std::endl;

  // YAML 输出得太丑了，我来自己写
  std::cerr << "Writing output file..." << std::flush;
  std::ofstream(argc > 3 ? argv[3] : argv[2]) << [&]
  {
    std::stringstream print;
    auto format = input.Filter.value_or(true) ? 3 : 10;
    print << "QPointData:\n";
    for (auto& qpoint: output.QPointData)
    {
      print << fmt::format("  - QPoint: [ {1:.{0}f}, {2:.{0}f}, {3:.{0}f} ]\n",
        format, qpoint.QPoint[0], qpoint.QPoint[1], qpoint.QPoint[2]);
      print << fmt::format("    Source: [ {1:.{0}f}, {2:.{0}f}, {3:.{0}f} ]\n",
        format, qpoint.Source[0], qpoint.Source[1], qpoint.Source[2]);
      print << "    ModeData:\n";
      for (auto& mode: qpoint.ModeData)
        print << fmt::format("      - {{ Frequency: {1:.{0}f}, Weight: {2:.{0}f} }}\n",
          format, mode.Frequency, mode.Weight);
    }
    return print.str();
  }();
  std::cerr << "Done." << std::endl;
}

// 从文件中读取输入, 文件中应当包含: (大多数据可以直接从 phonopy 的输出中复制)
// 单胞的格矢: lattice 单位为埃 直接从 phonopy 的输出中复制
// 超胞的倍数: SuperCellMultiplier 手动输入, 为一个包含三个整数的数组
// 平面波的基矢个数: PrimativeCellBasisNumber 手动输入, 为一个包含三个整数的数组
// 超胞中原子的坐标: points[*].coordinates 单位为超胞的格矢 直接从 phonopy 的输出中复制
// 各个 Q 点的坐标: phonon[*].q-position 单位为超胞的倒格子的格矢 直接从 phonopy 的输出中复制
// 各个模式的频率: phonon[*].band[*].frequency 单位为 THz 直接从 phonopy 的输出中复制
// 各个模式的原子运动状态: phonon[*].band[*].eigenvector 直接从 phonopy 的输出中复制
// 文件中可以有多余的项目, 多余的项目不管.
Input::Input(std::string yaml_file, std::optional<std::string> hdf5_file)
{
  auto node = YAML::LoadFile(yaml_file);
  for (unsigned i = 0; i < 3; i++)
    for (unsigned j = 0; j < 3; j++)
      PrimativeCell(i, j) = node["lattice"][i][j].as<double>();

  for (unsigned i = 0; i < 3; i++)
    SuperCellMultiplier(i) = node["SuperCellMultiplier"][i].as<int>();

  for (unsigned i = 0; i < 3; i++)
    for (unsigned j = 0; j < 3; j++)
      SuperCellDeformation(i, j) = node["SuperCellDeformation"][i][j].as<double>();

  for (unsigned i = 0; i < 3; i++)
    PrimativeCellBasisNumber(i) = node["PrimativeCellBasisNumber"][i].as<int>();

  if (auto value = node["Filter"])
    Filter = value.as<bool>();

  auto points = node["points"].as<std::vector<YAML::Node>>();
  auto atom_position_to_super_cell = Eigen::MatrixX3d(points.size(), 3);
  for (unsigned i = 0; i < points.size(); i++)
    for (unsigned j = 0; j < 3; j++)
      atom_position_to_super_cell(i, j) = points[i]["coordinates"][j].as<double>();
  AtomPosition = atom_position_to_super_cell
    * (SuperCellDeformation * SuperCellMultiplier.cast<double>().asDiagonal() * PrimativeCell);

  if (hdf5_file)
  {
    HighFive::File file(*hdf5_file, HighFive::File::ReadOnly);
    auto size = file.getDataSet("/frequency").getDimensions();
    auto frequency = file.getDataSet("/frequency")
      .read<std::vector<std::vector<std::vector<double>>>>();
    auto eigenvector_vector = file.getDataSet("/eigenvector")
      .read<std::vector<std::vector<std::vector<std::vector<PhonopyComplex>>>>>();
    auto path = file.getDataSet("/path")
      .read<std::vector<std::vector<std::vector<double>>>>();
    QPointData.resize(size[0] * size[1]);
    for (unsigned i = 0; i < size[0]; i++)
      for (unsigned j = 0; j < size[1]; j++)
      {
        QPointData[i * size[1] + j].QPoint = Eigen::Vector3d(path[i][j].data());
        QPointData[i * size[1] + j].ModeData.resize(size[2]);
        for (unsigned k = 0; k < size[2]; k++)
        {
          QPointData[i * size[1] + j].ModeData[k].Frequency = frequency[i][j][k];
          Eigen::MatrixX3cd eigenvectors(AtomPosition.rows(), 3);
          for (unsigned l = 0; l < AtomPosition.rows(); l++)
            for (unsigned m = 0; m < 3; m++)
              eigenvectors(l, m)
                = eigenvector_vector[i][j][k][l * 3 + m].r + eigenvector_vector[i][j][k][l * 3 + m].i * 1i;
          QPointData[i * size[1] + j].ModeData[k].AtomMovement = eigenvectors / eigenvectors.norm();
        }
      }
  }
  else
  {
    auto phonon = node["phonon"].as<std::vector<YAML::Node>>();
    QPointData.resize(phonon.size());
    for (unsigned i = 0; i < phonon.size(); i++)
    {
      QPointData[i].QPoint.resize(3);
      for (unsigned j = 0; j < 3; j++)
        QPointData[i].QPoint(j) = phonon[i]["q-position"][j].as<double>();
      auto band = phonon[i]["band"].as<std::vector<YAML::Node>>();
      QPointData[i].ModeData.resize(band.size());
      for (unsigned j = 0; j < band.size(); j++)
      {
        QPointData[i].ModeData[j].Frequency = band[j]["frequency"].as<double>();
        auto eigenvector_vectors = band[j]["eigenvector"]
          .as<std::vector<std::vector<std::vector<double>>>>();
        Eigen::MatrixX3cd eigenvectors(AtomPosition.rows(), 3);
        for (unsigned k = 0; k < AtomPosition.rows(); k++)
          for (unsigned l = 0; l < 3; l++)
            eigenvectors(k, l)
              = eigenvector_vectors[k][l][0] + 1i * eigenvector_vectors[k][l][1];
        // 需要对读入的原子运动状态作相位转换, 使得它们与我们的约定一致(对超胞周期性重复)
        // 这里还要需要做归一化处理 (指将数据简单地作为向量处理的归一化)
        auto& AtomMovement = QPointData[i].ModeData[j].AtomMovement;
        // AtomMovement = eigenvectors.array().colwise() * (-2 * std::numbers::pi_v<double> * 1i
        //   * (atom_position_to_super_cell * input.QPointData[i].QPoint)).array().exp();
        // AtomMovement /= AtomMovement.norm();
        // phonopy 似乎已经进行了相位的转换！为什么？
        AtomMovement = eigenvectors / eigenvectors.norm();
      }
    }
  }
}
