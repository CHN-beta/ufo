# include <iostream>
# include <array>
# include <numbers>
# include <numeric>
# include <fstream>
# include <yaml-cpp/yaml.h>
# include <eigen3/Eigen/Dense>

using namespace std::literals;

// 在相位中, 约定为使用 $\exp (2 \pi i \vec{q} \cdot \vec{r})$ 来表示原子的运动状态
//  (而不是 $\exp (-2 \pi i \vec{q} \cdot \vec{r})$)
// 一些书定义的倒格矢中包含了 $2 \pi$ 的部分, 我们这里约定不包含这部分.
//  也就是说, 正格子与倒格子的转置相乘, 得到单位矩阵.

struct Input
{
  // 单胞的三个格矢，每行表示一个格矢的坐标，单位为埃
  Eigen::Matrix3d PrimativeCell;
  // 超胞在各个方向上是单胞的多少倍，这是一个对角矩阵
  Eigen::Matrix3i SuperCellMultiplier;
  // 在单胞内取几个平面波的基矢
  // 在 debug 阶段, 先仅取一个
  Eigen::Vector3i PrimativeCellBasisNumber;
  // 超胞中原子的坐标，每行表示一个原子的坐标，单位为埃
  Eigen::MatrixX3d AtomPosition;

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
      // 这个数据在读入后会立即被归一化处理.
      Eigen::MatrixX3cd AtomMovement;
    };
    std::vector<ModeDataType_> ModeData;
  };
  std::vector<QPointDataType_> QPointData;

  // 从文件中读取输入, 文件中应当包含: (大多数据可以直接从 phonopy 的输出中复制)
  // 单胞的格矢: lattice 单位为埃 直接从 phonopy 的输出中复制
  // 超胞的倍数: SuperCellMultiplier 手动输入, 为一个包含三个整数的数组
  // 平面波的基矢个数: PrimativeCellBasisNumber 手动输入, 为一个包含三个整数的数组
  // 超胞中原子的坐标: points[*].coordinates 单位为超胞的格矢 直接从 phonopy 的输出中复制
  // 各个 Q 点的坐标: phonon[*].q-position 单位为超胞的倒格子的格矢 直接从 phonopy 的输出中复制
  // 各个模式的频率: phonon[*].band[*].frequency 单位为 THz 直接从 phonopy 的输出中复制
  // 各个模式的原子运动状态: phonon[*].band[*].eigenvector 直接从 phonopy 的输出中复制
  // 文件中可以有多余的项目, 多余的项目不管.
  Input(std::string filename) : Input(YAML::LoadFile(filename)) {}
  Input(YAML::Node node);
};

struct Output
{
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
      // 模式的权重
      double Weight;
    };
    std::vector<ModeDataType_> ModeData;
  };
  std::vector<QPointDataType_> QPointData;

  // 将数据写入文件
  void write(std::string filename) const;
};

int main(int argc, const char** argv)
{
  if (argc != 3)
    throw std::runtime_error("Usage: " + std::string(argv[0]) + " input.yaml output.yaml");

  Input input(argv[1]);

  // 反折叠的原理: 将超胞中的原子运动状态, 投影到一组平面波构成的基矢中.
  // 每一个平面波的波矢由两部分相加得到: 一部分是单胞倒格子的整数倍, 所取的个数有一定任意性, 论文中建议取大约单胞中原子个数那么多个;
  //  对于没有缺陷的情况, 取一个应该就足够了.
  // 另一部分是超胞倒格子的整数倍, 取 n 个, n 为超胞对应的单胞的倍数, 其实也就是倒空间中单胞对应倒格子中超胞的格点.
  // 只要第一部分取得足够多, 那么单胞中原子的状态就可以完全被这些平面波描述.
  // 将超胞中原子的运动状态投影到这些基矢上, 计算出投影的系数, 就可以将超胞的原子运动状态分解到单胞中的多个 q 点上.

  // 构建基
  // 外层下标对应超胞倒格子的整数倍那部分(第二部分), 也就是对应不同反折叠后的 q 点
  // 内层下标对应单胞倒格子的整数倍那部分(第一部分), 也就是对应同一个反折叠后的 q 点上的不同平面波
  std::vector<std::vector<Eigen::VectorXcd>> basis;
  basis.resize(input.SuperCellMultiplier.determinant());
  for (int i = 0; i < input.SuperCellMultiplier(0, 0); i++)
    for (int j = 0; j < input.SuperCellMultiplier(1, 1); j++)
      for (int k = 0; k < input.SuperCellMultiplier(2, 2); k++)
      {
        // 反折叠后的某个 q 点 对应的所有的基矢
        auto& basis_at_unfolded_qpoint
          = basis[i * input.SuperCellMultiplier(1, 1) * input.SuperCellMultiplier(2, 2)
            + j * input.SuperCellMultiplier(2, 2) + k];
        basis_at_unfolded_qpoint.resize(input.PrimativeCellBasisNumber.prod());
        for (int x = 0; x < input.PrimativeCellBasisNumber(0); x++)
          for (int y = 0; y < input.PrimativeCellBasisNumber(1); y++)
            for (int z = 0; z < input.PrimativeCellBasisNumber(2); z++)
            {
              // 计算 q 点的坐标, 单位为相对于超胞的倒格矢
              auto qpoint_relative_to_super_cell =
                Eigen::Vector3i({{i}, {j}, {k}})
                + input.SuperCellMultiplier * Eigen::Vector3i({{x}, {y}, {z}});
              // 将 q 点坐标转换为埃^-1
              auto qpoint = (qpoint_relative_to_super_cell.transpose().cast<double>()
                * (input.SuperCellMultiplier.cast<double>().inverse().transpose())).transpose();
              // 计算基矢
              auto& single_basis = basis_at_unfolded_qpoint
                [x * input.PrimativeCellBasisNumber(1) * input.PrimativeCellBasisNumber(2)
                  + y * input.PrimativeCellBasisNumber(2) + z];
              single_basis = (-2 * std::numbers::pi_v<double> * 1i * (input.AtomPosition * qpoint)).array().exp();
            }
      }

  // 计算投影的结果
  // 最外层下标对应反折叠前的 q 点, 第二层下标对应不同模式, 第三层下标对应这个模式在反折叠后的 q 点
  std::vector<std::vector<std::vector<double>>> projection_coefficient;
  projection_coefficient.resize(input.QPointData.size());
  for (int i_of_folded_qpoint = 0; i_of_folded_qpoint < input.QPointData.size(); i_of_folded_qpoint++)
  {
    projection_coefficient[i_of_folded_qpoint].resize(input.QPointData[i_of_folded_qpoint].ModeData.size());
    for (int i_of_mode = 0; i_of_mode < projection_coefficient[i_of_folded_qpoint].size(); i_of_mode++)
    {
      auto& coefficient_at_mode = projection_coefficient[i_of_folded_qpoint][i_of_mode];
      coefficient_at_mode.resize(input.SuperCellMultiplier.determinant());
      for
      (int i_of_unfolded_qpoint = 0; i_of_unfolded_qpoint < coefficient_at_mode.size(); i_of_unfolded_qpoint++)
        // 对于 basis 中, 对应于单胞倒格子的部分, 以及对应于不同方向的部分, 分别求内积, 然后求绝对值, 然后求和
        for
        (
          int i_of_basis_in_primary_cell = 0;
          i_of_basis_in_primary_cell < basis[i_of_unfolded_qpoint].size();
          i_of_basis_in_primary_cell++
        )
          coefficient_at_mode[i_of_unfolded_qpoint] +=
          (
            basis[i_of_unfolded_qpoint][i_of_basis_in_primary_cell].transpose()
            * input.QPointData[i_of_folded_qpoint].ModeData[i_of_mode].AtomMovement
          ).array().abs2().sum();

      // 归一化
      auto sum = std::accumulate(coefficient_at_mode.begin(), coefficient_at_mode.end(), 0.);
      for (auto& coefficient : coefficient_at_mode)
        coefficient /= sum;
    }
  }

  // 填充输出对象
  Output output;
  for (int i_of_folded_qpoint = 0; i_of_folded_qpoint < input.QPointData.size(); i_of_folded_qpoint++)
    // for each unfolded q point corresponding to this folded q point
    for (int x = 0; x < input.SuperCellMultiplier(0, 0); x++)
      for (int y = 0; y < input.SuperCellMultiplier(1, 1); y++)
        for (int z = 0; z < input.SuperCellMultiplier(2, 2); z++)
        {
          auto& unfolded_qpoint = output.QPointData.emplace_back();
          unfolded_qpoint.QPoint = input.SuperCellMultiplier.cast<double>().inverse() *
          (input.QPointData[i_of_folded_qpoint].QPoint
              + Eigen::Vector3i({{x}, {y}, {z}}).cast<double>());
          for (int i_of_mode = 0; i_of_mode < input.QPointData[i_of_folded_qpoint].ModeData.size(); i_of_mode++)
            if (projection_coefficient[i_of_folded_qpoint][i_of_mode]
              [x * input.SuperCellMultiplier(1, 1) * input.SuperCellMultiplier(2, 2)
                + y * input.SuperCellMultiplier(2, 2) + z] > 1e-3)
            {
              auto& mode = unfolded_qpoint.ModeData.emplace_back();
              mode.Frequency = input.QPointData[i_of_folded_qpoint].ModeData[i_of_mode].Frequency;
              mode.Weight = projection_coefficient[i_of_folded_qpoint][i_of_mode]
                [x * input.SuperCellMultiplier(1, 1) * input.SuperCellMultiplier(2, 2)
                  + y * input.SuperCellMultiplier(2, 2) + z];
            }
        }

  output.write(argv[2]);
}

Input::Input(YAML::Node root)
{
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      PrimativeCell(i, j) = root["lattice"][i][j].as<double>();

  SuperCellMultiplier.setZero();
  for (int i = 0; i < 3; i++)
    SuperCellMultiplier(i, i) = root["SuperCellMultiplier"][i].as<int>();

  for (int i = 0; i < 3; i++)
    PrimativeCellBasisNumber(i) = root["PrimativeCellBasisNumber"][i].as<int>();

  auto points = root["points"].as<std::vector<YAML::Node>>();
  auto atom_position_to_super_cell = Eigen::MatrixX3d(points.size(), 3);
  for (int i = 0; i < points.size(); i++)
    for (int j = 0; j < 3; j++)
      atom_position_to_super_cell(i, j) = points[i]["coordinates"][j].as<double>();
  AtomPosition = atom_position_to_super_cell * (SuperCellMultiplier.cast<double>() * PrimativeCell);

  auto phonon = root["phonon"].as<std::vector<YAML::Node>>();
  QPointData.resize(phonon.size());
  for (int i = 0; i < phonon.size(); i++)
  {
    QPointData[i].QPoint.resize(3);
    for (int j = 0; j < 3; j++)
      QPointData[i].QPoint(j) = phonon[i]["q-position"][j].as<double>();
    auto band = phonon[i]["band"].as<std::vector<YAML::Node>>();
    QPointData[i].ModeData.resize(band.size());
    for (int j = 0; j < band.size(); j++)
    {
      QPointData[i].ModeData[j].Frequency = band[j]["frequency"].as<double>();
      auto eigenvectors = Eigen::MatrixX3cd(AtomPosition.rows(), 3);
      auto eigenvector_vectors = band[j]["eigenvector"]
        .as<std::vector<std::vector<std::vector<double>>>>();
      for (int k = 0; k < AtomPosition.rows(); k++)
        for (int l = 0; l < 3; l++)
          eigenvectors(k, l)
            = eigenvector_vectors[k][l][0] + 1i * eigenvector_vectors[k][l][1];
      // 需要对读入的原子运动状态作相位转换, 使得它们与我们的约定一致(对超胞周期性重复)
      // QPointData[i].ModeData[j].AtomMovement
      //   = eigenvectors.array().colwise() * (-2 * std::numbers::pi_v<double> * 1i
      //       * (AtomPosition * QPointData[i].QPoint)).array().exp();
      QPointData[i].ModeData[j].AtomMovement.resize(AtomPosition.rows(), 3);
      for (int k = 0; k < AtomPosition.rows(); k++)
        QPointData[i].ModeData[j].AtomMovement.row(k) = eigenvectors.row(k)
          * std::exp(-2 * std::numbers::pi_v<double> * 1i
              * (atom_position_to_super_cell.row(k).dot(QPointData[i].QPoint)));

      // print AtomMovement
      // std::cout << "AtomMovement" << std::endl;
      // std::cout << QPointData[i].ModeData[j].AtomMovement << std::endl;


      // 这里还要需要做归一化处理
      // phonopy 的文档似乎没有保证一定做了归一化处理, 再来做一遍吧.
      auto sum = QPointData[i].ModeData[j].AtomMovement.cwiseAbs2().sum();
      QPointData[i].ModeData[j].AtomMovement /= std::sqrt(sum);
    }
  }
}

void Output::write(std::string filename) const
{
  YAML::Node root;
  root["QPointData"] = YAML::Node(YAML::NodeType::Sequence);
  for (int i = 0; i < QPointData.size(); i++)
  {
    root["QPointData"][i]["QPoint"] = YAML::Node(YAML::NodeType::Sequence);
    for (int j = 0; j < 3; j++)
      root["QPointData"][i]["QPoint"][j] = QPointData[i].QPoint(j);
    root["QPointData"][i]["ModeData"] = YAML::Node(YAML::NodeType::Sequence);
    for (int j = 0; j < QPointData[i].ModeData.size(); j++)
    {
      root["QPointData"][i]["ModeData"][j]["Frequency"] = QPointData[i].ModeData[j].Frequency;
      root["QPointData"][i]["ModeData"][j]["Weight"] = QPointData[i].ModeData[j].Weight;
    }
  }
  std::ofstream(filename) << root;
}
