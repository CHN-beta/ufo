# include <iostream>
# include <array>
# include <numbers>
# include <numeric>
# include <fstream>
# include <optional>
# include <array>
# include <utility>
# include <yaml-cpp/yaml.h>
# include <eigen3/Eigen/Dense>
# include <concurrencpp/concurrencpp.h>
# include <fmt/format.h>

using namespace std::literals;

// 计算单胞中的 q 点在超胞中的对应

struct Input
{
  // 单胞的三个格矢，每行表示一个格矢的坐标，单位为埃
  Eigen::Matrix3d PrimativeCell;
  // 单胞到超胞的格矢转换时用到的矩阵
  Eigen::Vector<unsigned, 3> SuperCellMultiplier;
  Eigen::Matrix<double, 3, 3> SuperCellDeformation;

  // Q 点的坐标，单位为单胞的倒格矢
  std::vector<Eigen::Vector3d> QPointData;
};
template<> struct YAML::convert<Input> { static bool decode(const Node& node, Input& input); };

struct Output
{
  // Q 点的坐标，单位为超胞的倒格矢
  std::vector<Eigen::Vector3d> QPointData;
};

int main(int argc, char** argv)
{
  if (argc != 3)
    throw std::runtime_error("Usage: " + std::string(argv[0]) + " input.yaml output.yaml");

  auto input = YAML::LoadFile(argv[1]).as<Input>();
  Output output;

  for (auto qpoint_by_reciprocal_primative_cell : input.QPointData)
  {
    // 计算出这个 q 点的绝对坐标, 再计算出它相对于超胞倒格子的相对坐标. 将这个结果取小数部分, 就得到了 meta qpoint 的坐标(相对于超胞倒格子)
    std::cout << "PrimativeCell:\n" << input.PrimativeCell << "\n";
    auto reciprocal_primative_cell = input.PrimativeCell.inverse().transpose();
    std::cout << "reciprocal_primative_cell:\n" << reciprocal_primative_cell << "\n";
    auto qpoint = (qpoint_by_reciprocal_primative_cell.transpose() * reciprocal_primative_cell).transpose();
    std::cout << "qpoint:\n" << qpoint << "\n";
    auto reciprocal_super_cell = (input.SuperCellDeformation * input.SuperCellMultiplier.cast<double>().asDiagonal()
      * input.PrimativeCell).inverse().transpose();
    std::cout << "reciprocal_super_cell:\n" << reciprocal_super_cell << "\n";
    auto qpoint_by_reciprocal_super_cell = (qpoint.transpose() * reciprocal_super_cell.inverse()).transpose();
    std::cout << "qpoint_by_reciprocal_super_cell:\n" << qpoint_by_reciprocal_super_cell << "\n";
    auto meta_qpoint_by_reciprocal_super_cell = [&]
      { auto _ = qpoint_by_reciprocal_super_cell.array(); return _ - _.floor(); }();
    std::cout << "meta_qpoint_by_reciprocal_super_cell:\n" << meta_qpoint_by_reciprocal_super_cell << "\n";
    output.QPointData.push_back(meta_qpoint_by_reciprocal_super_cell);
  }

  std::ofstream(argv[2]) << [&]
  {
    std::stringstream print;
    print << "QPointData:\n";
    for (auto& qpoint: output.QPointData)
      print << fmt::format("  - [ {:.3f}, {:.3f}, {:.3f} ]\n", qpoint(0), qpoint(1), qpoint(2));
    return print.str();
  }();
}

bool YAML::convert<Input>::decode(const Node& node, Input& input)
{
  for (unsigned i = 0; i < 3; i++)
    for (unsigned j = 0; j < 3; j++)
      input.PrimativeCell(i, j) = node["lattice"][i][j].as<double>();

  input.SuperCellMultiplier.setZero();
  for (unsigned i = 0; i < 3; i++)
    input.SuperCellMultiplier(i) = node["SuperCellMultiplier"][i].as<int>();

  for (unsigned i = 0; i < 3; i++)
    for (unsigned j = 0; j < 3; j++)
      input.SuperCellDeformation(i, j) = node["SuperCellDeformation"][i][j].as<double>();

  auto points = node["points"].as<std::vector<std::vector<double>>>();

  input.QPointData.resize(points.size());
  for (unsigned i = 0; i < points.size(); i++)
  {
    for (unsigned j = 0; j < 3; j++)
      input.QPointData[i](j) = points[i][j];
  }

  return true;
}
