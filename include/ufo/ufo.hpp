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
# include <any>
# include <map>
# include <vector>
# include <span>
# include <yaml-cpp/yaml.h>
# include <Eigen/Dense>
# include <concurrencpp/concurrencpp.h>
# include <fmt/format.h>
# include <fmt/std.h>
# include <fmt/ranges.h>
# include <highfive/H5File.hpp>
# include <zpp_bits.h>
# include <matplot/matplot.h>

using namespace std::literals;

struct PhonopyComplex { double r, i; };
HighFive::CompoundType create_compound_complex();
HIGHFIVE_REGISTER_TYPE(PhonopyComplex, create_compound_complex)

namespace Eigen
{
  constexpr inline auto serialize(auto & archive, Eigen::Matrix3d& matrix)
    { return archive(std::span(matrix.data(), matrix.size())); }
  constexpr inline auto serialize(auto & archive, const Eigen::Matrix3d& matrix)
    { return archive(std::span(matrix.data(), matrix.size())); }
  constexpr inline auto serialize(auto & archive, Eigen::Vector3d& vector)
    { return archive(std::span(vector.data(), vector.size())); }
  constexpr inline auto serialize(auto & archive, const Eigen::Vector3d& vector)
    { return archive(std::span(vector.data(), vector.size())); }
}

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
  std::optional<Eigen::Matrix<double, 3, 3>> SuperCellDeformation;
  // 在单胞内取几个平面波的基矢
  Eigen::Vector<unsigned, 3> PrimativeCellBasisNumber;

  struct InputOutputFile_
  {
    std::string FileName;
    std::string Format;
    std::map<std::string, std::any> ExtraParameters;
  };

  // 从哪个文件读入 AtomPosition, 以及这个文件的格式, 格式可选值包括 "yaml"
  InputOutputFile_ AtomPositionInputFile;
  // 从哪个文件读入 QPointData, 以及这个文件的格式, 格式可选值包括 "yaml" 和 "hdf5"
  InputOutputFile_ QPointDataInputFile;

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
      // 这个数据在读入后不会被立即归一化.
      Eigen::MatrixX3cd AtomMovement;
    };
    std::vector<ModeDataType_> ModeData;
  };
  std::vector<QPointDataType_> QPointData;

  // 输出到哪些文件, 以及使用怎样的格式, 格式可选值包括:
  // yaml: 使用 yaml 格式输出
  // yaml-human-readable: 使用 yaml 格式输出, 但是输出的结果更适合人类阅读, 包括合并相近的模式, 去除权重过小的模式, 限制输出的小数位数.
  // zpp: 使用 zpp-bits 序列化, 可以直接被 plot.cpp 读取
  std::vector<InputOutputFile_> QPointDataOutputFile;

  // 从文件中读取输入 (包括一个较小的配置文件, 和一个 hdf5 或者一个 yaml 文件), 文件中应当包含:
  // 单胞的格矢: PrimativeCell 单位为埃 直接从 phonopy 的输出中复制
  // 超胞的倍数: SuperCellMultiplier 手动输入, 为一个包含三个整数的数组
  // 超胞的变形: SuperCellDeformation 手动输入, 为一个三阶方阵
  // 平面波的基矢个数: PrimativeCellBasisNumber 手动输入, 为一个包含三个整数的数组
  // 另外还有一个文件, 直接将 phonopy 的输出复制过来即可, 如果是 yaml, 应该包含下面的内容:
  // 超胞中原子的坐标: points[*].coordinates 单位为超胞的格矢 直接从 phonopy 的输出中复制
  // 各个 Q 点的坐标: phonon[*].q-position 单位为超胞的倒格子的格矢 直接从 phonopy 的输出中复制
  // 各个模式的频率: phonon[*].band[*].frequency 单位为 THz 直接从 phonopy 的输出中复制
  // 各个模式的原子运动状态: phonon[*].band[*].eigenvector 直接从 phonopy 的输出中复制
  // 文件中可以有多余的项目, 多余的项目不管.
  Input(std::string filename);
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
    std::size_t SourceIndex_;

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

  void write(std::string filename, std::string format, unsigned percision = 10) const;
  Output() = default;
  Output(std::string filename);

  using serialize = zpp::bits::members<1>;
};

concurrencpp::generator<std::pair<Eigen::Vector<unsigned, 3>, unsigned>>
  triplet_sequence(Eigen::Vector<unsigned, 3> range);
