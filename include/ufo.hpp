# pragma once
# include <biu.hpp>

namespace ufo
{
  // 在相位中, 约定为使用 $\exp (2 \pi i \vec{q} \cdot \vec{r})$ 来表示原子的运动状态
  //  (而不是 $\exp (-2 \pi i \vec{q} \cdot \vec{r})$)
  // 一些书定义的倒格矢中包含了 $2 \pi$ 的部分, 我们这里约定不包含这部分.
  //  也就是说, 正格子与倒格子的转置相乘, 得到单位矩阵.

  using namespace biu::literals;
  using namespace biu::stream_operators;

  // 许多函数通用的一个结构体
  struct CommonData
  {
    // 原胞相关信息，其中存储的信息是直接读入的，而不是反折叠后的结果
    struct
    {
      // 晶胞，单位为埃，每行代表一个格矢
      Eigen::Matrix3d Cell;
      // 原子种类及其个数
      std::vector<std::pair<std::string, std::size_t>> AtomType;
      // 原子位置，单位为原胞格矢
      Eigen::MatrixX3d AtomPosition;
      // 单胞中每个 q 点的信息
      struct QpointType
      {
        // q 点的坐标，单位为倒格矢
        Eigen::Vector3d Qpoint;
        struct ModeType
        {
          // 振动频率，单位为 THz
          double Frequency;
          // 振动模式，单位为埃
          Eigen::MatrixX3cd EigenVector;
          // 拉曼张量
          // TODO: 单位是什么？
          std::optional<double> WeightOnRaman;
          // work around for fpr
          std::optional<std::array<std::array<double, 3>, 3>> RamanTensor;
          using serialize = zpp::bits::members<4>;
        };
        // 各个模式的的信息
        std::vector<ModeType> Mode;
        using serialize = zpp::bits::members<2>;
      };
      std::vector<QpointType> Qpoint;
      using serialize = zpp::bits::members<4>;
    } Primative;

    // 超胞相关信息，包括反折叠后的信息
    struct
    {
      // 超胞的晶胞，单位为埃，每行代表一个格矢
      Eigen::Matrix3d Cell;
      // 单胞与到超胞的变换
      Eigen::Vector3i CellMultiplier;
      Eigen::Matrix3d CellDeformation;
      std::vector<std::pair<std::string, std::size_t>> AtomType;
      Eigen::MatrixX3d AtomPosition;
      std::optional<std::array<double, 3>> AtomTranslation;
      struct QpointType
      {
        Eigen::Vector3d Qpoint;
        std::vector<Eigen::Vector3d> SubQpoint;
        struct ModeType
        {
          double Frequency;
          Eigen::MatrixX3cd EigenVector;
          std::vector<double> WeightOnUnfold;
          std::optional<double> WeightOnSelectedAtom;
          std::optional<double> WeightOnRaman;
          std::optional<std::array<std::array<double, 3>, 3>> RamanTensor;
          using serialize = zpp::bits::members<6>;
        };
        std::vector<ModeType> Mode;
        using serialize = zpp::bits::members<3>;
      };
      std::vector<QpointType> Qpoint;
      using serialize = zpp::bits::members<7>;
    } Super;

    // 选择的原子
    std::optional<std::set<std::size_t>> SelectedAtom;

    // 拉曼散射选择的偏振方向
    std::optional<std::array<Eigen::Vector3d, 2>> RamanPolarization;

    // 各种原子的质量
    std::map<std::string, double> AtomMass;
    using serialize = zpp::bits::members<5>;
  };

  void fold(std::string config_file);
  void unfold(std::string config_file);
  void project_to_atom(std::string config_file);
  void project_to_mode(std::string config_file);
  void plot_band(std::string config_file);
  void plot_point(std::string config_file);
  void raman_create_displacement(std::string config_file);
  void raman_extract(std::string path);
  void raman_apply_contribution(std::string config_file);
  void raman_bec_rotate(std::string config_file);
  void raman_bec_apply_contribution(std::string config_file);
}
