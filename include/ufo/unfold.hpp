# pragma once
# include <ufo/solver.hpp>

namespace ufo
{
  // 反折叠的原理: 将超胞中的原子运动状态, 投影到一组平面波构成的基矢中.
  // 每一个平面波的波矢由两部分相加得到: 一部分是单胞倒格子的整数倍, 所取的个数有一定任意性, 论文中建议取大约单胞中原子个数那么多个;
  //  对于没有缺陷的情况, 取一个应该就足够了.
  // 另一部分是超胞倒格子的整数倍, 取 n 个, n 为超胞对应的单胞的倍数, 其实也就是倒空间中单胞对应倒格子中超胞的格点.
  // 只要第一部分取得足够多, 那么单胞中原子的状态就可以完全被这些平面波描述.
  // 将超胞中原子的运动状态投影到这些基矢上, 计算出投影的系数, 就可以将超胞的原子运动状态分解到单胞中的多个 q 点上.
  class UnfoldSolver : public Solver
  {
    public:
      struct InputType
      {
        // 单胞的三个格矢，每行表示一个格矢的坐标，单位为埃
        Eigen::Matrix3d PrimativeCell;
        // 单胞到超胞的格矢转换时用到的矩阵
        // SuperCellMultiplier 是一个三维列向量且各个元素都是整数，表示单胞在各个方向扩大到多少倍之后，可以得到和超胞一样的体积
        // SuperCsolver.hpp>mation 是一个行列式为 1 的矩阵，它表示经过 SuperCellMultiplier 扩大后，还需要怎样的变换才能得到超胞
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

        struct InputOutputFile
        {
          std::string FileName;
          std::string Format;
          std::map<std::string, std::any> ExtraParameters;
        };

        // 从哪个文件读入 AtomPosition, 以及这个文件的格式, 格式可选值包括 "yaml"
        InputOutputFile AtomPositionInputFile;
        // 从哪个文件读入 QpointData, 以及这个文件的格式, 格式可选值包括 "yaml" 和 "hdf5"
        InputOutputFile QpointDataInputFile;

        // 超胞中原子的坐标，每行表示一个原子的坐标，单位为埃
        Eigen::MatrixX3d AtomPosition;
        // 关于各个 Q 点的数据
        struct QpointDataType
        {
          // Q 点的坐标，单位为超胞的倒格矢
          Eigen::Vector3d Qpoint;

          // 关于这个 Q 点上各个模式的数据
          struct ModeDataType
          {
            // 模式的频率，单位为 THz
            double Frequency;
            /*
              phonopy 中定义的动态矩阵的特征向量, 代表模式中各个原子的运动状态
                (也即原子位移中, 关于超胞有周期性的那一部分, 再乘以原子质量的开方).
              phonopy 的动态矩阵的定义与一般的定义不同, 因此这里的这个特征向量不是一般意义上的动态矩阵的特征向量.
              尽管在计算时, 确实需要将这个特征向量转化为标准的动态矩阵的特征向量, 但为了节省计算量, 在反折叠时一起计算, 而不是在输入时计算.
            */
            Eigen::MatrixX3cd PhonopyEigenvector;
          };
          std::vector<ModeDataType> ModeData;
        };
        std::vector<QpointDataType> QpointData;

        // 输出到哪些文件, 以及使用怎样的格式, 格式可选值包括:
        // yaml: 使用 yaml 格式输出
        // yaml-human-readable: 使用 yaml 格式输出, 但是输出的结果更适合人类阅读,
        //  包括合并相近的模式, 去除权重过小的模式, 限制输出的小数位数.
        // zpp: 使用 zpp-bits 序列化, 可以直接被 plot.cpp 读取
        std::vector<InputOutputFile> QpointDataOutputFile;

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
        InputType(std::string filename);
      };
      struct OutputType
      {
        // 关于各个 Q 点的数据
        struct QpointDataType
        {
          // Q 点的坐标，单位为单胞的倒格矢
          Eigen::Vector3d Qpoint;

          // 来源于哪个 Q 点, 单位为超胞的倒格矢
          Eigen::Vector3d Source;
          std::size_t SourceIndex_;

          // 关于这个 Q 点上各个模式的数据
          struct ModeDataType
          {
            // 模式的频率，单位为 THz
            double Frequency;
            // 模式的权重
            double Weight;
          };
          std::vector<ModeDataType> ModeData;
        };
        std::vector<QpointDataType> QpointData;

        void write(decltype(InputType::QpointDataOutputFile) output_files) const;
        void write(std::string filename, std::string format, unsigned percision = 10) const;

        using serialize = zpp::bits::members<1>;

        virtual ~OutputType() = default;
      };
    protected:
      InputType Input_;
      std::optional<OutputType> Output_;

    public:
      UnfoldSolver(std::string config_file);
      UnfoldSolver& operator()() override;

      // 计算单个 MetaQpoint 到所有 SubQpoint 的投影系数, 是反折叠的核心步骤
      struct ConstructProjectionCoefficientInputType
      {
        struct
        {
          const Eigen::Matrix3d& PrimativeCell;
          const Eigen::Vector<unsigned, 3>& SuperCellMultiplier;
          const std::optional<Eigen::Matrix<double, 3, 3>>& SuperCellDeformation;
          const Eigen::Vector<unsigned, 3>& PrimativeCellBasisNumber;
          const Eigen::MatrixX3d& AtomPosition;
        } GlobalData;
        struct
        {
          const Eigen::Vector3d& MetaQpointByReciprocalSuperCell;
          const std::vector<InputType::QpointDataType::ModeDataType>& ModeData;
        } MetaQpointData;
      };
      struct ConstructProjectionCoefficientOutputType
      {
        struct SubQpointDataType
        {
          Eigen::Vector3d SubQpointByReciprocalPrimativeCell;
          std::vector<OutputType::QpointDataType::ModeDataType> ModeData;
        };
        std::vector<SubQpointDataType> SubQpointData;
      };
      static ConstructProjectionCoefficientOutputType construct_projection_coefficient
        (ConstructProjectionCoefficientInputType Input);
  };
}
