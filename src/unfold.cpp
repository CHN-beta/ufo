# include <ufo.hpp>
# include <thread>
# include <syncstream>
# include <execution>
# include <ranges>

void ufo::unfold(std::string config_file)
{
  // 反折叠的原理: 将超胞中的原子运动状态, 投影到一组平面波构成的基矢中.
  // 每一个平面波的波矢由两部分相加得到: 一部分是单胞倒格子的整数倍, 所取的个数有一定任意性, 论文中建议取大约单胞中原子个数那么多个;
  //  对于没有缺陷的情况, 取一个应该就足够了.
  // 这些平面波以原胞为周期。
  // 另一部分是超胞倒格子的整数倍, 取 n 个, n 为超胞对应的单胞的倍数, 其实也就是倒空间中单胞对应倒格子中超胞的格点.
  // 只要第一部分取得足够多, 那么单胞中原子的状态就可以完全被这些平面波描述.
  // 将超胞中原子的运动状态投影到这些基矢上, 计算出投影的系数, 就可以将超胞的原子运动状态分解到单胞中的多个 q 点上.

  struct Config
  {
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
    Eigen::Matrix3d SuperCellDeformation;
    Eigen::Vector3i SuperCellMultiplier;

    // 在单胞内取几个平面波的基矢
    Eigen::Vector<std::size_t, 3> PrimativeCellBasisNumber;

    // 单胞的 phonopy 输出的 phonopy.yaml，用来读入单胞的晶格、原子坐标、原子类型、原子质量
    std::string PrimativePhonopy;
    // 单胞的 phonopy 输出的 band.hdf5 或 qpoint.hdf5，用来读入 q 点和振动模式
    std::string PrimativeQpoint;
    // 同上，但是是超胞里的结果
    std::string SuperPhonopy;
    std::string SuperQpoint;

    std::string OutputFile;
  };

  // 从文件中读取 q 点数据
  // data 为 CommonData::Primative 或 CommonData::Super
  // 返回值为原子类型和原子质量的对应关系
  auto read_qpoint = [](std::string phonopy_file, std::string qpoint_file, auto& data)
  {
    // phonopy 的输出有两种可能。
    // 直接指定计算的 q 点时，frequency 是 2 维，这时第一个维度是 q 点，第二个维度是不同模式。
    // 计算能带时，frequency 是 3 维，相比于二维的情况多了第一个维度，表示 q 点所在路径。
    // qpoint 或 path，以及 eigenvector 也有类似的变化。
    // 除此以外，eigenvector 后两个维度分别指示模式的特征向量的各个维度和各个模式（而不是各个模式和特征向量的各个维度），
    // 因为后两个维度的尺寸总是一样的（模式个数等于原子坐标个数），非常容易搞错。
    std::vector<std::array<double, 3>> qpoint;
    std::vector<std::vector<double>> frequency;
    std::vector<std::vector<std::vector<std::complex<double>>>> eigenvector_vector;
    auto file = biu::Hdf5file(qpoint_file);

    if (file.File.getDataSet("/frequency").getDimensions().size() == 2)
      file.read("/frequency", frequency)
        .read("/eigenvector", eigenvector_vector)
        .read("/qpoint", qpoint);
    else
    {
      std::vector<std::vector<std::array<double, 3>>> temp_path;
      std::vector<std::vector<std::vector<double>>> temp_frequency;
      std::vector<std::vector<std::vector<std::vector<std::complex<double>>>>> temp_eigenvector_vector;
      file.read("/frequency", temp_frequency)
        .read("/eigenvector", temp_eigenvector_vector)
        .read("/path", temp_path);
      frequency = temp_frequency | ranges::views::join | ranges::to_vector;
      qpoint = temp_path | ranges::views::join | ranges::to_vector;
      eigenvector_vector = temp_eigenvector_vector | ranges::views::join | ranges::to_vector;
    }

    // 整理并写入得到结果
    auto number_of_qpoints = frequency.size(), number_of_modes = frequency[0].size();
    data.Qpoint.resize(number_of_qpoints);
    for (auto i : std::views::iota(0u, number_of_qpoints)) 
    {
      data.Qpoint[i].Qpoint = qpoint[i] | biu::toEigen<>;
      data.Qpoint[i].Mode.resize(number_of_modes);
      for (auto j : std::views::iota(0u, number_of_modes)) 
      {
        data.Qpoint[i].Mode[j].Frequency = frequency[i][j];
        auto number_of_atoms = number_of_modes / 3;
        Eigen::MatrixX3cd eigenvector(number_of_atoms, 3);
        for (auto k : std::views::iota(0u, number_of_atoms))
          for (auto l : std::views::iota(0u, 3u))
            eigenvector(k, l) = eigenvector_vector[i][k * 3 + l][j];
        // 原则上讲，需要对读入的原子运动状态作相位转换, 使得它们与我们的约定一致(对超胞周期性重复)，但这个转换 phonopy 已经做了
        // 这里还要需要做归一化处理 (指将数据简单地作为向量处理的归一化)
        data.Qpoint[i].Mode[j].EigenVector = eigenvector / eigenvector.norm();
      }
    }

    // 读取并写入其它数据
    YAML::Node phonopy = YAML::LoadFile(phonopy_file);
    data.Cell = phonopy["unit_cell"]["lattice"].as<std::array<std::array<double, 3>, 3>>() | biu::toEigen<>;
    auto points = phonopy["points"].as<std::vector<YAML::Node>>();
    data.AtomType = points
      | ranges::views::transform([](const YAML::Node& point) { return point["symbol"].as<std::string>(); })
      | ranges::views::chunk_by(std::ranges::equal_to{})
      | ranges::views::transform([](const auto& chunk) { return std::pair{chunk[0], chunk.size()}; })
      | ranges::to_vector;
    data.AtomPosition = points
      | ranges::views::transform([](const YAML::Node& point)
        { return point["coordinates"].as<std::array<double, 3>>(); })
      | ranges::to_vector
      | biu::toEigen<>;
    return points
      | ranges::views::transform([](const YAML::Node& point)
        { return std::pair(point["symbol"].as<std::string>(), point["mass"].as<double>()); })
      | ranges::views::chunk_by(std::ranges::equal_to{})
      | ranges::views::transform([](const auto& chunk) { return chunk[0]; })
      | ranges::to<std::map<std::string, double>>;
  };

  // 构建基
  // 每个 q 点对应一组 sub qpoint。不同的 q 点所对应的 sub qpoint 是不一样的，但 sub qpoint 与 q 点的相对位移在不同 q 点之间是相同的。
  // 由于基只与这个相对位置有关（也就是说，不同 q 点的基是一样的），因此可以先计算出所有的基，这样降低计算量。
  // 外层下标对应超胞倒格子的整数倍那部分(第二部分), 也就是不同的 sub qpoint
  // 内层下标对应单胞倒格子的整数倍那部分(第一部分), 也就是 sub qpoint 上的不同平面波（取的数量越多，结果越精确）
  auto construct_basis = []
  (
    Eigen::Matrix3d primative_cell, Eigen::Vector3i super_cell_multiplier,
    Eigen::Vector<std::size_t, 3> primative_cell_basis_number, Eigen::MatrixX3d atom_position
  )
  {
    biu::Logger::Guard log;
    std::vector<std::vector<Eigen::VectorXcd>> basis(super_cell_multiplier.prod());
    // diff_of_sub_qpoint 表示 sub qpoint 与 qpoint 的相对位置，单位为超胞的倒格矢
    for (auto [diff_of_sub_qpoint_by_reciprocal_modified_super_cell, i_of_sub_qpoint]
      : biu::sequence(super_cell_multiplier))
    {
      basis[i_of_sub_qpoint].resize(primative_cell_basis_number.prod());
      for (auto [xyz_of_basis, i_of_basis]
        : biu::sequence(primative_cell_basis_number))
      {
        // 计算 q 点的坐标, 单位为单胞的倒格矢
        auto diff_of_sub_qpoint_by_reciprocal_primative_cell = xyz_of_basis.cast<double>()
          + super_cell_multiplier.cast<double>().cwiseInverse().asDiagonal()
          * diff_of_sub_qpoint_by_reciprocal_modified_super_cell.cast<double>();
        // 将单位转换为埃^-1
        auto diff_of_sub_qpoint = (diff_of_sub_qpoint_by_reciprocal_primative_cell.transpose()
          * (primative_cell.transpose().inverse())).transpose();
        // 计算基矢
        basis[i_of_sub_qpoint][i_of_basis]
          = (2i * std::numbers::pi_v<double> * (atom_position * diff_of_sub_qpoint)).array().exp();
      }
    }
    return basis;
  };

  // 计算从超胞到原胞的投影系数（不是分原子的投影系数），是反折叠的核心步骤
  // 返回的投影系数是一个三维数组，第一维对应不同的 q 点，第二维对应不同的模式，第三维对应不同的 sub qpoint
  auto construct_projection_coefficient = []
  (
    const std::vector<std::vector<Eigen::VectorXcd>>& basis,
    const std::vector<std::reference_wrapper<const Eigen::MatrixX3cd>>& modes,
    std::atomic<std::size_t>& number_of_finished_modes
  )
  {
    // 第一层下标对应不同模式, 第二层下标对应这个模式在反折叠后的 q 点(sub qpoint)
    std::vector<std::vector<double>> projection_coefficient(modes.size());
    // 对每个模式并行
    std::transform
    (
      std::execution::par, modes.begin(), modes.end(),
      projection_coefficient.begin(), [&](const auto& mode_data)
      {
        // 这里, mode_data 和 projection_coefficient 均指对应于一个模式的数据
        std::vector<double> projection_coefficient(basis.size());
        for (auto i_of_sub_qpoint : std::views::iota(0u, basis.size()))
          // 对于 basis 中, 对应于单胞倒格子的部分, 以及对应于不同方向的部分, 分别求内积, 然后求模方和
          for (auto i_of_basis : std::views::iota(0u, basis[i_of_sub_qpoint].size()))
            projection_coefficient[i_of_sub_qpoint] +=
              (basis[i_of_sub_qpoint][i_of_basis].transpose().conjugate() * mode_data.get()).array().abs2().sum();
        // 如果是严格地将向量分解到一组完备的基矢上, 那么不需要对计算得到的权重再做归一化处理
        // 但这里并不是这样一个严格的概念. 因此对分解到各个 sub qpoint 上的权重做归一化处理
        auto sum = ranges::accumulate(projection_coefficient, 0.);
        for (auto& _ : projection_coefficient) _ /= sum;
        number_of_finished_modes++;
        return projection_coefficient;
      }
    );
    return projection_coefficient;
  };

  biu::Logger::Guard log;

  log.info("Reading input file...");
  auto config = YAML::LoadFile(config_file).as<Config>();
  CommonData output;
  output.AtomMass = read_qpoint
    (config.PrimativePhonopy, config.PrimativeQpoint, output.Primative);
  output.AtomMass.merge(read_qpoint
    (config.SuperPhonopy, config.SuperQpoint, output.Super));
  output.Super.CellDeformation = config.SuperCellDeformation;
  output.Super.CellMultiplier = config.SuperCellMultiplier;
  log.info("Done.");

  log.info("Constructing basis...");
  auto basis = construct_basis
  (
    output.Primative.Cell, output.Super.CellMultiplier,
    config.PrimativeCellBasisNumber,
    output.Super.AtomPosition
      * (output.Super.CellDeformation * output.Super.CellMultiplier.cast<double>().asDiagonal() * output.Primative.Cell)
  );
  log.info("Done.");

  std::clog << "Calculating projection coefficient... " << std::flush;
  // 将所有模式放到一维来处理
  {
    auto modes = output.Super.Qpoint
      | ranges::views::transform([](const auto& qpoint)
        { return qpoint.Mode | ranges::views::transform([](const auto& mode)
          { return std::cref(mode.EigenVector); }); })
      | ranges::views::join
      | ranges::to_vector;
    std::atomic<std::size_t> number_of_finished_modes(0);
    std::thread print_thread([&]
    {
      while (true)
      {
        std::osyncstream(std::clog)
          << "\rCalculating projection coefficient... ({}/{})"_f(number_of_finished_modes, modes.size())
          << std::flush;
        std::this_thread::sleep_for(100ms);
        if (number_of_finished_modes == modes.size()) break;
      }
    });
    auto projection_coefficient = construct_projection_coefficient
      (basis, modes, number_of_finished_modes);
    std::size_t i_of_modes = 0;
    for (auto& qpoint : output.Super.Qpoint) for (auto& mode : qpoint.Mode)
      mode.WeightOnUnfold = projection_coefficient[i_of_modes++];
    print_thread.join();
  }
  std::clog << "\33[2K\rCalculating projection coefficient... Done." << std::endl;

  log.info("Writing data... ");
  std::ofstream(config.OutputFile, std::ios::binary) << biu::serialize<char>(output);
  log.info("Done.");

  log.info("Summary:");

  for (auto i_of_qpoint : std::views::iota(0u, output.Super.Qpoint.size()))
    for (auto i_of_mode : std::views::iota(0u, output.Super.Qpoint[i_of_qpoint].Mode.size()))
      for
      (
        auto i_of_sub_qpoint
          : std::views::iota(0u, output.Super.Qpoint[i_of_qpoint].SubQpoint.size())
      )
      if (output.Super.Qpoint[i_of_qpoint].Mode[i_of_mode].WeightOnUnfold[i_of_sub_qpoint] > 0.01)
        log.info("{}:{:.2f}:{}:{:.2f} -> {:.2f} {:.2f}"_f
        (
          i_of_qpoint, fmt::join(output.Super.Qpoint[i_of_qpoint].Qpoint, ", "),
          i_of_mode, output.Super.Qpoint[i_of_qpoint].Mode[i_of_mode].Frequency,
          fmt::join(output.Super.Qpoint[i_of_qpoint].SubQpoint[i_of_sub_qpoint], ", "),
          output.Super.Qpoint[i_of_qpoint].Mode[i_of_mode].WeightOnUnfold[i_of_sub_qpoint]
        ));
}
