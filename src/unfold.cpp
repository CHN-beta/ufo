# include <ufo/unfold.hpp>

namespace ufo
{
  UnfoldSolver::InputType::InputType(std::string filename)
  {
    // read main input file
    {
      auto node = YAML::LoadFile(filename);
        for (unsigned i = 0; i < 3; i++)
          for (unsigned j = 0; j < 3; j++)
            PrimativeCell(i, j) = node["PrimativeCell"][i][j].as<double>();

      for (unsigned i = 0; i < 3; i++)
        SuperCellMultiplier(i) = node["SuperCellMultiplier"][i].as<int>();

      if (auto value = node["SuperCellDeformation"])
      {
        SuperCellDeformation.emplace();
        for (unsigned i = 0; i < 3; i++)
          for (unsigned j = 0; j < 3; j++)
            (*SuperCellDeformation)(i, j) = value[i][j].as<double>();
      }

      for (unsigned i = 0; i < 3; i++)
        PrimativeCellBasisNumber(i) = node["PrimativeCellBasisNumber"][i].as<int>();

      auto read_file_config = [filename](YAML::Node source, InputOutputFile& config)
      {
        if (auto _ = source["SameAsConfigFile"])
        {
          auto __ = _.as<bool>();
          config.ExtraParameters["SameAsConfigFile"] = __;
          if (__)
          {
            config.FileName = filename;
            config.Format = "yaml";
            return;
          }
        }
        config.FileName = source["FileName"].as<std::string>();
        config.Format = source["Format"].as<std::string>();
        if (auto _ = source["RelativeToConfigFile"])
        {
          auto __ = _.as<bool>();
          config.ExtraParameters["RelativeToConfigFile"] = __;
          if (__)
            config.FileName = std::filesystem::path(filename).parent_path() / config.FileName;
        }
      };
      read_file_config(node["AtomPositionInputFile"], AtomPositionInputFile);
      if (!std::set<std::string>{"yaml"}.contains(AtomPositionInputFile.Format))
        throw std::runtime_error(fmt::format
          ("Unknown AtomPositionInputFile.Format: {}, should be \"yaml\".", AtomPositionInputFile.Format));
      read_file_config(node["QPointDataInputFile"], QPointDataInputFile);
      if (!std::set<std::string>{"yaml", "hdf5"}.contains(QPointDataInputFile.Format))
        throw std::runtime_error(fmt::format
          ("Unknown QPointDataInputFile.Format: {}, should be \"yaml\" or \"hdf5\".", QPointDataInputFile.Format));
      if (auto value = node["QPointDataOutputFile"])
      {
        QPointDataOutputFile.resize(value.size());
        for (unsigned i = 0; i < value.size(); i++)
        {
          read_file_config(value[i], QPointDataOutputFile[i]);
          if
          (
            QPointDataOutputFile[i].ExtraParameters.contains("SameAsConfigFile")
            && std::any_cast<bool>(QPointDataOutputFile[i].ExtraParameters["SameAsConfigFile"])
          )
            throw std::runtime_error("QPointDataOutputFile.SameAsConfigFile should not be set.");
          if
          (
            !std::set<std::string>{"yaml", "yaml-human-readable", "zpp"}
              .contains(QPointDataOutputFile[i].Format)
          )
            throw std::runtime_error(fmt::format
            (
              "Unknown QPointDataOutputFile[{}].Format: {}, should be \"yaml\", \"yaml-human-readable\" or \"zpp\".",
              i, QPointDataOutputFile[i].Format
            ));
        }
      }
    }

    if (AtomPositionInputFile.Format == "yaml")
    {
      auto node = YAML::LoadFile(AtomPositionInputFile.FileName);
      std::vector<YAML::Node> points;
      if (auto _ = node["points"])
        points = _.as<std::vector<YAML::Node>>();
      else
        points = node["unit_cell"]["points"].as<std::vector<YAML::Node>>();
      auto atom_position_to_super_cell = Eigen::MatrixX3d(points.size(), 3);
      for (unsigned i = 0; i < points.size(); i++)
        for (unsigned j = 0; j < 3; j++)
          atom_position_to_super_cell(i, j) = points[i]["coordinates"][j].as<double>();
      auto super_cell = (SuperCellDeformation.value_or(Eigen::Matrix3d::Identity())
        * SuperCellMultiplier.cast<double>().asDiagonal() * PrimativeCell).eval();
      AtomPosition = atom_position_to_super_cell * super_cell;
    }
    if (QPointDataInputFile.Format == "yaml")
    {
      auto node = YAML::LoadFile(QPointDataInputFile.FileName);
      auto phonon = node["phonon"].as<std::vector<YAML::Node>>();
      QPointData.resize(phonon.size());
      for (unsigned i = 0; i < phonon.size(); i++)
      {
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
    else if (QPointDataInputFile.Format == "hdf5")
    {
      HighFive::File file(QPointDataInputFile.FileName, HighFive::File::ReadOnly);
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
                  = eigenvector_vector[i][j][l * 3 + m][k].r + eigenvector_vector[i][j][l * 3 + m][k].i * 1i;
            QPointData[i * size[1] + j].ModeData[k].AtomMovement = eigenvectors / eigenvectors.norm();
          }
        }
    }
  }

  void UnfoldSolver::OutputType::write
    (decltype(InputType::QPointDataOutputFile) output_files) const
  {
    for (auto& output_file : output_files)
      write(output_file.FileName, output_file.Format);
  }
  void UnfoldSolver::OutputType::write(std::string filename, std::string format, unsigned percision) const
  {
    if (format == "yaml")
      std::ofstream(filename) << [&]
      {
        std::stringstream print;
        print << "QPointData:\n";
        for (auto& qpoint: QPointData)
        {
          print << fmt::format("  - QPoint: [ {1:.{0}f}, {2:.{0}f}, {3:.{0}f} ]\n",
            percision, qpoint.QPoint[0], qpoint.QPoint[1], qpoint.QPoint[2]);
          print << fmt::format("    Source: [ {1:.{0}f}, {2:.{0}f}, {3:.{0}f} ]\n",
            percision, qpoint.Source[0], qpoint.Source[1], qpoint.Source[2]);
          print << "    ModeData:\n";
          for (auto& mode: qpoint.ModeData)
            print << fmt::format("      - {{ Frequency: {1:.{0}f}, Weight: {2:.{0}f} }}\n",
              percision, mode.Frequency, mode.Weight);
        }
        return print.str();
      }();
    else if (format == "yaml-human-readable")
    {
      std::remove_cvref_t<decltype(*this)> output;
      std::map<unsigned, std::vector<decltype(QPointData)::const_iterator>>
        meta_qpoint_to_sub_qpoint_iterators;
      for (auto it = QPointData.begin(); it != QPointData.end(); it++)
        meta_qpoint_to_sub_qpoint_iterators[it->SourceIndex_].push_back(it);
      for (auto [meta_qpoint_index, sub_qpoint_iterators] : meta_qpoint_to_sub_qpoint_iterators)
        for (auto& qpoint : sub_qpoint_iterators)
        {
          std::map<double, double> frequency_to_weight;
          for (unsigned i_of_mode = 0; i_of_mode < qpoint->ModeData.size(); i_of_mode++)
          {
            auto frequency = qpoint->ModeData[i_of_mode].Frequency;
            auto weight = qpoint->ModeData[i_of_mode].Weight;
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
          auto& _ = output.QPointData.emplace_back();
          _.QPoint = qpoint->QPoint;
          _.Source = qpoint->Source;
          _.SourceIndex_ = qpoint->SourceIndex_;
          for (auto [frequency, weight] : frequency_to_weight)
            if (weight > 0.1)
            {
              auto& __ = _.ModeData.emplace_back();
              __.Frequency = frequency;
              __.Weight = weight;
            }
        }
      output.write(filename, "yaml", 3);
    }
    else if (format == "zpp")
    {
      auto [data, out] = zpp::bits::data_out();
      out(*this).or_throw();
      static_assert(sizeof(char) == sizeof(std::byte));
      std::ofstream file(filename, std::ios::binary | std::ios::out);
      file.exceptions(std::ios::badbit | std::ios::failbit);
      file.write(reinterpret_cast<const char*>(data.data()), data.size());
    }
  }

  UnfoldSolver::UnfoldSolver(std::string config_file) : Input_([&]
  {
    std::clog << "Reading input file... " << std::flush;
    return config_file;
  }())
  {
    std::clog << "Done." << std::endl;
  }

  UnfoldSolver& UnfoldSolver::operator()()
  {
    if (!Basis_)
    {
      std::clog << "Constructing basis... " << std::flush;
      Basis_ = construct_basis
      (
        Input_.PrimativeCell, Input_.SuperCellMultiplier,
        Input_.PrimativeCellBasisNumber, Input_.AtomPosition
      );
      std::clog << "Done." << std::endl;
    }
    if (!Output_)
    {
      std::clog << "Calculating projection coefficient... " << std::flush;
      std::vector<std::reference_wrapper<const decltype
        (InputType::QPointDataType::ModeDataType::AtomMovement)>> mode_data;
      for (auto& qpoint : Input_.QPointData)
        for (auto& mode : qpoint.ModeData)
          mode_data.emplace_back(mode.AtomMovement);
      std::atomic<unsigned> number_of_finished_modes(0);
      std::thread print_thread([&]
      {
        unsigned n;
        while ((n = number_of_finished_modes) < mode_data.size())
        {
          std::osyncstream(std::cerr) << fmt::format("\rCalculating projection coefficient... ({}/{})",
            number_of_finished_modes, mode_data.size()) << std::flush;
          std::this_thread::sleep_for(100ms);
          number_of_finished_modes.wait(n);
        }
      });
      auto projection_coefficient = construct_projection_coefficient
        (*Basis_, mode_data, number_of_finished_modes);
      number_of_finished_modes = mode_data.size();
      print_thread.join();
      std::clog << "\rCalculating projection coefficient... Done." << std::endl;

      std::clog << "Constructing output... " << std::flush;
      std::vector<std::reference_wrapper<const decltype(InputType::QPointDataType::QPoint)>> qpoint;
      std::vector<std::vector<std::reference_wrapper<const
        decltype(InputType::QPointDataType::ModeDataType::Frequency)>>> frequency;
      for (auto& qpoint_data : Input_.QPointData)
      {
        qpoint.emplace_back(qpoint_data.QPoint);
        frequency.emplace_back();
        for (auto& mode_data : qpoint_data.ModeData)
          frequency.back().emplace_back(mode_data.Frequency);
      }
      Output_ = construct_output
      (
        Input_.PrimativeCell, Input_.SuperCellMultiplier,
        Input_.SuperCellDeformation, qpoint, frequency, projection_coefficient
      );
      std::clog << "Done." << std::endl;
    }
    std::clog << "Writing output... " << std::flush;
    Output_->write(Input_.QPointDataOutputFile);
    std::clog << "Done." << std::endl;
    return *this;
  }

  UnfoldSolver::BasisType UnfoldSolver::construct_basis
  (
    const decltype(InputType::PrimativeCell)& primative_cell,
    const decltype(InputType::SuperCellMultiplier)& super_cell_multiplier,
    const decltype(InputType::PrimativeCellBasisNumber)& primative_cell_basis_number,
    const decltype(InputType::AtomPosition)& atom_position
  )
  {
    BasisType basis(super_cell_multiplier.prod());
    // 每个 q 点对应的一组 sub qpoint。不同的 q 点所对应的 sub qpoint 是不一样的，但 sub qpoint 与 q 点的相对位置一致。
    // 这里 xyz_of_diff_of_sub_qpoint 即表示这个相对位置，单位为超胞的倒格矢
    for (auto [xyz_of_diff_of_sub_qpoint_by_reciprocal_modified_super_cell, i_of_sub_qpoint]
      : triplet_sequence(super_cell_multiplier))
    {
      basis[i_of_sub_qpoint].resize(primative_cell_basis_number.prod());
      for (auto [xyz_of_basis, i_of_basis] : triplet_sequence(primative_cell_basis_number))
      {
        // 计算 q 点的坐标, 单位为单胞的倒格矢
        auto diff_of_sub_qpoint_by_reciprocal_primative_cell = xyz_of_basis.cast<double>()
          + super_cell_multiplier.cast<double>().cwiseInverse().asDiagonal()
          * xyz_of_diff_of_sub_qpoint_by_reciprocal_modified_super_cell.cast<double>();
        // 将 q 点坐标转换为埃^-1
        auto qpoint = (diff_of_sub_qpoint_by_reciprocal_primative_cell.transpose()
          * (primative_cell.transpose().inverse())).transpose();
        // 计算基矢
        basis[i_of_sub_qpoint][i_of_basis]
          = (2i * std::numbers::pi_v<double> * (atom_position * qpoint)).array().exp();
      }
    }
    return basis;
  }

  std::vector<std::vector<double>> UnfoldSolver::construct_projection_coefficient
  (
    const BasisType& basis,
    const std::vector<std::reference_wrapper<const decltype
      (InputType::QPointDataType::ModeDataType::AtomMovement)>>& mode_data,
    std::atomic<unsigned>& number_of_finished_modes
  )
  {
    // 第一层下标对应不同模式, 第二层下标对应这个模式在反折叠后的 q 点(sub qpoint)
    std::vector<std::vector<double>> projection_coefficient(mode_data.size());
    // 对每个模式并行
    std::transform
    (
      std::execution::par, mode_data.begin(), mode_data.end(), projection_coefficient.begin(),
      [&](const auto& mode_data)
      {
        // 这里, mode_data 和 projection_coefficient 均指对应于一个模式的数据
        std::vector<double> projection_coefficient(basis.size());
        for (unsigned i_of_sub_qpoint = 0; i_of_sub_qpoint < basis.size(); i_of_sub_qpoint++)
          // 对于 basis 中, 对应于单胞倒格子的部分, 以及对应于不同方向的部分, 分别求内积, 然后求模方和
          for (unsigned i_of_basis = 0; i_of_basis < basis[i_of_sub_qpoint].size(); i_of_basis++)
            projection_coefficient[i_of_sub_qpoint] +=
              (basis[i_of_sub_qpoint][i_of_basis].transpose().conjugate() * mode_data.get())
                .array().abs2().sum();
        // 如果是严格地将向量分解到一组完备的基矢上, 那么不需要对计算得到的权重再做归一化处理
        // 但这里并不是这样一个严格的概念. 因此对分解到各个 sub qpoint 上的权重做归一化处理
        auto sum = std::accumulate
          (projection_coefficient.begin(), projection_coefficient.end(), 0.);
        for (auto& _ : projection_coefficient)
          _ /= sum;
        number_of_finished_modes++;
        return projection_coefficient;
      }
    );
    return projection_coefficient;
  }

  UnfoldSolver::OutputType UnfoldSolver::construct_output
  (
    const decltype(InputType::PrimativeCell)& primative_cell,
    const decltype(InputType::SuperCellMultiplier)& super_cell_multiplier,
    const decltype(InputType::SuperCellDeformation)& super_cell_deformation,
    const std::vector<std::reference_wrapper<const decltype
      (InputType::QPointDataType::QPoint)>>& qpoint,
    const std::vector<std::vector<std::reference_wrapper<const decltype
      (InputType::QPointDataType::ModeDataType::Frequency)>>>& frequency,
    const ProjectionCoefficientType_& projection_coefficient
  )
  {
    OutputType output;
    for (unsigned i_of_qpoint = 0, num_of_mode_in_all_qpoint = 0; i_of_qpoint < qpoint.size(); i_of_qpoint++)
    {
      // 当 SuperCellDeformation 不是单位矩阵时, input.QPointData[i_of_qpoint].QPoint 不一定在 reciprocal_primative_cell 中
      // 需要首先将 q 点平移数个周期, 进入不包含 SuperCellDeformation 的超胞的倒格子中
      auto qpoint_by_reciprocal_modified_super_cell_in_modified_reciprocal_super_cell
        = !super_cell_deformation ? qpoint[i_of_qpoint].get() : [&]
        {
          auto current_qpoint = qpoint[i_of_qpoint].get();
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
              (super_cell_deformation->inverse() * qpoint_by_reciprocal_super_cell).eval();
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
                  auto this_qpoint = current_qpoint
                    + Eigen::Matrix<int, 3, 1>{{x}, {y}, {z}}.cast<double>();
                  auto this_score = score(this_qpoint);
                  if (this_score < min_score)
                  {
                    min_score = this_score;
                    min_score_qpoint = this_qpoint;
                  }
                }
            current_qpoint = min_score_qpoint;
          }
          return super_cell_deformation->inverse() * current_qpoint;
        }();
      for (auto [xyz_of_diff_of_sub_qpoint_by_reciprocal_modified_super_cell, i_of_sub_qpoint]
        : triplet_sequence(super_cell_multiplier))
      {
        auto& _ = output.QPointData.emplace_back();
        auto reciprocal_modified_super_cell =
          (super_cell_multiplier.cast<double>().asDiagonal() * primative_cell).inverse().transpose();
        // sub qpoint 的坐标，单位为埃^-1
        auto sub_qpoint = ((xyz_of_diff_of_sub_qpoint_by_reciprocal_modified_super_cell.cast<double>()
          + qpoint_by_reciprocal_modified_super_cell_in_modified_reciprocal_super_cell)
          .transpose() * reciprocal_modified_super_cell).transpose();
        // 将坐标转换为相对于单胞的倒格矢的坐标并写入
        // 由 sub_qpoint.transpose() = sub_qpoint_by_reciprocal_primative_cell.transpose()
        //  * PrimativeCell.transpose().inverse()
        // 得到 sub_qpoint_by_reciprocal_primative_cell = PrimativeCell * sub_qpoint
        _.QPoint = primative_cell * sub_qpoint;
        _.Source = qpoint[i_of_qpoint];
        _.SourceIndex_ = i_of_qpoint;
        for (unsigned i_of_mode = 0; i_of_mode < frequency[i_of_qpoint].size(); i_of_mode++, num_of_mode_in_all_qpoint++)
        {
          auto& __ = _.ModeData.emplace_back();
          __.Frequency = frequency[i_of_qpoint][i_of_mode];
          __.Weight = projection_coefficient[num_of_mode_in_all_qpoint][i_of_sub_qpoint];
        }
      }
    }
    return output;
  }
}

// inline Output::Output(std::string filename)
// {
//   auto input = std::ifstream(filename, std::ios::binary | std::ios::in);
//   input.exceptions(std::ios::badbit | std::ios::failbit);
//   std::vector<std::byte> data;
//   {
//     std::vector<char> string(std::istreambuf_iterator<char>(input), {});
//     data.assign
//     (
//       reinterpret_cast<std::byte*>(string.data()),
//       reinterpret_cast<std::byte*>(string.data() + string.size())
//     );
//   }
//   auto in = zpp::bits::in(data);
//   in(*this).or_throw();
// }

