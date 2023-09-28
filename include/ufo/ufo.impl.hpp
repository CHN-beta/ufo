# include <ufo/ufo.hpp>

inline HighFive::CompoundType create_compound_complex()
  { return {{"r", HighFive::AtomicType<double>{}}, {"i", HighFive::AtomicType<double>{}}}; }

inline Input::Input(std::string filename)
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

    auto read_file_config = [filename](YAML::Node source, InputOutputFile_& config)
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
    Eigen::Matrix3d super_cell = SuperCellDeformation.value_or(Eigen::Matrix3d::Identity())
      * SuperCellMultiplier.cast<double>().asDiagonal() * PrimativeCell;
    std::cout << "SuperCell:\n" << super_cell << std::endl;
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

inline void Output::write(std::string filename, std::string format, unsigned percision) const
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
    Output output;
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

inline Output::Output(std::string filename)
{
  auto input = std::ifstream(filename, std::ios::binary | std::ios::in);
  input.exceptions(std::ios::badbit | std::ios::failbit);
  std::vector<std::byte> data;
  {
    std::vector<char> string(std::istreambuf_iterator<char>(input), {});
    data.assign
    (
      reinterpret_cast<std::byte*>(string.data()),
      reinterpret_cast<std::byte*>(string.data() + string.size())
    );
  }
  auto in = zpp::bits::in(data);
  in(*this).or_throw();
}

inline concurrencpp::generator<std::pair<Eigen::Vector<unsigned, 3>, unsigned>>
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
