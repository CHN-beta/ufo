# pragma once
# include <ufo/solver.impl.hpp>
# include <ufo/unfold.hpp>

namespace ufo
{
  inline UnfoldSolver::InputType::InputType(std::string filename)
  {
    // read main input file
    {
      auto node = YAML::LoadFile(filename);
      PrimativeCell = node["PrimativeCell"].as<std::vector<std::vector<double>>>() | biu::toEigen<3, 3>;
      SuperCellMultiplier = node["SuperCellMultiplier"].as<std::vector<unsigned>>() | toEigenVector<3>;
      if (auto value = node["SuperCellDeformation"])
        SuperCellDeformation = value.as<std::vector<std::vector<double>>>() | toEigenMatrix<3, 3>;
      PrimativeCellBasisNumber = node["PrimativeCellBasisNumber"].as<std::vector<unsigned>>()
        | toEigenVector<3>;
      AtomPositionInputFile = DataFile
      (
        node["AtomPositionInputFile"], {"yaml"},
        filename, true
      );
      QpointDataInputFile = DataFile
      (
        node["QpointDataInputFile"], {"yaml", "hdf5"},
        filename, true
      );
      if (auto value = node["QpointDataOutputFile"])
        for (auto&& v : value)
          QpointDataOutputFile.push_back(DataFile
          (
            v, {"yaml", "yaml-human-readable", "zpp", "hdf5"},
            filename, false
          ));
    }

    if (AtomPositionInputFile.Format == "yaml")
    {
      auto node = YAML::LoadFile(AtomPositionInputFile.Filename);
      YAML::Node source = node["points"];
      if (!source)
        source = node["unit_cell"]["points"];
      auto atom_position_to_super_cell = source
        | std::views::transform([](YAML::Node&& v)
          { return v["coordinates"].as<std::vector<double>>(); })
        | std::
        | toEigenMatrix<std::dynamic_extent, 3>

      auto atom_position_to_super_cell = Eigen::MatrixX3d(source, 3);
      for (unsigned i = 0; i < source.size(); i++)
        atom_position_to_super_cell.row(i) = source[i]["coordinates"].as<std::vector<double>>()
          | toEigenVector<3>;
      auto super_cell = (SuperCellDeformation.value_or(Eigen::Matrix3d::Identity())
        * SuperCellMultiplier.cast<double>().asDiagonal() * PrimativeCell).eval();
      AtomPosition = atom_position_to_super_cell * super_cell;
    }
    if (QpointDataInputFile.Format == "yaml")
    {
      auto node = YAML::LoadFile(QpointDataInputFile.Filename);
      QpointData = node["phonon"].as<std::vector<YAML::Node>>()
        | std::views::transform([](YAML::Node&& v)
          {
            return QpointDataType
            {
              .Qpoint = v["q-position"].as<std::vector<double>>() | toEigenVector<3>,
              .ModeData = std::vector(
                v["band"].as<std::vector<YAML::Node>>()
                | std::views::transform([](YAML::Node&& v)
                  {
                    return QpointDataType::ModeDataType
                    {
                      .Frequency = v["frequency"].as<double>(),
                      .EigenVector = 
                      v["frequency"].as<double>(),
                      Eigen::MatrixX3cd(v["eigenvector"]
                        .as<std::vector<std::vector<std::vector<double>>>>()
                        | std::views::transform([](auto&& v)
                          {
                            return Eigen::Vector3cd(v | toEigenVector<2>);
                          })
                        | toEigenMatrix)
                    };
                  })
                | toStdVector
              )

              )

              v["q-position"].as<std::vector<double>>() | toEigenVector<3>,
              v["band"].as<std::vector<YAML::Node>>()
              | std::views::transform([](auto&& v)
                {
                  return ModeDataType
                  {
                    v["frequency"].as<double>(),
                    Eigen::MatrixX3cd(v["eigenvector"]
                      .as<std::vector<std::vector<std::vector<double>>>>()
                      | std::views::transform([](auto&& v)
                        {
                          return Eigen::Vector3cd(v | toEigenVector<2>);
                        })
                      | toEigenMatrix)
                  };
                })
              | toStdVector
            };

            QpointDataType output;
            for (unsigned i = 0; i < 3; i++)
              output.Qpoint(i) = v["q-position"][i].as<double>();
            output.ModeData = v["band"].as<std::vector<YAML::Node>>()
              | std::views::transform([](auto&& v)
                {
                  ModeDataType output;
                  output.Frequency = v["frequency"].as<double>();
                  auto eigenvector_vectors = v["eigenvector"]
                    .as<std::vector<std::vector<std::vector<double>>>>();
                  Eigen::MatrixX3cd eigenvectors(AtomPosition.rows(), 3);
                  for (unsigned i = 0; i < AtomPosition.rows(); i++)
                    for (unsigned j = 0; j < 3; j++)
                      eigenvectors(i, j)
                        = eigenvector_vectors[i][j][0] + 1i * eigenvector_vectors[i][j][1];
                  output.EigenVector = eigenvectors / eigenvectors.norm();
                  return output;
                })
              | toStdVector;
            return output;
          })


      auto phonon = node["phonon"].as<std::vector<YAML::Node>>();




      QpointData = std::vector<QpointDataType>(phonon.size());
      QpointData.resize(phonon.size());
      for (unsigned i = 0; i < phonon.size(); i++)
      {
        for (unsigned j = 0; j < 3; j++)
          QpointData[i].Qpoint(j) = phonon[i]["q-position"][j].as<double>();
        auto band = phonon[i]["band"].as<std::vector<YAML::Node>>();
        QpointData[i].ModeData.resize(band.size());
        for (unsigned j = 0; j < band.size(); j++)
        {
          QpointData[i].ModeData[j].Frequency = band[j]["frequency"].as<double>();
          auto eigenvector_vectors = band[j]["eigenvector"]
            .as<std::vector<std::vector<std::vector<double>>>>();
          Eigen::MatrixX3cd eigenvectors(AtomPosition.rows(), 3);
          for (unsigned k = 0; k < AtomPosition.rows(); k++)
            for (unsigned l = 0; l < 3; l++)
              eigenvectors(k, l)
                = eigenvector_vectors[k][l][0] + 1i * eigenvector_vectors[k][l][1];
          // 需要对读入的原子运动状态作相位转换, 使得它们与我们的约定一致(对超胞周期性重复)
          // 这里还要需要做归一化处理 (指将数据简单地作为向量处理的归一化)
          auto& AtomMovement = QpointData[i].ModeData[j].AtomMovement;
          // AtomMovement = eigenvectors.array().colwise() * (-2 * std::numbers::pi_v<double> * 1i
          //   * (atom_position_to_super_cell * input.QpointData[i].Qpoint)).array().exp();
          // AtomMovement /= AtomMovement.norm();
          // phonopy 似乎已经进行了相位的转换！为什么？
          AtomMovement = eigenvectors / eigenvectors.norm();
        }
      }
    }
    else if (QpointDataInputFile.Format == "hdf5")
    {
      std::vector<std::vector<std::vector<double>>> frequency, path;
      std::vector<std::vector<std::vector<std::vector<PhonopyComplex>>>> eigenvector_vector;
      Hdf5file{}.open_for_read(QpointDataInputFile.Filename).read(frequency, "/frequency")
        .read(eigenvector_vector, "/eigenvector")
        .read(path, "/path");
      std::vector size = { frequency.size(), frequency[0].size(), frequency[0][0].size() };
      QpointData.resize(size[0] * size[1]);
      for (unsigned i = 0; i < size[0]; i++)
        for (unsigned j = 0; j < size[1]; j++)
        {
          QpointData[i * size[1] + j].Qpoint = Eigen::Vector3d(path[i][j].data());
          QpointData[i * size[1] + j].ModeData.resize(size[2]);
          for (unsigned k = 0; k < size[2]; k++)
          {
            QpointData[i * size[1] + j].ModeData[k].Frequency = frequency[i][j][k];
            Eigen::MatrixX3cd eigenvectors(AtomPosition.rows(), 3);
            for (unsigned l = 0; l < AtomPosition.rows(); l++)
              for (unsigned m = 0; m < 3; m++)
                eigenvectors(l, m)
                  = eigenvector_vector[i][j][l * 3 + m][k].r + eigenvector_vector[i][j][l * 3 + m][k].i * 1i;
            QpointData[i * size[1] + j].ModeData[k].AtomMovement = eigenvectors / eigenvectors.norm();
          }
        }
    }
  }



}
