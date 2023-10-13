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
      PrimativeCell = [](YAML::Node node)
      {
        Eigen::Matrix3d matrix;
        for (unsigned i = 0; i < 3; i++)
          for (unsigned j = 0; j < 3; j++)
            matrix(i, j) = node[i][j].as<double>();
        return matrix;
      }(node["PrimativeCell"]);

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
      {
        QpointDataOutputFile.resize(value.size());
        for (unsigned i = 0; i < value.size(); i++)
          QpointDataOutputFile[i] = DataFile
          (
            value[i], {"yaml", "yaml-human-readable", "zpp", "hdf5"},
            filename, false
          );
      }
    }

    if (AtomPositionInputFile.Format == "yaml")
    {
      auto node = YAML::LoadFile(AtomPositionInputFile.Filename);
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
    if (QpointDataInputFile.Format == "yaml")
    {
      auto node = YAML::LoadFile(QpointDataInputFile.Filename);
      auto phonon = node["phonon"].as<std::vector<YAML::Node>>();
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
