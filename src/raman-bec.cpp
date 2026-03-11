# include <ufo.hpp>

namespace ufo
{
  void raman_bec_rotate(std::string config_file)
  {
    struct Config
    {
      Eigen::Matrix3d Cell;
    };

    biu::Logger::Guard log(config_file);
    auto config = YAML::LoadFile(config_file).as<Config>();
    std::map<std::string, Eigen::Matrix3d> rotation_matrix =
    {
      // 这里的xyz指的是旋转轴而不是施加电场的方向
      // 绕 z、x、y 轴旋转的分别只取旋转后 x、y、z 轴的电场结果
      // 我们认为原子和格矢不动、坐标轴顺时针（反向）旋转45度；因此，数值上来说，相当于格矢逆时针旋转45度
      {"original", Eigen::Matrix3d::Identity()},
      {
        "x",
        Eigen::Matrix3d
        {
          {1, 0, 0},
          {0, 1/std::sqrt(2), 1/std::sqrt(2)},
          {0, -1/std::sqrt(2), 1/std::sqrt(2)}
        }
      },
      {
        "y", 
        Eigen::Matrix3d
        {
          {1/std::sqrt(2), 0, -1/std::sqrt(2)},
          {0, 1, 0},
          {1/std::sqrt(2), 0, 1/std::sqrt(2)}
        }
      },
      {
        "z",
        Eigen::Matrix3d
        {
          {1/std::sqrt(2), 1/std::sqrt(2), 0},
          {-1/std::sqrt(2), 1/std::sqrt(2), 0},
          {0, 0, 1}
        }
      }
    };
    for (auto&& [name, rotation] : rotation_matrix)
    {
      auto rotated_cell = (rotation * config.Cell.transpose()).transpose().eval();
      log.info("rotation: {}"_f(name));
      for (std::size_t i = 0; i < 3; i++)
        log.info("{:0.16f} {:0.16f} {:0.16f}"_f
          (rotated_cell(i, 0), rotated_cell(i, 1), rotated_cell(i, 2)));
    }
  }

  void raman_bec_apply_contribution(std::string config_file)
  {
    struct Config
    {
      std::string InputDataFile;
      std::string OutputDataFile;
      std::set<std::size_t> SelectedQpoints;
      std::array<Eigen::Vector3d, 2> Polarization;
      // 4种旋转（原始、绕x、绕y、绕z）和2种电场方向（正、负）
      std::array<std::array<std::string, 2>, 4> VaspOutputFiles;
      double ElectricFieldStrength;
    };

    biu::Logger::Guard log(config_file);
    auto config = YAML::LoadFile(config_file).as<Config>();
    auto input = biu::deserialize<CommonData>
      (biu::read<std::byte>(config.InputDataFile));

    // read bec data
    // 下标：[旋转][电场正负][原子序号](位移方向, 电场方向)
    // Z_ij -> bec_data[][0][](j, i)
    // Z_-i,j -> bec_data[][1][](j, i)
    std::array<std::array<std::vector<Eigen::Matrix3d>, 2>, 4> bec_data;
    for (std::size_t i = 0; i < 4; i++)
      for (std::size_t j = 0; j < 2; j++)
      {
        biu::Hdf5file(config.VaspOutputFiles[i][j])
          .read("/results/linear_response/born_charges", bec_data[i][j]);
        for (std::size_t k = 0; k < bec_data[i][j].size(); k++)
          log.debug("read BEC data: {} {} {} {}"_f(i, j,  k, bec_data[i][j][k]));
        if (bec_data[i][j].size() != input.Super.AtomPosition.rows())
          throw std::runtime_error("Mismatch in BEC data size: {} vs {}"_f
            (bec_data[i][j].size(), input.Super.AtomPosition.rows()));
        // 将一个模型中所有原子的 BEC 加上一个偏移，使得它们求和为零
        auto offset = *std::ranges::fold_left_first
        (
          bec_data[i][j],
          [](auto&& a, auto&& b) { return (a + b).eval(); }
        ) / bec_data[i][j].size();
        for (std::size_t k = 0; k < bec_data[i][j].size(); k++) bec_data[i][j][k] -= offset;
      }

    // calculate atom raman
    // 下标：[原子序号][位移方向](电场方向，电场方向)
    // alpha_ijk -> atom_raman[][k](i, j)
    std::vector<std::array<Eigen::Matrix3d, 3>> atom_raman(input.Super.AtomPosition.rows());
    for (std::size_t i = 0; i < input.Super.AtomPosition.rows(); i++)
    {
      // 对角元素
      for (std::size_t j = 0; j < 3; j++) // 电场方向
        for (std::size_t k = 0; k < 3; k++) // 位移方向
          atom_raman[i][k](j, j) =
            (bec_data[0][0][i](k, j) - bec_data[0][1][i](k, j)) / config.ElectricFieldStrength;
      // 非对角元素
      for (std::size_t j = 0; j < 3; j++) // 旋转方向
      {
        atom_raman[i][j]((j+1)%3, (j+2)%3) = atom_raman[i][j]((j+2)%3, (j+1)%3) =
          1 / config.ElectricFieldStrength / std::sqrt(2)
            * (bec_data[j+1][0][i](j, (j+1)%3) - bec_data[j+1][1][i](j, (j+1)%3))
          - 0.5 * (atom_raman[i][j]((j+1)%3, (j+1)%3) + atom_raman[i][j]((j+2)%3, (j+2)%3));
        atom_raman[i][(j+1)%3]((j+1)%3, (j+2)%3) = atom_raman[i][(j+1)%3]((j+2)%3, (j+1)%3) =
          1 / config.ElectricFieldStrength / 2 *
          (
            bec_data[j+1][0][i]((j+1)%3, (j+1)%3) - bec_data[j+1][1][i]((j+1)%3, (j+1)%3)
              - bec_data[j+1][0][i]((j+2)%3, (j+1)%3) + bec_data[j+1][1][i]((j+2)%3, (j+1)%3)
          )
          - 0.5 *
            (atom_raman[i][(j+1)%3]((j+1)%3, (j+1)%3) + atom_raman[i][(j+1)%3]((j+2)%3, (j+2)%3));
        atom_raman[i][(j+2)%3]((j+1)%3, (j+2)%3) = atom_raman[i][(j+2)%3]((j+2)%3, (j+1)%3) =
          1 / config.ElectricFieldStrength / 2 *
          (
            bec_data[j+1][0][i]((j+1)%3, (j+1)%3) - bec_data[j+1][1][i]((j+1)%3, (j+1)%3)
              + bec_data[j+1][0][i]((j+2)%3, (j+1)%3) - bec_data[j+1][1][i]((j+2)%3, (j+1)%3)
          )
          - 0.5 *
            (atom_raman[i][(j+2)%3]((j+1)%3, (j+1)%3) + atom_raman[i][(j+2)%3]((j+2)%3, (j+2)%3));
      }
    }

    // 写入结果
    auto output = input;
    auto mass = input.Super.AtomType
      | ranges::views::transform([&](auto&& atom)
        { return ranges::views::repeat_n(input.AtomMass[atom.first], atom.second); })
      | ranges::views::join | ranges::to_vector;
    for (auto qpoint_index : config.SelectedQpoints)
      for (auto& mode : output.Super.Qpoint[qpoint_index].Mode)
      {
        Eigen::Matrix3d raman_tensor = Eigen::Matrix3d::Zero();
        for (std::size_t i = 0; i < output.Super.AtomPosition.rows(); i++)
          for (std::size_t j = 0; j < 3; j++)
            raman_tensor += mode.EigenVector(i, j).real() / std::sqrt(mass[i]) * atom_raman[i][j];
        mode.RamanTensor = raman_tensor | biu::fromEigen;
        mode.WeightOnRaman = config.Polarization[0].transpose() * raman_tensor * config.Polarization[1];
      }
    std::ofstream(config.OutputDataFile, std::ios::binary) << biu::serialize<char>(output);
  }
}
