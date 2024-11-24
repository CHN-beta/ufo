# include <ufo.hpp>

void ufo::raman_apply_contribution(std::string config_file)
{
  struct Input
  {
    std::string UnfoldedDataFile;
    std::string DisplacementDataFile;
    Eigen::Matrix3d OriginalSusceptibility;
    std::vector<Eigen::Matrix3d> Susceptibilities;
    std::string OutputDataFile;
    std::array<Eigen::Vector3d, 2> Polarization;
  };

  biu::Logger::Guard log(config_file);
  auto input = YAML::LoadFile(config_file).as<Input>();
  auto unfolded_data = biu::deserialize<UnfoldOutput>
    (biu::read<std::byte>(input.UnfoldedDataFile));
  auto displacement_data = biu::deserialize<DisplacementOutput>
    (biu::read<std::byte>(input.DisplacementDataFile));
  UnfoldOutput output;

  // 整理得到的数据，作放缩
  struct mode_t
  {
    std::size_t MetaQpointIndex;
    std::size_t ModeIndex;
    Eigen::Matrix3d RamanTensor;
    double Strength;
  };
  auto modes = ranges::views::zip(displacement_data.ModeData, input.Susceptibilities)
    | ranges::views::transform([&](const auto& data)
      {
        auto& [mode, susceptibility] = data;
        Eigen::Matrix3d raman_tensor = susceptibility / mode.Ratio / displacement_data.MaxDisplacement;
        double intensity = (input.Polarization[0].transpose() * raman_tensor * input.Polarization[1]).norm();
        log.debug("get raman tensor: {}"_f(raman_tensor));
        log.debug("get intensity: {}"_f(intensity));
        return mode_t{ mode.MetaQpointIndex, mode.ModeIndex, raman_tensor, intensity };
      })
    | ranges::to_vector;

  // 整理输出文件
  biu::for_each
  (
    [&](auto&& i) { output.*i = unfolded_data.*i; },
    std::tuple(&UnfoldOutput::PrimativeCell, &UnfoldOutput::SuperCellTransformation,
      &UnfoldOutput::SuperCellMultiplier, &UnfoldOutput::SuperCellDeformation, &UnfoldOutput::SelectedAtoms)
  );
  // 将以前的 meta qpoint 和 mode 的索引转换到新的
  std::map<std::size_t, std::size_t> meta_qpoint_map;
  std::vector<std::map<std::size_t, std::size_t>> mode_map;
  for (auto mode : displacement_data.ModeData)
  {
    auto old_meta_qpoint_index = mode.MetaQpointIndex;
    auto old_mode_index = mode.ModeIndex;
    auto new_meta_qpoint_index = meta_qpoint_map.contains(old_meta_qpoint_index)
      ? meta_qpoint_map[mode.MetaQpointIndex]
      : [&]
      {
        log.debug("map meta qpoint {} to {}"_f(old_meta_qpoint_index, meta_qpoint_map.size()));
        meta_qpoint_map[mode.MetaQpointIndex] = meta_qpoint_map.size();
        output.MetaQpointData.push_back({unfolded_data.MetaQpointData[old_mode_index].Qpoint, {}});
        return meta_qpoint_map.size() - 1;
      }();
    auto new_mode_index = mode_map[new_meta_qpoint_index].size();
    log.debug("map mode {} {} to {}"_f(old_meta_qpoint_index, old_mode_index, new_mode_index));
    mode_map[new_meta_qpoint_index][old_mode_index] = new_mode_index;
    output.MetaQpointData[new_meta_qpoint_index].ModeData.push_back
      (unfolded_data.MetaQpointData[old_meta_qpoint_index].ModeData[old_mode_index]);
  }
  for (auto old_qpoint_index : displacement_data.QpointIndices)
    if (meta_qpoint_map.contains(unfolded_data.QpointData[old_qpoint_index].SourceIndex))
    {
      log.debug("map qpoint {} {}"_f(old_qpoint_index, unfolded_data.QpointData[old_qpoint_index].Qpoint));
      auto new_meta_qpoint_index
        = meta_qpoint_map[unfolded_data.QpointData[old_qpoint_index].SourceIndex];
      output.QpointData.push_back
      ({
        unfolded_data.QpointData[old_qpoint_index].Qpoint,
        unfolded_data.QpointData[old_qpoint_index].Source,
        new_meta_qpoint_index,
        {}
      });
      for (std::size_t i = 0; i < unfolded_data.QpointData[old_qpoint_index].ModeData.size(); i++)
        if (mode_map[new_meta_qpoint_index].contains(i))
        {
          auto new_mode_index = mode_map[new_meta_qpoint_index][i];
          log.debug("map mode {} to {}"_f(i, new_mode_index));
          if (output.QpointData.back().ModeData.size() <= new_mode_index)
          {
            log.debug("resize to {}"_f(new_mode_index + 1));
            output.QpointData.back().ModeData.resize(new_mode_index + 1);
          }
          output.QpointData.back().ModeData[new_mode_index] =
          {
            unfolded_data.QpointData[old_qpoint_index].ModeData[i].Frequency,
            unfolded_data.QpointData[old_qpoint_index].ModeData[i].Weight * modes[new_mode_index].Strength
          };
        }
    }

  std::ofstream(input.OutputDataFile, std::ios::binary) << biu::serialize<char>(output);
}
