# include <ufo.hpp>
# include <ranges>

void ufo::project_to_atom(std::string config_file)
{
  struct Config
  {
    std::set<std::size_t> SelectedAtom;
    std::string InputFile;
    std::string OutputFile;
  };

  biu::Logger::Guard log(config_file);
  auto config = YAML::LoadFile(config_file).as<Config>();
  auto input = biu::deserialize<CommonData>(biu::read<std::byte>(config.InputFile));
  input.SelectedAtom = config.SelectedAtom;
  for (auto& qpoint : input.Super.Qpoint) for (auto& mode : qpoint.Mode)
    mode.WeightOnSelectedAtom = mode.EigenVector.rowwise().squaredNorm().cwiseProduct
    (
      ranges::views::iota(0ul, input.Super.AtomPosition.rows() * 1ul)
        | ranges::views::transform([&](auto i)
          { return config.SelectedAtom.contains(i) ? 1. : 0.; })
        | ranges::to_vector | biu::toEigen<>
    ).sum() * (input.Super.AtomPosition.rows() * 1. / config.SelectedAtom.size());
  std::ofstream(config.OutputFile, std::ios::binary) << biu::serialize<char>(input);
  log.info("Summary:");
  for (auto i_of_qpoint : std::views::iota(0ul, input.Super.Qpoint.size()))
    for (auto i_of_mode : std::views::iota(0ul, input.Super.Qpoint[i_of_qpoint].Mode.size()))
      if (*input.Super.Qpoint[i_of_qpoint].Mode[i_of_mode].WeightOnSelectedAtom > 0.1)
        log.info("{}:{:.2f}:{}:{:.2f}"_f
        (
          i_of_qpoint, fmt::join(input.Super.Qpoint[i_of_qpoint].Qpoint, ", "),
          i_of_mode, input.Super.Qpoint[i_of_qpoint].Mode[i_of_mode].Frequency,
          *input.Super.Qpoint[i_of_qpoint].Mode[i_of_mode].WeightOnSelectedAtom
        ));
}
