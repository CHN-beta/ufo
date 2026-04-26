# include <ufo.hpp>

void ufo::ipr(std::string config_file)
{
  struct Config { std::string InputDataFile, OutputFile; };
  struct Output { std::vector<std::vector<double>> Primative, Super; };

  biu::Logger::Guard log(config_file);
  auto config = YAML::LoadFile(config_file).as<Config>();
  auto input = biu::deserialize<CommonData>
    (biu::read<std::byte>(config.InputDataFile));
  auto calc_ipr = [](const auto& cell)
  {
    return cell.Qpoint
      | std::views::transform([](const auto& qpoint)
        {
          return qpoint.Mode
            | std::views::transform([](const auto& mode)
              { return mode.EigenVector.rowwise().squaredNorm().squaredNorm(); })
            | std::ranges::to<std::vector<double>>();
        })
      | std::ranges::to<std::vector<std::vector<double>>>();
  };
  Output output { .Primative = calc_ipr(input.Primative), .Super = calc_ipr(input.Super) };
  std::ofstream(config.OutputFile) << YAML::Node(output);
}
