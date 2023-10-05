# include <ufo/plot.hpp>

int main(int argc, const char** argv)
{
  if (argc != 2)
    throw std::runtime_error(fmt::format("Usage: {} task config.yaml", argv[0]));
  if (argv[1] == std::string("unfold"))
    ufo::UnfoldSolver{argv[1]}();
  else if (argv[1] == std::string("plot"))
    ufo::PlotSolver{argv[1]}();
  else
    throw std::runtime_error(fmt::format("Unknown task: {}", argv[1]));
}
