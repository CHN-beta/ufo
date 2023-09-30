# include <ufo/unfold.hpp>

int main(int argc, const char** argv)
{
  if (argc != 2)
    throw std::runtime_error(fmt::format("Usage: {} config.yaml", argv[0]));

  ufo::UnfoldSolver{argv[1]}();
}
