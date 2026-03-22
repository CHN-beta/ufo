# include <ufo.hpp>

int main(int argc, const char** argv)
{
  using namespace biu::literals;
  if (argv[1] == "fold"s) ufo::fold(argv[2]);
  else if (argv[1] == "unfold"s) ufo::unfold(argv[2]);
  else if (argv[1] == "zpp-to-yaml"s) std::cout << YAML::Node(biu::deserialize<ufo::CommonData>
    (biu::read<std::byte>(std::cin)));
  else if (argv[1] == "project-to-atom"s) ufo::project_to_atom(argv[2]);
  else if (argv[1] == "project-to-mode"s) ufo::project_to_mode(argv[2]);
  else if (argv[1] == "raman-create-displacement"s) ufo::raman_create_displacement(argv[2]);
  else if (argv[1] == "raman-extract"s) ufo::raman_extract(argv[2]);
  else if (argv[1] == "raman-apply-contribution"s) ufo::raman_apply_contribution(argv[2]);
  else if (argv[1] == "plot-band"s) ufo::plot_band(argv[2]);
  else if (argv[1] == "plot-point"s) ufo::plot_point(argv[2]);
  else if (argv[1] == "raman-bec-rotate"s) ufo::raman_bec_rotate(argv[2]);
  else if (argv[1] == "raman-bec-apply-contribution"s) ufo::raman_bec_apply_contribution(argv[2]);
  else if (argv[1] == "raman-bec-extract"s) ufo::raman_bec_extract(argv[2]);
  else throw std::runtime_error("Unknown task: {}"_f(argv[1]));
}
