# include <ufo/solver.hpp>

namespace ufo
{
  concurrencpp::generator<std::pair<Eigen::Vector<unsigned, 3>, unsigned>> Solver::triplet_sequence
    (Eigen::Vector<unsigned, 3> range)
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
}
