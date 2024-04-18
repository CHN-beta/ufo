# include <ufo/fold.hpp>
# include <ufo/unfold.hpp>
# include <ufo/plot.hpp>

int main(int argc, const char** argv)
{
  if (argc != 3)
    throw std::runtime_error(fmt::format("Usage: {} task config.yaml", argv[0]));
  if (argv[1] == std::string("fold"))
    ufo::FoldSolver{argv[2]}();
  else if (argv[1] == std::string("unfold"))
    ufo::UnfoldSolver{argv[2]}();
  else if (argv[1] == std::string("plot"))
    ufo::PlotSolver{argv[2]}();
  else if (argv[1] == std::string("special"))
  {
    std::vector<std::vector<double>>
      r(400, std::vector<double>(1024, 0)),
      g(400, std::vector<double>(1024, 0)),
      b(400, std::vector<double>(1024, 0)),
      a(400, std::vector<double>(1024, 0)),
      values1,  // 缺陷中心和附近原子的值
      values2;  // 远离缺陷的原子的值
    
    ufo::Solver::Hdf5file()
      .open_for_read("/home/chn/Documents/lammps-SiC/14/14.5/14.5.8/14.5.8.1/D141/plot.hdf5")
      .read(values1, "Values");
    ufo::Solver::Hdf5file()
      .open_for_read("/home/chn/Documents/lammps-SiC/14/14.5/14.5.8/14.5.8.2/D141/plot.hdf5")
      .read(values2, "Values");

    for (unsigned i = 0; i < 400; i++)
      for (unsigned j = 0; j < 1024; j++)
      {
        auto v1 = values1[j][i] / 72 * 256;
        auto v2 = values2[j][i] / 184 * 256;
        if (v1 < 0.05) v1 = 0;
        if (v2 < 0.05) v2 = 0;
        a[i][j] = v1 * 100 * 255 + v2 * 100 * 255;
        if (a[i][j] > 255)
          a[i][j] = 255;
        r[i][j] = 255 - v2 * 255;
        if (r[i][j] < 0)
          r[i][j] = 0;
        g[i][j] = 255 - v1 * 255 - v2 * 255;
        if (g[i][j] < 0)
          g[i][j] = 0;
        b[i][j] = 255 - v1 * 255;
      }
    auto f = matplot::figure<matplot::backend::gnuplot>(true);
    auto ax = f->current_axes();
    auto image = ax->image(std::tie(r, g, b));
    image->matrix_a(a);
    ax->y_axis().reverse(false);
    ax->x_axis().tick_values({ 178.75, 281.95, 488.36, 535.64, 714.40, 817.60 });
    ax->x_axis().tick_length(1);
    ax->y_axis().tick_values({ 50, 100, 150, 200, 250, 300, 350 });
    ax->y_axis().tick_length(1);
    f->save("output.png", "png");
  }
  else
    throw std::runtime_error(fmt::format("Unknown task: {}", argv[1]));
}
