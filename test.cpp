# include <eigen3/Eigen/Dense>
# include <zpp_bits.h>
# include <span>
# include <iostream>
# include <vector>

namespace Eigen
{
  constexpr auto serialize(auto & archive, Eigen::Matrix3d& matrix)
  {
    return archive(std::span<double, 9>(matrix.data(), 9));
  }
  constexpr auto serialize(auto & archive, const Eigen::Matrix3d& matrix)
  {
    return archive(std::span<const double, 9>(matrix.data(), 9));
  }
}

int main()
{
  auto [data, in, out] = zpp::bits::data_in_out();
  std::vector<Eigen::Matrix3d> PrimativeCell = { Eigen::Matrix3d::Identity() };
  out(PrimativeCell).or_throw();
  std::vector<Eigen::Matrix3d> PrimativeCell2;
  in(PrimativeCell2).or_throw();
  std::cout << PrimativeCell2[0] << std::endl;
}