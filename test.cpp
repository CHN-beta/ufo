# include <eigen3/Eigen/Dense>
# include <zpp_bits.h>
# include <span>
# include <iostream>

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
  Eigen::Matrix3d PrimativeCell = Eigen::Matrix3d::Identity();
  out(PrimativeCell);
  Eigen::Matrix3d PrimativeCell2;
  in(PrimativeCell2);
  std::cout << PrimativeCell2 << std::endl;
}