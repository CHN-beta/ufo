# pragma once
# include <iostream>
# include <array>
# include <numbers>
# include <numeric>
# include <fstream>
# include <optional>
# include <array>
# include <utility>
# include <execution>
# include <syncstream>
# include <any>
# include <map>
# include <vector>
# include <span>
# include <yaml-cpp/yaml.h>
# include <Eigen/Dense>
# include <concurrencpp/concurrencpp.h>
# include <fmt/format.h>
# include <fmt/std.h>
# include <fmt/ranges.h>
# include <highfive/H5File.hpp>
# include <zpp_bits.h>
# include <matplot/matplot.h>

// 在相位中, 约定为使用 $\exp (2 \pi i \vec{q} \cdot \vec{r})$ 来表示原子的运动状态
//  (而不是 $\exp (-2 \pi i \vec{q} \cdot \vec{r})$)
// 一些书定义的倒格矢中包含了 $2 \pi$ 的部分, 我们这里约定不包含这部分.
//  也就是说, 正格子与倒格子的转置相乘, 得到单位矩阵.

namespace Eigen
{
  constexpr inline auto serialize(auto & archive, Eigen::Matrix3d& matrix)
    { return archive(std::span(matrix.data(), matrix.size())); }
  constexpr inline auto serialize(auto & archive, const Eigen::Matrix3d& matrix)
    { return archive(std::span(matrix.data(), matrix.size())); }
  constexpr inline auto serialize(auto & archive, Eigen::Vector3d& vector)
    { return archive(std::span(vector.data(), vector.size())); }
  constexpr inline auto serialize(auto & archive, const Eigen::Vector3d& vector)
    { return archive(std::span(vector.data(), vector.size())); }
}

namespace ufo
{
  using namespace std::literals;
  struct PhonopyComplex { double r, i; };
  inline HighFive::CompoundType create_compound_complex()
    { return {{"r", HighFive::AtomicType<double>{}}, {"i", HighFive::AtomicType<double>{}}}; }
  class Solver
  {
    public:
      virtual Solver& operator()() = 0;
      ~Solver() = default;

      inline static concurrencpp::generator<std::pair<Eigen::Vector<unsigned, 3>, unsigned>>
        triplet_sequence(Eigen::Vector<unsigned, 3> range)
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
  };
}

HIGHFIVE_REGISTER_TYPE(ufo::PhonopyComplex, ufo::create_compound_complex)
