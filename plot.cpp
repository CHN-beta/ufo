# include <iostream>
# include <array>
# include <numbers>
# include <numeric>
# include <fstream>
# include <cassert>
# include <yaml-cpp/yaml.h>
# include <fmt/format.h>

using namespace std::literals;

// 一个临时的程序, 用于将数据导出画图

struct Output
{
  struct QPointDataType_
  {
    std::vector<double> QPoint;
    struct ModeDataType_
    {
      double Frequency;
      double Weight;
    };
    std::vector<ModeDataType_> ModeData;
  };
  std::vector<QPointDataType_> QPointData;
};

int main()
{
  YAML::Node root = YAML::LoadFile("out.yaml");
  Output output;
  output.QPointData.resize(root["QPointData"].size());
  for (int i = 0; i < root["QPointData"].size(); i++)
  {
    output.QPointData[i].QPoint = root["QPointData"][i]["QPoint"].as<std::vector<double>>();
    output.QPointData[i].ModeData.resize(root["QPointData"][i]["ModeData"].size());
    for (int j = 0; j < root["QPointData"][i]["ModeData"].size(); j++)
    {
      output.QPointData[i].ModeData[j].Frequency
        = root["QPointData"][i]["ModeData"][j]["Frequency"].as<double>();
      output.QPointData[i].ModeData[j].Weight
        = root["QPointData"][i]["ModeData"][j]["Weight"].as<double>();
    }
  }

  std::ofstream os("out.dat");
  // 13 个不同的 q 点
  // -5 \sim 30 THz, 分辨率为 0.01 THz
  // 先按照三角形进行涂抹, 虽然我觉得更合理的可能应该是正态分布
  // 当高度达到 1 时, 涂抹的宽度达到 0.4 THz

  // 构建未涂抹时的数据
  std::vector<std::vector<double>> original_data(13, std::vector<double>(3501));
  for (auto& data1: original_data)
    for (auto& data2: data1)
      assert(data2 == 0);

  for (auto qpoint: output.QPointData)
    for (auto mode: qpoint.ModeData)
      if (qpoint.QPoint[1] < 1e-3 && qpoint.QPoint[2] < 1e-3)
        original_data.at(static_cast<int>(qpoint.QPoint[0] * 12 + 0.5)).at(static_cast<int>(mode.Frequency / 0.01 + 500 + 0.5)) += mode.Weight;

  // 涂抹
  auto data = original_data;
  for (auto& data1: data)
    data1.assign(data1.size(), 0);

  for (int i = 0; i < 13; i++)
    for (int j = 0; j < 3501; j++)
    {
      data[i][j] += std::sqrt(original_data[i][j]);
      for (int k = 0; k < 3501; k++)
        if (std::sqrt(original_data[i][j]) - std::abs(k - j) * 0.05 > 0)
          data[i][k] += std::sqrt(original_data[i][j]) - std::abs(k - j) * 0.05;
    }

  // 输出
  for (int i = 0; i < 13; i++)
    os << fmt::format("\t{}", static_cast<double>(i) / 12);
  os << "\n";
  for (int i = 0; i < 3501; i++)
  {
    os << fmt::format("{}", static_cast<double>(i - 500) * 0.01);
    for (int j = 0; j < 13; j++)
      os << fmt::format("\t{}", data[j][i]);
    os << "\n";
  }
}
