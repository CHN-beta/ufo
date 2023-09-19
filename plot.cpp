# include <iostream>
# include <array>
# include <numbers>
# include <numeric>
# include <fstream>
# include <yaml-cpp/yaml.h>

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
  for (int i = 0; i < output.QPointData.size(); i++)
    for (int j = 0; j < output.QPointData[i].ModeData.size(); j++)
    {
      if (output.QPointData[i].QPoint[1] < 1e-3 && output.QPointData[i].QPoint[2] < 1e-3)
      {
        os << output.QPointData[i].QPoint[0] << "\t"
          << output.QPointData[i].ModeData[j].Frequency << "\t"
          << output.QPointData[i].ModeData[j].Weight << "\n";
      }
    }
}
