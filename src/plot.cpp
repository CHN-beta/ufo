# include <ufo/ufo.impl.hpp>

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

// BAND = 0.000000 0.000000 0.000000 0.500000 0.000000 0.000000 0.333333 0.333333 0.000000 0.000000 0.000000 0.000000 0 0 0.5 0.5 0 0.5 0.333333 0.333333 0.5 0 0 0.5
// BAND_LABELS = $\Gamma$ M K $\Gamma$ A L H A

// Gamma 0 0 0
// M 1/2 0 0
// K 1/3 1/3 0
// A 0 0 1/2
// L 1/2 0 1/2
// H 1/3 1/3 1/2

std::vector<std::vector<double>> Qpoints =
{
  {0, 0, 0},
  {0.025, 0, 0},
  {0.05, 0, 0},
  {0.075, 0, 0},
  {0.1, 0, 0},
  {0.125, 0, 0},
  {0.15, 0, 0},
  {0.175, 0, 0},
  {0.2, 0, 0},
  {0.225, 0, 0},
  {0.25, 0, 0},
  {0.275, 0, 0},
  {0.3, 0, 0},
  {0.325, 0, 0},
  {0.35, 0, 0},
  {0.375, 0, 0},
  {0.4, 0, 0},
  {0.425, 0, 0},
  {0.45, 0, 0},
  {0.475, 0, 0},
  {0.5, 0, 0}
};

int main(int argc, char** argv)
{
  YAML::Node root = YAML::LoadFile(argv[1]);
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

  // 外层表示 q 点坐标, 内层表示频率
  // 频率取 -5 到 30 THz, 每 0.1 THz 一个点
  std::vector<std::vector<double>> data(21);
  for (int i = 0; i < 21; i++)
  {
    data[i].resize(351);
    for (auto& qpoint : output.QPointData)
      if (std::abs(qpoint.QPoint[0] - Qpoints[i][0]) < 1e-3 &&
          std::abs(qpoint.QPoint[1] - Qpoints[i][1]) < 1e-3 &&
          std::abs(qpoint.QPoint[2] - Qpoints[i][2]) < 1e-3)
      {
        for (auto& mode : qpoint.ModeData)
          data[i][static_cast<int>((mode.Frequency + 5) * 10)] += mode.Weight;
        break;
      }
  }

  matplot::image(data, true);
  matplot::colorbar().limits({0, 0.5});
  matplot::show();
}
