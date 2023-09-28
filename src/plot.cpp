# include <ufo/ufo.impl.hpp>

// 要被用来画图的路径
std::vector<Eigen::Vector3d> Qpoints =
{
  { 0, 0, 0 },
  { 0.5, 0, 0 },
  { 1. / 3, 1. / 3, 0 },
  { 0, 0, 0 },
  { 0, 0, 0.5 },
  { 0.5, 0, 0.5 },
  { 1. / 3, 1. / 3, 0.5 },
  { 0, 0, 0.5 }
};
double Threshold = 0.001;

struct Point
{
  Eigen::Vector3d QPoint;
  Eigen::VectorXd Frequency, Weight;
  double Distance;
};

int main(int argc, char** argv)
{
  if (argc != 2)
    throw std::runtime_error("Usage: plot output.zpp");

  Output output(argv[1]);

  std::vector<Point> Points;

  double total_distance = 0;
  // 对于每一条路径进行搜索
  for (unsigned i = 0; i < Qpoints.size() - 1; i++)
  {
    std::vector<Point> point_of_this_path;
    // 对于 output 中的每一个点, 检查这个点是否在路径上. 如果在, 把它加入到 point_of_this_path 中
    for (auto& qpoint : output.QPointData)
    {
      // 计算三点围成的三角形的面积的两倍
      auto area = (Qpoints[i + 1] - Qpoints[i]).cross(qpoint.QPoint - Qpoints[i]).norm();
      // 计算这个点到前两个点所在直线的距离
      auto distance = area / (Qpoints[i + 1] - Qpoints[i]).norm();
      // 如果这个点到前两个点所在直线的距离小于阈值, 则认为这个点在路径上
      if (distance < Threshold)
      {
        // 计算这个点到前两个点的距离, 两个距离都应该小于两点之间的距离
        auto distance1 = (qpoint.QPoint - Qpoints[i]).norm();
        auto distance2 = (qpoint.QPoint - Qpoints[i + 1]).norm();
        auto distance3 = (Qpoints[i + 1] - Qpoints[i]).norm();
        if (distance1 < distance3 + Threshold && distance2 < distance3 + Threshold)
          // 如果这个点在终点处, 且这条路径不是最后一条, 则不加入
          if (distance2 > Threshold || i == Qpoints.size() - 2)
          {
            auto& _ = point_of_this_path.emplace_back();
            _.QPoint = qpoint.QPoint;
            _.Distance = distance1;
            _.Frequency.resize(qpoint.ModeData.size());
            _.Weight.resize(qpoint.ModeData.size());
            for (unsigned j = 0; j < qpoint.ModeData.size(); j++)
            {
              _.Frequency(j) = qpoint.ModeData[j].Frequency;
              _.Weight(j) = qpoint.ModeData[j].Weight;
            }
          }
      }
    }
    // 对筛选结果排序
    std::sort(point_of_this_path.begin(), point_of_this_path.end(),
        [](const Point& a, const Point& b) { return a.Distance < b.Distance; });
    // 去除非常接近的点
    for (unsigned j = 1; j < point_of_this_path.size(); j++)
      while
      (
        j < point_of_this_path.size()
        && point_of_this_path[j].Distance - point_of_this_path[j - 1].Distance < Threshold
      )
        point_of_this_path.erase(point_of_this_path.begin() + j);
    // 将结果加入
    for (auto& point : point_of_this_path)
      Points.emplace_back(point.QPoint, point.Frequency, point.Weight, point.Distance + total_distance);
    total_distance += (Qpoints[i + 1] - Qpoints[i]).norm();
  }

  // 打印结果看一下
  for (auto& point : Points)
    std::cout << point.Distance << " " << point.QPoint.transpose() << std::endl;

  // 对结果插值
  std::vector<Point> interpolated_points;
  for (unsigned i = 0; i < 1024; i++)
  {
    auto current_distance = i * total_distance / 1024;
    auto& _ = interpolated_points.emplace_back();
    _.Distance = current_distance;
    auto it = std::lower_bound(Points.begin(), Points.end(), current_distance,
      [](const Point& a, double b) { return a.Distance < b; });
    // 如果是开头或者结尾, 直接赋值, 否则插值
    if (it == Points.begin())
    {
      _.Frequency = Points.front().Frequency;
      _.Weight = Points.front().Weight;
    }
    else if (it == Points.end() - 1)
    {
      _.Frequency = Points.back().Frequency;
      _.Weight = Points.back().Weight;
    }
    else
    {
      _.Frequency = (it->Frequency * (it->Distance - current_distance)
        + (it - 1)->Frequency * (current_distance - (it - 1)->Distance)) / (it->Distance - (it - 1)->Distance);
      _.Weight = (it->Weight * (it->Distance - current_distance)
        + (it - 1)->Weight * (current_distance - (it - 1)->Distance)) / (it->Distance - (it - 1)->Distance);
    }
  }

  // 将结果对应到像素上的值
  std::vector<std::vector<double>> weight(1024, std::vector<double>(400, 0));
  for (auto& point : interpolated_points)
  {
    int x = point.Distance / total_distance * 1024;
    if (x < 0)
      x = 0;
    else if (x >= 1024)
      x = 1023;
    for (unsigned i = 0; i < point.Frequency.size(); i++)
    {
      auto y = (point.Frequency(i) + 5) * 10;
      if (y < 0)
        y = 0;
      else if (y >= 400)
        y = 399;
      weight[x][y] += point.Weight(i);
    }
  }

  std::ofstream fout("weight.txt");
  for (unsigned i = 0; i < 400; i++)
    fout << fmt::format(" {:.6f}", i * 0.1 - 5);
  fout << std::endl;
  for (unsigned i = 0; i < 1024; i++)
  {
    fout << fmt::format("{:.6f} ", total_distance / 1024 * i);
    for (unsigned j = 0; j < 400; j++)
      fout << fmt::format("{:.6f} ", weight[i][j]);
    fout << std::endl;
  }

  std::vector<std::vector<double>>
    r(400, std::vector<double>(1024, 0)),
    g(400, std::vector<double>(1024, 0)),
    b(400, std::vector<double>(1024, 0));
  for (unsigned i = 0; i < 400; i++)
    for (unsigned j = 0; j < 1024; j++)
    {
      r[i][j] = 255;
      g[i][j] = 255 - weight[j][i] * 2 * 255;
      if (g[i][j] < 0)
        g[i][j] = 0;
      b[i][j] = 255 - weight[j][i] * 2 * 255;
      if (b[i][j] < 0)
        b[i][j] = 0;
    }
  auto f = matplot::figure(true);
  auto ax = f->current_axes();
  ax->image(std::tie(r, g, b));
  ax->y_axis().reverse(false);
  f->show();
}
