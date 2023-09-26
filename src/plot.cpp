# include <ufo/ufo.impl.hpp>

class FloatVector
{
  protected:
    std::vector<double> Data_;
    double LowerBound_, UpperBound_, Step_;
  public:
    FloatVector(double LowerBound, double UpperBound, double Step)
      : LowerBound_(LowerBound), UpperBound_(UpperBound), Step_(Step), Data_((UpperBound - LowerBound) / Step + 1) {}
    FloatVector(const FloatVector&) = default;
    FloatVector(FloatVector&&) = default;
    FloatVector& operator=(const FloatVector&) = default;
    FloatVector& operator=(FloatVector&&) = default;
    double& operator[](double i) { return Data_[static_cast<int>((i - LowerBound_) / Step_)]; }
    double operator[](double i) const { return Data_[static_cast<int>((i - LowerBound_) / Step_)]; }
    double lower_bound() const { return LowerBound_; }
    double upper_bound() const { return UpperBound_; }
    double step() const { return Step_; }
    int size() const { return Data_.size(); }
    std::map<double, double> to_map() const
    {
      std::map<double, double> result;
      for (int i = 0; i < Data_.size(); i++)
        result[LowerBound_ + i * Step_] = Data_[i];
      return result;
    }
};

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

int main(int argc, char** argv)
{
  if (argc != 2)
    throw std::runtime_error("Usage: plot output.zpp");

  Output output(argv[1]);

  //

  struct Point
  {
    Eigen::Vector3d QPoint;
    FloatVector Weight;
    double Distance;
  };
  std::vector<Point> Points;

  double current_distance = 0;
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
            auto& _ = point_of_this_path.emplace_back
              (qpoint.QPoint, FloatVector{-5, 35, 0.1}, distance1);
            for (auto& mode : qpoint.ModeData)
              _.Weight[mode.Frequency] += mode.Weight;
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
      Points.emplace_back(point.QPoint, point.Weight, point.Distance + current_distance);
    current_distance += (Qpoints[i + 1] - Qpoints[i]).norm();
  }
}
