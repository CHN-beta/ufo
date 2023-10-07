# include <ufo/plot.hpp>

namespace ufo
{
  PlotSolver::InputType::SourceType::SourceType(std::string filename)
  {
    auto input = std::ifstream(filename, std::ios::binary | std::ios::in);
    input.exceptions(std::ios::badbit | std::ios::failbit);
    static_assert(sizeof(std::byte) == sizeof(char));
    std::vector<std::byte> data;
    {
      std::vector<char> string(std::istreambuf_iterator<char>(input), {});
      data.assign
      (
        reinterpret_cast<std::byte*>(string.data()),
        reinterpret_cast<std::byte*>(string.data() + string.size())
      );
    }
    auto in = zpp::bits::in(data);
    UnfoldSolver::OutputType output;
    in(output).or_throw();
    static_cast<UnfoldSolver::OutputType&>(*this) = std::move(output);
  }

  PlotSolver::InputType::InputType(std::string config_file)
  {
    auto input = YAML::LoadFile(config_file);
    for (unsigned i = 0; i < 3; i++)
      for (unsigned j = 0; j < 3; j++)
        PrimativeCell(i, j) = input["PrimativeCell"][i][j].as<double>();
    for (auto& figure : input["Figures"].as<std::vector<YAML::Node>>())
    {
      Figures.emplace_back();
      auto qpoints = figure["Qpoints"]
        .as<std::vector<std::vector<std::vector<double>>>>();
      for (auto& line : qpoints)
      {
        Figures.back().Qpoints.emplace_back();
        for (auto& point : line)
          Figures.back().Qpoints.back().emplace_back(point.at(0), point.at(1), point.at(2));
        if (Figures.back().Qpoints.back().size() < 2)
          throw std::runtime_error("Not enough points in a line");
      }
      if (Figures.back().Qpoints.size() < 1)
        throw std::runtime_error("Not enough lines in a figure");
      Figures.back().Resolution = figure["Resolution"].as<std::pair<unsigned, unsigned>>();
      Figures.back().Range = figure["Range"].as<std::pair<double, double>>();
      Figures.back().Filename = figure["Filename"].as<std::string>();
      if (figure["YTicks"])
        Figures.back().YTicks = figure["YTicks"].as<std::vector<double>>();
    }
    SourceFilename = input["SourceFilename"].as<std::string>();
    Source = SourceType(SourceFilename);
  }

  PlotSolver::PlotSolver(std::string config_file) : Input_(config_file) {}

  PlotSolver& PlotSolver::operator()()
  {
    for (auto& figure : Input_.Figures)
    {
      // 外层表示不同的线段的端点，内层表示这个线段上的 q 点
      std::vector<std::vector<std::reference_wrapper<const UnfoldSolver::OutputType::QpointDataType>>> qpoints;
      std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> lines;
      for (auto& path : figure.Qpoints)
        for (unsigned i = 0; i < path.size() - 1; i++)
        {
          lines.emplace_back(path[i], path[i + 1]);
          qpoints.push_back(search_qpoints
          (
            lines.back(), Input_.Source.QpointData,
            0.001, i != path.size() - 2
          ));
        }
      auto [values, x_ticks] = calculate_values
        (lines, qpoints, figure.Resolution, figure.Range);
      auto y_ticks = figure.YTicks.value_or(std::vector<double>{});
      for (auto& _ : y_ticks)
        _ = _ / (figure.Range.second - figure.Range.first) * figure.Resolution.second;
      plot(values, figure.Filename, x_ticks, y_ticks);
    }
    return *this;
  }

  std::vector<std::reference_wrapper<const UnfoldSolver::OutputType::QpointDataType>> PlotSolver::search_qpoints
  (
    const std::pair<Eigen::Vector3d, Eigen::Vector3d>& path,
    const decltype(InputType::SourceType::QpointData)& available_qpoints,
    double threshold, bool exclude_endpoint
  )
  {
    std::multimap<double, std::reference_wrapper<const UnfoldSolver::OutputType::QpointDataType>> selected_qpoints;
    // 对于 output 中的每一个点, 检查这个点是否在路径上. 如果在, 把它加入到 selected_qpoints 中
    for (auto& qpoint : available_qpoints)
    {
      // 计算三点围成的三角形的面积的两倍
      auto area = (path.second - path.first).cross(qpoint.Qpoint - path.first).norm();
      // 计算这个点到前两个点所在直线的距离
      auto distance = area / (path.second - path.first).norm();
      // 如果这个点到前两个点所在直线的距离小于阈值, 则认为这个点在路径上
      if (distance < threshold)
      {
        // 计算这个点到前两个点的距离, 两个距离都应该小于两点之间的距离
        auto distance1 = (qpoint.Qpoint - path.first).norm();
        auto distance2 = (qpoint.Qpoint - path.second).norm();
        auto distance3 = (path.second - path.first).norm();
        if (distance1 < distance3 + threshold && distance2 < distance3 + threshold)
          // 如果这个点不在终点处, 或者不排除终点, 则加入
          if (distance2 > threshold || !exclude_endpoint)
            selected_qpoints.emplace(distance1, std::ref(qpoint));
      }
    }
    // 去除非常接近的点
    for (auto it = selected_qpoints.begin(); it != selected_qpoints.end();)
    {
      auto next = std::next(it);
      if (next == selected_qpoints.end())
        break;
      else if (next->first - it->first < threshold)
        selected_qpoints.erase(next);
      else
        it = next;
    }
    if (selected_qpoints.empty())
      throw std::runtime_error("No q points found");
    std::vector<std::reference_wrapper<const UnfoldSolver::OutputType::QpointDataType>> result;
    for (auto& qpoint : selected_qpoints)
      result.push_back(qpoint.second);
    return result;
  }

  std::tuple<std::vector<std::vector<double>>, std::vector<double>> PlotSolver::calculate_values
  (
    const std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& path,
    const std::vector<std::vector<std::reference_wrapper<const UnfoldSolver::OutputType::QpointDataType>>>& qpoints,
    const decltype(InputType::FigureConfigType::Resolution)& resolution,
    const decltype(InputType::FigureConfigType::Range)& range
  )
  {
    // 整理输入
    std::map<double, std::reference_wrapper<const UnfoldSolver::OutputType::QpointDataType>> qpoints_with_distance;
    double total_distance = 0;
    std::vector<double> x_ticks;
    for (unsigned i = 0; i < path.size(); i++)
    {
      for (auto& _ : qpoints[i])
        qpoints_with_distance.emplace(total_distance + (_.get().Qpoint - path[i].first).norm(), _);
      total_distance += (path[i].second - path[i].first).norm();
      if (i != path.size() - 1)
        x_ticks.push_back(total_distance);
    }
    for (auto& _ : x_ticks)
      _ = _ / total_distance * resolution.first;

    // 插值
    std::vector<std::vector<double>> values;
    auto blend = []
    (
      const UnfoldSolver::OutputType::QpointDataType& a,
      const UnfoldSolver::OutputType::QpointDataType& b,
      double ratio, unsigned resolution, std::pair<double, double> range
    ) -> std::vector<double>
    {
      // 计算插值结果
      std::vector<double> frequency, weight;
      for (unsigned i = 0; i < a.ModeData.size(); i++)
      {
        frequency.push_back(a.ModeData[i].Frequency * ratio + b.ModeData[i].Frequency * (1 - ratio));
        weight.push_back(a.ModeData[i].Weight * ratio + b.ModeData[i].Weight * (1 - ratio));
      }
      std::vector<double> result(resolution);
      for (unsigned i = 0; i < frequency.size(); i++)
      {
        int index = (frequency[i] - range.first) / (range.second - range.first) * resolution;
        if (index >= 0 && index < static_cast<int>(resolution))
          result[index] += weight[i];
      }
      return result;
    };
    for (unsigned i = 0; i < resolution.first; i++)
    {
      auto current_distance = total_distance * i / resolution.first;
      auto it = qpoints_with_distance.lower_bound(current_distance);
      if (it == qpoints_with_distance.begin())
        values.push_back(blend(it->second.get(), it->second.get(), 1, resolution.second, range));
      else if (it == qpoints_with_distance.end())
        values.push_back(blend(std::prev(it)->second.get(), std::prev(it)->second.get(), 1, resolution.second,
          range));
      else
        values.push_back(blend
        (
          std::prev(it)->second.get(), it->second.get(),
          (it->first - current_distance) / (it->first - std::prev(it)->first),
          resolution.second, range)
        );
    }
    return {values, x_ticks};
  }
  void PlotSolver::plot
  (
    const std::vector<std::vector<double>>& values,
    const decltype(InputType::FigureConfigType::Filename)& filename,
    const std::vector<double>& x_ticks, const std::vector<double>& y_ticks
  )
  {
    std::vector<std::vector<double>>
      r(values[0].size(), std::vector<double>(values.size(), 0)),
      g(values[0].size(), std::vector<double>(values.size(), 0)),
      b(values[0].size(), std::vector<double>(values.size(), 0)),
      a(values[0].size(), std::vector<double>(values.size(), 0));
    for (unsigned i = 0; i < values[0].size(); i++)
      for (unsigned j = 0; j < values.size(); j++)
      {
        r[i][j] = 255;
        g[i][j] = 0;
        b[i][j] = 0;
        a[i][j] = values[j][i] * 2 * 255;
        if (a[i][j] > 255)
          a[i][j] = 255;
      }
    auto f = matplot::figure(true);
    auto ax = f->current_axes();
    auto image = ax->image(std::tie(r, g, b));
    image->matrix_a(a);
    ax->y_axis().reverse(false);
    ax->x_axis().tick_values(x_ticks);
    ax->x_axis().tick_length(1);
    ax->y_axis().tick_values(y_ticks);
    ax->y_axis().tick_length(1);
    f->save(filename);
  }
}
