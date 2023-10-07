# pragma once
# include <ufo/unfold.hpp>

namespace ufo
{
  class PlotSolver : public Solver
  {
    public:
      struct InputType
      {
        Eigen::Matrix3d PrimativeCell;

        struct FigureConfigType
        {
          std::vector<std::vector<Eigen::Vector3d>> Qpoints;
          std::pair<unsigned, unsigned> Resolution;
          std::pair<double, double> Range;
          std::string Filename;
        };
        std::vector<FigureConfigType> Figures;

        struct SourceType : public UnfoldSolver::OutputType
        {
          SourceType(std::string filename);
          SourceType() = default;
        };
        std::string SourceFilename;
        SourceType Source;

        InputType(std::string config_file);
      };
    protected:
      InputType Input_;
    public:
      PlotSolver(std::string config_file);
      PlotSolver& operator()() override;

      // 根据 q 点路径, 搜索要使用的 q 点
      static std::vector<std::reference_wrapper<const UnfoldSolver::OutputType::QpointDataType>> search_qpoints
      (
        const std::pair<Eigen::Vector3d, Eigen::Vector3d>& path,
        const decltype(InputType::SourceType::QpointData)& available_qpoints,
        double threshold, bool exclude_endpoint = false
      );
      // 根据搜索到的 q 点, 计算每个点的数值
      static std::vector<std::vector<double>> calculate_values
      (
        const std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& path,
        const std::vector<std::vector<std::reference_wrapper<const UnfoldSolver::OutputType::QpointDataType>>>& qpoints,
        const decltype(InputType::FigureConfigType::Resolution)& resolution,
        const decltype(InputType::FigureConfigType::Range)& range
      );
      // 根据数值, 画图
      static void plot
      (
        const std::vector<std::vector<double>>& values,
        const decltype(InputType::FigureConfigType::Filename)& filename
      );
  };
}
