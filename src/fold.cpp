# include <ufo.hpp>

void ufo::fold(std::string config_file)
{
  struct Config
  {
    Eigen::Matrix3d SuperCellDeformation;
    Eigen::Vector3i SuperCellMultiplier;
    std::vector<Eigen::Vector3d> Qpoints;
  };

  auto fold = []
  (
    Eigen::Matrix3d super_cell_deformation,
    Eigen::Vector3i super_cell_multiplier,
    Eigen::Vector3d qpoint_in_reciprocal_primitive_cell_by_reciprocal_primitive_cell
  ) -> Eigen::Vector3d
  {
    /*
      首先需要将 q 点坐标的单位转换为 ModifiedSuperCell 的格矢，可知：
      QpointByReciprocalModifiedSuperCell =
        SuperCellMultiplier.cast<double>().asDiagonal() * QpointByReciprocalPrimitiveCell;
      接下来考虑将 q 点坐标的单位转换为 SuperCell 的格矢
      ModifiedSuperCell = SuperCellMultiplier.transpose() * PrimativeCell;
      SuperCell = SuperCellDeformation * ModifiedSuperCell;
      ReciprocalModifiedSuperCell = ModifiedSuperCell.inverse().transpose();
      ReciprocalSuperCell = SuperCell.inverse().transpose();
      Qpoint = QpointByReciprocalModifiedSuperCell.transpose() * ReciprocalModifiedSuperCell;
      Qpoint = QpointByReciprocalSuperCell.transpose() * ReciprocalSuperCell;
      整理可以得到:
      QpointByReciprocalSuperCell = SuperCellDeformation * QpointByReciprocalModifiedSuperCell;
      两个式子结合，可以得到：
      QpointByReciprocalSuperCell =
        SuperCellDeformation * SuperCellMultiplier.cast<double>().asDiagonal() * QpointByReciprocalPrimitiveCell;
    */
    auto qpoint_in_reciprocal_primitive_cell_by_reciprocal_super_cell =
    (
      super_cell_deformation * super_cell_multiplier.cast<double>().asDiagonal()
        * qpoint_in_reciprocal_primitive_cell_by_reciprocal_primitive_cell
    ).eval();
    /*
      到目前为止，我们还没有移动过 q 点的坐标。现在，我们将它移动整数个 ReciprocalSuperCell，直到它落在超胞的倒格子中。
      这等价于直接取 QpointByReciprocalSuperCell - QpointByReciprocalSuperCell.floor()。
    */
    return (qpoint_in_reciprocal_primitive_cell_by_reciprocal_super_cell.array()
      - qpoint_in_reciprocal_primitive_cell_by_reciprocal_super_cell.array().floor()).matrix();
  };

  biu::Logger::Guard log(config_file);
  auto input = YAML::LoadFile(config_file).as<Config>();
  for (const auto& qpoint : input.Qpoints) log.info("{} -> {}"_f
    (
      qpoint,
      fold(input.SuperCellDeformation, input.SuperCellMultiplier, qpoint)
    ));
}
