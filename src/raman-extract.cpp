# include <ufo.hpp>

void ufo::raman_extract(std::vector<std::string> files)
{
  biu::Logger::Guard log(files);
  std::vector<Eigen::Matrix3d> electricity_tensors;

  for (const auto& file : files)
  {
    auto in_stream = std::ifstream(file);
    std::string line;
    // search line containing: MACROSCOPIC STATIC DIELECTRIC TENSOR
    while (std::getline(in_stream, line))
    {
      if (line.find("MACROSCOPIC STATIC DIELECTRIC TENSOR") != std::string::npos)
      {
        log.debug("find in {}"_f(file));
        // skip 1 line
        std::getline(in_stream, line);
        electricity_tensors.emplace_back();
        for (std::size_t i = 0; i < 3; i++) for (std::size_t j = 0; j < 3; j++)
          in_stream >> electricity_tensors.back()(i, j);
        log.debug("read tensor {}"_f(electricity_tensors.back()));
        break;
      }
    }
    // test if the file is read correctly
    if (!in_stream) throw std::runtime_error("Error reading file: {}"_f(file));
  }

  // output the result
  for (auto e: electricity_tensors) std::cout <<
R"(- - [ {}, {}, {} ]
  - [ {}, {}, {} ]
  - [ {}, {}, {} ]
)"_f
    (
      e(0, 0), e(0, 1), e(0, 2),
      e(1, 0), e(1, 1), e(1, 2),
      e(2, 0), e(2, 1), e(2, 2)
    );
}
