#include "io.hpp"

#include <iostream>
#include <fstream>

#include <boost/algorithm/string.hpp>

#include "vertex_impl.hpp"

namespace graph::io {

std::unique_ptr<RawVertices> read(const std::filesystem::path& inPath) {
  auto vertices = std::make_unique<RawVertices>();
  std::ifstream inFile{inPath};
  if (inFile.bad()) {
    std::cerr << "Bad input file" << std::endl;
    return {};
  }
  std::string line;
  std::getline(inFile, line);
  int verticesCount = std::atoi(line.c_str() + 4);
  vertices->reserve(verticesCount);
  while (std::getline(inFile, line)) {
    std::vector<std::string> coordinates(2);
    boost::split(coordinates, line, boost::is_any_of("\t"),
                 boost::token_compress_on);
    vertices->emplace_back(std::atoi(coordinates[0].c_str()),
                           std::atoi(coordinates[1].c_str()));
  }
  inFile.close();
  return vertices;
}

void write(const Mvdst& problem, const std::filesystem::path& outPath) {
  auto solution = problem.getSolution();
  std::filesystem::remove(outPath);
  std::ofstream outFile{outPath};
  outFile << "c Вес дерева = " << solution.getWeight()
          << ", число листьев = " << solution.getLeavesCount() << "\n";
  outFile << "p edge " << problem.getProblemSize() << " " << solution.size()
          << "\n";
  for (const auto& edge : solution.edges) {
    outFile << "e " << edge->m_source + 1 << " " << edge->m_target + 1 << "\n";
  }
  outFile.close();
}

}  // namespace graph::io
