#include <vertex_impl.hpp>

#include <cmath>
#include <iostream>

namespace graph {

  int manhattanDistance(const VertexImpl& point1, const VertexImpl& point2) {
    return std::abs(point2.x - point1.x) + std::abs(point2.y - point1.y);
  }

}

