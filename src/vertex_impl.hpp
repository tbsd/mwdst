#pragma once

#include <cmath>

namespace graph {

  struct VertexImpl {
    int x = 0;
    int y = 0;
    int maxDistance = 0;
    int degree = 0;
    float raito = 0;  // отношение длины текущего ребра к минимально возможному

    bool inCurrentSolution = false;

    VertexImpl() = default;
    VertexImpl(int x, int y) : x(x), y(y) {}
  };

int manhattanDistance(const VertexImpl& point1, const VertexImpl& point2);

}
