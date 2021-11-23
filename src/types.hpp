#pragma once

#include <memory>
#include <vector>
#include <boost/graph/adjacency_matrix.hpp>
#include <boost/graph/graph_mutability_traits.hpp>
#include <boost/graph/graph_selectors.hpp>
#include <boost/graph/graph_traits.hpp>

#include "vertex_impl.hpp"
#include "edge_impl.hpp"

namespace graph {

  using RawVertices = std::vector<VertexImpl>;
  using RawEdge = std::pair<int, int>;
  using GraphMatrix = boost::adjacency_matrix<boost::undirectedS, VertexImpl, EdgeImpl>;

}
