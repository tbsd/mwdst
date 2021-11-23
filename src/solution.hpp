#pragma once

#include <iostream>
#include <list>
#include <memory>
#include <unordered_set>

#include <boost/graph/filtered_graph.hpp>

#include "types.hpp"

namespace graph {

template <typename GraphMatrix>
class Solution {
  public:
  using Vertex = boost::graph_traits<GraphMatrix>::vertex_descriptor;
  using Edge = GraphMatrix::edge_descriptor;
  //  using EdgeIt = GraphMatrix::out_edge_iterator;
  using EdgeIt = GraphMatrix::edge_iterator;
  using VertexIt = GraphMatrix::vertex_iterator;

  Solution() = default;

  Solution(std::shared_ptr<GraphMatrix> matrix) : matrix(matrix) {}

  Solution(const Solution<GraphMatrix>& rhs)
      : edges(rhs.edges), matrix(rhs.matrix), weight(rhs.weight){};

  Solution(Solution<GraphMatrix>&& rhs)
      : edges(std::move(rhs.edges)), matrix(rhs.matrix), weight(rhs.weight){};

  Solution<GraphMatrix>& operator=(Solution<GraphMatrix>&& rhs) {
    edges = std::move(rhs.edges);
    matrix = rhs.matrix;
    weight = rhs.weight;
    return *this;
  }

  bool operator<(const Solution<GraphMatrix>& rhs) {
    return weight == rhs.weight ? getLeavesCount() < rhs.getLeavesCount()
                                : weight < rhs.weight;
  }

  void addEdge(EdgeIt e) {
    edges.push_back(e);
    getV(e->m_source).inCurrentSolution = true;
    getV(e->m_target).inCurrentSolution = true;
    weight += (*matrix)[*e].weight;
  }

  void removeLeaf(EdgeIt e, VertexIt v) {
    (*matrix)[*v].inCurrentSolution = false;
    edges.remove(e);
    weight -= (*matrix)[*e].weight;
  }

  bool isAddable(const EdgeIt& e) {
    return (!getV(e->m_source).inCurrentSolution ||
            !getV(e->m_target).inCurrentSolution) &&
           !(*matrix)[*e].isTried;
  }

  // Число рёбер в решении
  int size() const { return edges.size(); }

  int getWeight() const { return weight; }

  int getLeavesCount() const {
    std::unordered_set<int> vertices;
    for (auto& e : edges) {
      vertices.insert(e->m_source);
      vertices.insert(e->m_target);
    }
    int count = 0;
    for (auto& v : vertices)
      if (isLeaf(getVIt(v)))
        ++count;
    return count;
  }

  // Получить вершину-лист по ребру
  std::optional<VertexIt> getLeaf(EdgeIt e) {
    if (isLeaf(getVIt(e->m_source)))
      return getVIt(e->m_source);
    if (isLeaf(getVIt(e->m_target)))
      return getVIt(e->m_target);
    return std::nullopt;
  }

  auto findMaxLeaf() {
    // Надо бы учитывать длину ветвей

    std::optional<std::pair<VertexIt, EdgeIt>> maxLeaf;
    int maxW = -1;
    for (auto edge : edges) {
      auto curLeaf = getLeaf(edge);
      if (curLeaf) {
        int curW = (*matrix)[*edge].weight;
        if (curW > maxW) {
          maxW = curW;
          maxLeaf = {*curLeaf, edge};
        }
      }
    }
    return maxLeaf;
  }

  private:
  VertexIt getVIt(int index) const {
    auto [vert_begin, vert_end] = boost::vertices(*matrix);
    return vert_begin + index;
  }

  VertexImpl& getV(int index) const { return (*matrix)[*getVIt(index)]; }

  bool isLeaf(VertexIt v) const {
    int count = std::count_if(edges.begin(), edges.end(), [&](EdgeIt e) {
      return getVIt(e->m_source) == v || getVIt(e->m_target) == v;
    });
    return count == 1;
  }

  public:
  std::list<EdgeIt> edges;  // Лучше использовать std::unordered_set, но нужно
                            // написать нормальный hash

  private:
  std::shared_ptr<GraphMatrix> matrix = nullptr;
  int weight = 0;
};

}  // namespace graph
