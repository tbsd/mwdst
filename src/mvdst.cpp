#include "mvdst.hpp"

#include <iostream>
#include <set>
#include <algorithm>

#include <boost/graph/undirected_dfs.hpp>

#include "io.hpp"

namespace graph {

Mvdst::Mvdst(std::unique_ptr<RawVertices> vertices)
    : initSize(vertices->size()), targetDiameter(vertices->size() / 32 + 2) {
  matrix = std::make_shared<GraphMatrix>(vertices->size());
  solutionGraph = std::make_shared<GraphMatrix>(vertices->size());
  best = Solution(solutionGraph);
  auto j = vertices->begin();
  auto [vert_begin, vert_end] = boost::vertices(*matrix);
  for (VertexIt i = vert_begin; i != vert_end; ++i, ++j) {
    getVM(*i) = *j;
  }
  j = vertices->begin();
  auto [vert_begin_s, vert_end_s] = boost::vertices(*solutionGraph);
  for (VertexIt i = vert_begin_s; i != vert_end_s; ++i, ++j) {
    getVS(*i) = *j;
  }
  for (VertexIt i = vert_begin; i != vert_end; ++i) {
    for (VertexIt n = i + 1; n != vert_end; ++n) {
      boost::add_edge(*i, *n, *matrix);
    }
  }
  auto [edges_begin, edges_end] = boost::edges(*matrix);
  for (auto i = edges_begin; i != edges_end; ++i) {
    (*matrix)[*i].weight =
        manhattanDistance(getVM(i->m_source), getVM(i->m_target));
  }
  std::cout << "Vertices count: " << boost::num_vertices(*matrix)
            << ", edges count: " << boost::num_edges(*matrix)
            << ", target diameter (D): " << targetDiameter << std::endl;
  std::cout << "Initialization done" << std::endl;
}

/*
std::optional<Mvdst::EdgeIt> Mvdst::findMinEdge(
    Solution<GraphMatrix>& solution) {
  // Надо бы отдавать приоритет самому лёгкому ребру, образующему длиннейшую
  // цепочку, распологающемуся ближе к изначальной точке (раз уж такое
  // распределение). Может потом.

  using EdgeItRange = std::pair<EdgeIt, EdgeIt>;
  auto less = [](const EdgeItRange& lhs, const EdgeItRange& rhs) {
    return lhs.first->m_source < rhs.first->m_source;
  };
  std::set<EdgeItRange, decltype(less)> newEdges;
  for (auto e : solution.edges) {
    newEdges.insert(boost::out_edges(e->m_source, *matrix));
    newEdges.insert(boost::out_edges(e->m_target, *matrix));
  }
  std::optional<EdgeIt> minEdge;
  int minEdgeWeight = std::numeric_limits<int>::max();
  for (auto [begin, end] : newEdges) {
    for (auto i = begin; i != end; ++i) {
      int iWeight = (*matrix)[*i].weight;
      if (iWeight < minEdgeWeight && solution.isAddable(i)) {
        minEdge = i;
        minEdgeWeight = iWeight;
      }
    }
  }
  return minEdge;
}
*/

Mvdst::VertexIt Mvdst::getVItS(int index) {
  auto [vert_begin, vert_end] = boost::vertices(*solutionGraph);
  return vert_begin + index;
}

VertexImpl& Mvdst::getVS(int index) {
  return (*solutionGraph)[*getVItS(index)];
}

Mvdst::VertexIt Mvdst::getVItM(int index) {
  auto [vert_begin, vert_end] = boost::vertices(*matrix);
  return vert_begin + index;
}

VertexImpl& Mvdst::getVM(int index) {
  return (*matrix)[*getVItM(index)];
}

/*
void Mvdst::deleteHeavierEdges(int weight) {
  int initEdgesCount = boost::num_edges(*matrix);
  auto [all_edges_begin, all_edges_end] = boost::edges(*matrix);
  for (auto i = all_edges_begin; i != all_edges_end; ++i) {
    if ((*matrix)[*i].weight > weight - targetDiameter + 1) {
      boost::remove_edge(*i, *matrix);
    }
  }
  int deletedEdgesCount = initEdgesCount - boost::num_edges(*matrix);
  std::cout << "Deleted edges: " << deletedEdgesCount << "/" << initEdgesCount
            << " ("
            << static_cast<double>(deletedEdgesCount) / initEdgesCount * 100
            << "%)" << std::endl;
}
*/

/*
void Mvdst::removeBigEdges(int startShift) {
  std::cout << "Removing big edges" << std::endl;
  auto [vert_begin, vert_end] = boost::vertices(*matrix);
  auto initV = vert_begin + startShift;
  std::map<int, VertexIt>
      distances;  // расстояния от изначальной вершины до каждой другой.
  auto [initEdgesBegin, initEdgesEnd] = boost::out_edges(*initV, *matrix);
  for (auto i = initEdgesBegin; i != initEdgesEnd; ++i)
    distances.emplace((*matrix)[*i].weight, i->m_target);
  auto cur_v = initV;
  std::cout << "Start vertex: " << startShift << " (" << getVM(*initV).x << ", "
            << getVM(*initV).y << ")" << std::endl;
  Solution<GraphMatrix> curSolution(matrix);
  auto [edges_i, edges_end] = boost::out_edges(*cur_v, *matrix);

  auto min_e =
      std::min_element(edges_i, edges_end, [this](auto lhs, auto rhs) -> bool {
        return (*matrix)[lhs].weight < (*matrix)[rhs].weight;
      });
  curSolution.addEdge(min_e);
  while (curSolution.size() < targetDiameter) {
    auto newEdge = findMinEdge(curSolution);
    if (newEdge) {
      curSolution.addEdge(*newEdge);
      if (curSolution.size() % 10 == 0)
        std::cout << curSolution.size() << " " << std::flush;
    } else {
      std::cout << "\nCurrent solution is incomplete" << std::endl;
      return;
    }
  }
  std::cout << "\nCurrent solution complete. Total weight: "
            << curSolution.getWeight()
            << " Leaves count: " << curSolution.getLeavesCount() << std::endl;
  deleteHeavierEdges(curSolution.getWeight());

  std::cout << "Minimizeing solution" << std::endl;

  // находим минимальный лист и смотрим, можем ли добавить меньший где-то ещё
  while (true) {
    auto maxLeaf = curSolution.findMaxLeaf();
    if (maxLeaf) {
      auto [maxV, maxE] = *maxLeaf;
      curSolution.removeLeaf(maxE, maxV);
      auto minE = findMinEdge(curSolution);
      if (!minE) {
        std::cerr << "Something gone wrong. No unused edges found" << std::endl;
        return;
      }
      if ((*matrix)[*maxE].weight == (*matrix)[**minE].weight) {
        (*matrix)[*maxE].isTried = true;
        (*matrix)[**minE].isTried = true;
        curSolution.addEdge(*minE);
      } else if ((*matrix)[*maxE].weight < (*matrix)[**minE].weight) {
        curSolution.addEdge(maxE);
        break;
      } else {
        curSolution.addEdge(*minE);
        auto [edges_begin, edges_end] = boost::edges(*matrix);
        for (auto i = edges_begin; i != edges_end; ++i) {
          (*matrix)[*i].isTried = false;
        }
      }
    } else {
      break;
    }
  }
  std::cout << "Minimization complete. Total weight: "
            << curSolution.getWeight()
            << " Leaves count: " << curSolution.getLeavesCount() << std::endl;
  deleteHeavierEdges(curSolution.getWeight());

  if (curSolution < best || best.size() != targetDiameter)
    best = std::move(curSolution);
  clearMarks();
}
*/

int Mvdst::getStartVertex() {
  switch (boost::num_vertices(*matrix)) {
    case 64:
      return 1;
    case 128:
      return 42;
    case 512:
      return 0;
    case 2048:
      return 1080;
    case 4096:
      return 3015;
    default:
      return 0;
  }
}

void Mvdst::approximate() {
  // Раз уж работает так быстро, то в качестве начальных можно просто
  // перебрать все вершины из скопления
  // Или вообще все вершины. Для 2048 считалось 20 минут, для 4096 — 8 часов 47
  // минут. Но это вызывает какие-то непонятные краши на средних графах,
  // связаные с выделением/освобождением памяти. Подозреваю, из-за фрагментации
  // памяти. Как быстро побороть не знаю, пришлось перебирать кусками. А ещё
  // иногда запись в файл фейлится, возможно из-за этого же
  //  const int graphSz = boost::num_vertices(*matrix);
  //  for (int i = 0; i < graphSz; ++i) {
  //    removeBigEdges(i);
  //    std::cout << "Start points checked: " << i + 1 << "/" << graphSz
  //              << std::endl;
  //  }

  const int graphSz = boost::num_vertices(*matrix);
  int i = getStartVertex();
  //  for (int i = 0; i < graphSz; ++i) {
  createStar(i);
  replaceEdges();
  best = Solution(solutionGraph);
  auto [edges_begin, edges_end] = boost::edges(*solutionGraph);
  for (auto i = edges_begin; i != edges_end; ++i) {
    (*solutionGraph)[*i].weight =
        (*matrix)[boost::edge(i->m_source, i->m_target, *matrix).first].weight;
    best.addEdge(i);
  }
  best.diameter = this->curDiameter;
  std::cout << best.diameter;

  io::write(*this, "../solutions/" + std::to_string(graphSz) + "_" +
                       std::to_string(best.getWeight()) + "_" +
                       std::to_string(i) + ".txt");
  std::cout << "wheight: " << best.getWeight() << " size: " << best.size()
            << std::endl;
  //  }
}

const Solution<GraphMatrix> Mvdst::getSolution() const {
  return best;
}

int Mvdst::getProblemSize() const {
  return initSize;
}

bool Mvdst::hasCycle(unsigned long startV, int parentV) {
  getVM(startV).inCurrentSolution = true;
  for (auto& i : best.edges)
    if (i->m_source == startV || i->m_target == startV) {
      int newV = i->m_target == startV ? i->m_source : i->m_target;
      if (!getVM(newV).inCurrentSolution) {
        if (hasCycle(newV, startV))
          return true;
      } else if (newV != parentV)
        return true;
    }
  return false;
}

void Mvdst::clearMarks() {
  auto [vBegin, vEnd] = boost::vertices(*matrix);
  for (VertexIt i = vBegin; i != vEnd; ++i) {
    getVM(*i).inCurrentSolution = false;
  }
  auto [edges_begin, edges_end] = boost::edges(*matrix);
  for (auto i = edges_begin; i != edges_end; ++i) {
    (*matrix)[*i].isTried = false;
  }
}

bool Mvdst::check() {
  //  if (best.edges.size() !=
  //      static_cast<decltype(best.edges.size())>(targetDiameter))
  //    return false;
  if (hasCycle((*best.edges.begin())->m_source, -1))
    return false;
  //  for (auto& e : best.edges) {
  //    if (!getVM(e->m_source).inCurrentSolution ||
  //        !getVM(e->m_target).inCurrentSolution)
  //      return false;
  //  }
  //  clearMarks();
  return true;
}
void Mvdst::createStar(unsigned long centerInex) {
  auto [vertBegin, vertEnd] = boost::vertices(*solutionGraph);
  for (auto i = vertBegin; i != vertEnd; ++i)
    boost::clear_vertex(*i, *solutionGraph);
  for (auto i = vertBegin; i != vertEnd; ++i)
    if (*i != centerInex) {
      boost::add_edge(centerInex, *i, *solutionGraph);
      getVS(centerInex).degree++;
      getVS(*i).degree++;
      getVS(*i).maxDistance = 2;
    }
  auto [edgesBegin, edgesEnd] = boost::edges(*solutionGraph);
  for (auto i = edgesBegin; i != edgesEnd; ++i) {
    // получить вес из полного графа
    (*solutionGraph)[*i].weight =
        //        (*matrix)[boost::edge(i->m_source, i->m_target,
        //        *matrix).first].weight;
        manhattanDistance(getVM(i->m_source), getVM(i->m_target));
  }
  getVS(centerInex).maxDistance = 1;
  curDiameter = 2;
}

void Mvdst::replaceEdges() {
  auto [vertBegin, vertEnd] = boost::vertices(*solutionGraph);
  //  std::set<unsigned long> notToAdd;
  //  std::set<unsigned long> prevNotToAdd;
  unsigned long start = 0;
  while (true) {
    float maxRatio = 0;
    unsigned long maxVert = *vertBegin;
    unsigned long newEdgeVert = *vertBegin;
    for (auto i = vertBegin + start; i != vertEnd; ++i) {
      //      std::cout << *i << std::endl;
      auto vert = getVS(*i);
      // replace with bfs and sum check HEHMDA
      // чекать, что к этому вершине можно добавить
      if (vert.degree == 1 && vert.maxDistance <= targetDiameter) {
        for (auto j = vertBegin; j != vertEnd; ++j) {
          //          std::cout << "i: " << *i << " j: " << *j << std::endl;
          if (*i != *j) {
            if (getVS(*j).maxDistance > targetDiameter)
              std::cout << "HEHMDA wtf!!!!!!!!!!!!!!!!!!!" << std::endl;
            if (getVS(*j).maxDistance == targetDiameter)
              continue;
            auto [curE, dummy] = boost::out_edges(*i, *solutionGraph);
            float raito = static_cast<float>((*solutionGraph)[*curE].weight) /
                          (*matrix)[boost::edge(*i, *j, *matrix).first].weight;
            //                          manhattanDistance(getVS(*i), getVS(*j));
            //            std::cout << "curE: " << *curE << " curE weight: "
            //                      <<
            //                      static_cast<float>((*solutionGraph)[*curE].weight)
            //                      << " ij: " << boost::edge(*i, *j,
            //                      *matrix).first
            //                      << " ij weight: "
            //                      << (*matrix)[boost::edge(*i, *j,
            //                      *matrix).first].weight
            //                      << " raito: " << raito << " maxVert: " <<
            //                      maxVert
            //                      << " maxRaito: " << maxRatio << '\n';
            if (raito > maxRatio) {
              //              std::cout << "+" << '\n';
              maxRatio = raito;
              maxVert = *i;
              newEdgeVert = *j;
            }
          }
        }
        /*
        auto [curE, dummy] = boost::out_edges(*i, *solutionGraph);
        auto [outEBegin, outEEnd] = boost::out_edges(*i, *matrix);
        unsigned long maxEIndex = 0;
        for (auto j = outEBegin; j != outEEnd; ++j) {
          ++maxEIndex;
          if (*curE != *j) {
            auto otherVert = j->m_source == *i ? j->m_target : j->m_source;
            if (getVS(otherVert).maxDistance > targetDiameter)
              std::cout << "HEHMDA wtf!!!!!!!!!!!!!!!!!!!" << std::endl;
            if (getVS(otherVert).maxDistance == targetDiameter)
              continue;
            float raito = static_cast<float>((*solutionGraph)[*curE].weight) /
                          (*matrix)[*j].weight;
            if (raito > maxRatio) {
              maxRatio = raito;
              maxVert = *i;
              maxVertNewEIndex = maxEIndex;
            }
          }
        }
*/
      }
      // HEHMDA not searching for max of all vertices
      start = *i;
      if (maxRatio > 1.0)
        break;
    }
    std::cout << "max ratio: " << maxRatio
              << " edges count: " << boost::num_edges(*solutionGraph)
              << std::endl;
    if (maxRatio <= 1.0) {
      //      if (notToAdd == prevNotToAdd)
      return;
      //      prevNotToAdd = notToAdd;
      //      notToAdd.clear();
    }
    auto [curE, dummy] = boost::out_edges(maxVert, *solutionGraph);
    //        auto [outEBegin, outEEnd] = boost::out_edges(maxVert, *matrix);
    //        auto newE = boost::edge(curVert, maxVert, *solutionGraph);
    //        std::advance(newE, maxVertNewEIndex);
    //        auto newEdgeVert =
    //            curE->m_source == maxVert ? curE->m_target : curE->m_source;
    int oldWeight = (*solutionGraph)[*curE].weight;
    boost::remove_edge(curE->m_source, curE->m_target, *solutionGraph);
    boost::add_edge(maxVert, newEdgeVert, *solutionGraph);
    (*solutionGraph)[boost::edge(maxVert, newEdgeVert, *solutionGraph).first]
        .weight = manhattanDistance(getVM(maxVert), getVM(newEdgeVert));
    ++(getVS(newEdgeVert).degree);
    int maxDist = setMaxDistances(maxVert);
    std::cout << "curE: " << *curE << " newE: "
              << boost::edge(maxVert, newEdgeVert, *solutionGraph).first
              << " curE->m_source: " << curE->m_source
              << " curE->m_target: " << curE->m_target
              << " maxVert: " << maxVert << " newEdgeVert: " << newEdgeVert
              << " edges count: " << boost::num_edges(*solutionGraph)
              << " curE->m_source degree: " << getVS(curE->m_source).degree
              << " curE->m_target degree: " << getVS(curE->m_target).degree
              << " maxVert degree: " << getVS(maxVert).degree
              << " newEdgeVert degree: " << getVS(newEdgeVert).degree
              << " maxDist: " << maxDist << std::endl;
    if (maxDist > targetDiameter) {
      boost::add_edge(curE->m_source, curE->m_target, *solutionGraph);
      boost::remove_edge(maxVert, newEdgeVert, *solutionGraph);
      (*solutionGraph)
          [boost::edge(curE->m_source, curE->m_target, *solutionGraph).first]
              .weight = oldWeight;
      --(getVS(newEdgeVert).degree);

      //      notToAdd.insert(maxVert);

      return;
      //      // HEHMDA
      //      static int i = 0;
      //      i++;
      //      if (i > 5)
      //        return;

      //      std::cout << "HEHMDA never gonna get here !!!!!!!" << std::endl;
      //      boost::remove_edge(newE->m_source, newE->m_target,
      //      *solutionGraph); boost::add_edge(curE->m_source, curE->m_target,
      //      *solutionGraph); std::cout << "curE->m_source: "
      //                << curE->m_source
      //                << " curE->m_target: " << curE->m_target
      //                << " newE->m_source: " << newE->m_source
      //                << " newE->m_target: " << newE->m_target
      //                << " edges count: " <<
      //                boost::num_edges(*solutionGraph)
      //                << std::endl;
      //      std::cout << "curE->m_source: " << curE->m_source
      //                << " curE->m_target: " << curE->m_target
      //                << " newE->m_source: " << newE->m_source
      //                << " newE->m_target: " << newE->m_target
      //                << " edges count: " <<
      //                boost::num_edges(*solutionGraph)
      //                << std::endl;
      //      return;
    }
    if (maxDist > this->curDiameter)
      this->curDiameter = maxDist;

    // HEHMDA
    //    const int graphSz = boost::num_vertices(*matrix);
    //    best = Solution(solutionGraph);
    //    auto [edges_begin, edges_end] = boost::edges(*solutionGraph);
    //    for (auto i = edges_begin; i != edges_end; ++i) {
    //      (*solutionGraph)[*i].weight =
    //          (*matrix)[boost::edge(i->m_source, i->m_target, *matrix).first]
    //              .weight;
    //      best.addEdge(i);
    //    }
    //    best.diameter = curDiameter;
    //    static int t = 0;
    //    io::write(*this, "../tmpPics/" + std::to_string(graphSz) + "_" +
    //                         std::to_string(t) + ".txt");
    //    ++t;

    /*
    // добавляем к ней впервые
    if (++(getVS(maxVert).degree) == 2) {
      ++curDiameter;
      for (auto i = vertBegin; i != vertEnd; ++i) {
        if (*i != maxVert) {
          ++(getVS(*i).maxDistance);
        }
      }
    }
*/
  }
}

int Mvdst::setMaxDistances(unsigned long v) {
  return setMaxDistancesImpl(v, std::numeric_limits<unsigned long>::max(), 0);
}

int Mvdst::setMaxDistancesImpl(unsigned long v, unsigned long prev, int dist) {
  auto vertex = getVS(v);
  int max = dist;
  auto [vBegin, vEnd] = boost::adjacent_vertices(v, *solutionGraph);
  for (auto i = vBegin; i != vEnd; ++i) {
    if (*i != prev) {
      int cur = setMaxDistancesImpl(*i, v, dist + 1);
      if (cur > max)
        max = cur;
    }
  }
  if (vertex.maxDistance < dist)
    vertex.maxDistance = dist;
  if (vertex.maxDistance < max - dist)
    vertex.maxDistance = max - dist;
  return max;
}

}  // namespace graph
