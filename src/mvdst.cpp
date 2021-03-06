#include "mvdst.hpp"

#include <algorithm>
#include <exception>
#include <iostream>
#include <set>

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
    int weight = manhattanDistance(getVM(i->m_source), getVM(i->m_target));
    (*matrix)[*i].weight = weight;  // HEHMDA do i need it anymore?
    getVS(i->m_source).weights.insert({weight, i->m_target});
    getVS(i->m_target).weights.insert({weight, i->m_source});
  }
  std::cout << "Vertices count: " << boost::num_vertices(*matrix)
            << ", edges count: " << boost::num_edges(*matrix)
            << ", target diameter (D): " << targetDiameter << std::endl;
  std::cout << "Initialization done" << std::endl;
}


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

int initVertNum = 0;

void Mvdst::approximate() {
  const int graphSz = boost::num_vertices(*matrix);
  //  int i = getStartVertex();
  for (int i = getStartVertex(); i < graphSz; ++i) {
    initVertNum = i;
    try {
      createStar(i);
      replaceEdges();
    } catch (...) {
      // Если что-то пошло не так, просто скипаем решение
      continue;
    }

    best = Solution(solutionGraph);
    auto [edges_begin, edges_end] = boost::edges(*solutionGraph);
    for (auto i = edges_begin; i != edges_end; ++i) {
      (*solutionGraph)[*i].weight =
          (*matrix)[boost::edge(i->m_source, i->m_target, *matrix).first]
              .weight;
      best.addEdge(i);
    }
    best.diameter = this->curDiameter;

    io::write(*this, "../solutions/" + std::to_string(graphSz) + "_" +
                         std::to_string(best.getWeight()) + "_" +
                         std::to_string(curDiameter) + "_" + std::to_string(i) +
                         ".txt");
    std::cout << "wheight: " << best.getWeight() << " size: " << best.size()
              << " diameter: " << best.diameter << std::endl;
  }
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
  auto [vBegin, vEnd] = boost::vertices(*solutionGraph);
  for (VertexIt i = vBegin; i != vEnd; ++i) {
    getVS(*i).inCurrentSolution = false;
  }
}

bool Mvdst::check() {
  if (hasCycle((*best.edges.begin())->m_source, -1))
    return false;
  return true;
}
void Mvdst::createStar(unsigned long centerInex) {
  std::cout << "Creating a star" << std::endl;
  std::set<unsigned long> initVertices;
  initVertices.insert(centerInex);
  while (initVertices.size() < targetDiameter - 1)
    initVertices.insert(rand() % initSize);
  auto [vertBegin, vertEnd] = boost::vertices(*solutionGraph);
  for (auto i = vertBegin; i != vertEnd; ++i) {
    boost::clear_vertex(*i, *solutionGraph);
    getVS(*i).inCurrentSolution = false;
    getVS(*i).maxDistance = 0;
  }

  auto prev = initVertices.begin();
  for (auto i = ++initVertices.begin(); i != initVertices.end(); ++i) {
    boost::add_edge(*prev, *i, *solutionGraph);
    prev = i;
  }

  for (auto i = vertBegin; i != vertEnd; ++i)
    if (!initVertices.contains(*i)) {
      auto initBegin = initVertices.begin();
      unsigned long minVert = *initBegin;
      int minWeight =
          (*matrix)[boost::edge(*initBegin, *i, *matrix).first].weight;
      for (auto j = ++(initVertices.begin()); j != initVertices.end(); ++j) {
        int curWeight = (*matrix)[boost::edge(*j, *i, *matrix).first].weight;
        if (curWeight < minWeight) {
          minWeight = curWeight;
          minVert = *j;
        }
      }
      boost::add_edge(minVert, *i, *solutionGraph);
    }
  auto [edgesBegin, edgesEnd] = boost::edges(*solutionGraph);
  for (auto i = edgesBegin; i != edgesEnd; ++i) {
    (*solutionGraph)[*i].weight =
        manhattanDistance(getVM(i->m_source), getVM(i->m_target));
  }

  for (auto i = vertBegin; i != vertEnd; ++i) {
    clearMarks();
    int d = setMaxDistances(*i);
    if (d > curDiameter)
      curDiameter = d;
  }

  std::cout << "Done with the star" << std::endl;
}

void Mvdst::replaceEdges() {
  auto [vertBegin, vertEnd] = boost::vertices(*solutionGraph);
  unsigned long start = 0;
  while (true) {
    float maxRatio = 0;
    unsigned long maxVert = *vertBegin;
    unsigned long oldEdgeVert = *vertBegin;
    unsigned long newEdgeVert = *vertBegin;
    unsigned long maxBranchDist = 0;
    for (auto i = vertBegin + start; i != vertEnd; ++i) {
      //      std::cout << *i << std::endl;
      auto [eBegin, eEnd] = boost::out_edges(*i, *solutionGraph);
      auto curE = eEnd;
      for (auto k = eBegin; k != eEnd; ++k) {
        if (getVS(k->m_target).maxDistance < getVS(*i).maxDistance) {
          curE = k;
          break;
        }
      }
      if (curE == eEnd)
        continue;
      //      std::cout << "i:" << *i << " curE" << *curE
      //                << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
      unsigned long targetV = curE->m_target;
      boost::remove_edge(curE->m_source, curE->m_target, *solutionGraph);
      clearMarks();
      int branchDist = setMaxDistancesImpl(*i, curE->m_target, 0, false);
      int curWeight = (*matrix)[*curE].weight;
      auto vBegin = getVS(*i).weights.begin();
      auto vEnd = getVS(*i).weights.end();
      for (auto j = vBegin; j != vEnd; ++j) {
        // HEHMDA add sorted list
        if (*i != j->second && j->second != targetV &&
            !boost::edge(*i, j->second, *solutionGraph).second) {
          if (getVS(j->second).maxDistance > targetDiameter) {
            std::cerr << "distance overflow\n";
            throw std::runtime_error("distance overflow");
          }
          // Считаем длину, если заменить ребро
          clearMarks();
          int mainDist = setMaxDistancesImpl(j->second, *i, 0, false);
          int newDist = branchDist + mainDist + 1;
          //          std::cout << "branchDist:" << branchDist << " mainDist:"
          //          << mainDist
          //                    << " newDist:" << newDist
          //                    << " targetDiameter:" << targetDiameter << '\n';
          if (newDist > targetDiameter)
            continue;

          float raito =
              static_cast<float>(curWeight) /  // HEHMDA *solutionGraph???
              (*matrix)[boost::edge(*i, j->second, *matrix).first].weight;
          //          std::cout << "curE: " << *curE << " curE weight: "
          //                    <<
          //                    static_cast<float>((*solutionGraph)[*curE].weight)
          //                    << " ij: " << boost::edge(*i, *j, *matrix).first
          //                    << " ij weight: "
          //                    << (*matrix)[boost::edge(*i, *j,
          //                    *matrix).first].weight
          //          std::cout << "raito: " << raito << " maxRaito: " <<
          //          maxRatio << '\n';
          //          if (raito > maxRatio) {
          if (raito > maxRatio &&
              //          if (raito > maxRatio &&
              !isReachable(*i, j->second,
                           std::numeric_limits<unsigned long>::max())) {
            //            std::cout << "+" << '\n';
            //            std::cout << "i:" << *i << " j:" << *j << "\n";
            maxRatio = raito;
            maxVert = *i;
            newEdgeVert = j->second;
            maxBranchDist = branchDist;
            oldEdgeVert = curE->m_target;
            break;
          }
        }
      }
      boost::add_edge(curE->m_source, curE->m_target, *solutionGraph);

      // not searching for max of all vertices
      //            start = *i;
      if (maxRatio > 1.0)
        break;
    }
    //    std::cout << "max ratio: " << maxRatio
    //              << " edges count: " << boost::num_edges(*solutionGraph)
    //              << std::endl;
    if (maxRatio <= 1.0) {
      return;
    }
    boost::remove_edge(maxVert, oldEdgeVert, *solutionGraph);
    boost::add_edge(maxVert, newEdgeVert, *solutionGraph);
    //    std::cout
    //        << "curE: " << boost::edge(maxVert, oldEdgeVert,
    //        *solutionGraph).first
    //        << " newE: " << boost::edge(maxVert, newEdgeVert,
    //        *solutionGraph).first
    //        << " " << maxVert << ":" << getVS(maxVert).maxDistance << " "
    //        << oldEdgeVert << ":" << getVS(oldEdgeVert).maxDistance << " "
    //        << maxVert << ":" << getVS(maxVert).maxDistance << " " <<
    //        newEdgeVert
    //        << ":" << getVS(newEdgeVert).maxDistance << " edges count: "
    //        << boost::num_edges(*solutionGraph)
    //              << " curE->m_source degree: " <<
    //              getVS(curE->m_source).degree
    //              << " curE->m_target degree: " <<
    //              getVS(curE->m_target).degree
    //              << " maxVert degree: " << getVS(maxVert).degree
    //              << " newEdgeVert degree: " << getVS(newEdgeVert).degree
    //        << std::endl;

    clearMarks();
    int mainDist =
        setMaxDistancesImpl(newEdgeVert, maxVert, maxBranchDist + 1, true);
    clearMarks();
    int branchDist = setMaxDistancesImpl(maxVert, newEdgeVert, mainDist, true);
    int maxDist = branchDist;  // учитывает вес остального дерева + ветви
    std::cout
        << "curE: " << boost::edge(maxVert, oldEdgeVert, *solutionGraph).first
        << " newE: " << boost::edge(maxVert, newEdgeVert, *solutionGraph).first
        << " " << maxVert << ":" << getVS(maxVert).maxDistance << " "
        << oldEdgeVert << ":" << getVS(oldEdgeVert).maxDistance << " "
        << " " << maxVert << ":" << getVS(maxVert).maxDistance << " "
        << newEdgeVert << ":" << getVS(newEdgeVert).maxDistance
        << " maxBranchDist: " << maxBranchDist << " mainDist: " << mainDist
        << " branchDist: " << branchDist
        << " edges count: " << boost::num_edges(*solutionGraph) << " maxRaito: "
        << maxRatio
        //              << " curE->m_source degree: " <<
        //              getVS(curE->m_source).degree
        //              << " curE->m_target degree: " <<
        //              getVS(curE->m_target).degree
        //              << " maxVert degree: " << getVS(maxVert).degree
        //              << " newEdgeVert degree: " << getVS(newEdgeVert).degree
        << " maxDist: " << maxDist << std::endl;
    if (maxDist > targetDiameter) {
      boost::add_edge(maxVert, oldEdgeVert, *solutionGraph);
      boost::remove_edge(maxVert, newEdgeVert, *solutionGraph);
      return;
    }
    if (maxDist > this->curDiameter)
      this->curDiameter = maxDist;


    static int t = 0;
    const int period = initSize == 512 ? 50 : initSize == 2048 ? 20 : 1;
    if (t % period == 0) {
      const int graphSz = boost::num_vertices(*matrix);
      best = Solution(solutionGraph);
      auto [edges_begin, edges_end] = boost::edges(*solutionGraph);
      for (auto i = edges_begin; i != edges_end; ++i) {
        (*solutionGraph)[*i].weight =
            (*matrix)[boost::edge(i->m_source, i->m_target, *matrix).first]
                .weight;
        best.addEdge(i);
      }
      best.diameter = curDiameter;
      io::write(*this, "../tmpPics/" + std::to_string(graphSz) + "_" +
                           std::to_string(best.getWeight()) + "_" +
                           std::to_string(curDiameter) + "_" +
                           std::to_string(initVertNum) + "_" +
                           std::to_string(t) + ".txt");
    }
    ++t;

  }
}

int Mvdst::setMaxDistances(unsigned long v, bool isWriting) {
  return setMaxDistancesImpl(v, std::numeric_limits<unsigned long>::max(), 0,
                             isWriting);
}

int Mvdst::setMaxDistancesImpl(unsigned long v,
                               unsigned long prev,
                               int dist,
                               bool isWriting) {
  auto& vertex = getVS(v);
  //  std::cout << "v:" << v << " prev:" << prev << " dist:" << dist
  //            << " isWriting:" << isWriting << " v.m:" << vertex.maxDistance
  //            << '\n';
  int max = dist;
  auto [vBegin, vEnd] = boost::adjacent_vertices(v, *solutionGraph);
  vertex.inCurrentSolution = true;
  for (auto i = vBegin; i != vEnd; ++i) {
    if (*i != prev) {
      if (getVS(*i).inCurrentSolution) {
        std::cerr << "cycled vertex i: " << *i << std::endl;
        throw std::runtime_error("cycle found!! ");
      }
      int cur = setMaxDistancesImpl(*i, v, dist + 1, isWriting);
      if (cur > max)
        max = cur;
    }
  }
  if (isWriting) {
    if (vertex.maxDistance < dist)
      vertex.maxDistance = dist;
    if (vertex.maxDistance < max - dist)
      vertex.maxDistance = max - dist;
  }
  //  std::cout << "v:" << v << " prev:" << prev << " dist:" << dist
  //            << " isWriting:" << isWriting << " v.m:" << vertex.maxDistance
  //            << " max:" << max << '\n';
  return max;
}

bool Mvdst::isReachable(unsigned long source,
                        unsigned long target,
                        unsigned long prev) {
  if (source == target)
    return true;
  auto [vBegin, vEnd] = boost::adjacent_vertices(source, *solutionGraph);
  for (auto i = vBegin; i != vEnd; ++i) {
    if (*i != prev && isReachable(*i, target, source))
      return true;
  }
  return false;
}

}  // namespace graph
