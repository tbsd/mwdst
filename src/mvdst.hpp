#pragma once

#include <filesystem>
#include <fstream>
#include <utility>

#include "solution.hpp"
#include "types.hpp"

namespace graph {
  // k-minimum spanning tree
  class Mvdst {
    using Vertex = boost::graph_traits<GraphMatrix>::vertex_descriptor;
    using Edge = GraphMatrix::edge_descriptor;
    using EdgeIt = GraphMatrix::out_edge_iterator;
    using VertexIt = GraphMatrix::vertex_iterator;

    public:
    Mvdst(std::unique_ptr<RawVertices> vertices);

    // Найти приближённое решение
    void approximate();

    const Solution<GraphMatrix> getSolution() const;

    int getProblemSize() const;

    // Соответствует ли решение условиям задачи
    bool check();

    private:
    //    void removeBigEdges(int startShift);

    /* Получить индекс вершины, оптимальной для начала поиска решения.
     * Вершины изначально выбраны руками на основе полученных картинок (потом
     * сделал перебор по всем вершинам).
     * Для общего случая имеет смысл искать
     * скопления точек и в них выбирать центральную или с наименьшим ребом
     */
    int getStartVertex();

    // Вершина по индексу
    VertexImpl& getVS(int index);

    // Итератор вершины по индексу
    VertexIt getVItS(int index);

    // Вершина по индексу
    VertexImpl& getVM(int index);

    // Итератор вершины по индексу
    VertexIt getVItM(int index);

    //    std::optional<EdgeIt> findMinEdge(Solution<GraphMatrix>& solution);

    //    void deleteHeavierEdges(int weight);

    bool hasCycle(unsigned long startV, int parentV);

    void clearMarks();

    void createStar(unsigned long centerInex);

    void replaceEdges();

    /*
     * Пересчитать максимальные расстояния относительно вершины v:
     * если расстояние от v до u больше u.maxDistance, то записываем его, иначе
     * ничего не делаем Возвращает максимальное расстояние от v до любой другой
     * вершины.
     */
    int setMaxDistances(unsigned long v, bool isWriting = true);

    int setMaxDistancesImpl(unsigned long v,
                            unsigned long prev,
                            int dist,
                            bool isWriting);

    bool isReachable(unsigned long source,
                     unsigned long target,
                     unsigned long prev);

    const int initSize;
    const int targetDiameter;
    int curDiameter = 0;
    std::shared_ptr<GraphMatrix> matrix;
    std::shared_ptr<GraphMatrix> solutionGraph;
    int bestWeight = std::numeric_limits<int>::max();
    int bestLeavesCount = std::numeric_limits<int>::max();
    Solution<GraphMatrix> best;
  };

  }  // namespace graph
