#include <cassert>
#include <iostream>
#include <limits>

#include "../src/Graph.hpp"

void testGraph() {
    // Test initialization with larger graph
    Graph graph(25);

    // Add edges to create a more complex graph
    for (int i = 0; i < 24; ++i) {
        graph.addEdge(i, i + 1, 1.0 + (i % 3), true);
    }
    graph.addEdge(0, 5, 2.5, true);
    graph.addEdge(4, 10, 1.5, true);
    graph.addEdge(6, 15, 3.0, true);
    graph.addEdge(8, 20, 2.0, true);

    // Dijkstra basic test with larger graph
    auto result = graph.dijkstra(0, 20);
    assert(result.second > 0);
    assert(!result.first.empty());

    // Dijkstra unreachable test
    Graph isolatedGraph(30);
    isolatedGraph.addEdge(0, 1, 1.0, false);
    result = isolatedGraph.dijkstra(29, 0);
    assert(result.second == std::numeric_limits<double>::infinity());
    assert(result.first.empty());

    // Multi-source multi-target test
    std::vector<int> sources = {0, 5, 10};
    std::vector<int> targets = {15, 20, 24};
    auto multi_result = graph.multiSourceMultiTargetDijkstra(sources, targets, true);
    assert(multi_result.size() == 3);

    // Test longestShortestPath
    std::vector<Point> points = {Point(0, 0, 0), Point(1, 1, 10), Point(2, 2, 20)};
    double longest_shortest = graph.longestShortestPath(points);
    assert(longest_shortest > 0);

    // Test WSPD Check with well-separated points
    std::vector<Point> set1 = {Point(0, 0, 0), Point(0, 1, 5), Point(1, 0, 10)};
    std::vector<Point> set2 = {Point(2, 2, 15), Point(3, 3, 20), Point(4, 4, 24)};
    double d = -1.0;
    auto wspd_result = graph.wspdCheck(set1, set2, d, d, 1.0);
    assert(!wspd_result.empty());

    std::cout << "All tests passed!" << std::endl;
}

int main() {
    testGraph();
    return 0;
}