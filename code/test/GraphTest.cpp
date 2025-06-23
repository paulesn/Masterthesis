#include <cassert>
#include <iostream>
#include <limits>

#include "../src/Graph.hpp"
#include "../src/Dataloader.hpp"

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
    std::vector<Point> points = {Point(0, 0, 0, 0), Point(1, 1, 10, 0), Point(2, 2, 20, 0)};
    //double longest_shortest = graph.longestShortestPath(points);
    //assert(longest_shortest > 0);

    // Test WSPD Check with well-separated points
    std::vector<Point> set1 = {Point(0, 0, 0, 0), Point(0, 1, 5, 0), Point(1, 0, 10, 0)};
    std::vector<Point> set2 = {Point(2, 2, 15, 0), Point(3, 3, 20, 0), Point(4, 4, 24, 0)};
    double d = -1.0;
    //auto wspd_result = graph.wspdCheck(set1, set2, d, d, 1.0);
    //assert(!std::get<0>(wspd_result).empty());

    std::cout << "All tests passed!" << std::endl;
}

void testHub_labels() {
    std::string path = "../../data/mini-ch.fmi";

    ///////////////////////////////////////////////////////////////////////////////////
    /// LOAD THE ORIGINAL GRAPH
    ///////////////////////////////////////////////////////////////////////////////////
    auto tup = load_fmi_ch(path,  -1, true);
    auto systems = std::get<0>(tup);
    Graph graph = std::get<1>(tup);
    int number_of_nodes = graph.adj.size();
    // init the hub labels for faster shortest path distance calculation
    graph.init_hub_labels();
    bool test = true;

    for (int i = 0; i < 100; i++) {
        int random_source = rand() % number_of_nodes;
        int random_target = rand() % number_of_nodes;
        auto original_dist = graph.dijkstra(random_source, random_target).second;
        auto hl_dist = graph.hl_distance(random_source, random_target);
        if (hl_dist != original_dist) {
            std::cout << "Hub label distance does not match Dijkstra distance for nodes "
                      << random_source << " and " << random_target << ": "
                      << hl_dist << " != " << original_dist << std::endl;
            test = false;
        }
        assert(test);
    }
}

int main(){
    testGraph();
    testHub_labels();
    return 0
    ;
}