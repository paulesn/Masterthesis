#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <algorithm>
#include <vector>
#include <utility>
#include <iostream>
#include <unordered_set>


struct pair_hash {
    size_t operator()(const std::pair<int, int>& p) const {
        return std::hash<int>()(p.first) ^ (std::hash<int>()(p.second) << 1);
    }
};

// Define Point structure
struct Point {
    double x, y;
    int id; // Unique identifier for the point, -1 if not set

    Point(double x_, double y_);
    Point(double x_, double y_, int id_, int level_=-1);

    bool operator==(const Point& other) const;
    bool operator<(const Point& other) const;

    friend std::ostream& operator<<(std::ostream& os, const Point& p);
};

struct Edge {
    int source;
    int target;
    double weight;
    Edge(const int source_, const int target_, const double weight_) : source(source_), target(target_), weight(weight_) {}
    Edge() : source(-1), target(-1), weight(0.0) {}
    bool operator<(const Edge& other) const {return weight < other.weight;}
};

// Define Graph class
struct Graph {
public:
    explicit Graph(int nodes);

    explicit Graph(const std::vector<Point>& points);

    explicit Graph(const Graph* graph);

    std::vector<Point> id_point_map;
    std::unordered_set<std::pair<int, int>, pair_hash> existance; // to check if an edge exists

    int number_of_edges = 0;
    std::vector<std::vector<std::tuple<int,int>>> forward_hub_labels; // forward hub labels for each node
    std::vector<std::vector<std::tuple<int,int>>> backward_hub_labels; // backward hub labels for each node

    void addEdge(int u, int v, double w);
    std::pair<std::vector<int>, double> dijkstra(int src, int dest, double maximum=-1);

    /**
     * This is a variant of the dijkstra algorithm designed by spira et al. in the paper
     * https://link.springer.com/article/10.1007/BF01190847?utm_source=chatgpt.com
     * @param src the node id of the source node
     * @param dest the node id of the target node
     * @param maximum the maximum distance explored. currently not used. only incuded to allow drop in replacement with the dijkstra implementation
     * @return
     */
    double spira_sp(int src, int dest, double maximum=-1);
    std::vector<std::pair<std::vector<int>, double>> multiSourceMultiTargetDijkstra(
        const std::vector<int>& sources,
        const std::vector<int>& targets,
        bool all = false
    );

    double hl_distance(int source, int target);



    void sort_adj() {
        for (int i = 0; i < n; ++i) {
            std::sort(adj[i].begin(), adj[i].end(), [](Edge& a, Edge& b) {return a.weight < b.weight;});
        }
    }

    int n;
    std::vector<std::vector<Edge>> adj;

    void attach_lone_points();

    void store_to_disk(const std::string& path) const;

    void draw();
};

#endif // GRAPH_HPP
