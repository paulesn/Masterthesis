#ifndef GRAPH_HPP
#define GRAPH_HPP

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
struct Pointc {
    double x, y;
    int id; // Unique identifier for the point, -1 if not set
    int meta; // additional metadata if needed

    Pointc(double x_, double y_);
    Pointc(double x_, double y_, int id_, int meta_ = -1);

    bool operator==(const Pointc& other) const;
    bool operator<(const Pointc& other) const;
    int level; // level of CH

    friend std::ostream& operator<<(std::ostream& os, const Pointc& p);
};

struct Edge {
    int source;
    int target;
    double weight;
    double t_value = 0.0; // only used in some algorithms
    Edge(const int source_, const int target_, const double weight_) : source(source_), target(target_), weight(weight_) {}
    Edge() : source(-1), target(-1), weight(0.0) {}
    bool operator<(const Edge& other) const {return weight < other.weight;}
};

// Define Graph class
struct Graph {
public:
    explicit Graph(int nodes);

    explicit Graph(const std::vector<Pointc>& points);

    explicit Graph(const Graph* graph);

    std::vector<Pointc> id_point_map;
    std::unordered_set<std::pair<int, int>, pair_hash> existance; // to check if an edge exists

    size_t number_of_edges = 0;

    bool addEdge(int u, int v, double w, bool safe = true);
    bool hasEdge(int u, int v) {
        if (u < 0 || u >= n || v < 0 || v >= n) return false;
        return existance.find({u, v}) != existance.end();
    }
    std::pair<std::vector<int>, double> dijkstra(int src, int dest, double maximum=-1);
    std::vector<std::pair<std::vector<int>, double>> multiSourceMultiTargetDijkstra(
        const std::vector<int>& sources,
        const std::vector<int>& targets,
        bool all = false
    );
    std::vector<std::pair<std::vector<int>, double>> multiSourceMultiTargetEdgeFrontierDijkstra(
        const std::vector<int>& sources,
        const std::vector<int>& targets,
        bool all = false
    );

    double longestShortestPath(const std::vector<Pointc>& set, int source = -1);
    double getDistance(int source, int target);
    void sort_edges();

    int n;
    std::vector<std::vector<Edge>> adj;
};

#endif // GRAPH_HPP
