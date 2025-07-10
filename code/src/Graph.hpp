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
struct Point {
    double x, y;
    int id; // Unique identifier for the point, -1 if not set

    Point(double x_, double y_, int level_=-1);
    Point(double x_, double y_, int id_, int level_=-1);

    bool operator==(const Point& other) const;
    bool operator<(const Point& other) const;
    int level; // level of CH

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

    void addEdge(int u, int v, double w, bool undirected = true);
    std::pair<std::vector<int>, double> dijkstra(int src, int dest, double maximum=-1);
    std::vector<std::pair<std::vector<int>, double>> multiSourceMultiTargetDijkstra(
        const std::vector<int>& sources,
        const std::vector<int>& targets,
        bool all = false
    );

    double hl_distance(int source, int target);

    void init_hub_labels();

    double longestShortestPath(const std::vector<Point>& set, int source = -1);

    /**
     *
     * @param set1 a set of points in the graph
     * @param set1_prec_dist the inner radius of the set in shortest paht distance -1 if not known
     * @param set2 another set of points in the graph
     * @param set2_prec_dist the inner radius of the set in shortest path distance -1 if not known
     * @param s the seperation constant, the two sets are well-separated if the distance between them is at least s times the maximum of the inner radii
     * @return
     */
    std::tuple<std::vector<int>,double> wspdCheck(std::vector<Point> &set1, std::vector<Point> &set2, double &set1_prec_dist, double &set2_prec_dist, double s=2);

    int n;
    std::vector<std::vector<Edge>> adj;
};

#endif // GRAPH_HPP
