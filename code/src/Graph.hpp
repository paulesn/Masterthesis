#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <vector>
#include <utility>
#include <iostream>

// Define Point structure
struct Point {
    double x, y;
    int id; // Unique identifier for the point, -1 if not set

    Point(double x_, double y_);
    Point(double x_, double y_, int id_);

    bool operator==(const Point& other) const;
    bool operator<(const Point& other) const;

    friend std::ostream& operator<<(std::ostream& os, const Point& p);
};

struct Edge {
    int source;
    int target;
    double weight;
    Edge(const int source_, const int target_, const double weight_) : source(source_), target(target_), weight(weight_) {}
    bool operator<(const Edge& other) const {return weight < other.weight;}
};

// Define Graph class
struct Graph {
public:
    explicit Graph(int nodes);

    explicit Graph(const std::vector<Point>& points);

    std::vector<Point> id_point_map;
    std::vector<std::vector<bool>> existance; // edge existance matrix for faster edge checks

    int number_of_edges = 0;

    void addEdge(int u, int v, double w, bool undirected = true);
    std::pair<std::vector<int>, double> dijkstra(int src, int dest);
    std::vector<std::pair<std::vector<int>, double>> multiSourceMultiTargetDijkstra(
        const std::vector<int>& sources,
        const std::vector<int>& targets,
        bool all = false
    );

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
    std::vector<int> wspdCheck(std::vector<Point> &set1, std::vector<Point> &set2, double &set1_prec_dist, double &set2_prec_dist, double s=2);

    int n;
    std::vector<std::vector<Edge>> adj;
};

#endif // GRAPH_HPP
