#include "Graph.hpp"
#include <vector>
#include <queue>
#include <limits>
#include <utility>
#include <algorithm>
#include <chrono>
#include <csignal>
#include <fstream>
#include <iostream>
#include <unordered_set>
#include <oneapi/tbb/task_arena.h>

#include "HubLabels.hpp"
#include "Quadtree.hpp"

using namespace std;



Point::Point(const double x_, const double y_): x(x_), y(y_), id(-1){}

Point::Point(const double x_, const double y_, const int id_): x(x_), y(y_), id(id_) {}

std::ostream& operator<<(std::ostream& os, const Point& p) {
    os << "(" << p.x << ", " << p.y << ")";
    return os;
}


bool Point::operator==(const Point &other) const {
    return (x == other.x && y == other.y);
}

bool Point::operator<(const Point &other) const {
    return (x < other.x || (x == other.x && y < other.y));
}

Graph::Graph(int nodes) : n(nodes), adj(nodes){
    for (int i = 0; i < n; ++i) {
        id_point_map.emplace_back(0.0, 0.0, i, -1); // Initialize points with default coordinates and unique IDs
    }
}
Graph::Graph(const vector<Point>& points) : n(points.size()), adj(n) {
    id_point_map = points;
    for (int i = 0; i < n; ++i) {
        id_point_map[i].id = i; // Assign unique IDs to points
    }
}

Graph::Graph(const Graph *graph) {
    // copy constructor
    n = graph->n;
    adj = graph->adj;
    id_point_map = graph->id_point_map;
    existance = graph->existance;
    number_of_edges = graph->number_of_edges;
    number_of_edges = graph->number_of_edges;
}

/**
 * dd edge from u to v with weight w. If undirected is true, also add edge from v to u.
 * @param u source node index
 * @param v target node index
 * @param w weight of the edge
 * @param undirected
 */
void Graph::addEdge(int u, int v, double w) {
    if (u > v) std::swap(u, v); // Ensure u is always less than v for consistent edge representation
    if (u < 0 || u >= n || v < 0 || v >= n) {
        cout << "Invalid node index: " << u << " or " << v << endl;
        raise(SIGINT);
    }
    // Check if the edge already exists // TODO remove
    if (existance.find({u, v}) != existance.end()) {
        return; // Edge already exists, do not add again
    }

    Edge e1 {u, v, w};
    auto &edges_u = adj[u];
    edges_u.emplace_back(e1); // this is not very efficant TODO
    existance.emplace(u, v);// Mark the edge as existing

    // As the graph is undirected, we also add the inverse direction
    Edge e2{v, u, w};
    auto &edges_v = adj[v];
    edges_v.emplace_back(e2);

    number_of_edges++;
}

/**
 * Dijkstra's algorithm: returns pair of (shortest path as list of nodes, total distance)
 * @param src source node index
 * @param dest destination node index
 * @return the path as a vector of node indices and the total distance as a double.
 */
pair<vector<int>, double> Graph::dijkstra(int src, int dest, double maximum) {
    const double INF = numeric_limits<double>::infinity();
    vector<double> dist(n, INF);
    vector<int> prev(n, -1);
    dist[src] = 0.0;

    // Min-heap priority queue: (distance, node)
    using NodeDist = pair<double, int>;
    priority_queue<NodeDist, vector<NodeDist>, greater<>> pq;
    pq.emplace(0.0, src);

    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();
        if (d > dist[u]) continue;
        if (u == dest) break;  // Stop early if we reached destination

        if (maximum >= 0 && d > maximum) {
            // If we have a maximum distance and the current distance exceeds it, stop
            return {{}, INF};
        }


        for (auto &edge : adj[u]) {
            int v = edge.target;
            double w = edge.weight;
            double nd = d + w;
            if (nd < dist[v]) {
                dist[v] = nd;
                prev[v] = u;
                pq.emplace(nd, v);
            }
        }
    }

    vector<int> path;
    if (dist[dest] == INF) {
        // No path found
        return {path, INF};
    }

    // Reconstruct path from dest to src
    for (int at = dest; at != -1; at = prev[at]) {
        path.push_back(at);
    }
    reverse(path.begin(), path.end());

    return {path, dist[dest]};
}

using NodeDist = pair<double, Edge>;
bool q_comp(NodeDist a, NodeDist b) {
    return a.first > b.first;
}

double Graph::spira_sp(int src, int dest, double maximum) {
    const double INF = numeric_limits<double>::infinity();
    vector<double> dist(n, INF);
    vector<int> idx(n, 0);
    vector<int> prev(n, -1);
    dist[src] = 0.0;

    // Min-heap priority queue: (distance, node)

priority_queue<NodeDist, vector<NodeDist>, decltype(&q_comp)> pq(q_comp);

    // add the smallest edge in the list
    pq.emplace(adj[src][0].weight, adj[src][0]);
    idx[src] = 0;
    // iterate
    while (!pq.empty()) {

        auto [d, e] = pq.top();
        pq.pop();

        int c_src = e.source;
        int c_target = e.target;
        //std::cout << "----------------------------------------------------------------------" << endl;
        //std::cout << "exploring edge: (" << c_src << " - " << c_target <<")" << std::endl;
        //std::cout << c_src << "has " << adj[c_src].size() << " edges. This is edge Nr. " << idx[c_src] << std::endl;

        // update distance if shorter
        if (d < dist[c_target]) {
            dist[c_target] = d;
            prev[c_target] = c_src;
        }

        if (c_target == dest) {
            return d;
        }

        // add the next edge if exists
        if (idx[c_src] < adj[c_src].size()) {
            Edge n_edge = adj[c_src][idx[c_src]++];
            pq.emplace(dist[n_edge.source] + n_edge.weight, n_edge);
            //std::cout << "Adding edge (" << n_edge.source << " - " << n_edge.target << "| " << d+ n_edge.weight <<")" << std::endl;
        }

        if (idx[c_target] < adj[c_target].size()) {
            Edge n_edge = adj[c_target][idx[c_target]++];
            pq.emplace(dist[n_edge.source] + n_edge.weight, n_edge);
            //std::cout << "Adding edge(" << n_edge.source << " - " << n_edge.target << "| " << d+ n_edge.weight <<")" << std::endl;
        }
    }
    return dist[dest];
}

vector<pair<vector<int>, double>> Graph::multiSourceMultiTargetDijkstra(
    const vector<int>& sources,
    const vector<int>& targets,
    bool all)
{
    const double INF = numeric_limits<double>::infinity();
    vector<double> dist(n, INF);
    vector<int> prev(n, -1);
    vector<bool> isTarget(n, false);
    std::unordered_set<int> remainingTargets;

    // Mark target nodes
    for (int t : targets) {
        if (t >= 0 && t < n) {
            isTarget[t] = true;
            remainingTargets.insert(t);
        }
    }

    // Min-heap priority queue: (distance, node)
    using NodeDist = pair<double, int>;
    priority_queue<NodeDist, vector<NodeDist>, greater<>> pq;

    for (int src : sources) {
        if (src >= 0 && src < n) {
            dist[src] = 0.0;
            pq.emplace(0.0, src);
        }
    }

    vector<pair<vector<int>, double>> results;

    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();

        if (d > dist[u]) continue;

        if (isTarget[u]) {
            // Reconstruct path
            vector<int> path;
            for (int at = u; at != -1; at = prev[at]) {
                path.push_back(at);
            }
            reverse(path.begin(), path.end());
            results.emplace_back(path, dist[u]);
            remainingTargets.erase(u);

            if (!all) break;
            if (remainingTargets.empty()) break;
        }

        for (const auto &e : adj[u]) {
            int v = e.target;
            double w = e.weight;
            double nd = d + w;
            if (nd < dist[v]) {
                dist[v] = nd;
                prev[v] = u;
                pq.emplace(nd, v);
            }
        }
    }

    return results;
}

double Graph::hl_distance(int source, int target) {
    auto fl = forward_hub_labels[source];
    auto bl = backward_hub_labels[target];
    auto dist =  optimized_query(fl, bl);
    return get<1>(dist);
}


void Graph::attach_lone_points() {
    for (int i=0; i<this->adj.size();i++){
        auto list = this->adj[i];
        if (list.empty()){
            Edge edge = Edge(i, i-1, euklidian_distance(this->id_point_map[i], this->id_point_map[i-1]));
            this->adj[i].emplace_back(edge);
        }
    }
}

void Graph::store_to_disk(const std::string &path) const {
    auto out = std::ofstream(path);
    out << "#\n#\n#\n#\n \n" << n << "\n" << number_of_edges << "\n";
    for (const auto &point : id_point_map) {
        out << point.id << " 0 " <<  point.x << " " << point.y << " 0 \n";
    }
    for (auto edge_list : adj) {
        for (const auto &edge : edge_list) {
            out << edge.source << " " << edge.target << " " << edge.weight << " 0 0 .\n";
        }
    }
    out.close();
}
