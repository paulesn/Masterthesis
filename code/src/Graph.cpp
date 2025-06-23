#include "Graph.hpp"
#include <vector>
#include <queue>
#include <limits>
#include <utility>
#include <algorithm>
#include <chrono>
#include <csignal>
#include <iostream>
#include <unordered_set>

#include "HubLabels.hpp"

using namespace std;



Point::Point(const double x_, const double y_, int level_): x(x_), y(y_), id(-1), level(level_) {}

Point::Point(const double x_, const double y_, const int id_, int level_): x(x_), y(y_), id(id_), level(level_) {}

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
void Graph::addEdge(int u, int v, double w, bool undirected) {
    undirected = true; // TODO remove
    if (u > v) std::swap(u, v); // Ensure u is always less than v for consistent edge representation
    if (u < 0 || u >= n || v < 0 || v >= n) {
        cout << "Invalid node index: " << u << " or " << v << endl;
        raise(SIGINT);
    }
    // Check if the edge already exists
    if (existance.find({u, v}) != existance.end()) {
        return; // Edge already exists, do not add again
    }

    Edge e1 {u, v, w};
    auto &edges_u = adj[u];
    auto it = lower_bound(edges_u.begin(), edges_u.end(), e1);
    edges_u.insert(it, e1);
    existance.emplace(u, v);// Mark the edge as existing

    if (undirected) {
        Edge e2{v, u, w};
        auto &edges_v = adj[v];
        it = lower_bound(edges_v.begin(), edges_v.end(), e2);
        edges_v.insert(it, e2);
        number_of_edges++;
    }
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
    priority_queue<NodeDist, vector<NodeDist>, greater<NodeDist>> pq;
    pq.push({0.0, src});

    while (!pq.empty()) {
        auto [d, u] = pq.top();
        pq.pop();
        if (d > dist[u]) continue;
        if (u == dest) break;  // Stop early if we reached destination
        if (maximum >= 0 && d > maximum) {
            // If we have a maximum distance and the current distance exceeds it, stop
            return {{}, maximum};
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
    priority_queue<NodeDist, vector<NodeDist>, greater<NodeDist>> pq;

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

void Graph::init_hub_labels() {
    vector<vector<int>> levels;
    // generate upwards edges
    std::vector<std::vector<std::tuple<int,int>>> upwards_edges (n);
    for (int u = 0; u < n; ++u) {
        // add to the levels vector
        Point* u_point = &id_point_map[u];
        if (u_point->level < 0) {
            cout << "Node " << u << " has no level assigned, cannot generate hub labels." << endl;
            raise(SIGINT); // Raise SIGINT to terminate the program
        }
        if (levels.size() <= u_point->level) {
            levels.resize(u_point->level + 1, {});
        }
        levels[u_point->level].push_back(u);
        for (const auto &edge : adj[u]) {
            int v = edge.target;
            Point* v_point = &id_point_map[v];
            if (v_point->level > u_point->level) { // Only consider edges going upwards in the order
                upwards_edges[u].emplace_back(v, edge.weight);
            }
        }
    }
    // generate downwards edges
    std::vector<std::vector<std::tuple<int,int>>> downwards_edges (n);
    for (int u = 0; u < n; ++u) {
        for (const auto &edge : adj[u]) {
            int v = edge.target;
            Point* u_point = &id_point_map[u];
            Point* v_point = &id_point_map[v];
            if (v_point->level < u_point->level) { // Only consider edges going downwards in the order
                downwards_edges[u].emplace_back(v, edge.weight);
            }
        }
    }
    // call hublabel generator
    auto hub_labels = generate_hub_labels(levels, upwards_edges, downwards_edges, 1);
    // store hub labels in a suitable data structure
    forward_hub_labels = std::get<0>(hub_labels);
    backward_hub_labels = std::get<1>(hub_labels);

}

double Graph::hl_distance(int source, int target) {
    auto fl = forward_hub_labels[source];
    auto bl = backward_hub_labels[target];
    auto dist =  optimized_query(fl, bl);
    return get<1>(dist);
}


double ::Graph::longestShortestPath(const vector<Point>& set, int source) {
    // check trivial cases
    if (set.empty()) {
        //cout << "Set is empty, cannot calculate longest shortest path." << endl;
        return 0.0;
    }
    if (set.size() == 1) {
        //cout << "Set has only one point, cannot calculate longest shortest path." << endl;
        return 0.0; // No valid path to calculate
    }
    // set default source
    if (source==-1) source = set[0].id; // Use the first point as source if not specified


    double max_dist = 0;
    int new_source = -1;
    for (const auto &p : set) {
        double dist = hl_distance(source, p.id);
        if (max_dist < dist) max_dist = dist;
        new_source = p.id;
    }
    if (new_source == -1) {
        cout << "No valid source found in the set." << endl;
        raise(SIGINT); // Raise SIGINT to terminate the program
    }

    // now we have the longest path distance, we can use it to find the maximum distance from the new source to all targets
    // which is the radius of the set
    for (const auto &p : set) {
        double dist = hl_distance(new_source, p.id);
        if (max_dist < dist) max_dist = dist;
        new_source = p.id;
    }
    if (new_source == -1) {
        cout << "No valid source found in the set." << endl;
        raise(SIGINT); // Raise SIGINT to terminate the program
    }
    return max_dist; // Return the maximum distance found
}

/**
 *
 * @return if two sets of nodes are well-separated by using the shortest path distance.
 */
std::tuple<std::vector<int>,double> Graph::wspdCheck(std::vector<Point> &set1, std::vector<Point> &set2, double& set1_prec_dist, double& set2_prec_dist, double s) {

    if (set1.empty() || set2.empty()) {
        //cout << "One of the sets is empty, cannot calculate WSPD." << endl;
        return {}; // Return empty path if one of the sets is empty
    }

    // check if all indicies are in the graph if not, return 0.0
    for (const auto p : set1) {
        if (p.id < 0 || p.id >= n) {
            cout << "Invalid node index in set1: " << p.id << endl;
            raise (SIGINT); // Raise SIGINT to terminate the program
        }
    }
    for (const auto p : set2) {
        if (p.id < 0 || p.id >= n) {
            cout << "Invalid node index in set1: " << p.id << endl;
            raise (SIGINT); // Raise SIGINT to terminate the program
        }
    }

    // calculate the shortest longest path for both sets
    double rad1 = (set1_prec_dist==-1) ? longestShortestPath(set1) : set1_prec_dist;
    set1_prec_dist = rad1;
    double rad2 = (set2_prec_dist==-1) ? longestShortestPath(set2) : set2_prec_dist;
    set2_prec_dist = rad2;

    // calculate the distance between the two sets using multi-source multi-target Dijkstra
    vector<int> sources;
    vector<int> targets;
    for (const auto &p : set1) {
        if (p.id >= 0 && p.id < n) {
            sources.push_back(p.id);
        } else {
            cout << "Invalid node index in set1: " << p.id << endl;
            raise(SIGINT); // Raise SIGINT to terminate the program
        }
    }
    for (const auto &p : set2) {
        if (p.id >= 0 && p.id < n) {
            targets.push_back(p.id);
        } else {
            cout << "Invalid node index in set2: " << p.id << endl;
            raise(SIGINT); // Raise SIGINT to terminate the program
        }
    }
    auto distances = multiSourceMultiTargetDijkstra(
        sources, targets, true
    );
    if (distances.empty()) {
        cout << "No path found between the two sets." << endl;
        return std::make_tuple(std::vector<int>{}, std::numeric_limits<double>::infinity()); // No path found
    }
    double dist = distances[0].second; // Get the distance from the first source to the first target
    vector<int> ret_path = distances[0].first;
    for (auto &result : distances) {
        auto &path = result.first;
        double d = result.second;
        if (d < dist) {
            dist = d; // Find the minimum distance
            ret_path = result.first; // Update the path
        }
    }

    if (dist >= s * max(rad1, rad2)) {
        return {ret_path, dist}; // TODO maybe this is uneccesary repacking
    }
    return std::make_tuple(std::vector<int>{}, std::numeric_limits<double>::infinity()); // No path found
}
