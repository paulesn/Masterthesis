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
#include "Quadtree.hpp"

using namespace std;



Pointc::Pointc(const double x_, const double y_): x(x_), y(y_), id(-1) {}

Pointc::Pointc(const double x_, const double y_, const int id_): x(x_), y(y_), id(id_) {}

std::ostream& operator<<(std::ostream& os, const Pointc& p) {
    os << "(" << p.x << ", " << p.y << ")";
    return os;
}


bool Pointc::operator==(const Pointc &other) const {
    return (x == other.x && y == other.y);
}

bool Pointc::operator<(const Pointc &other) const {
    return (x < other.x || (x == other.x && y < other.y));
}

Graph::Graph(int nodes) : n(nodes), adj(nodes){
    for (int i = 0; i < n; ++i) {
        id_point_map.emplace_back(0.0, 0.0, i); // Initialize points with default coordinates and unique IDs
    }
}
Graph::Graph(const vector<Pointc>& points) : n(points.size()), adj(n) {
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
bool Graph::addEdge(int u, int v, double w, bool safe) {
    bool undirected = true; // TODO remove
    if (u > v) std::swap(u, v); // Ensure u is always less than v for consistent edge representation
    if (u < 0 || u >= n || v < 0 || v >= n) {
        cout << "Invalid node index: " << u << " or " << v << endl;
        return false;
    }
    // Check if the edge already exists // TODO remove
    if (safe) {
        if (existance.find({u, v}) != existance.end()) {
            return false; // Edge already exists, do not add again
        }
    }

    Edge e1 {u, v, w};
    auto &edges_u = adj[u];
    edges_u.emplace_back(e1); // this is not very efficant TODO
    if (safe) existance.emplace(u, v);// Mark the edge as existing

    if (undirected) {
        Edge e2{v, u, w};
        auto &edges_v = adj[v];
        edges_v.emplace_back(e2);
        number_of_edges++;
    }
    number_of_edges++;
    return true; // Edge added successfully
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

    vector<int> path = {};
    if (dist[dest] == INF) {
        // No path found
        return {path, INF};
    }

    // Reconstruct path from dest to src
    if (prev.size() == 0) {
        return {path, INF};
    }
    for (int at = dest; at != -1; at = prev[at]) {
        if (at < 0 || at >= n) {
            // Invalid node index encountered, return empty path
            std::cerr << "Invalid node index encountered during path reconstruction: " << at << std::endl;
        }
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
    vector<double> dist(adj.size(), INF);
    vector<int> prev(adj.size(), -1);
    vector<bool> isTarget(adj.size(), false);
    std::unordered_set<int> remainingTargets;

    // Mark target nodes
    for (int t : targets) {
        if (t >= 0 && t < adj.size()) {
            isTarget[t] = true;
            remainingTargets.insert(t);
        }
    }

    // Min-heap priority queue: (distance, node)
    using NodeDist = pair<double, int>;
    priority_queue<NodeDist, vector<NodeDist>, greater<NodeDist>> pq;

    for (int src : sources) {
        if (src >= 0 && src < adj.size()) {
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

vector<pair<vector<int>, double>> Graph::multiSourceMultiTargetEdgeFrontierDijkstra(
    const vector<int>& sources,
    const vector<int>& targets,
    bool all)
{
    //std::cout << "Using multiSourceMultiTargetEdgeFrontierDijkstra This means the edges have to be sorted according to weight" << std::endl;
    const double INF = numeric_limits<double>::infinity();
    vector<double> dist(adj.size(), INF);
    vector<int> prev(adj.size(), -1);
    vector<int> edgeCounter(adj.size(), 0);
    vector<bool> isTarget(adj.size(), false);
    std::unordered_set<int> remainingTargets;

    // Mark target nodes
    for (int t : targets) {
        if (t >= 0 && t < adj.size()) {
            isTarget[t] = true;
            remainingTargets.insert(t);
        }
    }

    // Min-heap priority queue: (distance, source, target)
    using NodeDist = tuple<double, int, int>;
    priority_queue<NodeDist, vector<NodeDist>, greater<NodeDist>> pq;

    for (int src : sources) {
        if (src >= 0 && src < adj.size()) {
            dist[src] = 0.0;
            if (!adj[src].empty()){
                pq.emplace(adj[src][0].weight, src, adj[src][0].target);
                //std::cout << "Adding initial edge from source " << src << " to target " << adj[src][0].target << " with weight " << adj[src][0].weight << std::endl;
                edgeCounter[src]++;
            }
        }
    }

    vector<pair<vector<int>, double>> results;

    while (!pq.empty()) {
        auto [d, v, u] = pq.top();
        //std::cout << "Popping edge from " << v << " to " << u << " with distance " << d << std::endl;
        pq.pop();


        // 2. We've processed an edge from 'v', so queue its *next* edge.
        //    (This block is correct and belongs here).
        if (edgeCounter[v] < adj[v].size()) {
            double nd = dist[v] + adj[v][edgeCounter[v]].weight;
            if (edgeCounter[v] + 1 < adj[v].size() && adj[v][edgeCounter[v]].weight > adj[v][edgeCounter[v]+1].weight) {
                //std::cerr << "Warning: Edges are not sorted by weight for node " << v << std::endl;
                return {};
            }
            pq.emplace(nd, v, adj[v][edgeCounter[v]].target);
            //std::cout << "Expanding edge from " << v << " to " << adj[v][edgeCounter[v]].target << " with weight " << adj[v][edgeCounter[v]].weight << std::endl;
            edgeCounter[v]++;
        }

        // 1. Prune stale paths.
        if (d > dist[u]) continue;

        // 3. Check if this is a *new* shortest path to 'u'.
        //    If d == dist[u], we've been here before. Do nothing.
        if (d < dist[u]) {
            // --- THIS IS THE NEW SHORTEST PATH ---

            // 4. Update distance and predecessor *first*.
            dist[u] = d;
            prev[u] = v;

            // 5. *Now* check if this newly-found node is a target.
            if (isTarget[u]) {
                vector<int> path;
                for (int at = u; at != -1; at = prev[at]) {
                    path.push_back(at);
                }
                reverse(path.begin(), path.end());
                // Now dist[u] correctly holds the value 'd'.
                results.emplace_back(path, dist[u]);
                remainingTargets.erase(u);

                if (!all) break;
                if (remainingTargets.empty()) break;
            }

            // 6. Expand from 'u' *only* on its first finalization.
            //    (Block 1 belongs *inside* here, with a safety check).
            if (edgeCounter[u] == 0 && !adj[u].empty()) {
                double nd = dist[u] + adj[u][0].weight; // 'd' is the new dist[u]
                pq.emplace(nd, u, adj[u][0].target);
                //std::cout << "Adding first edge from node " << u << " to target " << adj[u][0].target << " with weight " << adj[u][0].weight << std::endl;
                edgeCounter[u]++;
            }
        }
    }

    return results;
}






double ::Graph::longestShortestPath(const vector<Pointc>& set, int source) {
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
        double dist = get<1>(dijkstra(source, p.id));
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
        double dist = get<1>(dijkstra(new_source, p.id));
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
std::tuple<std::vector<int>,double> Graph::wspdCheck(std::vector<Pointc> &set1, std::vector<Pointc> &set2, double& set1_prec_dist, double& set2_prec_dist, double s) {

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

void Graph::attach_lone_points() {
    for (int i=0; i<this->adj.size();i++){
        auto list = this->adj[i];
        if (list.empty()){
            //Edge edge = Edge(i, i-1, euklidian_distance(this->id_point_map[i], this->id_point_map[i-1]));
            Edge edge = Edge(i, i-1, 1);
            this->adj[i].emplace_back(edge);
        }
    }
}

void Graph::sort_edges() {
    for (int i=0; i<this->adj.size();i++){
        auto &list = this->adj[i];
        std::sort(list.begin(), list.end(), [](const Edge& a, const Edge& b) {
            return a.weight < b.weight;
        });
    }
}
