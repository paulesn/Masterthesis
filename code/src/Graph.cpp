#include <vector>
#include <queue>
#include <limits>
#include <utility>
#include <algorithm>

using namespace std;

struct Graph {
    int n;
    vector<vector<pair<int,double>>> adj;

    // Constructor: initialize graph with given number of nodes
    Graph(int nodes) : n(nodes), adj(nodes) {}

    // Add edge from u to v with weight w. If undirected is true, also add edge from v to u.
    void addEdge(int u, int v, double w, bool undirected = true) {
        adj[u].push_back({v, w});
        if (undirected) {
            adj[v].push_back({u, w});
        }
    }

    // Dijkstra's algorithm: returns pair of (shortest path as list of nodes, total distance)
    pair<vector<int>, double> dijkstra(int src, int dest) {
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

            for (auto &edge : adj[u]) {
                int v = edge.first;
                double w = edge.second;
                double nd = d + w;
                if (nd < dist[v]) {
                    dist[v] = nd;
                    prev[v] = u;
                    pq.push({nd, v});
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
};
