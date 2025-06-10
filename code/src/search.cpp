
#include <vector>
#include <tuple>
#include <unordered_map>
#include <queue>
#include <cmath>
#include <limits>
#include <algorithm>
namespace tree {
    #include "Octtree.cpp"
}
using namespace std;

vector<int> aStar(
    const vector<tree::Star>& nodes,
    int startID,
    int targetID,
    double jumpRange
) {
    int n = nodes.size();

    auto otree = new tree::Octtree(10, 1000000);
    for (auto& system : nodes) {
        otree->insert(&system);
    }

    // map ids to indices
    unordered_map<int,int> id2idx;
    for (int i = 0; i < n; ++i) {
        id2idx[nodes[i].id] = i;
    }
    if (!id2idx.count(startID) || !id2idx.count(targetID))
        return {};  // start or goal not in nodes

    int start = id2idx[startID];
    int target  = id2idx[targetID];

    const double INF = numeric_limits<double>::infinity();
    vector<double> gScore(n, INF), fScore(n, INF);
    vector<int> cameFrom(n, -1);
    vector<bool> closed(n, false);

    // min-heap by fScore
    struct Entry {
        int idx;
        double f;
        bool operator<(Entry const& o) const { return f > o.f; }
    };
    priority_queue<Entry> open;

    auto heuristic = [&](int i) {
        double dx = nodes[i].x - nodes[target].x;
        double dy = nodes[i].y - nodes[target].y;
        double dz = nodes[i].z - nodes[target].z;
        return sqrt(dx*dx + dy*dy + dz*dz);
    };

    gScore[start] = 0.0;
    fScore[start] = heuristic(start);
    open.push({start, fScore[start]});

    while (!open.empty()) {
        int current = open.top().idx;
        open.pop();
        if (closed[current]) continue;
        if (current == target) {
            // reconstruct path of ids
            vector<int> path;
            for (int u = target; u != -1; u = cameFrom[u])
                path.push_back(nodes[u].id);
            reverse(path.begin(), path.end());
            return path;
        }
        closed[current] = true;

        // precompute threshold squared
        double mult = nodes[current].jump_multiplicator;
        double thresh = jumpRange * mult;
        double thresh2 = thresh * thresh; // to prevent multiple root calculations

        // coords of current
        double cx = nodes[current].x;
        double cy = nodes[current].y;
        double cz = nodes[current].z;

        // examine all other nodes
        // TODO use a better data structure for this
        for (auto s: otree->find(cx, cy, cz, thresh)) {
            auto nei = id2idx[s->id];
            if (nei == current || closed[nei]) continue;
            double dx = nodes[nei].x - cx;
            double dy = nodes[nei].y - cy;
            double dz = nodes[nei].z - cz;
            double dist2 = dx*dx + dy*dy + dz*dz;
            if (dist2 >= thresh2) continue;  // no edge

            double tentativeG = gScore[current] + 1;
            if (tentativeG < gScore[nei]) {
                cameFrom[nei] = current;
                gScore[nei] = tentativeG;
                fScore[nei] = tentativeG + heuristic(nei);
                open.push({nei, fScore[nei]});
            }
        }
    }

    // no path found
    return {};
}
