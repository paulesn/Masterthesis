//
// Created by Sebastian Paule on 9/4/25.
//

#include "../io/Dataloader.hpp"
#include <omp.h>

#include "MultiThetaAnalysis.h"
#include "../spanner/ThetaSpanner.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include "../daniel/theta-graph/headers/Triangulation.h"

using namespace std;

// Function to load a CSV file into a 2D vector of strings
std::vector<std::vector<std::string>> loadCSV(const std::string& filename) {
    std::vector<std::vector<std::string>> data;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return data;
    }

    std::string line;
    // Read the file line by line
    while (std::getline(file, line)) {
        std::vector<std::string> row;
        std::stringstream lineStream(line);
        std::string cell;

        // Split the line by commas
        while (std::getline(lineStream, cell, ',')) {
            row.push_back(cell);
        }

        // Special handling if the line ends with a comma (adds an empty cell)
        if (!line.empty() && line.back() == ',') {
            row.push_back("");
        }

        data.push_back(row);
    }

    file.close();
    return data;
}

double triang_dijkstra(int src, int dest, Triangulation triang) {
    const double INF = numeric_limits<double>::infinity();
    int n = triang.getPoints()->size();
    vector<double> dist(n, INF);
    vector<int> edge_count(n, 0);
    dist[src] = 0.0;

    // Min-heap priority queue: (distance, node)
    using NodeDist = pair<double, int>;
    priority_queue<NodeDist, vector<NodeDist>, greater<NodeDist>> pq;
    pq.push({0.0, src});

    while (!pq.empty()) {
        auto [d, source] = pq.top();
        pq.pop();
        if (d > dist[source]) continue;
        if (source == dest) break;  // Stop early if we reached destination

        std::vector<GlobalID> temp = triang.oneToAllVisibility(source, false);
        for (auto t : temp) {
            int target = (int) t;
            Point u = triang.getPoint(source);
            Point v = triang.getPoint(target);
            double w = u.distance(v);
            double nd = d + w;
            if (nd < dist[target]) {
                dist[target] = nd;
                pq.emplace(nd, target);
            }
        }
    }
    return dist[dest];
}

int main(int argc, char* argv[]) {


    Triangulation triang;
    triang.readFromGraph("../../data/coastlines25.txt.graph");
    std::cout << "triangulation loaded" << endl;

    /*
     *This code tests the triangulation dijkstra against the graph dijkstra for 100 random pairs of nodes. It counts how many times the distances do not match and prints the results.
    Graph graph = get<1>(load_fmi("../../data/Archipelago.vis.fmi"));
    int miss = 0;
    for (int i = 0; i < 100; i++) {
        int s = (int)(random() % graph.adj.size());
        int t = (int)(random() % graph.adj.size());
        auto res = graph.dijkstra(s,t);
        double dist = triang_dijkstra(s,t,triang);
        std::cout << "Distance from " << s << " to " << t << ": " << res.second << " and " << dist << std::endl;
        if (res.second != dist) {
            miss++;
        }
    }
    cout << "There where " << miss << " misses" << endl;

    */

    // load csv data
    int counter = 0;
    auto csv = loadCSV("../../data/resultsBS.csv");
    for (int i = 1; i < csv.size(); i++) {
        // s,t,triangle,theta

        int src = stoi(csv[i][0]);
        int trg = stoi(csv[i][1]);

        if (src >= triang.getPoints()->size() || trg >= triang.getPoints()->size()) {
            std::cerr << "Node ID out of range: " << src << " or " << trg << std::endl;
            continue;
        }

        //auto results = graph.multiSourceMultiTargetDijkstra({src}, {trg});
        //double true_dist = results.empty() ? -1.0 : results[0].second;
        double true_dist;
        {
            auto t0 = std::chrono::high_resolution_clock::now();
            true_dist = triang_dijkstra(src, trg, triang);
            auto t1 = std::chrono::high_resolution_clock::now();
            auto dur = std::chrono::duration_cast<std::chrono::microseconds>(t1 - t0).count();

            std::cout << counter << ",";
            for (int j = 0; j < csv[i].size(); j++) {
                std::cout << csv[i][j] << ",";
            }
            std::cout << true_dist << ", " << dur << "us" << std::endl;
        }
        counter++;
    }
}
