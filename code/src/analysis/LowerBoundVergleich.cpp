//
// Created by Sebastian Paule on 9/4/25.
//

#include "../io/Dataloader.hpp"
#include "../spanner/ThetaSpanner.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

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

int main(int argc, char* argv[]) {


    Graph graph = get<1>(load_fmi("../../data/0025daniel.fmi"));
    std::cout << "graph loaded" << endl;

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
    auto csv = loadCSV("../../data/resultsBS1.csv");
    for (int i = 1; i < csv.size(); i++) {
        // s,t,triangle,theta

        int src = stoi(csv[i][0]);
        int trg = stoi(csv[i][1]);

        if (src >= graph.adj.size() || trg >= graph.adj.size()) {
            std::cerr << "Node ID out of range: " << src << " or " << trg << std::endl;
            continue;
        }

        //auto results = graph.multiSourceMultiTargetDijkstra({src}, {trg});
        //double true_dist = results.empty() ? -1.0 : results[0].second;
        double true_dist;
        {
            auto t0 = std::chrono::high_resolution_clock::now();
            true_dist = graph.dijkstra(src, trg, false).second;
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
