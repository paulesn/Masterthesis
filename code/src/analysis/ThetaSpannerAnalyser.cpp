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

    string base_graph_path;
    string spanner_path;
    string csv_path;


    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "-sg" && i + 1 < argc) {
            spanner_path = argv[++i];
            std::cout << "Spanner path: " << spanner_path << std::endl;
        } else if (arg == "-bg" && i + 1 < argc) {
            base_graph_path = argv[++i];
            std::cout << "Base graph path: " << base_graph_path << std::endl;
        } else if (arg == "-csv" && i + 1 < argc) {
            csv_path = argv[++i];
            std::cout << "CSV output path: " << csv_path << std::endl;
        }
    }

    if (spanner_path.empty() || base_graph_path.empty()) {
        std::cerr << "Usage: " << argv[0] << " -sg <string> -bg <string> -csv <string>\n";
        for (int i = 0; i < argc; i++) std::cerr << " " << argv[i];
        return 1;
    }

    ///////////////////////////////////////////////////////////////////////////////////
    /// LOAD THE GRAPHS
    ///////////////////////////////////////////////////////////////////////////////////

    auto res = load_fmi("../../data/0025daniel.fmi");
    Graph graph = std::get<1>(res);
    std::cout << "Loaded  graph\n" << std::endl;

    // load csv data
    int counter = 0;
    auto csv = loadCSV("../../data/resultsBS.csv");
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
        true_dist = graph.getDistance(src, trg);
        std::cout << counter << ",";
        for (int j = 0; j < csv[i].size(); j++) {
            std::cout << csv[i][j] << ",";
        }
        std::cout << true_dist << std::endl;
        counter++;
    }
    //get true dist for each s-t pair in the csv
}

// TODO auch die worst case edge ausrechnen

