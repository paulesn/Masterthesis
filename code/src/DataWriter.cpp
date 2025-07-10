//
// Created by sebastian on 04.07.25.
//

#include "DataWriter.h"
#include "Graph.hpp"

#include <fstream>
#include <iostream>

int write_gf(const Graph& g, std::string path) {
    // Open file with default mode: out | trunc (creates file if it doesn't exist, overwrites if it does)
    std::ofstream outFile(path);

    // calculate coordinate transformation to project graph on a globe
    // the current point coordinates are xy coordinates
    // we need to transform them into degree coordinates
    double min_x = 0;
    double max_x = 0;
    double min_y = 0;
    double max_y = 0;
    for (int n = 0; n < g.n; n++) {
        auto node = g.id_point_map[n];
        if (node.x > max_x) max_x = node.x;
        if (node.y > max_y) max_y = node.y;
        if (node.x < min_x) min_x = node.x;
        if (node.y < min_y) min_y = node.y;
    }

    double mod_x = (max_x - min_x)/360;
    double mod_y = (max_y - min_y)/360;


    if (!outFile) {
        std::cerr << "Failed to create or open the file." << std::endl;
        return 1;
    }

    outFile << g.n << std::endl;
    outFile << (g.number_of_edges/2)-1 << std::endl;
    for (const auto& point : g.id_point_map) {

        double coord_x = (point.x+(-1*min_x)) * mod_x;
        double coord_y = (point.y+(-1*min_y)) * mod_y;

        outFile << point.id << " " << coord_x << " " << coord_y  << std::endl;
    }
    int counter = 0;
    for (int u = 0; u < g.n; ++u) {
        for (const auto& edge : g.adj[u]) {
            if (edge.target > u) { // Avoid duplicate edges in undirected graph
                outFile << counter << " " << u << " " << edge.target << " " << "1 1" << std::endl; // the last two numbers are line thikness and color
                counter++;
            }
        }
    }

    outFile.close();
    std::cout << counter << " vs " << (g.number_of_edges/2)-1;
    return 0;
}