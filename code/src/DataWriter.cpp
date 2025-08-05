//
// Created by sebastian on 04.07.25.
//

#include "DataWriter.h"
#include "Graph.hpp"

#include <fstream>
#include <iostream>
#include <math.h>

WGS84Coordinates pointInWGS84(Point p) {
    double lon = (p.x / 20037508.34) * 180;
    double lat = (p.y / 20037508.34) * 180;

    lat = 180 / M_PI * (2 * atan(exp(lat * M_PI / 180)) - M_PI / 2);
    return std::make_pair(lat, lon);
}

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
    outFile << (g.number_of_edges) << std::endl;
    int counter_nodes = 0;
    for (const auto& point : g.id_point_map) {
        counter_nodes++;

        auto wgs84 = pointInWGS84(point);

        double coord_x = wgs84.first;
        double coord_y = wgs84.second;

        outFile << coord_x << " " << coord_y  << std::endl;
    }
    int counter = 0;
    int max_idx = 0;
    for (int u = 0; u < g.n; ++u) {
        for (const auto& edge : g.adj[u]) {
            if (edge.target > u) { // Avoid duplicate edges in undirected graph
                outFile << u << " " << edge.target << " " << "1 1" << std::endl; // the last two numbers are line thikness and color
                counter++;
                if (edge.target > max_idx) {
                    max_idx = edge.target;
                }
                if (u > max_idx) {
                    max_idx = u;
                }
            }
        }
    }
    std::cout << counter_nodes << " " << max_idx << std::endl;
    outFile.close();
    return 0;
}

int write_gf_with_highlight(const Graph& g, std::string path, std::vector<int> highlight_nodes) {
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
    outFile << (g.number_of_edges) << std::endl;
    int counter_nodes = 0;
    for (const auto& point : g.id_point_map) {
        counter_nodes++;

        auto wgs84 = pointInWGS84(point);

        double coord_x = wgs84.first;
        double coord_y = wgs84.second;

        outFile << coord_x << " " << coord_y  << std::endl;
    }
    int counter = 0;
    int max_idx = 0;
    // Write highlighted nodes
    for (int u = 0; u < highlight_nodes.size()-1; u++) {
        int i = highlight_nodes[u];
        for (const auto& edge : g.adj[i]) {
            if (edge.target > i) { // Avoid duplicate edges in undirected graph
                outFile << i << " " << edge.target << " " << "1 1" << std::endl; // the last two numbers are line thikness and color
                counter++;
                if (edge.target > max_idx) {
                    max_idx = edge.target;
                }
                if (i > max_idx) {
                    max_idx = i;
                }
            }
        }
        if (highlight_nodes.size()-1 == i) {
            break;
        }
        int s = highlight_nodes[i];
        int t = highlight_nodes[i+1];
        outFile << s << " " << t << " " << "5 5" << std::endl; // the last two numbers are line thikness and color
        std::cout << "Highlighting edge from " << s << " to " << t << std::endl;
    }
    std::cout << counter_nodes << " " << max_idx << std::endl;
    outFile.close();
    return 0;
}