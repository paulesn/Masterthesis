//
// Created by sebastian on 04.07.25.
//

#include "DataWriter.h"
#include "Graph.hpp"

#include <fstream>
#include <iostream>
#include <cmath>

#include "Dataloader.hpp"
#include "Quadtree.hpp"

typedef std::pair<double, double> WGS84Coordinates;

WGS84Coordinates pointInWGS84(double x, double y) {
    double lon = (x / 20037508.34) * 180;
    double lat = (y / 20037508.34) * 180;

    lat = 180 / M_PI * (2 * atan(exp(lat * M_PI / 180)) - M_PI / 2);
    return std::make_pair(lat, lon);
}

int write_gf_with_graph(std::string path, int number_of_cones, std::vector<Point> path_to_draw, int cone_base, Graph graph) {

    if (cone_base > path_to_draw.size()) {
        std::cerr << "Cone base index is out of bounds." << std::endl;
        return -1;
    }

    int number_of_nodes = 0;
    int number_of_edges = 0;
    std::vector<Point> points;
    std::vector<std::tuple<int, int, int, int>> edges; // source, target, color, thickness

    // draw the path to draw
    for (int i = 0; i < path_to_draw.size(); ++i) {
        points.emplace_back(path_to_draw[i].x, path_to_draw[i].y, number_of_nodes);
        number_of_nodes++;
        if (i > 0) {
            edges.emplace_back(number_of_nodes - 2, number_of_nodes - 1, 1, 3);
            number_of_edges++;
        }
    }
    // draw the cones around the first point of the path
    for (int i = 0; i < number_of_cones; ++i) {
        int last_node_id = cone_base; // center node
        for (int length: {500000}) { // draw multiple lengths
            double angle = i * (2 * M_PI / number_of_cones);
            double x = path_to_draw[0].x + length * cos(angle);
            double y = path_to_draw[0].y + length * sin(angle);
            points.emplace_back(x, y, number_of_nodes);
            number_of_nodes++;
            edges.emplace_back(last_node_id, number_of_nodes - 1, 1, 2); // color 2 for cones
            last_node_id = number_of_nodes - 1;
            number_of_edges++;
        }
    }
    // draw the coastlines
    int offset = number_of_nodes;
    Graph coastline = load_coastline("../../data/coastlines-mercator.txt.pruned.wc.txt.shrunk.0025POI.32.txt");
    for (const auto &point : coastline.id_point_map) {
        points.emplace_back(point.x, point.y, number_of_nodes);
        number_of_nodes++;
    }
    for (const auto &edge_list : coastline.adj) {
        for (const auto &edge : edge_list) {
            if (edge.source < edge.target) // to avoid double edges in undirected graph
            {
                edges.emplace_back(edge.source + offset, edge.target + offset, 1, 1); // color 3 for coastline
                number_of_edges++;
            }
        }
    }
    if (graph.n > 1) {
        offset = number_of_nodes;
        for (const auto &point : graph.id_point_map) {
            points.emplace_back(point.x, point.y, number_of_nodes);
            number_of_nodes++;
        }
        for (const auto &edge_list : graph.adj) {
            for (const auto &edge : edge_list) {
                if (edge.source < edge.target) // to avoid double edges in undirected graph
                {
                    // check if edge source or edge target is sufficiantly close to path start
                    //if(euklidian_distance(graph.id_point_map[edge.source], path_to_draw[0]) > 2000 &&
                    //   euklidian_distance(graph.id_point_map[edge.target], path_to_draw[0]) > 2000){
                    //    continue;
                    //}
                    edges.emplace_back(edge.source + offset, edge.target + offset, 1, 4); // color 4 for spanner
                    number_of_edges++;
                }
            }
        }
    }

    // commit the data to the file
    auto out = std::ofstream(path);
    if (!out.is_open()) {
        std::cerr << "Failed to open file: " << path << std::endl;
        return -1;
    }
    out << number_of_nodes << std::endl;
    out << number_of_edges << std::endl;
    for (const auto &point : points) {
        auto wgs84 = pointInWGS84(point.x, point.y);
        out << wgs84.first << " " << wgs84.second << std::endl;
    }
    for (const auto &edge : edges) {
        out << std::get<0>(edge) << " " << std::get<1>(edge) << " " << std::get<2>(edge) << " " << std::get<3>(edge) << std::endl;
    }

    out.close();
    return 1;
}

int write_gf(std::string path, int number_of_cones, std::vector<Point> path_to_draw, int cone_base) {
    write_gf_with_graph(path, number_of_cones, path_to_draw, cone_base, Graph(1));
}
