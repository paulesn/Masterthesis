//
// Created by sebastian on 04.07.25.
//

#ifndef DATAWRITER_H
#define DATAWRITER_H

#include "Graph.hpp"
#include <fstream>
#include <iostream>

typedef std::pair<double, double> WGS84Coordinates;
/**
 * Writes the graph to a file in a format that can be read by the simple graph renderer.
 * The file will contain the number of nodes, the number of edges, and the edges themselves.
 * The edges are written in the format: source target weight
 * @param g The graph to write to the file
 * @return 0 on success, 1 on failure
 */
int write_gf(const Graph& g, std::string path);

/**
 * This function does the same as write_gf, but it also highlights certain nodes in the graph.
 * The highlighted nodes will be written to the file in a special format that can be read by the rednder to render them in a different color.
 * (The edges between the highlighted path will be rendered in a different color as well.)
 * @param g
 * @param path
 * @param highlight_nodes
 * @return
 */
int write_gf_with_highlight(const Graph& g, std::string path, std::vector<int> highlight_nodes);

WGS84Coordinates pointInWGS84(Point p);

#endif //DATAWRITER_H
