//
// Created by sebastian on 04.07.25.
//

#ifndef DATAWRITER_H
#define DATAWRITER_H

#include "../structure/Graph.hpp"
#include <fstream>
#include <iostream>
/**
 * Writes the graph to a file in a format that can be read by the simple graph renderer.
 * The file will contain the number of nodes, the number of edges, and the edges themselves.
 * The edges are written in the format: source target weight
 * @param g The graph to write to the file
 * @return 0 on success, 1 on failure
 */
int write_gf(std::string path, int number_of_cones, std::vector<Pointc> path_to_draw, int cone_base=0);
int write_gf_with_graph(std::string path, int number_of_cones, std::vector<Pointc> path_to_draw, int cone_base=0, Graph graph = Graph(1));
int write_fmi(std::string path, const Graph& g);
#endif //DATAWRITER_H
