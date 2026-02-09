//
// Created by Sebastian Paule on 12/16/25.
//

#ifndef MASTERTHESIS_LOWERBOUND_GRAPH_H
#define MASTERTHESIS_LOWERBOUND_GRAPH_H

#include <string>
#include <fstream>
#include "../structure/Graph.hpp"

/**
 * This class is designed to be selfcontained so that it can be used in other projects more easily
 */
class Lowerbound_Graph {

    Graph graph;

public:

    Lowerbound_Graph() : graph(1) {
    };

    bool load_graph_from_fmi_file(std::string filepath);
    double get_lower_bound_path(int source, int target, int k);

    double identify_closest_cone_neighbor_offset(Pointc p_t);

    double identify_cone_offset(Pointc p_t, int k);
};


#endif //MASTERTHESIS_LOWERBOUND_GRAPH_H