//
// Created by sebastian on 13.06.25.
//

#ifndef THETASPANNER_H
#define THETASPANNER_H
#include "Graph.hpp"

/**
 *
 * @param graph the original graph
 * @param theta the number of edges to keep per node (number of zones)
 * @return
 */
Graph create_theta_spanner_graph(Graph* graph, int theta);

#endif //THETASPANNER_H
