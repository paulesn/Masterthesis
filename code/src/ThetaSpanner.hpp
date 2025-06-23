//
// Created by sebastian on 13.06.25.
//

#ifndef THETASPANNER_H
#define THETASPANNER_H
#include "Graph.hpp"
#include "Quadtree.hpp"

/**
 *
 * @param graph the original graph
 * @param theta the number of edges to keep per node (number of zones)
 * @return
 */
Graph create_theta_spanner_graph(Graph* graph, const int theta, bool spd_fallback=false, Quadtree* quadtree= nullptr);

#endif //THETASPANNER_H
