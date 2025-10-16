//
// Created by sebastian on 13.06.25.
//

#ifndef THETASPANNER_H
#define THETASPANNER_H
#include "../structure/Graph.hpp"
#include "../structure/Quadtree.hpp"

/**
 * This function creates a theta graph based on a input graph. The resulting graph only contains edges in the oroiginal graph
 * @param graph the original graph the theta spanner should be based on
 * @param theta the number of zones around each node for edges
 * @param early_break if true, the algorithm stops adding edges for a node if all zones are filled only use if edges are sorted by length
 * @return a new graph
 */
Graph create_theta_spanner_graph(Graph* graph, const int theta, bool early_break = false);


void dynamic_theta_update(Graph* graph, Graph* spanner, double t);

/**
 * This function calculates the relative angle between the x axis and the vector from point a to point b.
 * @param a
 * @param b
 * @return the angle in radians between the x-axis and the vector from a to b.
 */
double angle_between(Pointc a, Pointc b);

#endif //THETASPANNER_H
