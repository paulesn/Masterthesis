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

/**
 * this variant of the theta spanner creation always adds the shortest edge in each cone.
 * If the edge has a to large angle to any of the cone borders, the edge and cone border are another cone.

 * @param graph
 * @param k
 * @param early_break
 * @param max_an the maximum angle between the edge and the furthest cone border
 * @return
 */
Graph create_theta_spanner_graph_with_max_angle(Graph* graph, const int k, double max_an);

/**
 * this function iterates through all adjacent angles on each node and for each angle between edge uv and uw it adds the edge vw if the angle is to small
 * this version uses a triangulation to check for visibility
 * @param base_graph_path the path to the coastlines file for the triangulation
 * @param spanner
 * @param min_an
 * @return
 */
void enforce_small_angle_constraint(std::string base_graph_path, Graph* spanner, double min_an);

/**
 * this function iterates through all adjacent angles on each node and for each angle between edge uv and uw it adds the edge vw if the angle is to small
 * this variant uses a vis graph to check for visibility
 * @param base_graph the base graph containing all visibility edges
 * @param spanner
 * @param min_an
 * @return
 */
void enforce_small_angle_constraint(Graph* base_graph, Graph* spanner, double min_an);

#endif //THETASPANNER_H
