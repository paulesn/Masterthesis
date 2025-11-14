//
// Created by Sebastian Paule on 9/5/25.
//

#ifndef MULTITHETAANALYSIS_H
#define MULTITHETAANALYSIS_H

double normalize(double value);

/**
 *
 * @param base_graph a graph that contains all nodes but only the polygon edges, not the visibility edges
 * @param spanner_graph the spanner graph to be analysed
 * @param csv_out_path  the path to the csv file where the histogramm of t-values will be stored
 * @return
 */
int analyse_spanner(Graph base_graph, Graph spanner_graph, std::string csv_out_path, std::string graph_path, int max=-1, bool multi_threading=true);
std::vector<Edge> analyse_spanner_with_vis_graph(Graph base_graph, Graph spanner_graph, std::string csv_out_path, std::string all_edges_path, double percent);

/**
 *
 * @param base_graph the orignial graph the spanner is based on
 * @param spanner_graph the spanner
 * @param csv_out_path the place to store analyzed data that is to big for the terminal
 * @param max the number of random paths to check
 * @param theta the theta value of the spanner for the file paths
 * @return
 */
void analyse_random_paths_with_vis_graph(Graph base_graph, Graph spanner_graph, std::string csv_out_path, int max, int theta);

/**
 * ensures that the nodes are at least > 1/2 the furthest distance in the graph away from each other
 * @param base_graph the orignial graph the spanner is based on
 * @param spanner_graph the spanner
 * @param csv_out_path the place to store analyzed data that is to big for the terminal
 * @param max the number of random paths to check
 * @param theta the theta value of the spanner for the file paths
 * @param plimit the percent of the furthest distance two nodes have to be away from each other to be considered
 * @return
 */
void analyse_long_random_paths_with_vis_graph(Graph base_graph, Graph spanner_graph, std::string csv_out_path, int max, int theta, double plimit=0.8);

/**
 * this function analyses all edges in the base_graph_coastline. it is very inefficiant as it checks all edges.
 * @param base_graph_coastline the coastline graph, the spanner is based on
 * @param spanner_graph the spanner
 * @param csv_out_path the place to store analyzed data that is to big for the terminal
 * @param all_edges_path the place to store all edges with their t-values
 * @param t_cutoff the cutoff to consider an edge as "worst" to return it
 * @return all edges with a t-value higher than t_cutoff
 */
std::vector<Edge> analyse_spanner_with_coastline_graph(std::string base_graph_coastline_path, Graph spanner_graph, double t_cutoff);


void dijkstra_debugging(Graph graph);
#endif //MULTITHETAANALYSIS_H


