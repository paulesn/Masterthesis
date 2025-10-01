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
int analyse_spanner(Graph base_graph, Graph spanner_graph, std::string csv_out_path, std::string graph_path);
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

#endif //MULTITHETAANALYSIS_H


