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
std::vector<Edge> analyse_spanner_with_vis_graph(Graph base_graph, Graph spanner_graph, std::string csv_out_path, std::string graph_path, double percent);

#endif //MULTITHETAANALYSIS_H
