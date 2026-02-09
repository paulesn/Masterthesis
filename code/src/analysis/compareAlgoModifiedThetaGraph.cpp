//
// Created by Sebastian Paule on 12/18/25.
//
/**
 * The  idea of this analysis is to compare the ways to modify the theta grapph
 *  - in parallel
 *  - from short edge to long edge
 *  - from long edge to short edge
 *
 */

#include "../io/Dataloader.hpp"
#include "MultiThetaAnalysis.h"
#include "../daniel/theta-graph/headers/Progressbar.h"
#include "../spanner/ThetaSpanner.hpp"
#include <iostream>
#include <string>
using namespace std;



int main(int argc, char* argv[]) {

    string base_graph_path;
    string spanner_path;
    double cutoff = -1.0;
    int k_in = 75;


    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "-vg" && i + 1 < argc) {
            base_graph_path = argv[++i];
            std::cout << "VisGraph graph path: " << base_graph_path << std::endl;
        }else if (arg == "-p" && i + 1 < argc) {
            cutoff = stod(argv[++i]);
            std::cout << "cutoff for worst edges to add to spanner: " << cutoff << std::endl;
        }else if (arg == "-k" && i + 1 < argc) {
            k_in  = stoi(argv[++i]);
            std::cout << "k cones for the spanner (ehem. theta): " << k_in << std::endl;
        }
    }

    if (base_graph_path.empty() || cutoff < 0) {
        std::cerr << "Usage: " << argv[0] << " -vg <string> -csv <string> -p <double> -k <int>\n";
        for (int i = 0; i < argc; i++) std::cerr << " " << argv[i];
        return 1;
    }

    ///////////////////////////////////////////////////////////////////////////////////
    /// LOAD THE GRAPHS
    ///////////////////////////////////////////////////////////////////////////////////
    auto base_graph = get<1>(load_fmi(base_graph_path));
    std::cout << "Graph loaded. Sorting Edges" << std::endl;
    base_graph.sort_edges();
    std::cout << "Graph sorted. Analysing Spanner" << std::endl;

    int k = k_in;
    std::cout << "k:p->" << k << " : " << cutoff << std::endl;


    auto spanner_graph = create_theta_spanner_graph(&base_graph, k, false);
    update_spanner_with_too_long_edges(&base_graph, &spanner_graph, cutoff, 0);
    analyse_spanner_with_vis_graph(base_graph, spanner_graph, "../data/temp.csv", "../data/temp.csv", 0.0);

    spanner_graph = create_theta_spanner_graph(&base_graph, k, false);
    update_spanner_with_too_long_edges(&base_graph, &spanner_graph, cutoff, 1);
    analyse_spanner_with_vis_graph(base_graph, spanner_graph, "../data/temp.csv", "../data/temp.csv", 0.0);

    spanner_graph = create_theta_spanner_graph(&base_graph, k, false);
    update_spanner_with_too_long_edges(&base_graph, &spanner_graph, cutoff, 2);
    analyse_spanner_with_vis_graph(base_graph, spanner_graph, "../data/temp.csv", "../data/temp.csv", 0.0);

}
