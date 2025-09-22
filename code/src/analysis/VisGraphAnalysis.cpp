//
// Created by Sebastian Paule on 9/15/25.
//
/**
* This analysis file starts a analysis, then adds the X% worst edges and analyses it again
*/


//
// Created by Sebastian Paule on 9/4/25.
//

#include "../io/Dataloader.hpp"
#include <omp.h>

#include "MultiThetaAnalysis.h"
#include "../daniel/theta-graph/headers/Progressbar.h"
#include "../spanner/ThetaSpanner.hpp"


using namespace std;



int main(int argc, char* argv[]) {

    string base_graph_path;
    string spanner_path;
    string csv_path;
    double percent = -1.0;


    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "-sg" && i + 1 < argc) {
            spanner_path = argv[++i];
            std::cout << "Spanner path: " << spanner_path << std::endl;
        } else if (arg == "-vg" && i + 1 < argc) {
            base_graph_path = argv[++i];
            std::cout << "VisGraph graph path: " << base_graph_path << std::endl;
        } else if (arg == "-csv" && i + 1 < argc) {
            csv_path = argv[++i];
            std::cout << "CSV output path: " << csv_path << std::endl;
        }else if (arg == "-p" && i + 1 < argc) {
            percent = stod(argv[++i]);
            std::cout << "percent of worst edges to add to spanner: " << percent << std::endl;
        }
    }

    if (spanner_path.empty() || base_graph_path.empty() || percent < 0) {
        std::cerr << "Usage: " << argv[0] << " -sg <string> -vg <string> -csv <string> -p <double>\n";
        for (int i = 0; i < argc; i++) std::cerr << " " << argv[i];
        return 1;
    }

    ///////////////////////////////////////////////////////////////////////////////////
    /// LOAD THE GRAPHS
    ///////////////////////////////////////////////////////////////////////////////////
    auto base_graph = get<1>(load_fmi(base_graph_path, -1));
    auto spanner_graph = create_theta_spanner_graph(&base_graph, 24); //get<1>(load_fmi(spanner_path,  -1));

    vector<Edge> edges_to_add = analyse_spanner_with_vis_graph(base_graph, spanner_graph, csv_path, base_graph_path, percent);
    for (Edge edge: edges_to_add) {
        spanner_graph.addEdge(edge.source, edge.target, edge.weight, true);
    }
    cout << "Added " << edges_to_add.size() << " edges to the spanner." << endl;
    string new_csv_path = csv_path.substr(0, csv_path.find_last_of('.')) + "_with_worst_edges.csv";
    analyse_spanner_with_vis_graph(base_graph, spanner_graph, new_csv_path, base_graph_path, 0.0);
}


