//
// Created by Sebastian Paule on 11/4/25.
//
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

    string coastline_path;
    string spanner_path;
    string csv_path;
    double cutoff = 1.1;
    int k = 75;


    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "-sg" && i + 1 < argc) {
            spanner_path = argv[++i];
            std::cout << "Spanner path: " << spanner_path << std::endl;
        } else if (arg == "-cg" && i + 1 < argc) {
            coastline_path = argv[++i];
            std::cout << "Coastline graph path: " << coastline_path << std::endl;
        }else if (arg == "-c" && i + 1 < argc) {
            cutoff = stod(argv[++i]);
            std::cout << "percent of worst edges to add to spanner: " << cutoff << std::endl;
        }else if (arg == "-k" && i + 1 < argc) {
            k = stoi(argv[++i]);
            std::cout << "k cones for the spanner (ehem. theta): " << k << std::endl;
        }
    }

    if (spanner_path.empty() || coastline_path.empty() || cutoff < 0) {
        std::cerr << "Usage: " << argv[0] << " -sg <string> -vg <string> -csv <string> -p <double> -k <int>\n";
        for (int i = 0; i < argc; i++) std::cerr << " " << argv[i];
        return 1;
    }

    ///////////////////////////////////////////////////////////////////////////////////
    /// LOAD THE GRAPHS
    ///////////////////////////////////////////////////////////////////////////////////
    auto spanner_graph = get<1>(load_fmi(spanner_path,  -1));

    analyse_spanner_with_coastline_graph(coastline_path, spanner_graph, cutoff);


    /*vector<Edge> edges_to_add = analyse_spanner_with_vis_graph(base_graph, spanner_graph, csv_path+"-"+to_string(theta)+".csv", "../../data/all_edges-"+to_string(theta)+".csv", percent);
    for (Edge edge: edges_to_add) {
        spanner_graph.addEdge(edge.source, edge.target, edge.weight, true);
    }
    cout << "Added " << edges_to_add.size() << " edges to the spanner." << endl;
    string new_csv_path = csv_path.substr(0, csv_path.find_last_of('.')) + "_with_worst_edges.csv";
    */
    //analyse_spanner_with_vis_graph(base_graph, spanner_graph, new_csv_path, "../../data/all_edges-"+to_string(theta)+"-"+to_string(percent)+".csv", 0.0);
}


