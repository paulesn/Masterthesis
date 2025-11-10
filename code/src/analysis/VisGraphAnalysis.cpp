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
    int k_in = 75;


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
        }else if (arg == "-k" && i + 1 < argc) {
            k_in  = stoi(argv[++i]);
            std::cout << "k cones for the spanner (ehem. theta): " << k_in << std::endl;
        }
    }

    if (base_graph_path.empty() || percent < 0) {
        std::cerr << "Usage: " << argv[0] << " -sg <string> -vg <string> -csv <string> -p <double> -k <int>\n";
        for (int i = 0; i < argc; i++) std::cerr << " " << argv[i];
        return 1;
    }

    ///////////////////////////////////////////////////////////////////////////////////
    /// LOAD THE GRAPHS
    ///////////////////////////////////////////////////////////////////////////////////
    auto base_graph = get<1>(load_fmi(base_graph_path));
    //auto spanner_graph = create_theta_spanner_graph(&base_graph, theta); //get<1>(load_fmi(spanner_path,  -1));
    std::cout << "Graph loaded. Sorting Edges" << std::endl;
    base_graph.sort_edges();
    std::cout << "Graph sorted. Analysing Spanner" << std::endl;
    for (int k :{24,32,50,64,75,100,128}){
        for (double p :{1.0, 0.9, 0.8, 0.7, 0.6, 0.5}){
            std::cout << "k:p->" << k << " : " << p << std::endl;
            auto spanner_graph = create_theta_spanner_graph_with_max_angle(&base_graph, k, ((2*M_PI)/k)*p);

            /*vector<Edge> edges_to_add = analyse_spanner_with_vis_graph(base_graph, spanner_graph, csv_path+"-"+to_string(theta)+".csv", "../../data/all_edges-"+to_string(theta)+".csv", percent);
            for (Edge edge: edges_to_add) {
                spanner_graph.addEdge(edge.source, edge.target, edge.weight, true);
            }
            cout << "Added " << edges_to_add.size() << " edges to the spanner." << endl;
            string new_csv_path = csv_path.substr(0, csv_path.find_last_of('.')) + "_with_worst_edges.csv";
            */
            analyse_spanner_with_vis_graph(base_graph, spanner_graph, csv_path, "../data/temp.csv", 0.0);
            analyse_random_paths_with_vis_graph(base_graph, spanner_graph, csv_path, 10000, k);
            analyse_long_random_paths_with_vis_graph(base_graph, spanner_graph, csv_path, 10000, k, 0.9);
}}}


