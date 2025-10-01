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
    int paths = -1;
    int theta = 24;


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
            paths = stoi(argv[++i]);
            std::cout << "number of paths to check: " << paths << std::endl;
        }else if (arg == "-t" && i + 1 < argc) {
            theta = stoi(argv[++i]);
            std::cout << "theta for spanner: " << theta << std::endl;
        }
    }

    if (spanner_path.empty() || base_graph_path.empty() || paths < 0) {
        std::cerr << "Usage: " << argv[0] << " -sg <string> -vg <string> -csv <string> -p <double>\n";
        for (int i = 0; i < argc; i++) std::cerr << " " << argv[i];
        return 1;
    }

    ///////////////////////////////////////////////////////////////////////////////////
    /// LOAD THE GRAPHS
    ///////////////////////////////////////////////////////////////////////////////////
    auto base_graph = get<1>(load_fmi(base_graph_path, -1));
    auto spanner_graph = create_theta_spanner_graph(&base_graph, theta); //get<1>(load_fmi(spanner_path,  -1));

    analyse_random_paths_with_vis_graph(base_graph, spanner_graph, csv_path, paths, theta);

}


