//
// Created by Sebastian Paule on 9/4/25.
//

#include "../src/Dataloader.hpp"
#include <omp.h>

#include "MultiThetaAnalysis.h"
#include "Progressbar.h"


using namespace std;



int main(int argc, char* argv[]) {

    string base_graph_path;
    string spanner_path;
    string csv_path;


    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "-sg" && i + 1 < argc) {
            spanner_path = argv[++i];
            std::cout << "Spanner path: " << spanner_path << std::endl;
        } else if (arg == "-bg" && i + 1 < argc) {
            base_graph_path = argv[++i];
            std::cout << "Base graph path: " << base_graph_path << std::endl;
        } else if (arg == "-csv" && i + 1 < argc) {
            csv_path = argv[++i];
            std::cout << "CSV output path: " << csv_path << std::endl;
        }
    }

    if (spanner_path.empty() || base_graph_path.empty()) {
        std::cerr << "Usage: " << argv[0] << " -sg <string> -bg <string> -csv <string>\n";
        for (int i = 0; i < argc; i++) std::cerr << " " << argv[i];
        return 1;
    }

    ///////////////////////////////////////////////////////////////////////////////////
    /// LOAD THE GRAPHS
    ///////////////////////////////////////////////////////////////////////////////////
    auto base_graph = get<1>(load_fmi(base_graph_path,  -1));
    auto spanner_graph = get<1>(load_fmi(spanner_path,  -1));

    analyse_spanner(base_graph, spanner_graph, csv_path);
}

// TODO auch die worst case edge ausrechnen

