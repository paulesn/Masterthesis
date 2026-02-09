//
// Created by Sebastian Paule on 9/4/25.
//

#include "../io/Dataloader.hpp"
#include <omp.h>

#include "MultiThetaAnalysis.h"
#include "../daniel/theta-graph/headers/Progressbar.h"
#include "../spanner/ThetaSpanner.hpp"
#include "../io/DataWriter.h"


using namespace std;



int main(int argc, char* argv[]) {

    string coastline_path;
    string csv_path;
    string temp_path = "../../data/temp.fmi";


    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "-c" && i + 1 < argc) {
            coastline_path = argv[++i];
            std::cout << "coastline path: " << coastline_path << std::endl;
        } else if (arg == "-csv" && i + 1 < argc) {
            csv_path = argv[++i];
            std::cout << "CSV output path: " << csv_path << std::endl;
        }
    }

    if (coastline_path.empty()) {
        std::cerr << "Usage: " << argv[0] << " -c <string> -csv <string>\n";
        for (int i = 0; i < argc; i++) std::cerr << " " << argv[i];
        return 1;
    }

    ///////////////////////////////////////////////////////////////////////////////////
    /// LOAD THE GRAPHS
    ///////////////////////////////////////////////////////////////////////////////////
    auto base_graph = load_coastline(coastline_path);
    auto spanner_graph = create_theta_spanner_graph(&base_graph,24, false);
    cout << "Writing coastline file to temp file as graph" << endl;
    //write_fmi(temp_path, base_graph);

    analyse_spanner(base_graph, spanner_graph, csv_path, temp_path);
}

// TODO auch die worst case edge ausrechnen

