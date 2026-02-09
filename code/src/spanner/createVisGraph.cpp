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

#include "../daniel/theta-graph/headers/Progressbar.h"
#include "../spanner/ThetaSpanner.hpp"
#include "../daniel/theta-graph/headers/Triangulation.h"
#include "../io/DataWriter.h"


using namespace std;


int main(int argc, char* argv[]) {

    string coastline_path;
    string output_path;


    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "-o" && i + 1 < argc) {
            output_path = argv[++i];
            std::cout << "Output path: " << output_path << std::endl;
        } else if (arg == "-c" && i + 1 < argc) {
            coastline_path = argv[++i];
            std::cout << "Coastline graph path: " << coastline_path << std::endl;
        }
    }

    if (output_path.empty() || coastline_path.empty()) {
        std::cerr << "Usage: " << argv[0] << " -c <string> -o <string> \n";
        for (int i = 0; i < argc; i++) std::cerr << " " << argv[i];
        return 1;
    }

    ///////////////////////////////////////////////////////////////////////////////////
    /// LOAD THE GRAPHS
    ///////////////////////////////////////////////////////////////////////////////////
    auto coastline_graph = load_coastline(coastline_path);

    Triangulation tri = Triangulation();
    tri.readFromGraph(coastline_path);
    vector<std::vector<GlobalID>> visEdges;
    vector<vector<double>> visCosts;
    tri.createVisibilityGraph(visEdges, visCosts, false);
    tri.saveVisibilityGraph(visEdges, visCosts, output_path);
    ///////////////////////////////////////////////////////////////////////////////////
    ///
    ////////////////////////////////////////////////////////////////////////////////

    }


