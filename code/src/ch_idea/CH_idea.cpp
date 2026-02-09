//
// Created by Sebastian Paule on 9/18/25.
//

#include "../io/Dataloader.hpp"
#include "../structure/Graph.hpp"
#include <fstream>
#include <omp.h>
#include "../daniel/theta-graph/headers/Progressbar.h"
#include "../daniel/theta-graph/headers/Triangulation.h"
#include "../daniel/theta-graph/headers/Timer.h"

using namespace std;

/**
 * An attempt to order nodes based on a page rank algorithm
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char* argv[]) {
    auto path = "../../data/coastlines-mercator.txt.pruned.wc.txt.graph";
    int number_of_nodes = 15355611;
    ///////////////////////////////////////////////////////////////////////////////////
    /// TEST random queries on the spanners
    /// analyse t-value
    ///////////////////////////////////////////////////////////////////////////////////

    vector<int> edges_v = vector<int>(number_of_nodes, 0);

    Triangulation triangulation;
    triangulation.readFromGraph(path);

    int num_threads = std::max(1,omp_get_num_procs()-1);  // Get the number of available processors
    //num_threads = 10; // TODO remove limit


    int max = number_of_nodes; // TODO remove limit
    int counter = 0;

    #pragma omp parallel for num_threads(num_threads) shared(counter)
    for (int source = 0; source < max; source++) {

        #pragma omp critical
        {
            if (max > 100 && counter % (max/100) == 0) {
                std::cout << "|";
                cout.flush();
            }
            counter++;
        }

        auto targets = triangulation.oneToAllVisibility(source, false);
        edges_v[source] = targets.size();
    }
    int total_edges = 0;
    for (auto e:edges_v) total_edges += e;
    cout << "Total edges in visibility graph: " << total_edges << endl;
}