#include <cassert>
#include <chrono>
#include <iostream>
#include <limits>

#include "../src/Graph.hpp"
#include "../src/Dataloader.hpp"

using namespace std;

void testGraph() {

    constexpr int number_of_tests = 400;
    string path = "../../data/0200.32.fmi";

    ///////////////////////////////////////////////////////////////////////////////////
    /// LOAD THE ORIGINAL GRAPH
    ///////////////////////////////////////////////////////////////////////////////////
    auto tup = load_fmi(path,  -1);
    auto systems = get<0>(tup);
    auto used_points = vector<int>();
    Graph graph = get<1>(tup);
    graph.sort_adj();

    // add all points to used_points
    for (int p = 0; p < graph.adj.size(); p++) {
        used_points.push_back(p);
    }


    ///////////////////////////////////////////////////////////////////////////////////
    /// TEST random queries on the spanners
    /// analyse t-value
    ///////////////////////////////////////////////////////////////////////////////////

    // RANDOM PATH TESTING
    double time_sum_dijkstra = 0;
    double time_sum_spira = 0;
    int skipped_tests = 0;

    for (int i = 0; i < number_of_tests; i++) {
        int random_source = used_points[rand() % used_points.size()];
        int num_targets = graph.adj[random_source].size();
        if (num_targets == 0) {
            std::cout << "No targets for source " << random_source << ". Skipping test." << std::endl;
            skipped_tests++;
            continue; // Skip if no targets are available
        }
        int random_target = used_points[rand() % used_points.size()];

        auto start_dijkstra = std::chrono::high_resolution_clock::now();
        auto dikstra_dist = graph.dijkstra(random_source, random_target).second;
        auto end_dijkstra = std::chrono::high_resolution_clock::now();

        auto start_spira = std::chrono::high_resolution_clock::now();
        auto spira_dist = graph.spira_sp(random_source, random_target);
        auto end_spira = std::chrono::high_resolution_clock::now();

        if (dikstra_dist != spira_dist) {
            cout << "Distane ERROR !!!! <--------------------------------------------" << endl;
            cout << "Source: " << random_source << ", Target: " << random_target << endl;
            cout << "Dijkstra: " << dikstra_dist << endl;
            cout << "Spira: " << spira_dist << endl;
        }


        time_sum_dijkstra += std::chrono::duration_cast<std::chrono::milliseconds>(end_dijkstra - start_dijkstra).count();
        time_sum_spira += std::chrono::duration_cast<std::chrono::milliseconds>(end_spira - start_spira).count();

    }
    cout << "average Dijkstra time: " << time_sum_dijkstra/(number_of_tests-skipped_tests) << endl;
    cout << "average Spira time: " << time_sum_spira/(number_of_tests-skipped_tests) << endl;
}



int main(){
    testGraph();
    return 0
    ;
}