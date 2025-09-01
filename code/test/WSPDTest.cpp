
#include "../src/Dataloader.hpp"
#include "../src/Quadtree.hpp"
#include "../src/ThetaSpanner.hpp"
#include "../src/DataWriter.h"
#include <ostream>
#include <iostream>
#include <set>
#include <chrono>
#include <csignal>
#include <omp.h>


using namespace std;

/**
 * floating point normalization to prevent floating point error missmatches
 * @param value
 * @return
 */
double normalize(double value) {
    return std::round(value * 10000.0) / 10000.0;
}



int main() {
    int theta = 24; // Number of zones for theta spanner
    double s = 4.5; // Separation factor for WSPD
    bool using_theta = true;

    string path = "../../data/0025.32.fmi";
    string out_path = "../../data/theta_spanner.gl";

    ///////////////////////////////////////////////////////////////////////////////////
    /// LOAD THE ORIGINAL GRAPH
    ///////////////////////////////////////////////////////////////////////////////////
    auto tup = load_fmi(path,  -1);
    auto systems = get<0>(tup);
    auto used_points = vector<int>();
    Graph graph = get<1>(tup);
    // init the hub labels for faster shortest path distance calculation
    //graph.init_hub_labels();
    std::cout << "Loaded " << systems.size() << " points." << std::endl;

    Graph spanner_theta = Graph(graph.n);

    long time_t =-1;





    auto start3 = std::chrono::high_resolution_clock::now();
    spanner_theta = create_theta_spanner_graph(&graph, theta);
    std::cout << "Created theta spanner graph with " << spanner_theta.number_of_edges << " edges." << std::endl;
    //dynamic_theta_update(&graph, &spanner_theta, 1.1);
    std::cout << "Updated theta spanner graph with 1.1 zones. and has " << spanner_theta.number_of_edges << std::endl;
    auto end3 = std::chrono::high_resolution_clock::now();

    time_t = std::chrono::duration_cast<std::chrono::milliseconds>(end3 - start3).count();

    std::cout << "theta("<< theta << ") took "
          << time_t
          << " ms. And created " << spanner_theta.number_of_edges << " edges." << std::endl;
    std::cout << "For   " << systems.size() << " nodes." << std::endl;



    ///////////////////////////////////////////////////////////////////////////////////
    /// TEST random queries on the spanners
    /// analyse t-value
    ///////////////////////////////////////////////////////////////////////////////////

    constexpr int number_of_tests = 50000;
    bool t_has_inf = false;
    double theta_max_t = 0.0;
    double theta_mean_t = 0.0;
    bool theta_has_neg_t = false;
    vector <int> max_path;


    // RANDOM PATH TESTING
    int num_threads = std::max(1,omp_get_num_procs()-1);  // Get the number of available processors

    #pragma omp parallel for num_threads(num_threads) shared(spanner_theta, t_has_inf, theta_max_t, theta_mean_t, theta_has_neg_t, max_path) schedule(dynamic)
    for (int i = 0; i < number_of_tests; i++) {
        int random_source = rand() % spanner_theta.adj.size();
        int num_targets = graph.adj[random_source].size();
        if (num_targets == 0) {
            std::cout << "No targets for source " << random_source << ". Skipping test." << std::endl;
            i--;
            continue; // Skip if no targets are available
        }
        int random_target = graph.adj[random_source][rand() % num_targets].target;
        auto original_dist = graph.dijkstra(random_source, random_target).second;
        auto spanner_theta_dist = numeric_limits<double>::infinity();
        auto s_dist = spanner_theta.dijkstra(random_source, random_target);

        #pragma omp critical
        {
            if (s_dist.second < numeric_limits<double>::infinity()) {
                spanner_theta_dist = s_dist.second;
            }
            if (spanner_theta_dist == numeric_limits<double>::infinity()) {t_has_inf = true;}
            if (normalize(spanner_theta_dist) < normalize(original_dist)) {theta_has_neg_t = true;}
            if (spanner_theta_dist/original_dist > theta_max_t && spanner_theta_dist/original_dist < numeric_limits<double>::infinity()) {
                theta_max_t = spanner_theta_dist/original_dist;
                max_path = s_dist.first;
            }
            theta_mean_t += spanner_theta_dist/original_dist;
        }
    }

    int c2 = 0;
    double spanner_theta_max_t = 0.0;
    cout << ">";
    // ALL EDGES TEST


    cout << std::endl;

    cout << "----------------------------------------------------------------------------------" << std::endl;
    std::
    cout << std::endl << "Results after " << number_of_tests << " tests:" << std::endl;
    cout << "----------------------------------------------------------------------------------" << std::endl;
    std::cout << "Spanner Theta has " << (using_theta ? "" : "not ") << "been used." << std::endl;
    std::cout << "Spanner Theta has been calculated in " << time_t << " ms." << std::endl;
    std::cout << "Spanner Theta has infinity: " << (t_has_inf ? "Yes" : "No") << std::endl;
    std::cout << "Spanner Theta has negative t-value: " << (theta_has_neg_t ? "Yes" : "No") << std::endl;
    std::cout << "Spanner Theta max t-value: " << theta_max_t << std::endl;
    std::cout << "Spanner Theta mean t-value: " << (theta_mean_t / number_of_tests) << std::endl;
    std::cout << "All tests completed." << std::endl;

    vector<Point> path_to_draw;
    for (auto &id : max_path) {
        path_to_draw.emplace_back(graph.id_point_map[id].x, graph.id_point_map[id].y);
    }

    write_gf(out_path, theta, path_to_draw);
    write_flat(out_path + ".flat", theta, path_to_draw);
}

// TODO auch die worst case edge ausrechnen

