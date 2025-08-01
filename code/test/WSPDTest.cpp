
#include "../src/Dataloader.hpp"
#include "../src/Quadtree.hpp"
#include "../src/ThetaSpanner.hpp"
#include "../src/DataWriter.h"
#include <ostream>
#include <iostream>
#include <set>
#include <chrono>
#include <csignal>

#include "../lib/json.hpp"


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
    double t_target = 1.1;
    vector<double> early_stops = {0.01, 0.05, 0.1, 0.25, 0.4, 0.5}; // Early stop condition for dynamic theta update, 1.0 means no early stop
    bool using_wspd_e = true;
    bool using_wspd_spd = false;
    bool using_theta = true;

    string path = "../../data/0200.32.fmi";
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


    ///////////////////////////////////////////////////////////////////////////////////
    /// INSERT SYSTEMS INTO THE QUADTREE
    ///////////////////////////////////////////////////////////////////////////////////
    auto tree = new Quadtree(32);
    int counter = 0;
    std::cout << "Inserting systems into the quadtree..." << std::endl;
    for (int p = 0; p < graph.adj.size(); p++) {
        used_points.push_back(p);
    }

    auto all_points = tree->get_all_points();
    if (all_points.size() != counter) {
        std::cout << "Error: Quadtree contains " << all_points.size() << " points, but expected " << systems.size() << "." << std::endl;
        raise(SIGINT); // Raise SIGINT to terminate the program
    }

    for (auto early_stop : early_stops){
        if (using_theta) {
            auto start3 = std::chrono::high_resolution_clock::now();
            spanner_theta = create_theta_spanner_graph(&graph, theta);
            dynamic_theta_update(&graph, &spanner_theta, t_target, early_stop);
            auto end3 = std::chrono::high_resolution_clock::now();

            time_t = std::chrono::duration_cast<std::chrono::milliseconds>(end3 - start3).count();

            std::cout << "theta("<< theta << ") took "
                  << time_t
                  << " ms. And created " << spanner_theta.number_of_edges << " edges." << std::endl;
            std::cout << "For   " << systems.size() << " nodes." << std::endl;
        }



        ///////////////////////////////////////////////////////////////////////////////////
        /// TEST random queries on the spanners
        /// analyse t-value
        ///////////////////////////////////////////////////////////////////////////////////

        constexpr int number_of_tests = 5000;
        bool t_has_inf = false;
        double theta_max_t = 0.0;
        double theta_mean_t = 0.0;
        bool theta_has_neg_t = false;


        // RANDOM PATH TESTING
        for (int i = 0; i < number_of_tests; i++) {
            int random_source = used_points[rand() % used_points.size()];
            int num_targets = graph.adj[random_source].size();
            if (num_targets == 0) {
                //std::cout << "No targets for source " << random_source << ". Skipping test." << std::endl;
                i-=1; // Decrement i to repeat this iteration
                continue; // Skip if no targets are available
            }
            int random_target = graph.adj[random_source][rand() % num_targets].target;
            auto original_dist = graph.dijkstra(random_source, random_target).second;

            auto spanner_theta_dist = numeric_limits<double>::infinity();
            if (using_theta) spanner_theta_dist = spanner_theta.dijkstra(random_source, random_target).second;

            if (spanner_theta_dist == numeric_limits<double>::infinity()) {t_has_inf = true;}

            if (spanner_theta_dist/original_dist > theta_max_t && spanner_theta_dist/original_dist < numeric_limits<double>::infinity()) {theta_max_t = spanner_theta_dist/original_dist;}

            theta_mean_t += spanner_theta_dist/original_dist;

            if (spanner_theta_dist == numeric_limits<double>::infinity()) {
                std::cout << "Test " << i + 1
                          << " (" << random_source << " -> " << random_target << "): "
                          << "Original distance: " << original_dist << ", "
                          << "Spanner Theta distance: " << spanner_theta_dist << std::endl;
            }
        }

        int noe = 0;
        for (int i = 0; i < graph.adj.size(); i++) {
            noe += spanner_theta.adj[i].size();
        }

        cout << std::endl;

        cout << "----------------------------------------------------------------------------------" << std::endl;
        std::
        cout << std::endl << "Results after " << number_of_tests << " tests:" << std::endl;
        cout << "----------------------------------------------------------------------------------" << std::endl;
        std::cout << "Spanner Theta has " << (using_theta ? "" : "not ") << "been used." << std::endl;
        std::cout << "Spanner Theta has " << noe << " Edges." << std::endl;
        std::cout << "Spanner Theta has been calculated in " << time_t << " ms." << std::endl;
        std::cout << "Spanner Theta has infinity: " << (t_has_inf ? "Yes" : "No") << std::endl;
        std::cout << "Spanner Theta has negative t-value: " << (theta_has_neg_t ? "Yes" : "No") << std::endl;
        std::cout << "Spanner Theta max t-value: " << theta_max_t << std::endl;
        std::cout << "Spanner Theta mean t-value: " << (theta_mean_t / number_of_tests) << std::endl;
        std::cout << "All tests completed." << std::endl;

        spanner_theta.store_to_disk("../../data/theta_spanner.fmi");
        write_gf(spanner_theta, out_path);
        //spanner_theta.draw();
    }
}

// TODO auch die worst case edge ausrechnen

