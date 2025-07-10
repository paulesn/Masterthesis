//
// Created by sebastian on 09.07.25.
//


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
    vector<int> theta_values = {24,30,50,75,126,150,200}; // Number of zones for theta spanner
    double s = 4.5; // Separation factor for WSPD
    constexpr int number_of_tests =5000; // Number of random tests to perform
    bool using_theta = true;

    vector<double> mean_vector;
    vector<double> max_vector;
    vector<double> time_vector;
    vector<double> inf_vector;
    vector<int> edges_vector;

    string path = "../../data/0100.32.fmi";

    ///////////////////////////////////////////////////////////////////////////////////
    /// LOAD THE ORIGINAL GRAPH
    ///////////////////////////////////////////////////////////////////////////////////
    auto tup = load_fmi(path,  -1);
    auto systems = get<0>(tup);
    auto used_points = vector<int>();
    Graph graph = get<1>(tup);
    std::cout << "Loaded " << systems.size() << " points." << std::endl;

    ///////////////////////////////////////////////////////////////////////////////////
    /// INSERT SYSTEMS INTO THE QUADTREE
    ///////////////////////////////////////////////////////////////////////////////////
    auto tree = new Quadtree(32);
    int counter = 0;
    std::cout << "Inserting systems into the quadtree..." << std::endl;
    for (int p = 0; p < graph.adj.size(); p++) {
        if (tree->insert(graph.id_point_map[p])) {
            counter++;
            used_points.push_back(p);
        }
    }

    for (auto theta:theta_values){
        Graph spanner_theta = Graph(graph.n);

        long time_t =-1;



        auto all_points = tree->get_all_points();
        if (all_points.size() != counter) {
            std::cout << "Error: Quadtree contains " << all_points.size() << " points, but expected " << systems.size() << "." << std::endl;
            raise(SIGINT); // Raise SIGINT to terminate the program
        }


        std::cout << "Inserted all systems into the quadtree." << std::endl;
        ///////////////////////////////////////////////////////////////////////////////////
        /// CALCULATE WSPD (euklidian and shortest path distance)
        ///////////////////////////////////////////////////////////////////////////////////


        auto start3 = std::chrono::high_resolution_clock::now();
        spanner_theta = create_theta_spanner_graph(&graph, theta);
        auto end3 = std::chrono::high_resolution_clock::now();

        time_t = std::chrono::duration_cast<std::chrono::milliseconds>(end3 - start3).count();

        std::cout << "theta(" << theta << ") took "
              << time_t
              << " ms. And created " << spanner_theta.number_of_edges << " edges." << std::endl;
        std::cout << "For   " << systems.size() << " nodes." << std::endl;




        ///////////////////////////////////////////////////////////////////////////////////
        /// TEST random queries on the spanners
        /// analyse t-value
        ///////////////////////////////////////////////////////////////////////////////////
        double theta_max_t = 0.0;
        double theta_mean_t = 0.0;
        int t_has_inf = 0;


        // RANDOM PATH TESTING
        for (int i = 0; i < number_of_tests; i++) {
            int random_source = used_points[rand() % used_points.size()];
            int random_target = used_points[rand() % used_points.size()];
            auto original_dist = graph.dijkstra(random_source, random_target).second;
            auto spanner_theta_dist = numeric_limits<double>::infinity();
            spanner_theta_dist = spanner_theta.dijkstra(random_source, random_target).second;
            if (spanner_theta_dist == numeric_limits<double>::infinity()) {
                if (original_dist != numeric_limits<double>::infinity()) {
                    t_has_inf++;
                }
                continue;
            }

            if (spanner_theta_dist/original_dist > theta_max_t && spanner_theta_dist/original_dist < numeric_limits<double>::infinity()) {
                theta_max_t = spanner_theta_dist/original_dist;
            }

            theta_mean_t += spanner_theta_dist/original_dist;
        }

        max_vector.emplace_back(theta_max_t);
        mean_vector.emplace_back(theta_mean_t / number_of_tests);
        time_vector.emplace_back(time_t);
        edges_vector.emplace_back(spanner_theta.number_of_edges);
        inf_vector.emplace_back(t_has_inf);
    }
    ///////////////////////////////////////////////////////////////////////////////////
    /// OUTPUT
    ///////////////////////////////////////////////////////////////////////////////////

    cout << "Theta, max, mean, time, edges, inf" << endl;
    int c = 0;
    for (auto theta: theta_values) {
        std::cout << theta
                  << ", " << max_vector[c]
                  << ", " << mean_vector[c]
                  << ", " << time_vector[c]
                  << ", " << edges_vector[c]
                  << ", " << inf_vector[c]
                  << std::endl;
        c++;
    }
}
