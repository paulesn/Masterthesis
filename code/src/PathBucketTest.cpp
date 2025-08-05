
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
    double early_stop = 0.01; // Early stop condition for dynamic theta update, 1.0 means no early stop
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

    int counter = 0;
    std::cout << "Inserting systems into the quadtree..." << std::endl;
    for (int p = 0; p < graph.adj.size(); p++) {
        used_points.push_back(p);
    }

    if (using_theta) {
        auto start3 = std::chrono::high_resolution_clock::now();
        spanner_theta = create_theta_spanner_graph(&graph, theta);
        //dynamic_theta_update(&graph, &spanner_theta, t_target, early_stop);
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

    constexpr int number_of_tests = 500000;
    bool t_has_inf = false;
    double theta_max_t = 0.0;
    double theta_mean_t = 0.0;
    bool theta_has_neg_t = false;
    vector<int> max_path;

    vector<int> path_buckets_original_length = {};
    vector<vector<double>> path_buckets_spanner_t = {};

    counter = 0;

    std::cout << "Checking spanning tree..." << std::endl << ">";

    // RANDOM PATH TESTING
    for (auto edge_list : graph.adj) {

        // print progress
        if (graph.adj.size() > 100 && counter % (graph.adj.size()/100) == 0) {
            std::cout << "|";
            std::cout.flush();
        }
        counter++;

        for (Edge edge : edge_list) {
            int random_source = edge.source;
            int random_target = edge.target;
            auto result_d = graph.dijkstra(random_source, random_target);
            double original_dist = result_d.second;

            int bucket = static_cast<int>(original_dist/100);
            if (bucket >= path_buckets_original_length.size()) {
                path_buckets_original_length.resize(bucket + 1, 0);
                path_buckets_spanner_t.resize(bucket + 1, {});
            }
            path_buckets_original_length[bucket]++;

            auto spanner_theta_dist = numeric_limits<double>::infinity();

            auto result_t = spanner_theta.dijkstra(random_source, random_target);
            spanner_theta_dist = result_t.second;

            double t = spanner_theta_dist/original_dist;
            path_buckets_spanner_t[bucket].emplace_back(t);

            if (spanner_theta_dist == numeric_limits<double>::infinity()) {t_has_inf = true;}

            if (spanner_theta_dist/original_dist > theta_max_t && spanner_theta_dist/original_dist < numeric_limits<double>::infinity()) {
                theta_max_t = spanner_theta_dist/original_dist;
                max_path = result_t.first;
            }

            theta_mean_t += spanner_theta_dist/original_dist;

            if (spanner_theta_dist == numeric_limits<double>::infinity()) {
                std::cout << "Test :"
                          << " (" << random_source << " -> " << random_target << "): "
                          << "Original distance: " << original_dist << ", "
                          << "Spanner Theta distance: " << spanner_theta_dist << std::endl;
            }
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
    //spanner_theta.draw();

    // PRINT BUCKET LOG into custom out stream to file
    auto out_stream = std::ofstream(out_path);
    cout << "----------------------------------------------------------------------------------" << std::endl;
    out_stream << "Path Buckets Original Length:" << std::endl;
    for (int i = 0; i < path_buckets_original_length.size(); i++) {
        out_stream << "Bucket " << i << ": " << path_buckets_original_length[i] << std::endl;
    }
    out_stream << "Path Buckets Spanner T-values:" << std::endl;
    for (int i = 0; i < path_buckets_spanner_t.size(); i++) {
        if (path_buckets_spanner_t[i].empty()) continue; // Skip empty buckets
        out_stream << "Bucket " << i << ", ";
        for (double t : path_buckets_spanner_t[i]) {
            out_stream << t << ", ";
        }
        out_stream << std::endl;
    }
}

// TODO auch die worst case edge ausrechnen

