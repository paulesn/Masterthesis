
#include "../src/Dataloader.hpp"
#include "../src/Quadtree.hpp"
#include "../src/ThetaSpanner.hpp"
#include "../src/HubLabels.hpp"
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
    int theta = 8; // Number of zones
    double s = 2;
    bool using_wspd_e = true;
    bool using_wspd_spd = true;
    bool using_theta = true;

    string path = "../../data/mini-ch.fmi";

    ///////////////////////////////////////////////////////////////////////////////////
    /// LOAD THE ORIGINAL GRAPH
    ///////////////////////////////////////////////////////////////////////////////////
    auto tup = load_fmi(path,  -1);
    auto systems = get<0>(tup);
    auto used_points = vector<int>();
    Graph graph = get<1>(tup);
    // init the hub labels for faster shortest path distance calculation
    graph.init_hub_labels();
    std::cout << "Loaded " << systems.size() << " points." << std::endl;

    Graph spanner_e = Graph(graph.n);
    Graph spanner_sp = Graph(graph.n);
    Graph spanner_theta = Graph(graph.n);

    long time_t =-1;
    long time_sp =-1;
    long time_e =-1;



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

    auto all_points = tree->get_all_points();
    if (all_points.size() != counter) {
        std::cout << "Error: Quadtree contains " << all_points.size() << " points, but expected " << systems.size() << "." << std::endl;
        raise(SIGINT); // Raise SIGINT to terminate the program
    }


    std::cout << "Inserted all systems into the quadtree." << std::endl;
    ///////////////////////////////////////////////////////////////////////////////////
    /// CALCULATE WSPD (euklidian and shortest path distance)
    ///////////////////////////////////////////////////////////////////////////////////

    if (using_wspd_e){
        auto start1 = std::chrono::high_resolution_clock::now();
        spanner_e = wspd(tree, s, &graph);
        auto end1 = std::chrono::high_resolution_clock::now();

        time_e = std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1).count();
        std::cout << "wspd_e("<< s << ") took "
              << time_e
              << " ms. And created "<< spanner_e.number_of_edges
              << " edges." << std::endl;
    }

    if (using_wspd_spd) {
        auto start2 = std::chrono::high_resolution_clock::now();
        spanner_sp = wspd_spd(tree, s, graph);
        auto end2 = std::chrono::high_resolution_clock::now();

        time_sp = std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count();

        std::cout << "wspd_sp("<< s << ") took "
              << time_sp
              << " ms. And created " << spanner_sp.number_of_edges << " edges." << std::endl;
        std::cout << "For   " << systems.size() << " nodes." << std::endl;
    }

    if (using_theta) {
        auto start3 = std::chrono::high_resolution_clock::now();
        spanner_theta = create_theta_spanner_graph(&graph, theta);
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

    constexpr int number_of_tests = 50;
    bool e_has_inf = false;
    bool sp_has_inf = false;
    bool t_has_inf = false;
    double e_max_t = 0.0;
    double sp_max_t = 0.0;
    double theta_max_t = 0.0;
    double e_mean_t = 0.0;
    double sp_mean_t = 0.0;
    double theta_mean_t = 0.0;
    bool e_has_neg_t = false;
    bool sp_has_neg_t = false;
    bool theta_has_neg_t = false;


    for (int i = 0; i < number_of_tests; i++) {
        int random_source = used_points[rand() % used_points.size()];
        int random_target = used_points[rand() % used_points.size()];
        auto original_dist = graph.dijkstra(random_source, random_target).second;

        auto spanner_sp_dist = numeric_limits<double>::infinity();
        if (using_wspd_spd) spanner_sp_dist = spanner_sp.dijkstra(random_source, random_target).second;

        auto spanner_e_dist = numeric_limits<double>::infinity();
        if (using_wspd_e) spanner_e_dist = spanner_e.dijkstra(random_source, random_target).second;

        auto spanner_theta_dist = numeric_limits<double>::infinity();
        if (using_theta) spanner_theta_dist = spanner_theta.dijkstra(random_source, random_target).second;

        if (spanner_e_dist == numeric_limits<double>::infinity()) {e_has_inf = true;}
        if (spanner_sp_dist == numeric_limits<double>::infinity()) {sp_has_inf = true;}
        if (spanner_theta_dist == numeric_limits<double>::infinity()) {t_has_inf = true;}
        if (normalize(spanner_e_dist) < normalize(original_dist)) {
            e_has_neg_t = true;
        }
        if (normalize(spanner_sp_dist) < normalize(original_dist)) {
            sp_has_neg_t = true;
        }
        if (normalize(spanner_theta_dist) < normalize(original_dist)) {theta_has_neg_t = true;}
        if (spanner_e_dist/original_dist > e_max_t && spanner_e_dist/original_dist < numeric_limits<double>::infinity()) {e_max_t = spanner_e_dist/original_dist;}
        if (spanner_sp_dist/original_dist > sp_max_t && spanner_sp_dist/original_dist < numeric_limits<double>::infinity()) {
            sp_max_t = spanner_sp_dist/original_dist;
        }
        if (spanner_theta_dist/original_dist > theta_max_t && spanner_theta_dist/original_dist < numeric_limits<double>::infinity()) {theta_max_t = spanner_theta_dist/original_dist;}

        sp_mean_t += spanner_sp_dist/original_dist;
        e_mean_t += spanner_e_dist/original_dist;
        theta_mean_t += spanner_theta_dist/original_dist;

        std::cout << "Test " << i + 1
                  << " (" << random_source << " -> " << random_target << "): "
                  << "Original distance: " << original_dist << ", "
                  << "Spanner SPD distance: " << spanner_sp_dist << ", "
                  << "Spanner E distance: " << spanner_e_dist << std::endl
                  << "Spanner Theta distance: " << spanner_theta_dist << std::endl;

    }
    cout << "----------------------------------------------------------------------------------" << std::endl;
    std::cout << std::endl << "Results after " << number_of_tests << " tests:" << std::endl;
    cout << "----------------------------------------------------------------------------------" << std::endl;
    std::cout << "Spanner SPD has " << (using_wspd_spd ? "" : "not ") << "been used." << std::endl;
    std::cout << "Spanner SPD has " << spanner_sp.number_of_edges << " Edges." << std::endl;
    std::cout << "Spanner SPD has been calculated in " << time_sp << " ms." << std::endl;
    std::cout << "Spanner SPD has infinity: " << (sp_has_inf ? "Yes" : "No") << std::endl;
    std::cout << "Spanner SPD has negative t-value: " << (sp_has_neg_t ? "Yes" : "No") << std::endl;
    std::cout << "Spanner SPD max t-value: " << sp_max_t << std::endl;
    std::cout << "Spanner SPD mean t-value: " << (sp_mean_t / number_of_tests) << std::endl << std::endl;
    std::cout << "Spanner E has " << (using_wspd_e ? "" : "not ") << "been used." << std::endl;
    std::cout << "Spanner E has " << spanner_e.number_of_edges << " Edges." << std::endl;
    std::cout << "Spanner E has been calculated in " << time_e << " ms." << std::endl;
    std::cout << "Spanner E has infinity: " << (e_has_inf ? "Yes" : "No") << std::endl;
    std::cout << "Spanner E has negative t-value: " << (e_has_neg_t ? "Yes" : "No") << std::endl;
    std::cout << "Spanner E max t-value: " << e_max_t << std::endl;
    std::cout << "Spanner E mean t-value: " << (e_mean_t/ number_of_tests)<< std::endl << std::endl;
    std::cout << "Spanner Theta has " << (using_theta ? "" : "not ") << "been used." << std::endl;
    std::cout << "Spanner Theta has " << spanner_theta.number_of_edges << " Edges." << std::endl;
    std::cout << "Spanner Theta has been calculated in " << time_t << " ms." << std::endl;
    std::cout << "Spanner Theta has infinity: " << (t_has_inf ? "Yes" : "No") << std::endl;
    std::cout << "Spanner Theta has negative t-value: " << (theta_has_neg_t ? "Yes" : "No") << std::endl;
    std::cout << "Spanner Theta max t-value: " << theta_max_t << std::endl;
    std::cout << "Spanner Theta mean t-value: " << (theta_mean_t / number_of_tests) << std::endl;
    std::cout << "All tests completed." << std::endl;

}


