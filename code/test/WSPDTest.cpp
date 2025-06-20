
#include "../src/Dataloader.hpp"
#include "../src/Quadtree.hpp"
#include "../src/ThetaSpanner.hpp"
#include <ostream>
#include <iostream>
#include <set>
#include <chrono>
#include <csignal>

#include "../lib/json.hpp"


using namespace std;



int main() {
    int theta = 4; // Number of zones
    double s = 2;
    bool using_wspd_e = true;
    bool using_wspd_spd = false;
    bool using_theta = true;

    string path = "../../data/0025.32.fmi";

    ///////////////////////////////////////////////////////////////////////////////////
    /// LOAD THE ORIGINAL GRAPH
    ///////////////////////////////////////////////////////////////////////////////////
    auto tup = load_fmi(path,  -1);
    auto systems = get<0>(tup);
    Graph graph = get<1>(tup);
    std::cout << "Loaded " << systems.size() << " systems." << std::endl;

    Graph spanner_e = Graph(graph.n);
    Graph spanner_sp = Graph(graph.n);
    Graph spanner_theta = Graph(graph.n);

    ///////////////////////////////////////////////////////////////////////////////////
    /// INSERT SYSTEMS INTO THE QUADTREE
    ///////////////////////////////////////////////////////////////////////////////////
    auto tree = new Quadtree(32);
    int counter = 0;
    std::cout << "Inserting systems into the quadtree..." << std::endl;
    for (auto& system : systems) {
        tree->insert(Point(system.x, system.y, counter));
        counter++;
    }
    auto all_points = tree->get_all_points();
    if (all_points.size() != systems.size()) {
        std::cout << "Error: Quadtree contains " << all_points.size() << " points, but expected " << systems.size() << "." << std::endl;
        raise(SIGINT); // Raise SIGINT to terminate the program
    }


    std::cout << "Inserted all systems into the quadtree." << std::endl;
    ///////////////////////////////////////////////////////////////////////////////////
    /// CALCULATE WSPD (euklidian and shortest path distance)
    ///////////////////////////////////////////////////////////////////////////////////

    if (using_wspd_e){
        auto start1 = std::chrono::high_resolution_clock::now();
        auto result_e = wspd(tree, s);
        auto end1 = std::chrono::high_resolution_clock::now();

        std::cout << "wspd_e("<< s << ") took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end1 - start1).count()
              << " ms. And created "<< result_e.size()
              << " edges." << std::endl;
        // Create a spanner graph from the result_e
        spanner_e = Graph(graph.n);
        for (const auto& pair : result_e) {
            if (get<0>(pair) == get<1>(pair)) continue; // skip self-pairs
            auto p1 = get<0>(pair);
            auto p2 = get<1>(pair);
            spanner_e.addEdge(p1.id, p2.id, euklidian_distance(p1,p2), true);
        }
    }

    if (using_wspd_spd) {
        auto start2 = std::chrono::high_resolution_clock::now();
        auto result_sp = wspd_spd(tree, s, graph);
        auto end2 = std::chrono::high_resolution_clock::now();

        std::cout << "wspd_sp("<< s << ") took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end2 - start2).count()
              << " ms. And created " << result_sp.size() << " edges." << std::endl;
        std::cout << "For   " << systems.size() << " nodes." << std::endl;

        // Create a spanner graph from the result_sp
        spanner_sp = Graph(graph.n);
        for (const auto& pair : result_sp) {
            if (get<0>(pair) == get<1>(pair)) continue; // skip self-pairs
            auto p1 = get<0>(pair);
            auto p2 = get<1>(pair);
            spanner_sp.addEdge(p1.id, p2.id, euklidian_distance(p1,p2));
        }
    }

    if (using_theta) {
        auto start3 = std::chrono::high_resolution_clock::now();
        spanner_theta = create_theta_spanner_graph(&graph, theta);
        auto end3 = std::chrono::high_resolution_clock::now();

        std::cout << "theta("<< theta << ") took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end3 - start3).count()
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
        int random_source = rand() % graph.n;
        int random_target = rand() % graph.n;
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
        if (spanner_e_dist < original_dist) {
            e_has_neg_t = true;
        }
        if (spanner_sp_dist < original_dist) {
            sp_has_neg_t = true;
        }
        if (spanner_theta_dist < original_dist) {theta_has_neg_t = true;}
        if (spanner_e_dist/original_dist > e_max_t) {e_max_t = spanner_e_dist/original_dist;}
        if (spanner_sp_dist/original_dist > sp_max_t) {sp_max_t = spanner_sp_dist/original_dist;}
        if (spanner_theta_dist/original_dist > theta_max_t) {theta_max_t = spanner_theta_dist/original_dist;}

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
    std::cout << std::endl << "Results after " << number_of_tests << " tests:" << std::endl;
    std::cout << "Spanner SPD has infinity: " << (sp_has_inf ? "Yes" : "No") << std::endl;
    std::cout << "Spanner SPD has negative t-value: " << (sp_has_neg_t ? "Yes" : "No") << std::endl;
    std::cout << "Spanner SPD max t-value: " << sp_max_t << std::endl;
    std::cout << "Spanner SPD mean t-value: " << (sp_mean_t / number_of_tests) << std::endl;
    std::cout << "Spanner E has infinity: " << (e_has_inf ? "Yes" : "No") << std::endl;
    std::cout << "Spanner E has negative t-value: " << (e_has_neg_t ? "Yes" : "No") << std::endl;
    std::cout << "Spanner E max t-value: " << e_max_t << std::endl;
    std::cout << "Spanner E mean t-value: " << (e_mean_t/ number_of_tests)<< std::endl;
    std::cout << "Spanner Theta has infinity: " << (t_has_inf ? "Yes" : "No") << std::endl;
    std::cout << "Spanner Theta has negative t-value: " << (theta_has_neg_t ? "Yes" : "No") << std::endl;
    std::cout << "Spanner Theta max t-value: " << theta_max_t << std::endl;
    std::cout << "Spanner Theta mean t-value: " << (theta_mean_t / number_of_tests) << std::endl;
    std::cout << "All tests completed." << std::endl;

}


