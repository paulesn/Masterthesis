//
// Created by Sebastian Paule on 9/12/25.
//
#include <iostream>
#include "../io/Dataloader.hpp"
#include "../io/DataWriter.h"
#include "../spanner/ThetaSpanner.hpp"
#include "../structure/Quadtree.hpp"

using namespace std;

int main(int argc, char* argv[]) {

    string base_graph_path;
    string spanner_graph_path;
    string out_path;
    int theta;
    double seperation_factor;


    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "-o" && i + 1 < argc) {
            out_path = argv[++i];
            std::cout << "output path: " << out_path << std::endl;
        } else if (arg == "-bg" && i + 1 < argc) {
            base_graph_path = argv[++i];
            std::cout << "Base vis graph path: " << base_graph_path << std::endl;
        } else if (arg == "-t" && i + 1 < argc) {
            theta = stoi(argv[++i]);
            std::cout << "theta value: " << theta << std::endl;
        }else if (arg == "-s" && i + 1 < argc) {
            seperation_factor = stod(argv[++i]);
            std::cout << "seperation factor: " << seperation_factor << std::endl;
        }else if (arg == "-sg" && i + 1 < argc) {
            spanner_graph_path = argv[++i];
            std::cout << "spanner graph path: " << spanner_graph_path << std::endl;
    }
    }

    ///////////////////////////////////////////////////////////////////////////////////
    /// LOAD THE ORIGINAL GRAPH
    ///////////////////////////////////////////////////////////////////////////////////
    auto tup = load_fmi(base_graph_path,  -1);
    auto points = get<0>(tup);
    auto used_points = vector<int>();
    Graph graph = get<1>(tup);
    // init the hub labels for faster shortest path distance calculation
    //graph.init_hub_labels();
    std::cout << "Loaded " << points.size() << " points." << std::endl;

    Graph spanner_theta = std::get<1>(load_fmi(spanner_graph_path,  -1));
    Graph spanner_wspd = Graph(graph.n);
    Quadtree tree = Quadtree();

    for (auto &p : points) {
        tree.insert(p);
    }
    std::cout << "Created quadtree." << std::endl;

    spanner_wspd = wspd_for_visibility(&tree, seperation_factor, &graph);
    std::cout << "Created WSPD spanner graph with " << spanner_wspd.number_of_edges << " edges." << std::endl;

    // merge the two spanners
    Graph spanner_mix = Graph(graph.n);
    spanner_mix = spanner_theta;
    for (int source = 0; source < spanner_wspd.adj.size(); source++) {
        for (const auto &edge : spanner_wspd.adj[source]) {
            spanner_mix.addEdge(edge.source, edge.target, edge.weight);
        }
    }
    std::cout << "Created mixed spanner graph with " << spanner_mix.number_of_edges << " edges." << std::endl;

    // store graph to disk
    write_fmi(out_path, spanner_mix);
    write_gf(out_path + ".gf", theta, {{0,0}}, 0);
}