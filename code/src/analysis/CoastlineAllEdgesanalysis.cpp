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

#include "MultiThetaAnalysis.h"
#include "../daniel/theta-graph/headers/Progressbar.h"
#include "../spanner/ThetaSpanner.hpp"
#include "../daniel/theta-graph/headers/Triangulation.h"
#include "../io/DataWriter.h"


using namespace std;

vector<unsigned long> bucket_edge_counter(std::string base_graph_path, Graph* spanner_graph, int number_of_buckets, double longest_edge_distance = -1.0, bool log_dist=false) {
    Triangulation triangulation;
    triangulation.readFromGraph(base_graph_path);

    if (longest_edge_distance == -1) {
        double longest_edge_distance = 0.0;
        int max_dist_node = -1;
        // randomly choose 3 nodes
        int node1 = rand() % spanner_graph->adj.size();
        // identify maximum distance between them and all nodes
        auto res1 = spanner_graph->multiSourceMultiTargetDijkstra({node1}, vector<int>(spanner_graph->adj.size()), true);
        // use that as max distance
        for (auto pair : res1) {
            if (pair.second != numeric_limits<double>::infinity() && pair.second > longest_edge_distance) {
                longest_edge_distance = pair.second;
                max_dist_node = pair.first[pair.first.size() - 1];
            }
        }
        auto res2 = spanner_graph->multiSourceMultiTargetDijkstra({max_dist_node}, vector<int>(spanner_graph->adj.size()), true);
        for (auto pair : res2) {
            if (pair.second != numeric_limits<double>::infinity() && pair.second > longest_edge_distance) {
                longest_edge_distance = pair.second;
            }
        }
    }
    double longest_found_edge_distance = 0.0;

    std::cout << "Triangulation loaded" << std::endl;

    int num_threads = std::max(1,omp_get_num_procs()-1);  // Get the number of available processors
    // num_threads = 1;
    std::vector<std::vector<unsigned long>> buckets_per_thread = std::vector<std::vector<unsigned long>>();
    for (int i = 0; i < num_threads; i++) {
        buckets_per_thread.push_back(vector<unsigned long>(number_of_buckets, 0));
    }

    std::cout << "Attention! total elements in the next lines are presented per thread! The true total is this number times the number of threads" << std::endl;

    int target_size = spanner_graph->adj.size();
    //target_size = 10000; // TODO remove limit

    #pragma omp parallel for num_threads(num_threads)
    for (int u = 0; u < target_size; u++) {
        // progress
        auto thread_numt = omp_get_thread_num();
        if (thread_numt == 0) {
            if (u%1000==0) {
                std::cout << "Thread #" << thread_numt << ": " <<u << " of " << (spanner_graph->adj.size()/num_threads) << std::endl;
            }
        }

        auto targets = triangulation.oneToAllVisibility(u, false);
        for (auto v_id : targets) {
            int v = static_cast<int>(v_id);
            if (u == v) continue;
            auto a = spanner_graph->id_point_map[u];
            auto b = spanner_graph->id_point_map[v];
            auto dx = abs(a.x - b.x);
            auto dy = abs(a.y - b.y);
            double dist = sqrt(dx*dx + dy*dy);
            if (dist > longest_found_edge_distance) {
                longest_found_edge_distance = dist;
            }

            // Calculate bucket
            double rel_dist;
            if (!log_dist) rel_dist = (double) dist / (double) longest_edge_distance;
            else rel_dist = (double) log(dist+1) / (double) log(longest_edge_distance); // log scale
            int bucket_id = (int)(rel_dist * number_of_buckets);

            // Clamp the bucket_id to be in the valid range [0, number_of_buckets - 1]
            if (bucket_id >= number_of_buckets) {
                bucket_id = number_of_buckets - 1;
            } else if (bucket_id < 0) {
                // This shouldn't happen if dist is always positive, but good to have
                std::cout << "WARNING: Triangulation returned negative number: " << dist << std::endl;
                bucket_id = 0;
            }

            // Now, increment the correct counter using the thread number
            buckets_per_thread[thread_numt][bucket_id] += 1;
        }
    }
    std::vector<unsigned long> buckets;


    for (int i = 0; i < number_of_buckets; i++) {
        int sum = 0;
        for (int t = 0; t < num_threads; t++) {
            sum += buckets_per_thread[t][i];
        }
        buckets.push_back(sum);
        std::cout << "Bucket " << i << ": "<< "up to :" << (longest_edge_distance/number_of_buckets)*(i+1) << ":"<< sum << ": edges." << std::endl;
    }
    std::cout << "Number of buckets: " << buckets.size() << std::endl;
    std::cout << "longest found edge: " << longest_found_edge_distance << std::endl;
    return buckets;
}



int main(int argc, char* argv[]) {

    string coastline_path;
    string spanner_path;
    string csv_path;
    double cutoff = 1.1;
    int k = 75;


    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "-sg" && i + 1 < argc) {
            spanner_path = argv[++i];
            std::cout << "Spanner path: " << spanner_path << std::endl;
        } else if (arg == "-cg" && i + 1 < argc) {
            coastline_path = argv[++i];
            std::cout << "Coastline graph path: " << coastline_path << std::endl;
        }else if (arg == "-c" && i + 1 < argc) {
            cutoff = stod(argv[++i]);
            std::cout << "percent of worst edges to add to spanner: " << cutoff << std::endl;
        }else if (arg == "-k" && i + 1 < argc) {
            k = stoi(argv[++i]);
            std::cout << "k cones for the spanner (ehem. theta): " << k << std::endl;
        }
    }

    if (spanner_path.empty() || coastline_path.empty() || cutoff < 0) {
        std::cerr << "Usage: " << argv[0] << " -sg <string> -vg <string> -csv <string> -p <double> -k <int>\n";
        for (int i = 0; i < argc; i++) std::cerr << " " << argv[i];
        return 1;
    }

    ///////////////////////////////////////////////////////////////////////////////////
    /// LOAD THE GRAPHS
    ///////////////////////////////////////////////////////////////////////////////////
    auto spanner_graph = get<1>(load_fmi(spanner_path,  -1));
    auto coastline_graph = load_coastline(coastline_path);

    //sanity check that this is the spanner for this coastline
    if (spanner_graph.adj.size() != coastline_graph.adj.size()) {
        std::cerr << "Error: Spanner graph and coastline graph have different number of nodes!" << std::endl;
        return 1;
    }

    bucket_edge_counter(coastline_path, &spanner_graph, 100, 4.00751e+07, false);

    //auto add_edges = identify_bad_edges_sorted_by_length(coastline_path, &spanner_graph, cutoff, 100);

    //spanner_graph.sort_edges();
    cout << "";
    //dijkstra_debugging(spanner_graph);
    //auto edges = analyse_spanner_with_coastline_graph(coastline_path, spanner_graph, cutoff);


    //for (Edge edge: edges) {
    //    spanner_graph.addEdge(edge.source, edge.target, edge.weight, true);
    //}
    //cout << "Added " << edges.size() << " edges to the spanner." << endl;
    //write_fmi("../../data/spanner_with_worst_edges-"+to_string(k)+"-"+to_string(cutoff)+".fmi", spanner_graph);
    //string new_csv_path = csv_path.substr(0, csv_path.find_last_of('.')) + "_with_worst_edges.csv";
    //analyse_spanner_with_vis_graph(base_graph, spanner_graph, new_csv_path, "../../data/all_edges-"+to_string(theta)+"-"+to_string(percent)+".csv", 0.0);
}


