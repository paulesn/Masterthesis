//
// Created by Sebastian Paule on 9/15/25.
//
/**
* This analysis file starts a analysis, then adds the X% worst edges and analyses it again
*/


//
// Created by Sebastian Paule on 9/4/25.
//

#include <fstream>

#include "../io/Dataloader.hpp"
#include "../io/DataWriter.h"
#include <omp.h>

#include "MultiThetaAnalysis.h"
#include "Triangulation.h"
#include "../daniel/theta-graph/headers/Progressbar.h"
#include "../spanner/ThetaSpanner.hpp"
#include "../structure/Graph.hpp"


using namespace std;

double greedy_optimize_path(Graph &spanner, string coastline_path, vector<int> &path) {
    /*Triangulation tri = Triangulation();
    tri.readFromGraph(coastline_path);
    double path_length;
    for (int i = 0; i < path.size(); i++) {
        for (int j = i + 2; j < path.size(); j++) {
            if (tri.directVisibility(path[i], path[j])) {
                // remove all nodes between i and j
                path.erase(path.begin() + i + 1, path.begin() + j);
                // recalculate path length
                path_length = 0.0;
                for (int k = 0; k < path.size() - 1; k++) {
                    path_length += euk_dist(spanner.id_point_map[k], spanner.id_point_map[k + 1]);
                }
                // restart from beginning
                i = -1;
                break;
            }
        }
    }*/
    return 1.0;
}


vector<tuple<int,int,double>> analysis_with_distance_file(Graph &vis_graph, Graph &spanner, vector<pair<int,int>> &node_pairs) {
    vector<vector<tuple<int, int,double>>> distances = vector<vector<tuple<int,int,double>>>(omp_get_max_threads());
    int num_threads = std::max(1,omp_get_num_procs()-1);  // Get the number of available processors
    int count = 0;
    int max = node_pairs.size();
    #pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (int i = 0; i < node_pairs.size(); i++) {
        auto pair = node_pairs[i];

        if (omp_get_thread_num() == 0) {
            count++;
            if (count % 100 == 0) {
                std::cout << "Analyzing pair " << count << "/" << max/num_threads << std::endl;
            }
        }
        auto res = spanner.dijkstra(pair.first, pair.second);
        double spanner_dist = res.second;


        distances[omp_get_thread_num()].push_back({pair.first, pair.second,spanner_dist});
    }
    vector<tuple<int,int,double>> all_distances = vector<tuple<int,int, double>>();
    for (auto &thread_distances : distances) {
        all_distances.insert(all_distances.end(), thread_distances.begin(), thread_distances.end());
    }
    return all_distances;
}

vector<tuple<int,int,double>> analysis_with_distance_file_greedy_optimized(Graph &vis_graph, Graph &spanner, vector<pair<int,int>> &node_pairs, string coastline_path) {
    vector<vector<tuple<int, int,double>>> distances = vector<vector<tuple<int,int,double>>>(omp_get_max_threads());
    int num_threads = 1; //std::max(1,omp_get_num_procs()-1);  // Get the number of available processors
    int count = 0;
    int max = node_pairs.size();
#pragma omp parallel for num_threads(num_threads) schedule(dynamic)
    for (int i = 0; i < node_pairs.size(); i++) {
        auto pair = node_pairs[i];

        if (omp_get_thread_num() == 0) {
            count++;
            if (count % 100 == 0) {
                std::cout << "Analyzing pair " << count << "/" << max/num_threads << std::endl;
            }
        }
        auto res = spanner.dijkstra(pair.first, pair.second);
        double spanner_dist = greedy_optimize_path(spanner, coastline_path, res.first);


        distances[omp_get_thread_num()].push_back({pair.first, pair.second,spanner_dist});
    }
    vector<tuple<int,int,double>> all_distances = vector<tuple<int,int, double>>();
    for (auto &thread_distances : distances) {
        all_distances.insert(all_distances.end(), thread_distances.begin(), thread_distances.end());
    }
    return all_distances;
}

tuple<vector<pair<int,int>>,vector<double>,vector<double>> load_distance_data(string dist_data_path) {
    vector<pair<int,int>> node_pairs;
    vector<double> node_distances;
    vector<double> original_distances;

    std::ifstream infile(dist_data_path);
    string line;
    while (getline(infile, line)) {
        // Assuming each line in the file is formatted as: u,v,distance, org_dist
        stringstream ss(line);
        string u_str, v_str, dist_str, org_dist_str;
        getline(ss, u_str, ',');
        getline(ss, v_str, ',');
        getline(ss, dist_str, ',');
        getline(ss,org_dist_str, ',');
        int u = stoi(u_str);
        int v = stoi(v_str);
        double distance = stod(dist_str);
        double org_dist = stod(org_dist_str);
        node_pairs.push_back({u, v});
        node_distances.push_back(distance);
        original_distances.push_back(org_dist);
    }
    return {node_pairs, node_distances, original_distances};
}



int main(int argc, char* argv[]) {

    string vis_graph_path;
    //string dist_data_path;
    string daniel_graph_path;
    string out_path;
    double error_limit = 1.1;
    int k = 24;


    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        if (arg == "-v" && i + 1 < argc) {
            vis_graph_path = argv[++i];
            std::cout << "vis graph path: " << vis_graph_path << std::endl;
        } else if (arg == "-d" && i + 1 < argc) {
            daniel_graph_path = argv[++i];
            std::cout << "daniel graph path: " << daniel_graph_path << std::endl;
        }else if (arg == "-e" && i + 1 < argc) {
            error_limit = stod(argv[++i]);
            std::cout << "error limit for modified theta spanner: " << error_limit << std::endl;
        }else if (arg == "-k" && i + 1 < argc) {
            k  = stoi(argv[++i]);
            std::cout << "k cones for the spanner (ehem. theta): " << k << std::endl;
        }else if (arg == "-o" && i + 1 < argc) {
            out_path  = argv[++i];
            std::cout << "path to store data: " << out_path << std::endl;
        }
    }

    if (daniel_graph_path.empty() || out_path.empty() || vis_graph_path.empty()) {
        std::cerr << "Usage: " << argv[0] << " -v <string> -d <string> -o <string> -e <double> -k <int>\n";
        for (int i = 0; i < argc; i++) std::cerr << " " << argv[i];
        return 1;
    }

    ///////////////////////////////////////////////////////////////////////////////////
    /// LOAD THE GRAPHS
    ///////////////////////////////////////////////////////////////////////////////////
    auto base_graph = get<1>(load_fmi(vis_graph_path));
    //auto spanner_graph = create_theta_spanner_graph(&base_graph, theta); //get<1>(load_fmi(spanner_path,  -1));
    std::cout << "Graph loaded. Sorting Edges" << std::endl;
    base_graph.sort_edges();
    std::cout << "Graph sorted. Loading Daniel" << std::endl;
    Graph daniel_graph = get<1>(load_fmi(daniel_graph_path));
    std::cout << "Graph loaded. Creating Spanner with k= " << k << std::endl;

    Graph spanner_graph = create_theta_spanner_graph(&base_graph, k, false);
    write_fmi(out_path, spanner_graph);

    std::cout << "Spanner has " << spanner_graph.number_of_edges << " edges" << std::endl;

    std::cout << "Spanner Created. Analysing Graphs" << std::endl;

    vector<pair<int,int>> node_pairs = {};

    for (int i = 0; i < 100; i++) {
        // select two random nodes from the spanner
        int u = rand() % spanner_graph.n;
        int v = rand() % spanner_graph.n;
        node_pairs.push_back({u, v});

        // sanity check if it is the same point in both graphs
        auto p_u_s = spanner_graph.id_point_map[u];
        auto p_u_d = daniel_graph.id_point_map[u];
        auto p_v_s = spanner_graph.id_point_map[v];
        auto p_v_d = daniel_graph.id_point_map[v];
        if (normalize(euk_dist(p_u_s, p_u_d)) > 0.0001 || normalize(euk_dist(p_v_s, p_v_d)) > 0.0001) {
            std::cerr << "ERROR: Node " << u << " or " << v << " differ between spanner and daniel graph." << std::endl;
        }
    }

    string c_path = "../data/world-shrunk-coastlines-0500.32.txt";

    auto results_spanner = analysis_with_distance_file(base_graph, spanner_graph, node_pairs);
    //sort results_spanner by first and second
    sort(results_spanner.begin(), results_spanner.end(), [](const tuple<int,int,double> &a, const tuple<int,int,double> &b) {
        if (get<0>(a) == get<0>(b)) {
            return get<1>(a) < get<1>(b);
        }
        return get<0>(a) < get<0>(b);
    });
    auto results_daniel = analysis_with_distance_file(base_graph, daniel_graph, node_pairs);
    //sort results_daniel by first and second
    sort(results_daniel.begin(), results_daniel.end(), [](const tuple<int,int,double> &a, const tuple<int,int,double> &b) {
        if (get<0>(a) == get<0>(b)) {
            return get<1>(a) < get<1>(b);
        }
        return get<0>(a) < get<0>(b);
    });

    double sum_vis_dist = 0.0;
    double sum_spanner_dist = 0.0;
    double sum_daniel_dist = 0.0;

    int err_count = 0;
    for (int i = 0; i < node_pairs.size(); i++) {

        auto [u_s,v_s, spanner_dist] = results_spanner[i];
        auto [u_d, v_d,daniel_dist] = results_daniel[i];


        if (u_s != u_d || v_s != v_d) {
            std::cerr << "ERROR: Node pairs do not match between spanner and daniel results: "
                      << u_s << "," << v_s << " vs " << u_d << "," << v_d << std::endl;
            err_count++;
            continue;
        }

        double vis_dist = base_graph.dijkstra(u_s, v_s).second;
        if (vis_dist == numeric_limits<double>::infinity()) {
            std::cerr << "ERROR: No path found in vis graph between nodes " << u_s << " and " << v_s << std::endl;
            err_count++;
            continue;
        }

        if (vis_dist > spanner_dist || vis_dist > daniel_dist) {
            std::cerr << "ERROR: Vis graph distance is greater than spanner or daniel distance for nodes " << u_s << " and " << v_s << ": " << vis_dist << " vs " << spanner_dist << " and " << daniel_dist << std::endl;
            err_count++;
            continue;
        }

        sum_spanner_dist += spanner_dist;
        sum_daniel_dist += daniel_dist;
        sum_vis_dist += vis_dist;

        double t_spanner = spanner_dist / vis_dist;
        double t_daniel = daniel_dist / vis_dist;

        cout << "Nodes: " << u_s << " - " << v_s << " | Vis Dist: " << vis_dist
             << " | Spanner Dist: " << spanner_dist << " (t=" << t_spanner << ")"
             << " | Daniel Dist: " << daniel_dist << " (t=" << t_daniel << ")" << endl;
    }
    cout << "----------------------------------------" << endl;
    int total = node_pairs.size()-err_count;
    cout << "Average Vis Dist: " << sum_vis_dist / total << endl;
    cout << "Average Spanner Dist: " << sum_spanner_dist / total << " (t=" << (sum_spanner_dist / sum_vis_dist) << ")" << endl;
    cout << "Average Daniel Dist: " << sum_daniel_dist / total << " (t=" << (sum_daniel_dist / sum_vis_dist) << ")" << endl;
}


