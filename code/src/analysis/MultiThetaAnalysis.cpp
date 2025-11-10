#include <fstream>
#include <cmath>
#include <stdint.h>
#include <atomic>
#include <random>
#include <chrono>

#include "../io/Dataloader.hpp"
#include <omp.h>
#include "../daniel/theta-graph/headers/Progressbar.h"
#include "../daniel/theta-graph/headers/Triangulation.h"
#include "../spanner/ThetaSpanner.hpp"
#include "MultiThetaAnalysis.h"
#include "../daniel/theta-graph/headers/Timer.h"
#include "../io/DataWriter.h"


using namespace std;

/**
 * floating point normalization to prevent floating point error missmatches
 * @param value
 * @return
 */
double normalize(double value) {
    return (value * 10000.0)/ 10000.0;
}

double euk_dist(Pointc p1, Pointc p2) {
    return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

int analyse_spanner(Graph base_graph, Graph spanner_graph, string csv_out_path, string graph_path, int max, bool multi_threading) {
    ///////////////////////////////////////////////////////////////////////////////////
    /// TEST random queries on the spanners
    /// analyse t-value
    ///////////////////////////////////////////////////////////////////////////////////

    if (max == -1) {
        max = base_graph.adj.size();
    }

    vector<double> t_values_sum = vector<double>(base_graph.adj.size(), 0);
    vector<double> t_values_max = vector<double>(base_graph.adj.size(), 0);
    vector<int> edges_v = vector<int>(base_graph.adj.size(), 0);
    vector<vector<int>> t_histogram = vector<vector<int>>(); // histogram with 20 bins
    for (int i = 0; i < base_graph.adj.size(); i++) {
        t_histogram.push_back(vector<int>(100));
    }

    Triangulation triangulation;
    triangulation.readFromGraph(graph_path);

    cout << "Triangulation loaded";

    int num_threads = std::max(1,omp_get_num_procs()-1);  // Get the number of available processors
    if (!multi_threading) {
        num_threads = 1;
    }
    cout << "Using " << num_threads << " threads for analysis." << endl;

    Timer timer;
    timer.start();
    ProgressBar progressBar;
    progressBar.start(max);

#pragma omp parallel for num_threads(num_threads) shared(spanner_graph, base_graph, t_values_sum, t_values_max, triangulation) schedule(dynamic)
    for (int source = 0; source < max; source++) {
        progressBar.update(1);
        if (base_graph.adj[source].size() == 0) {
            //std::cout << "Node " << source << " has no edges in the base graph. Skipping." << std::endl;
            continue;
        }

        std::vector<GlobalID> temp = triangulation.oneToAllVisibility(source, false);
  std::vector<int> targets;
  targets.reserve(temp.size());
  for (auto id : temp) targets.push_back(static_cast<int>(id));

        if (targets.size() == 0) {
            std::cerr << "Node " << source << " has no visible targets in the triangulation. Skipping." << std::endl;
            continue;
        }
        if (spanner_graph.adj[source].size() == 0) {
            //std::cerr << "Node " << source << " has no edges in the spanner graph. Skipping." << std::endl;
            continue;
        }

        auto res = spanner_graph.multiSourceMultiTargetDijkstra({source}, targets, true);

        for (int i = 0; i < targets.size(); i++) {

            int target = targets[i];
            if (source > target) {
                continue; // to avoid double checks in undirected graph
            }

            double original_dist = euk_dist(spanner_graph.id_point_map.at(source), spanner_graph.id_point_map.at(target));
            auto spanner_theta_dist = numeric_limits<double>::infinity();

            auto path = res[i].first;
            // santiy check for path
            int path_target = path[path.size() - 1];
            if (path_target != target) {
                //std::cout << "ERROR: Path target does not match expected target for source " << source << " and target " << target << std::endl;
                continue;
            }
            double spanner_dist  = res[i].second;

            double t_value = spanner_dist/original_dist;
            vector<Pointc> path_points;
            for (int node_id : path) {
                path_points.push_back(base_graph.id_point_map.at(node_id));
            }
            if (t_value > 2.0) write_gf("../../data/daniel/"+to_string(t_value)+".gf", 0, path_points, 0);


            t_values_sum[source] += t_value;
            edges_v[source] += 1;
            if (t_value > t_values_max[source]) {
                t_values_max[source] = t_value;
            }
            int histogram_index = std::min(static_cast<int>((t_value-1) * 1000), 99); // cap at 10.0
            t_histogram[i][histogram_index] += 1;



                    // sanity checks
                    if (normalize(spanner_dist) < normalize(original_dist)) {
                        //std::cout << "ERROR: Spanner distance is less than original distance for source " << source << " and target " << target << std::endl;
                        //std::cout << "Original distance: " << original_dist << ", Spanner distance: " << spanner_dist << std::endl;
                    }
                    if (spanner_dist == numeric_limits<double>::infinity()) {
                        //std::cout << "ERROR: Spanner distance is infinity for source " << source << " and target " << target << std::endl;
                    }
                }
            }
            progressBar.stop();

            ////////////////////////////////////////////////////////////////////////////////////
            /// WRITE THE RESULTS
            ////////////////////////////////////////////////////////////////////////////////////

            double t_values_all_sum = 0.0;
            int edges = 0;
            double max_t_all = 0.0;
            vector<int> full_histogramm = vector<int>(100);

            for (int i = 0; i < t_values_sum.size(); i++) {
                if (base_graph.adj[i].size() > 0) {
                    t_values_all_sum += t_values_sum[i];
                    edges += edges_v[i];
                    if (t_values_max[i] > max_t_all) {
                        max_t_all = t_values_max[i];
                    }
                    for (int j = 0; j < t_histogram[i].size(); j++) {
                        full_histogramm[j] += t_histogram[i][j];
                    }
                }
            }

            cout << "----------------------------------------------------------------------------------" << endl;
            cout << endl << "Results after analysing all edges:" << endl;
            cout << "----------------------------------------------------------------------------------" << endl;
            cout << "Average t-value: " << normalize(t_values_all_sum/edges) << endl;
            cout << "Max t-value: " << normalize(max_t_all) << endl;
            cout << "Number of edges: " << edges << endl;
            cout << "----------------------------------------------------------------------------------" << endl;

            auto log_stream = ofstream("./../data/analysis_log.csv", ios_base::app);
            log_stream  << "i," << normalize(t_values_all_sum/edges) << "," << normalize(max_t_all) << endl;
            cout  << "i," << normalize(t_values_all_sum/edges) << "," << normalize(max_t_all) << endl;

            // TODO store full histogramm in csv
            auto csv_out_stream = ofstream(csv_out_path);
            csv_out_stream << "t_value,frequency" << endl;
            for (int i = 0; i < full_histogramm.size(); i++) {
                csv_out_stream << normalize(1.0 + i * 0.001) << "," << full_histogramm[i] << endl;
            }
            csv_out_stream.close();
            return 0;
        }

std::vector<Edge> analyse_spanner_with_vis_graph(Graph base_graph, Graph spanner_graph, std::string csv_out_path, std::string all_edges_path, double percent) {
    ///////////////////////////////////////////////////////////////////////////////////
    /// TEST random queries on the spanners
    /// analyse t-value
    ///////////////////////////////////////////////////////////////////////////////////

    vector<double> t_values_sum = vector<double>(base_graph.adj.size(), 0);
    vector<double> t_values_max = vector<double>(base_graph.adj.size(), 0);
    vector<int> edges_v = vector<int>(base_graph.adj.size(), 0);
    vector<vector<int>> t_histogram = vector<vector<int>>(); // histogram with 20 bins
    vector<vector<Edge>> worst_edges = vector<vector<Edge>>(base_graph.adj.size());
    vector<string> all_edges = vector<string>(base_graph.adj.size());
    for (int i = 0; i < base_graph.adj.size(); i++) {
        t_histogram.push_back(vector<int>(100));
        worst_edges.push_back(vector<Edge>());
    }

    int num_threads = std::max(1,omp_get_num_procs()-1);  // Get the number of available processors
    //num_threads = 10; // TODO remove limit

    Timer timer;
    timer.start();
    ProgressBar progressBar;
    int max = base_graph.adj.size(); // TODO remove limit
    progressBar.start(max);
    int worst_edges_vector_size = (int)(base_graph.adj.size()*percent); // keep track of the 100 worst edges per node

#pragma omp parallel for num_threads(num_threads) shared(spanner_graph, base_graph, t_values_sum, t_values_max) schedule(dynamic)
    for (int source = 0; source < max; source++) {
        progressBar.update(1);
        if (base_graph.adj[source].size() == 0) {
            //std::cout << "Node " << source << " has no edges in the base graph. Skipping." << std::endl;
            continue;
        }

        std::vector<Edge> edges =  base_graph.adj[source];

        if (edges.size() == 0) {
            std::cerr << "Node " << source << " has no visible targets in the triangulation. Skipping." << std::endl;
            continue;
        }
        if (spanner_graph.adj[source].size() == 0) {
            //std::cerr << "Node " << source << " has no edges in the spanner graph. Skipping." << std::endl;
            continue;
        }
        // iterate over each target
        for (int i = 0; i < edges.size(); i++) {

            int target = edges[i].target;
            Edge edge = edges[i];
            if (source > target) {
                continue; // to avoid double checks in undirected graph
            }

            double original_dist = edge.weight;
            auto spanner_theta_dist = numeric_limits<double>::infinity();
            double spanner_dist  = spanner_graph.dijkstra(source, target).second;
            if (normalize(spanner_dist) < normalize(original_dist)) {
                std::cerr << "Spanner distance is less than original distance for edge " << source << " - " << target
                    << ": " << spanner_dist << " vs " << original_dist << std::endl;
            }

            double t_value = spanner_dist/original_dist;
            worst_edges[source].push_back(Edge(source, target, edge.weight));
            worst_edges[source].back().t_value = t_value;
            all_edges[source] += to_string(source) + "," + to_string(target) + "," + to_string(normalize(t_value)) + "," + to_string(original_dist) + "\n";

            t_values_sum[source] += t_value;
            edges_v[source] += 1;
            if (t_value > t_values_max[source]) {
                t_values_max[source] = t_value;
            }
            int histogram_index = std::min(static_cast<int>((t_value-1) * 1000), 99); // cap at 10.0
            t_histogram[i][histogram_index] += 1;

                    // sanity checks
                    if (normalize(spanner_dist) < normalize(original_dist)) {
                        //std::cout << "ERROR: Spanner distance is less than original distance for source " << source << " and target " << target << std::endl;
                        //std::cout << "Original distance: " << original_dist << ", Spanner distance: " << spanner_dist << std::endl;
                    }
                    if (spanner_dist == numeric_limits<double>::infinity()) {
                        //std::cout << "ERROR: Spanner distance is infinity for source " << source << " and target " << target << std::endl;
                    }
                }
            }
            progressBar.stop();

            ////////////////////////////////////////////////////////////////////////////////////
            /// WRITE THE RESULTS
            ////////////////////////////////////////////////////////////////////////////////////

            double t_values_all_sum = 0.0;
            int edges = 0;
            double max_t_all = 0.0;
            vector<int> full_histogramm = vector<int>(100);
            vector<Edge> full_worst_edges;
            auto edge_out_stream = ofstream(all_edges_path, ios_base::app);
            //edge_out_stream << "source,target,t_value,original_distance\n";

            for (int i = 0; i < t_values_sum.size(); i++) {
                if (base_graph.adj[i].size() > 0) {
                    t_values_all_sum += t_values_sum[i];
                    edges += edges_v[i];
                    if (t_values_max[i] > max_t_all) {
                        max_t_all = t_values_max[i];
                    }
                    for (int j = 0; j < t_histogram[i].size(); j++) {
                        full_histogramm[j] += t_histogram[i][j];
                    }
                    // merge worst edges results
                    for(Edge edge:worst_edges[i]) {
                        full_worst_edges.push_back(edge);
                    }
                    // write all edges to file
                    edge_out_stream << all_edges[i];
                }
            }
            edge_out_stream.close();
            // sort worst edges
            std::sort(full_worst_edges.begin(), full_worst_edges.end(), [](const Edge &a, const Edge &b) {
                return a.t_value > b.t_value;
            });
            if (full_worst_edges.size() > worst_edges_vector_size) {
                full_worst_edges.resize(worst_edges_vector_size);
            }

            cout << "----------------------------------------------------------------------------------" << endl;
            cout << endl << "Results after analysing all edges:" << endl;
            cout << "----------------------------------------------------------------------------------" << endl;
            std::cout << "The Spanner has now " << spanner_graph.number_of_edges << " edges." << std::endl;
            cout << "Average t-value: " << normalize(t_values_all_sum/edges) << endl;
            cout << "Max t-value: " << normalize(max_t_all) << endl;
            cout << "Number of paths checked: " << edges << endl;
            cout << "----------------------------------------------------------------------------------" << endl;

            // TODO store full histogramm in csv
            auto csv_out_stream = ofstream(csv_out_path);
            csv_out_stream << "t_value,frequency" << endl;
            for (int i = 0; i < full_histogramm.size(); i++) {
                csv_out_stream << normalize(1.0 + i * 0.001) << "," << full_histogramm[i] << endl;
            }
            csv_out_stream.close();
            return full_worst_edges;
        }

double get_furthest_distance_in_graph(Graph* graph) {
    // identify the furthest distance in the graph
    std::vector<int> v;
    v.reserve(graph->n);
    for (int i = 1; i <= graph->n; ++i) v.push_back(i);
    vector<pair<vector<int>, double>> dists = graph->multiSourceMultiTargetDijkstra({0}, v, true);
    // get the node with the maximum distance
    double max_dist = 0.0;
    int max_node_src = -1;
    for (const auto& dist_pair : dists) {
        double dist = dist_pair.second;
        int node = dist_pair.first[dist_pair.first.size() - 1];
        if (dist > max_dist && dist != numeric_limits<double>::infinity()) {
            max_dist = dist;
            max_node_src = node;
        }
    }
    // repeat the process to get the true max distance
    dists = graph->multiSourceMultiTargetDijkstra({max_node_src}, v, true);
    max_dist = 0.0;
    int max_node_trgt = -1;
    for (const auto& dist_pair : dists) {
        double dist = dist_pair.second;
        int node = dist_pair.first[dist_pair.first.size() - 1];
        if (dist > max_dist && dist != numeric_limits<double>::infinity()) {
            max_dist = dist;
            max_node_trgt = node;
        }
    }
    return euk_dist(graph->id_point_map[max_node_src], graph->id_point_map[max_node_trgt]);
}


#include <atomic>
#include <random>

// ...

void analyse_long_random_paths_with_vis_graph(Graph base_graph, Graph spanner_graph,
    std::string csv_out_path, int max, int theta, double plimit) {

    const double limit = get_furthest_distance_in_graph(&base_graph) * plimit;
    std::cout << "Path search is limited to nodes further apart than: " << limit << std::endl;

    // Pre-size per-sample storage
    std::vector<double> t_values_sum(max, 0.0);
    std::vector<double> t_values_max(max, 0.0);
    std::vector<int>    edges_v(max, 0);
    std::vector<std::vector<std::pair<int64_t,int64_t>>> edge_log(max);
    std::vector<std::vector<int64_t>> node_log(max);
    std::vector<std::string> all_edges(max);

    const int64_t N = static_cast<int64_t>(base_graph.adj.size());
    const int num_threads = std::max(1, omp_get_num_procs() - 1);

    std::atomic<int> produced{0};

    #pragma omp parallel num_threads(num_threads) shared(produced)
    {
        // Thread-local RNG
        std::mt19937_64 rng(std::random_device{}() + std::hash<int>{}(omp_get_thread_num()));
        std::uniform_int_distribution<int64_t> pick(0, N - 1);

        while (true) {
            // Early check: if already reached, stop
            if (produced.load(std::memory_order_relaxed) >= max) break;

            // Draw a candidate pair
            int64_t source = pick(rng);
            int64_t target = pick(rng);
            if (source == target) continue;

            if (euk_dist(base_graph.id_point_map[source],
                         base_graph.id_point_map[target]) < limit) {
                continue; // too close -> retry
            }

            // Compute paths
            auto original = base_graph.dijkstra(source, target);
            const double original_dist = original.second;
            const auto& original_path  = original.first;

            auto spanner_res   = spanner_graph.dijkstra(source, target);
            const double spanner_dist = spanner_res.second;
            const auto& spanner_path  = spanner_res.first;

            if (spanner_dist == std::numeric_limits<double>::infinity()) continue;
            if (normalize(spanner_dist) < normalize(original_dist)) {
                std::cout << "ERROR: Spanner distance < original distance for (" << source
                          << "," << target << ")\n";
                continue;
            }

            double t_value = spanner_dist / original_dist;
            if (std::isnan(t_value)) {
                std::cout << "ERROR: t-value is NaN for (" << source << "," << target << ")\n";
                continue;
            }

            // Reserve a unique sample slot
            int slot = produced.fetch_add(1, std::memory_order_relaxed);
            if (slot >= max) break; // we raced past the cap

            // --- Write results to "slot" (unique; no data race) ---
            all_edges[slot] += std::to_string(source) + "," + std::to_string(target) + "," +
                               std::to_string(normalize(t_value)) + "," +
                               std::to_string(original_dist) + "\n";

            t_values_sum[slot] += t_value;
            edges_v[slot] += 1;
            t_values_max[slot] = std::max(t_values_max[slot], t_value);

            int64_t last = source;
            for (size_t k = 1; k < original_path.size(); ++k) {
                int64_t p = original_path[k];
                if (p < 0 || p >= N) continue;
                edge_log[slot].push_back({last, p});
                node_log[slot].push_back(p);
                last = p;
            }
            // -----------------------------------------------
        }
    } // end parallel region

    // Aggregate strictly over the slots we actually produced
    const int produced_count = std::min(max, (int)std::count_if(edges_v.begin(), edges_v.end(),
                                                                [](int e){ return e > 0; }));

    double t_values_all_sum = 0.0;
    int edges = 0;
    double max_t_all = 0.0;
    std::vector<std::pair<int,int>> all_logged_edges;
    std::vector<int> all_logged_nodes;

    for (int i = 0; i < max; ++i) {
        if (edges_v[i] == 0) continue;          // was this slot filled?
        t_values_all_sum += t_values_sum[i];
        edges += edges_v[i];
        max_t_all = std::max(max_t_all, t_values_max[i]);

        for (auto &e : edge_log[i])  all_logged_edges.push_back(e);
        for (auto &n : node_log[i])  all_logged_nodes.push_back(n);
        // if you want: edge_out_stream << all_edges[i];
    }

    if (csv_out_path != "-1") {
        auto edge_log_stream = std::ofstream("../../data/logged_edges-" + std::to_string(theta) + ".csv");
        for (auto &e : all_logged_edges) edge_log_stream << e.first << "," << e.second << "\n";
        edge_log_stream.close();

        auto node_log_stream = std::ofstream("../../data/logged_nodes-" + std::to_string(theta) + ".csv");
        for (auto n : all_logged_nodes) node_log_stream << n << "\n";
        node_log_stream.close();
    }

    std::cout << "----------------------------------------------------------------------------------\n";
    std::cout << "\nResults after analysing up to " << max
              << " paths of euclidean distance > " << limit << ":\n";
    std::cout << "plimit was set to: " << plimit << "\n";
    std::cout << "----------------------------------------------------------------------------------\n";
    std::cout << "Average t-value: " << normalize(edges ? (t_values_all_sum / edges) : 0.0) << "\n";
    std::cout << "Max t-value: " << normalize(max_t_all) << "\n";
    std::cout << "Number of valid paths tested: " << edges << " (samples: " << produced_count << ")\n";
    std::cout << "size of the spanner: " << spanner_graph.number_of_edges << "\n";
    std::cout << "----------------------------------------------------------------------------------\n";
}


void analyse_random_paths_with_vis_graph(Graph base_graph, Graph spanner_graph, std::string csv_out_path, int max, int theta) {
    ///////////////////////////////////////////////////////////////////////////////////
    /// TEST random queries on the spanners
    /// analyse t-value
    ///////////////////////////////////////////////////////////////////////////////////

    vector<double> t_values_sum = vector<double>(max, 0);
    vector<double> t_values_max = vector<double>(max, 0);
    vector<vector<pair<int64_t, int64_t>>> edge_log = vector<vector<pair<int64_t, int64_t>>>(); // log of (source, target) pairs per thread
    vector<vector<int64_t>> node_log = vector<vector<int64_t>>(); // log of nodes used per thread
    vector<int> edges_v = vector<int>(max, 0);
    vector<vector<int>> t_histogram = vector<vector<int>>(); // histogram with 20 bins
    vector<string> all_edges = vector<string>(max);
    for (int i = 0; i < max; i++) {
        t_histogram.push_back(vector<int>(100));
        edge_log.push_back(vector<pair<int64_t, int64_t>>());
        node_log.push_back(vector<int64_t>());
    }

    int num_threads = std::max(1,omp_get_num_procs()-1);  // Get the number of available processors
    //num_threads = 10; // TODO remove limit

    ProgressBar progressBar;

#pragma omp parallel for num_threads(num_threads) shared(spanner_graph, base_graph, t_values_sum, t_values_max) schedule(dynamic)
    for (int c= 0; c < max; c++) {

        int source = rand() % base_graph.adj.size();
        int target = rand() % base_graph.adj.size();
        if (source == target) {c--; continue;}

        auto original = base_graph.dijkstra(source, target);
        double original_dist = original.second;
        auto original_path = original.first;
        auto spanner_theta_dist = numeric_limits<double>::infinity();
        auto spanner_res  = spanner_graph.dijkstra(source, target);
        double spanner_dist  = spanner_res.second;
        auto spanner_path = spanner_res.first;

        // sanity checks
        if (normalize(spanner_dist) < normalize(original_dist)) {
            std::cout << "ERROR: Spanner distance is less than original distance for source " << source << " and target " << target << std::endl;
            std::cout << "Original distance: " << original_dist << ", Spanner distance: " << spanner_dist << std::endl;
        }
        if (spanner_dist == numeric_limits<double>::infinity()) {
            //std::cout << "ERROR: Spanner distance is infinity for source " << source << " and target " << target << std::endl;
            //std::cout << "Original distance: " << original_dist << ", Spanner distance: " << spanner_dist << std::endl;
            c--;
            continue;
        }

        double t_value = spanner_dist/original_dist;
        if (t_value != t_value) { // check for nan
            std::cout << "ERROR: t-value is greater than nan for source " << source << " and target " << target << std::endl;
            std::cout << "Original distance: " << original_dist << ", Spanner distance: " << spanner_dist << ", t-value: " << t_value << std::endl;
            continue;
        }
        all_edges[c] += to_string(source) + "," + to_string(target) + "," + to_string(normalize(t_value)) + "," + to_string(original_dist) + "\n";

        t_values_sum[c] += t_value;
        edges_v[c] += 1;
        if (t_value > t_values_max[source]) {
            t_values_max[c] = t_value;
        }

        int last = source;
        for (int count = 1; count<original_path.size(); count++) {
            if (original_path[count] < 0 || original_path[count] >= base_graph.adj.size()) {
                std::cerr << "ERROR: Spanner path contains invalid node index: " << original_path[count] << std::endl;
                continue;
            }
            auto p = original_path[count];
            edge_log[c].push_back({last, p});
            node_log[c].push_back(p);
            last = p;
        }


        }

            ////////////////////////////////////////////////////////////////////////////////////
            /// WRITE THE RESULTS
            ////////////////////////////////////////////////////////////////////////////////////

            double t_values_all_sum = 0.0;
            int edges = 0;
            double max_t_all = 0.0;
            vector<int> full_histogramm = vector<int>(100);
            auto edge_out_stream = ofstream("../../data/all_edges-" + to_string(theta)+".csv", ios_base::app);
            vector<pair<int, int>> all_logged_edges;
            vector<int> all_logged_nodes;

                for (int i = 0; i < t_values_sum.size(); i++) {
                    if (base_graph.adj[i].size() > 0) {
                        t_values_all_sum += t_values_sum[i];
                        edges += edges_v[i];
                        if (t_values_max[i] > max_t_all) {
                            max_t_all = t_values_max[i];
                        }
                        for (auto e: edge_log[i]) {
                            all_logged_edges.push_back(e);
                        }
                        for (auto n: node_log[i]) {
                            all_logged_nodes.push_back(n);
                        }
                        // write all edges to file
                        //edge_out_stream << all_edges[i];
                    }
                }

                if (csv_out_path != "-1") {
                auto edge_log_stream = ofstream("../../data/logged_edges-" + to_string(theta)+".csv");
                for (auto e: all_logged_edges) {
                    if (e.first < 0 || e.second < 0) {
                        //std::cerr << "ERROR: Logged edge contains node that is less than 0: " << e.first << ", " << e.second << std::endl;
                    }
                    edge_log_stream << e.first << "," << e.second << "\n";
                }
                edge_log_stream.close();
                auto node_log_stream = ofstream("../../data/logged_nodes-" + to_string(theta)+".csv");
                for (auto n: all_logged_nodes) {
                    if (n < 0) {
                        //std::cerr << "ERROR: Logged node is less than 0: " << n << std::endl;
                    }
                    node_log_stream << n << "\n";
                }
                node_log_stream.close();
            }


            cout << "----------------------------------------------------------------------------------" << endl;
            cout << endl << "Results after analysing " << max << " paths:" << endl;
            cout << "----------------------------------------------------------------------------------" << endl;
            cout << "Average t-value: " << normalize(t_values_all_sum/edges) << endl;
            cout << "Max t-value: " << normalize(max_t_all) << endl;
            cout << "Number of paths tested: " << edges << endl;
            cout << "size of the spanner: " << spanner_graph.number_of_edges << endl;
            cout << "----------------------------------------------------------------------------------" << endl;

        }


std::vector<Edge> analyse_spanner_with_coastline_graph(string base_graph_coastline_path, Graph spanner_graph, double t_cutoff) {
    ///////////////////////////////////////////////////////////////////////////////////
    /// TEST random queries on the spanners
    /// analyse t-value
    ///////////////////////////////////////////////////////////////////////////////////
    int num_threads = std::max(1,omp_get_num_procs()-1);  // Get the number of available processors
    // num_threads = 1;


    vector<double> t_values_sum = vector<double>(num_threads+1, 0);
    vector<double> t_values_max = vector<double>(num_threads+1, 0);
    vector<int> number_of_edges = vector<int>(num_threads+1, 0);
    vector<vector<Edge>> worst_edges = vector<vector<Edge>>(num_threads+1);
    vector<vector<double>> triang_times = vector<vector<double>>(num_threads+1);
    vector<vector<double>> dijkstr_times = vector<vector<double>>(num_threads+1);
    for (int i = 0; i < num_threads; i++) {
        worst_edges.push_back(vector<Edge>());
    }

    int max = 10000; //spanner_graph.adj.size();
    cout << "Prepare Triangulation" << endl;
    Triangulation triangulation;
    triangulation.readFromGraph(base_graph_coastline_path);
    cout << "Triangulation loaded" << endl;

    int percent_count = 0;
    std::cout << ">";
    std::cout.flush();

#pragma omp parallel for num_threads(num_threads) shared(spanner_graph, base_graph_coastline_path, t_values_sum, t_values_max, triangulation) schedule(dynamic)
    for (int source = 0; source < max; source++) {
        // get id of current threat
        auto thread_num = omp_get_thread_num();
        if (thread_num == 0) {
            if ((source/max)*100 > percent_count) {
            //if (source % 100 == 0) {
                percent_count++;
                //std::cout << "|";
                // 1. Get the current time
                auto now = std::chrono::system_clock::now();

                // 2. Convert it to a C-style time_t
                // cpp
                using namespace std::chrono;
                auto my_now = system_clock::now();
                auto ms  = duration_cast<milliseconds>(my_now.time_since_epoch()).count(); // milliseconds since epoch
                std::cout << "Progress: " << percent_count << "%: " << source <<"of " << max << "\t(" << ms << ")"<< std::endl;
                //std::cout.flush();
            }
        }

        auto time_tri_start = std::chrono::system_clock::now();
        std::vector<GlobalID> temp = triangulation.oneToAllVisibility(source, false);
        auto time_tri_end = std::chrono::system_clock::now();
        std::chrono::duration<double> elapsed_seconds = time_tri_end - time_tri_start;
        triang_times[thread_num].push_back(elapsed_seconds.count());
        std::vector<int> targets;
        targets.reserve(temp.size());
        for (auto id : temp) targets.push_back(static_cast<int>(id));

        if (targets.size() == 0) {
            //std::cerr << "Node " << source << " has no visible targets in the triangulation. Skipping." << std::endl;
            continue;
        }
        if (spanner_graph.adj[source].size() == 0) {
            //std::cerr << "Node " << source << " has no edges in the spanner graph. Skipping." << std::endl;
            continue;
        }

        // single source multiple targets dijkstra
        auto time_dijk_start = std::chrono::system_clock::now();
        auto res = spanner_graph.multiSourceMultiTargetDijkstra({source}, targets, true);
        auto time_dijk_end = std::chrono::system_clock::now();
        elapsed_seconds = time_dijk_end - time_dijk_start;
        dijkstr_times[thread_num].push_back(elapsed_seconds.count());

        // iterate over each target
        for (int i = 0; i < targets.size(); i++) {


            int target = targets[i];
            if (source > target) {
                continue; // to avoid double checks in undirected graph
            }

            double original_dist = euk_dist(spanner_graph.id_point_map.at(source), spanner_graph.id_point_map.at(target));
            auto spanner_theta_dist = numeric_limits<double>::infinity();

            auto path = res[i].first;
            // santiy check for path
            int path_target = path[path.size() - 1];
            if (path_target != target) {
                //std::cout << "ERROR: Path target does not match expected target for source " << source << " and target " << target << std::endl;
                continue;
            }
            double spanner_dist  = res[i].second;

            double t_value = spanner_dist/original_dist;

            t_values_sum[thread_num] += t_value;
            if (t_value > t_values_max[thread_num]) t_values_max[thread_num] = t_value;

            if (t_value > t_cutoff) {
                //worst_edges[thread_num].push_back(Edge(source, target, original_dist));
                //worst_edges[thread_num].back().t_value = t_value;
            }
            number_of_edges[thread_num]++;

        }
    }

            ////////////////////////////////////////////////////////////////////////////////////
            /// WRITE THE RESULTS
            ////////////////////////////////////////////////////////////////////////////////////

            double t_values_all_sum = 0.0;
            double t_values_max_all = 0.0;
            int full_number_of_edges = 0;
            int full_number_of_worst_edges = 0;

            double full_total_time_triang = 0.0;
            double full_total_time_dijkstra = 0.0;

            vector<Edge> all_worst_edges = vector<Edge>();

            for (int i=0; i<num_threads; i++) {
                t_values_all_sum += t_values_sum[i];
                if (t_values_max[i] > t_values_max_all) t_values_max_all = t_values_max[i];
                full_number_of_worst_edges += worst_edges[i].size();
                full_number_of_edges += number_of_edges[i];
                for (auto edge: worst_edges[i]) {
                    all_worst_edges.push_back(edge);
                }
                for (auto t: triang_times[i]) {
                    full_total_time_triang += t;
                }
                for (auto t: dijkstr_times[i]) {
                    full_total_time_dijkstra += t;
                }
            }

            cout << "AVG time for triangulation queries: " << full_total_time_triang/max << " seconds." << endl;
            cout << "AVG time for dijkstra queries: " << full_total_time_dijkstra/max << " seconds." << endl;

            cout << "----------------------------------------------------------------------------------" << endl;
            cout << endl << "Results after analysing all edges:" << endl;
            cout << "----------------------------------------------------------------------------------" << endl;
            std::cout << "The Spanner has now " << spanner_graph.number_of_edges << " edges." << std::endl;
            cout << "Average t-value: " << normalize(t_values_all_sum/full_number_of_edges) << endl;
            cout << "Max t-value: " << normalize(t_values_max_all) << endl;
            cout << "Number of edges checked: " << full_number_of_edges << endl;
            cout << " To reach a t-value cutoff of " << t_cutoff << ",one has to add additional " << full_number_of_worst_edges << " edges." << endl;
            cout << "----------------------------------------------------------------------------------" << endl;


            return all_worst_edges;
        }