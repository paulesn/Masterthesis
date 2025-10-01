#include <fstream>
#include <cmath>
#include <stdint.h>

#include "../io/Dataloader.hpp"
#include <omp.h>
#include "../daniel/theta-graph/headers/Progressbar.h"
#include "../daniel/theta-graph/headers/Triangulation.h"
#include "../spanner/ThetaSpanner.hpp"
#include "MultiThetaAnalysis.h"
#include "../daniel/theta-graph/headers/Timer.h"


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

int analyse_spanner(Graph base_graph, Graph spanner_graph, string csv_out_path, string graph_path) {
    ///////////////////////////////////////////////////////////////////////////////////
    /// TEST random queries on the spanners
    /// analyse t-value
    ///////////////////////////////////////////////////////////////////////////////////

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
    //num_threads = 10; // TODO remove limit

    Timer timer;
    timer.start();
    ProgressBar progressBar;
    int max = 1000; //base_graph.adj.size(); // TODO remove limit
    progressBar.start(max);

#pragma omp parallel for num_threads(num_threads) shared(spanner_graph, base_graph, t_values_sum, t_values_max, triangulation) schedule(dynamic)
    for (int source = 0; source < max; source++) {
        progressBar.update(1);
        if (base_graph.adj[source].size() == 0) {
            std::cout << "Node " << source << " has no edges in the base graph. Skipping." << std::endl;
            continue;
        }

        std::vector<GlobalID> targets = triangulation.oneToAllVisibility(source, false);

        if (targets.size() == 0) {
            std::cerr << "Node " << source << " has no visible targets in the triangulation. Skipping." << std::endl;
            continue;
        }
        if (spanner_graph.adj[source].size() == 0) {
            std::cerr << "Node " << source << " has no edges in the spanner graph. Skipping." << std::endl;
            continue;
        }
        // iterate over each target
        for (int i = 0; i < targets.size(); i++) {

            int target = targets[i];
            if (source > target) {
                continue; // to avoid double checks in undirected graph
            }

            double original_dist = euk_dist(base_graph.id_point_map.at(source), base_graph.id_point_map.at(target));
            auto spanner_theta_dist = numeric_limits<double>::infinity();
            double spanner_dist  = spanner_graph.dijkstra(source, target).second;

            double t_value = spanner_dist/original_dist;


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
            std::cout << "Node " << source << " has no edges in the base graph. Skipping." << std::endl;
            continue;
        }

        std::vector<Edge> edges =  base_graph.adj[source];

        if (edges.size() == 0) {
            std::cerr << "Node " << source << " has no visible targets in the triangulation. Skipping." << std::endl;
            continue;
        }
        if (spanner_graph.adj[source].size() == 0) {
            std::cerr << "Node " << source << " has no edges in the spanner graph. Skipping." << std::endl;
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
                std::cerr << "Spanner distance is less than original distance for edge " << source << " - " << target << std::endl;
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
            cout << endl << "Results after analysing all paths:" << endl;
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
            //edge_out_stream << "source,target,t_value,original_distance\n";

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
            edge_out_stream.close();
            // edge and node log to file

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


            cout << "----------------------------------------------------------------------------------" << endl;
            cout << endl << "Results after analysing all edges:" << endl;
            cout << "----------------------------------------------------------------------------------" << endl;
            cout << "Average t-value: " << normalize(t_values_all_sum/edges) << endl;
            cout << "Max t-value: " << normalize(max_t_all) << endl;
            cout << "Number of edges: " << edges << endl;
            cout << "----------------------------------------------------------------------------------" << endl;

            // TODO store full histogramm in csv
            //auto csv_out_stream = ofstream(csv_out_path);
            //csv_out_stream << "t_value,frequency" << endl;
            //for (int i = 0; i < full_histogramm.size(); i++) {
            //    csv_out_stream << normalize(1.0 + i * 0.001) << "," << full_histogramm[i] << endl;
            //}
            //csv_out_stream.close();
        }