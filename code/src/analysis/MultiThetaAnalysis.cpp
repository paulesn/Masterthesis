#include <fstream>

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
    int max = base_graph.adj.size(); // TODO remove limit
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



/**
 *int main(int argc, char* argv[]) {

    string base_graph_path = "../../data/0025.32.fmi";
    string base_csv_path = "../../data/0025.32";
    vector<int> thetas = {24,50,64,75,100,124};

    ///////////////////////////////////////////////////////////////////////////////////
    /// LOAD THE GRAPH
    ///////////////////////////////////////////////////////////////////////////////////
    auto base_graph = get<1>(load_fmi(base_graph_path,  -1));

    for (int theta: thetas) {
        cout << "Analysing theta = " << theta << endl;
        auto spanner_graph = create_theta_spanner_graph(&base_graph, theta);
        analyse_spanner(base_graph, spanner_graph, base_csv_path + "_theta_" + to_string(theta) + ".csv");
    }


}*/