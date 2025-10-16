//
// Created by Sebastian Paule on 10/1/25.
//
/**
 *
 *This file was created to analyse the difference between the t value of a theta graph for a point cloud and a visibilit graph
 *
 *
 *
 */

#include "../structure/Graph.hpp"
#include "../io/Dataloader.hpp"
#include "../analysis/MultiThetaAnalysis.h"
#include "../spanner/ThetaSpanner.hpp"

using namespace std;

int main() {

    string base_graph_path = "../../data/0025.32.fmi";
    string csv_path = "-1";
    int paths = 50000;
    vector<int> thetas = {12,24,32,50,64,75,100,124,128};

    ///////////////////////////////////////////////////////////////////////////////////
    /// LOAD THE GRAPHS
    ///////////////////////////////////////////////////////////////////////////////////
    auto base_graph = get<1>(load_fmi(base_graph_path, -1));
    auto base_point_cloud = get<1>(load_fmi(base_graph_path, -1));


    for (int s = 0; s < base_point_cloud.adj.size(); s++) {
        if (s % 1000 == 0) cout << "Checking point " << s << endl;
        base_point_cloud.adj[s].reserve(base_point_cloud.adj.size());
        for (int t = 0; t < s; t++) {
            Pointc sp = base_point_cloud.id_point_map[s];
            Pointc tp = base_point_cloud.id_point_map[t];
            base_point_cloud.addEdge(s,t,(sp.x-tp.x)*(sp.x-tp.x)+(sp.y-tp.y)*(sp.y-tp.y), false);
        }
    }
    cout << "Created Point Cloud" << endl;
    base_point_cloud.sort_edges();
    cout << "Sorted Point Cloud" << endl;


    for (auto theta: thetas) {
        cout << "Analysing theta: " << theta << endl;
        auto base_point_spanner = create_theta_spanner_graph(&base_graph, theta, true); //get<1>(load_fmi(spanner_path,  -1));
        analyse_random_paths_with_vis_graph(base_point_cloud, base_point_spanner, csv_path, paths, theta);

        auto spanner_graph = create_theta_spanner_graph(&base_graph, theta); //get<1>(load_fmi(spanner_path,  -1));
        analyse_random_paths_with_vis_graph(base_graph, spanner_graph, csv_path, paths, theta);
    }


}

