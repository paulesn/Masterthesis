//
// Created by sebastian on 13.06.25.
//
#include "ThetaSpanner.hpp"
#include "Dataloader.hpp"
#include <cmath>
#include <csignal>

#include "Graph.hpp"
#include "Quadtree.hpp"
# include <omp.h>
#include <queue>
#include <vector>


double angle_between(Point a, Point b) {
    double x_a = a.x;
    double y_a = a.y;
    double x_b = b.x;
    double y_b = b.y;
    // using the point (x_b,y_a) as the  right angle corner, we can use the arctan (y_b - y_a, x_b - x_a) to calculate the angle
    double dx = x_b - x_a;
    double dy = y_b - y_a;
    double angle = std::atan2(dy, dx); // angle in radians
    if (angle < 0) {
        angle += 2 * M_PI; // normalize angle to [0, 2π]
    }
    return angle; // angle in radians
}

struct EdgeWeightCompare {
    bool operator()(const Edge& a, const Edge& b) const {
        return a.weight > b.weight; // min heap: smallest weight first
    }
};

Graph create_theta_spanner_graph(Graph* graph, const int theta) {
    Graph spanner (graph->id_point_map);
    int counter = 0;
    int fallback_counter = 0;
    std::cout << "Creating theta spanner graph with " << theta << " zones per node." << std::endl << ">";

    for (int node_id= 0; node_id < graph->adj.size(); node_id++) {
        counter++;

        // print progress
        if (graph->adj.size() > 100 && counter % (graph->adj.size()/100) == 0) {
            std::cout << "|";
            std::cout.flush();
        }

        Point source = graph->id_point_map[node_id];
        std::vector<bool> edges = std::vector(theta, false);
        std::vector<Edge> spanner_edges = std::vector<Edge>(theta);
        int edge_count = 0;


        for (Edge edge : graph->adj[node_id]) {
            Point target = graph->id_point_map[edge.target];
            // calculate angle between the edge and the x-axis
            double radians = angle_between(source, target);
            int zone = static_cast<int>(std::floor(radians/(2*M_PI/theta)));
            if (zone >= theta) zone = theta - 1;

            // if no edge is already in the zone, add it to the spanner
            if (!edges[zone]) {
                edges[zone] = true;
                spanner_edges[zone] = edge;
            }
            else if (edge.weight < spanner_edges[zone].weight) {
                // if the edge is shorter than the existing edge in the zone, replace it
                spanner_edges[zone] = edge;
            }
        }
        for (int i = 0; i < theta; i++) { // TODO warum hat jeder Knoten nur eine Edge
            if (edges[i]) {
                // add the edge to the spanner graph
                spanner.addEdge(node_id, spanner_edges[i].target, spanner_edges[i].weight);
                edge_count++;
            }
        }
    }
    // finalize progress output
    std::cout << std::endl;
    std::cout << "Created theta spanner graph with " << fallback_counter << " fallbacks." << std::endl;

    return spanner;
}

void dynamic_theta_update(Graph *graph, Graph* spanner, double t) {

    /*
     *In this variant of the dynamic theta update, we iterate over nodes instead of edges.
     *this potentially increases the number of edges in the spanner graph, but it is more efficient, because we can use a multiple target djikstra
     *
     */

    // min heap for the edges in the graph
    std::priority_queue<Edge, std::vector<Edge>, EdgeWeightCompare> minHeap;
    for (int i = 0; i < graph->adj.size(); i++) {
        for (const auto &edge : graph->adj[i]) {
            minHeap.push(edge);
        }
    }

    int counter = 0;
    std::cout << "Updating theta spanner graph with " << t << " zones." << std::endl <<">";
    int number_of_added_edges = 0;

    for (auto edge_list: graph->adj) {
        // print progress
        if (graph->adj.size() > 100 && counter % (graph->adj.size()/100) == 0) {
            std::cout << "|";
            std::cout.flush();
        }
        counter++;
        std::vector<int> targets = {};
        std::vector<int> source = {edge_list[0].source}; // get the source node id
        for (const auto &edge : edge_list) {
            targets.push_back(edge.target);
        }
        // calculate the shortest path from the source to all targets
        auto distances_graph = graph->multiSourceMultiTargetDijkstra(source, targets, true);
        auto distances_spanner = spanner->multiSourceMultiTargetDijkstra(source, targets, true);

        //validate the distances
        for (int i = 0; i< distances_graph.size(); i++) {
            double true_dist = distances_graph[i].second;
            double spanner_dist = std::numeric_limits<double>::infinity();
            int target = distances_graph[i].first.back();
            // TODO fix this is quadratic
            for (auto dist : distances_spanner) {
                if (dist.first.back() == target) {
                    spanner_dist = dist.second;
                    break;
                }
            }


            if (spanner_dist > t * true_dist) {
                // add the edge to the spanner graph
                spanner->addEdge(edge_list[0].source, distances_graph[i].first.back(), true_dist);
                number_of_added_edges++;
            }
        }
    }
    std::cout << std::endl;
    std::cout << "Updated theta spanner graph with " << number_of_added_edges << " new edges." << std::endl;
}

