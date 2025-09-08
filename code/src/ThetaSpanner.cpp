//
// Created by sebastian on 13.06.25.
//
#include "ThetaSpanner.hpp"
#include "Dataloader.hpp"
#include <cmath>
#include <csignal>
#include <algorithm>
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

void dynamic_theta_update(Graph *graph, Graph* spanner, const double t) {

    /*
     *In this variant of the dynamic theta update, we iterate over nodes instead of edges.
     *this potentially increases the number of edges in the spanner graph, but it is more efficient, because we can use a multiple target djikstra
     *
     */


    int counter = 0;
    std::cout << "Updating theta spanner graph with " << t << " zones." << std::endl <<">";
    int number_of_added_edges = 0;
    int edge_existing = 0;
    int edge_not_existing = 0;

    int num_threads = std::max(1,omp_get_num_procs()-1);  // Get the number of available processors
    num_threads = 10; // TODO remove

    #pragma omp parallel for num_threads(num_threads) shared(number_of_added_edges, edge_existing, edge_not_existing, counter, graph, spanner, t) schedule(dynamic)
    for (auto edge_list: spanner->adj) {
        // print progress
        if (graph->adj.size() > 100 && counter % (graph->adj.size()/100) == 0) {
            std::cout << "|";
            std::cout.flush();
        }
        counter++;


        // Assuming you have access to the source Point for the edges
        std::sort(edge_list.begin(), edge_list.end(), [&](const Edge& a, const Edge& b) {
            Point source = spanner->id_point_map[a.source];
            double angle_a = angle_between(source, graph->id_point_map[a.target]);
            double angle_b = angle_between(source, graph->id_point_map[b.target]);
            return angle_a < angle_b;
        });

        // for each edge in the edge list, we identify the edge in the next zone
        for (int i = 0; i < edge_list.size(); i++) {
            auto edge = edge_list[i];
            auto edge2 = edge_list[(i + 1) % edge_list.size()]; // next edge in the list, wraps around
            if (edge.target == edge2.target) continue; // skip self-comparison

            // then we identify the distance between the ends of the two edges
            // check if the distance between the two edges is smaller than the two edges times t
            auto p1 = graph->id_point_map[edge.target];
            auto p2 = graph->id_point_map[edge2.target];
            if (euklidian_distance(p1, p2) > t*(edge.weight + edge2.weight)) {
                // if this distance is smaller than the two edges times the value t, we add the direct edge to the spanner graph
                if (!spanner->addEdge(edge.source, edge2.target, edge.weight + edge2.weight)) {
                    // if the edge already exists, we skip it
                    edge_existing++;
                } else {
                    edge_not_existing++;
                }
                number_of_added_edges++;
            }
        }
    }
    std::cout << std::endl;
    std::cout << "Updated theta spanner graph with " << number_of_added_edges << " new edges." << std::endl;
    std::cout << "Existing edges: " << edge_existing << ", Not existing edges: " << edge_not_existing << std::endl;
    std::cout << "Total edges in spanner graph: " << spanner->number_of_edges << std::endl;
}

