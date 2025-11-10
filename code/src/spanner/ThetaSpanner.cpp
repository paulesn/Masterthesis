//
// Created by sebastian on 13.06.25.
//
#include "ThetaSpanner.hpp"
#include "../io/Dataloader.hpp"
#include <cmath>
#include <csignal>
#include <algorithm>
#include "../structure/Graph.hpp"
#include "../structure/Quadtree.hpp"
# include <omp.h>
#include <queue>
#include <vector>
#include <vector>
#include <boost/fusion/container/vector/vector.hpp>

#include "../daniel/theta-graph/headers/Triangulation.h"


double angle_between(Pointc a, Pointc b) {
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

Graph create_theta_spanner_graph(Graph* graph, const int theta, bool early_break) {
    Graph spanner (graph->id_point_map);
    int counter = 0;
    int fallback_counter = 0;
    std::cout << "Creating theta spanner graph with " << theta << " zones per node." << std::endl << ">";
    if (early_break) std::cerr << "ATTENTION: early_break is active. This fails to create spanner if the adj list is not sorted by edge length!" << std::endl;

    for (int node_id= 0; node_id < graph->adj.size(); node_id++) {
        counter++;

        // print progress
        if (graph->adj.size() > 100 && counter % (graph->adj.size()/100) == 0) {
            std::cout << "|";
            std::cout.flush();
        }

        Pointc source = graph->id_point_map[node_id];
        std::vector<bool> edges = std::vector(theta, false);
        std::vector<Edge> spanner_edges = std::vector<Edge>(theta);
        int edge_count = 0;
        int filled_zones = 0;


        for (Edge edge : graph->adj[node_id]) {
            Pointc target = graph->id_point_map[edge.target];
            // calculate angle between the edge and the x-axis
            double radians = angle_between(source, target);
            int zone = static_cast<int>(std::floor(radians/(2*M_PI/theta)));
            if (zone >= theta) {
                zone = theta - 1;
            }

            // if no edge is already in the zone, add it to the spanner
            if (!edges[zone]) {
                edges[zone] = true;
                spanner_edges[zone] = edge;
                filled_zones++;
                if (early_break && filled_zones == theta) break; // if all zones are filled, we can stop
            }
            else if (edge.weight < spanner_edges[zone].weight) {
                // if the edge is shorter than the existing edge in the zone, replace it
                spanner_edges[zone] = edge;
            }
        }
        for (int i = 0; i < theta; i++) {
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
            Pointc source = spanner->id_point_map[a.source];
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
            std::cerr << "I REMOVED STUFF HERE. THE CODE DOES NOT WORK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
            /**if (euklidian_distance(p1, p2) > t*(edge.weight + edge2.weight)) {
                // if this distance is smaller than the two edges times the value t, we add the direct edge to the spanner graph
                if (!spanner->addEdge(edge.source, edge2.target, edge.weight + edge2.weight)) {
                    // if the edge already exists, we skip it
                    edge_existing++;
                } else {
                    edge_not_existing++;
                }
                number_of_added_edges++;
            }*/
        }
    }
    std::cout << std::endl;
    std::cout << "The Spanner has now " << spanner->number_of_edges << " edges." << std::endl;
    std::cout << "Updated theta spanner graph with " << number_of_added_edges << " new edges." << std::endl;
    std::cout << "Existing edges: " << edge_existing << ", Not existing edges: " << edge_not_existing << std::endl;
    std::cout << "Total edges in spanner graph: " << spanner->number_of_edges << std::endl;
}

std::vector<Edge> get_nearest_zone_between_angles(Graph* graph, int source_node, double start_angle, double end_angle, double max_an, bool early_stop) {

    //std::cout << "Searching for edge between angles " << start_angle << " and " << end_angle << " for node " << source_node << std::endl;

    if (start_angle< 0 || end_angle < 0 ) return {};

    std::vector<Edge> result = {};

    bool angle_found=false;
    double dist = -1;
    double best_angle = -1;
    Edge best_edge = {};

    for (Edge edge : graph->adj[source_node]) {
        // calculate the angle between the edge and the x-axis
        Pointc source = graph->id_point_map[source_node];
        Pointc target = graph->id_point_map[edge.target];
        double radians = angle_between(source, target);
        if (radians > 2 * M_PI) {
            radians -= 2 * M_PI;
        }
        // check if the angle is between the two angles
        if (radians > start_angle && radians < end_angle) {
            // check if this is the first edge found or if it is shorter than the current best edge
            if (!angle_found || edge.weight < dist) {
                dist = edge.weight;
                best_edge = edge;
                angle_found = true;
                double best_angle = radians;
                if (early_stop) break;
            }
        }
    }
    if (!angle_found) {
        //std::cerr << "WARNING: No edge found between angles " << start_angle << " and " << end_angle << " for node " << source_node << std::endl;
        return {};
    } else {
        //std::cout << "Found edge to node " << best_edge.target << " with angle " << best_angle << " and weight " << best_edge.weight << std::endl;
    }

    // check if the angle needs to be reduced recursivly
    double angle_to_start = std::abs(best_angle - start_angle);
    double angle_to_end = std::abs(best_angle - end_angle);
    if (angle_to_start > max_an) {
        auto a = get_nearest_zone_between_angles(graph, source_node, start_angle, best_angle, max_an, early_stop);
        for (auto e : a) {
            result.push_back(e);
        }
    }
    if (angle_to_end > max_an) {
        auto a = get_nearest_zone_between_angles(graph, source_node, best_angle, start_angle, max_an, early_stop);
        for (auto e : a) {
            result.push_back(e);
        }
    }
    return result;
}

Graph create_theta_spanner_graph_with_max_angle(Graph* graph, const int k, double max_an) {
    Graph spanner (graph->id_point_map);
    int counter = 0;
    int fallback_counter = 0;
    std::cout << "Creating theta spanner graph with " << k << " zones per node." << std::endl << ">";
    std::cout << "ATTENTION: early_break is active. This fails to create spanner if the adj list is not sorted by edge length!" << std::endl;

    for (int node_id= 0; node_id < graph->adj.size(); node_id++) {
        counter++;

        // print progress
        if (graph->adj.size() > 100 && counter % (graph->adj.size()/100) == 0) {
            std::cout << "|";
            std::cout.flush();
        }

        Pointc source = graph->id_point_map[node_id];
        std::vector<std::tuple<double, double>> empty_zones = {}; // start_angle, end_angle
        empty_zones.reserve(k);
        for (int i = 0; i < k; i++) {
            double zone_start_angle = i * (2*M_PI/k);
            double zone_end_angle = (i + 1) * (2*M_PI/k);
            empty_zones.push_back(std::make_tuple(zone_start_angle, zone_end_angle));
        }
        std::vector<Edge> spanner_edges = std::vector<Edge>();
        spanner_edges.reserve(k);
        int edge_count = 0;


        for (Edge edge : graph->adj[node_id]) {
            Pointc target = graph->id_point_map[edge.target];
            // calculate angle between the edge and the x-axis
            double radians = angle_between(source, target);
            for (int i = 0; i < empty_zones.size(); i++) {
                auto [start_angle, end_angle] = empty_zones[i];
                if (radians >= start_angle && radians < end_angle) {
                    empty_zones.erase(empty_zones.begin() + i);
                    //std::cout << edge.source << " -> " << edge.target << std::endl;
                    Edge temp = {edge.source+0, edge.target+0, edge.weight+0};
                    spanner_edges.push_back(temp);
                    // check if the zone requires an aditional edge due to the max angle constraint
                    double angle_to_start = std::abs(radians - start_angle);
                    double angle_to_end = std::abs(end_angle - radians);
                    if (angle_to_start > max_an) {

                    }
                    if (angle_to_end > max_an) {
                        empty_zones.push_back({radians, end_angle});
                    }
                    break;
                }
            }
        }
        for (int i = 0; i < spanner_edges.size(); i++) {
                // add the edge to the spanner graph
                //std::cout << spanner_edges[i].source << " -> " << spanner_edges[i].target << std::endl;
                spanner.addEdge(node_id, spanner_edges[i].target, spanner_edges[i].weight);
                edge_count++;
        }
    }
    // finalize progress output
    std::cout << std::endl;
    std::cout << "Created theta spanner graph with " << fallback_counter << " fallbacks." << std::endl;

    return spanner;
}

void enforce_small_angle_constraint(std::string base_graph_path, Graph* spanner, double min_an) {
    int num_threads = std::max(1,omp_get_num_procs()-1);  // Get the number of available processors
    std::vector<std::vector<Edge>> added_edges = std::vector<std::vector<Edge>>(num_threads);

    Triangulation triangulation;
    triangulation.readFromGraph(base_graph_path);


    for (int node_id = 0; node_id < spanner->adj.size(); node_id++) {
        std::vector<std::tuple<Edge*, double>> local_edges;
        for (int i = 0; i < spanner->adj[node_id].size(); i++) {
            auto edge = spanner->adj[node_id][i];
            Pointc source = spanner->id_point_map[edge.source];
            Pointc target = spanner->id_point_map[edge.target];
            double angle = angle_between(source, target);
            local_edges.push_back({&spanner->adj[node_id][i], angle});
        }
        // sort edges by angle
        std::sort(local_edges.begin(), local_edges.end(), [](const std::tuple<Edge*, double>& a, const std::tuple<Edge*, double>& b) {
            return std::get<1>(a) < std::get<1>(b);
        });
        // check angles between edges
        for (int i = 0; i < local_edges.size(); i++) {
            auto [edge_a, angle_a] = local_edges[i];
            auto [edge_b, angle_b] = local_edges[(i + 1) % local_edges.size()]; // next edge, wraps around
            double angle_diff = angle_b - angle_a;
            if (angle_diff < min_an) {
                InternalID ip1 = edge_a->target;
                InternalID ip2 = edge_b->target;
                if (triangulation.directVisibility(ip1, ip2)) {
                    // need to add an edge between edge_a and edge_b
                    Pointc p1 = spanner->id_point_map[edge_a->target];
                    Pointc p2 = spanner->id_point_map[edge_b->target];
                    Edge new_edge = {
                        edge_a->target,
                        edge_b->target,
                        sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y))
                    };
                    int thread_id = omp_get_thread_num();
                    added_edges[thread_id].push_back(new_edge);
                }
            }
        }
    }

    int c = 0;
    for (int i = 0; i < added_edges.size()-1; i++) {
        auto e = added_edges[i];
        for (auto edge: e) {
            if (spanner->addEdge(edge.source, edge.target, edge.weight)) c++;
        }
    }
    std::cout << "Enforcing angle constraint added " << c << " Edges to the spanner" << std::endl;
}


void enforce_small_angle_constraint(Graph* base_graph, Graph* spanner, double min_an) {
    int num_threads = std::max(1,omp_get_num_procs()-1);  // Get the number of available processors
    std::vector<std::vector<Edge>> added_edges = std::vector<std::vector<Edge>>(num_threads);

    for (int node_id = 0; node_id < spanner->adj.size(); node_id++) {
        std::vector<std::tuple<Edge*, double>> local_edges;
        for (int i = 0; i < spanner->adj[node_id].size(); i++) {
            auto edge = spanner->adj[node_id][i];
            Pointc source = spanner->id_point_map[edge.source];
            Pointc target = spanner->id_point_map[edge.target];
            double angle = angle_between(source, target);
            local_edges.push_back({&spanner->adj[node_id][i], angle});
        }
        // sort edges by angle
        std::sort(local_edges.begin(), local_edges.end(), [](const std::tuple<Edge*, double>& a, const std::tuple<Edge*, double>& b) {
            return std::get<1>(a) < std::get<1>(b);
        });
        // check angles between edges
        for (int i = 0; i < local_edges.size(); i++) {
            auto [edge_a, angle_a] = local_edges[i];
            auto [edge_b, angle_b] = local_edges[(i + 1) % local_edges.size()]; // next edge, wraps around
            double angle_diff = angle_b - angle_a;
            if (i == local_edges.size() - 1) {
                angle_diff += 2*M_PI; // wrap around case
            }
            if (angle_diff < min_an) {
                InternalID ip1 = edge_a->target;
                InternalID ip2 = edge_b->target;
                bool has_edge = false;
                double correct_weight = 0.0; // Create a variable to store the weight

                // Loop through the base_graph to find the edge AND its weight
                for (int i=0; i < base_graph->adj[ip1].size(); i++) {
                    if (base_graph->adj[ip1][i].target == ip2) {
                        has_edge = true;
                        correct_weight = base_graph->adj[ip1][i].weight; // <-- Get the correct weight
                        break;
                    }
                }
                if (has_edge) {
                    // need to add an edge between edge_a and edge_b
                    Pointc p1 = spanner->id_point_map[edge_a->target];
                    Pointc p2 = spanner->id_point_map[edge_b->target];
                    Edge new_edge = {
                        edge_a->target,
                        edge_b->target,
                        correct_weight
                    };
                    int thread_id = omp_get_thread_num();
                    added_edges[thread_id].push_back(new_edge);
                }
            }
        }
    }

    int c = 0;
    for (int i = 0; i < added_edges.size()-1; i++) {
        auto e = added_edges[i];
        for (auto edge: e) {
            if (spanner->addEdge(edge.source, edge.target, edge.weight)) c++;
        }
    }
    std::cout << "Enforcing angle constraint added " << c << " Edges to the spanner" << std::endl;
}




