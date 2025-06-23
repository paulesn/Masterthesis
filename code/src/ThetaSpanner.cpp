//
// Created by sebastian on 13.06.25.
//
#include "ThetaSpanner.hpp"
#include "Dataloader.hpp"
#include <cmath>
#include "Graph.hpp"
#include "Quadtree.hpp"

/**
 * This function creates a theta graph based on a input graph. The resulting graph only contains edges in the oroiginal graph
 * @param graph the original graph the theta spanner should be based on
 * @param theta the number of zones around each node for edges
 * @param spd_fallback if true, zones that contain no nodes identify the closest neighbor and add the dijkstra route to them
 * @param quadtree a quadtree that contains all points of the graph. Only used if spd_fallback is true.
 * @return a graph
 */
Graph create_theta_spanner_graph(Graph* graph, const int theta, bool spd_fallback, Quadtree* quadtree) {
    Graph spanner (graph->n);
    for (int node_id= 0; node_id < graph->adj.size(); node_id++) {
        std::vector<Edge*> edges = {};
        for (int i = 0; i < theta; ++i) {
            edges.emplace(edges.begin()+i,nullptr); // Initialize the map with -1 to indicate no edge in that zone
        }

        double sx = graph->id_point_map[node_id].x;
        double sy = graph->id_point_map[node_id].y;

        for (auto edge : graph->adj[node_id]) {
            // calculate angle between the edge and the x-axis
            //assuming the x asis is angle=0
            double tx = graph->id_point_map[edge.target].x;
            double ty = graph->id_point_map[edge.target].y;

            double angle = std::atan2(ty - sy, tx - sx);

            // identify the zone of the edge
            // Vector from (tx, ty) to (sx, sy)
            double dx = sx - tx;
            double dy = sy - ty;
            // Angle in radians
            double angleRadians = std::atan2(dy, dx);
            // Optionally convert to degrees
            double angleDegrees = angleRadians * (180.0 / M_PI);
            int zone = static_cast<int>(std::floor((angleDegrees + 180.0) / (360.0 / theta))) % theta;

            // if no edge is already in the zone, add it to the spanner
            if (edges[zone] == nullptr) {
                edges[zone] = &edge;
                spanner.addEdge(edge.source, edge.target, edge.weight, false); // Add the edge to the spanner graph
                bool full =  true;
                for (int i = 0; i < theta; ++i) {
                    if (edges[i] == nullptr) {
                        full = false;
                        break;
                    }
                }
                if (full) {
                    // if all zones are filled, break the loop
                    break;
                }
            }
        }
        if (!spd_fallback) return spanner;
        ////////////////////////////////////////////////////////////////////////////////////
        /// Shortest Path Distance Fallback
        ////////////////////////////////////////////////////////////////////////////////////
        // it is possible that the theta graph creates multiple connected components if it only relies on existing edges

        for (int e = 0; e < theta; ++e) {
            if (edges[e] == nullptr) {
                // this zone has no edge. We identify the closest node and calculate the dijkstra path. then we add a direct edge to that node
                int min_node = -1;
                double min_distance = std::numeric_limits<double>::infinity();
                //auto potentials = quadtree->angle_intersect(sx, sy, (e * 2 * M_PI / theta), ((e + 1) * 2 * M_PI / theta));
                // TODO brute force implementation REPLACE
                std::vector<int> potentials = {};
                for (int i = 0; i < graph->id_point_map.size(); i++) {
                    double angle = std::atan2(graph->id_point_map[i].y - sy, graph->id_point_map[i].x - sx);
                    int this_zone = static_cast<int>(std::floor((angle+ 180.0) / (360.0 / theta))) % theta;
                    if (this_zone == e) {
                        potentials.push_back(i);
                    }
                }
                for (auto node: potentials) {
                    // check if the node is closer than the current minimum
                    double distance = euklidian_distance(graph->id_point_map[node], graph->id_point_map[node_id]);
                    if (distance < min_distance) {
                        min_distance = distance;
                        min_node = node;
                    }
                }
                if (min_node != -1) {
                    // we have found a node in the zone, now we calculate the dijkstra path to it
                    auto dijkstra_result = graph->dijkstra(node_id, min_node);
                    if (dijkstra_result.first.empty()) {
                        // no path found, skip this zone
                        continue;
                    }
                    double dijkstra_distance = dijkstra_result.second;
                    // add the edge to the spanner graph
                    spanner.addEdge(node_id, min_node, dijkstra_distance, true);

                    // TODO REMOVE sanity check
                    double tx = graph->id_point_map[min_node].x;
                    double ty = graph->id_point_map[min_node].y;
                    double dx = sx - tx;
                    double dy = sy - ty;
                    // Angle in radians
                    double angleRadians = std::atan2(dy, dx);
                    // Optionally convert to degrees
                    double angleDegrees = angleRadians * (180.0 / M_PI);
                    int zone = static_cast<int>(std::floor((angleDegrees + 180.0) / (360.0 / theta))) % theta;
                    if (zone != e) {
                        std::cout << "Error: Zone mismatch! Expected zone " << e << ", but found zone " << zone << " for node " << min_node << std::endl;
                    }

                } else {
                    // no node found in the zone, skip this zone
                    continue;
                }
            }
        }
    }
    return spanner;
}

