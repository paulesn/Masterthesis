//
// Created by sebastian on 13.06.25.
//
#include "ThetaSpanner.hpp"
#include "Dataloader.hpp"
#include <cmath>
#include "Graph.hpp"

Graph create_theta_spanner_graph(Graph* graph, const int theta) {
    Graph spanner (graph->n);
    for (int node_id= 0; node_id < graph->adj.size(); node_id++) {
        std::vector<Edge*> edges = {};
        for (int i = 0; i < theta; ++i) {
            edges.emplace(edges.begin()+i,nullptr); // Initialize the map with -1 to indicate no edge in that zone
        }

        for (auto edge : graph->adj[node_id]) {
            // calculate angle between the edge and the x-axis
            //assuming the x asis is angle=0
            double sx = graph->id_point_map[node_id].x;
            double sy = graph->id_point_map[node_id].y;
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
    }
    return spanner;
}

