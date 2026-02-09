//
// Created by Sebastian Paule on 12/16/25.
//

#include "Lowerbound_Graph.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <limits>

/**
 * A utility function to split a string by a delimiter.
 * @param s
 * @param delimiter
 * @return
 */
std::vector<std::string> split(const std::string& s, const std::string& delimiter) {
    std::vector<std::string> tokens;
    size_t start = 0;
    size_t end = 0;
    while ((end = s.find(delimiter, start)) != std::string::npos) {
        tokens.push_back(s.substr(start, end - start));
        start = end + delimiter.length();
    }
    tokens.push_back(s.substr(start));
    return tokens;
}

double euk_dist(Pointc p1, Pointc p2) {
    return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

/**
 * This function calculates the angle between the line connecting points a and b and the x-axis.
 * The angle is returned in radians in the range [0, 2π].
 * @param a
 * @param b
 * @return angle in radians
 */
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

/**
* This function loads a graph in the fmi format.
* It is mostly copied from the Dataloader.cpp file to make this class selfcontained
* @param filepath the path to the fmi file
* @return
*/
bool Lowerbound_Graph::load_graph_from_fmi_file(const std::string filepath) {
    std::cout << "Loading fmi format from " << filepath << std::endl;
    std::vector<Pointc> systems; // id, x, y, z, jump multiplicator

    // print current working directory
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filepath << std::endl;
        return false;
    }

    int counter = 0;

    int number_of_nodes = 0;
    long number_of_edges = 0;
    int technical_counter = -2;
    Graph g = Graph(1);
    bool graph_created = false;

    for (std::string line; std::getline(file, line);) {
        // update counter
        counter++;

        // Skip empty lines and comments
        if (line.empty() || line[0] == '#') {
            continue; // Skip empty lines
        }
        if (technical_counter == -2) {
            // First line contains the number of nodes and edges
            number_of_nodes = std::stoi(line);
            technical_counter++;
            continue;
        }
        if (technical_counter == -1) {
            // Second line contains the number of edges
            number_of_edges = std::stol(line);
            technical_counter++;
            std::cout << "Loading Nodes \n>";
            continue;
        }
        if (technical_counter < number_of_nodes) {
            if (number_of_nodes > 100 && technical_counter != 0 && counter % (number_of_nodes/100) == 0) {
                std::cout << "|";
                std::cout.flush();
            }

            // Read node data
            technical_counter++;
            std::string delimiter = " ";

            auto arr = split(line, delimiter);
            // this creates the point. the data is in the format: id nonesense x y ...
            systems.emplace_back(std::stoi(arr[2]), std::stod(arr[3]), std::stod(arr[0]), std::stod(arr[4]));
        } else {
            if (technical_counter == number_of_nodes) {
                // We have read all nodes
                std::cout << std::endl << "Loading Edges\n>";
            }
            technical_counter++;
            if (technical_counter >= number_of_nodes + number_of_edges) {
                // We have read all nodes and edges
                break;
            }
            if (number_of_nodes > 100 && number_of_edges > 100 && (counter-number_of_nodes) % (number_of_edges/100) == 0) {
                std::cout << "|";
                std::cout.flush();
            }
            // Read edge data
            std::string delimiter = " ";
            auto arr = split(line, delimiter);
            if (!graph_created) {
                g = Graph(systems);
                graph_created = true;
            }
            g.addEdge(std::stoi(arr[0]), std::stoi(arr[1]), std::stod(arr[2]), false);
        }
    }
    graph = g;
    std::cout << std::endl << "Graph loaded with " << graph.n << " nodes and " << graph.number_of_edges << " edges." << std::endl;
    return true;
}

double Lowerbound_Graph::get_lower_bound_path(int source, int target, int k) {
    auto result = graph.dijkstra(source, target);
    double dist = result.second;
    auto path = result.first;
    auto p_s = graph.id_point_map[source];
    auto p_t = graph.id_point_map[target];
    if (p_s.meta != 1) {
        // the meta value indicates a spanner graph point
        // reduce the total distance by 2x the first edge
        double cone_offset = identify_cone_offset(p_s, k);
        double closest_neighbor_offset = identify_closest_cone_neighbor_offset(p_s);
        dist -= std::min(cone_offset, closest_neighbor_offset);
    }
    if (p_t.meta != 1) {
        // the meta value indicates a spanner graph point
        // reduce the total distance by 2x the first edge
        double cone_offset = identify_cone_offset(p_t, k);
        double closest_neighbor_offset = identify_closest_cone_neighbor_offset(p_t);
        dist -= std::min(cone_offset, closest_neighbor_offset);
    }
    return std::max(dist,0.0);
}

double Lowerbound_Graph::identify_closest_cone_neighbor_offset(Pointc p_t) {
    double shortest_edge_dist = std::numeric_limits<double>::infinity();
    Edge shortest_edge;
    for (Edge edge : graph.adj[p_t.id]) {
        if (edge.weight < shortest_edge_dist) {
            shortest_edge_dist = edge.weight;
            shortest_edge = edge;
        }
    }
    return 2*shortest_edge_dist;
}

double Lowerbound_Graph::identify_cone_offset(Pointc p, int k) {
    double cone_size = 2 * M_PI / k;
    double offset = 0.0001;
    double max_ret = 0.0;
    for (Edge edge : graph.adj[p.id]) {
        double angle = angle_between(p, graph.id_point_map[edge.target])+offset;
        int cone_low = std::floor(angle/cone_size);
        int cone_high = (cone_low + 1) % k;
        double larger_cone_angle = std::max(angle-cone_low, cone_high-angle);

        double length = edge.weight;
        double ret = length * std::sqrt(2*(1-std::cos(larger_cone_angle)));
        if (ret > max_ret) {
            max_ret = ret;
        }
    }
    return max_ret;
}

// trim helpers (optional)
static inline std::string trim(const std::string &s) {
    auto start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    auto end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

std::vector<std::vector<std::string>> load_csv(const std::string &filepath,
                                               char delimiter = ',',
                                               bool skip_empty_lines = true,
                                               bool trim_fields = true)
{
    std::ifstream in(filepath);
    if (!in.is_open()) {
        throw std::runtime_error("Could not open file: " + filepath);
    }

    std::vector<std::vector<std::string>> table;
    std::string line;
    while (std::getline(in, line)) {
        // remove possible trailing CR (for Windows-created files)
        if (!line.empty() && line.back() == '\r') line.pop_back();

        std::vector<std::string> row;
        std::string field;
        bool in_quotes = false;

        for (size_t i = 0; i < line.size(); ++i) {
            char c = line[i];
            if (c == '"') {
                if (in_quotes && i + 1 < line.size() && line[i + 1] == '"') {
                    // Escaped quote inside quoted field
                    field.push_back('"');
                    ++i; // skip next quote
                } else {
                    in_quotes = !in_quotes; // toggle quote state
                }
            } else if (c == delimiter && !in_quotes) {
                row.push_back(trim_fields ? trim(field) : field);
                field.clear();
            } else {
                field.push_back(c);
            }
        }
        // push last field
        row.push_back(trim_fields ? trim(field) : field);

        if (skip_empty_lines) {
            bool all_empty = true;
            for (auto &f : row) {
                if (!f.empty()) { all_empty = false; break; }
            }
            if (all_empty) continue;
        }
        table.push_back(std::move(row));
    }

    return table;
}

int main(int argc, char *argv[]) {
    // For Testing
    std::cout << "Running Lowerbound Graph" << std::endl;

    // TODO load lower bound graph
    auto lb = Lowerbound_Graph();
    lb.load_graph_from_fmi_file("../../data/thetamix.fmi");

    // load distances
    auto table = load_csv("../../data/results_lowerbound.csv");
    for (int i = 1; i < table.size(); i++) {
        int s = std::stoi(table[i][0]);
        int t = std::stoi(table[i][1]);
        double dist_daniel = std::stod(table[i][2]);
        double dist_theta = std::stod(table[i][3]);
        //calculate distances
        if (dist_theta < std::numeric_limits<double>::infinity()) {
            double dist = lb.get_lower_bound_path(s,	t, 12);
            std::cout << "theta Lowerbound distance from "<< s <<" to " << t <<": " << dist << std::endl;
        }
    }


}


