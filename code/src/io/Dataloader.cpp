//
// Created by sebastian on 05.05.25.
//

#include "Dataloader.hpp"
#include <csignal>
#include <iostream>
#include <fstream>

#include "../structure/Quadtree.hpp"

using namespace std;

Star::Star(int id, double x, double y, double z, double jump_multiplicator)
    : id(id), x(x), y(y), z(z), jump_multiplicator(jump_multiplicator) {}



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

/**
 * This function loads a graph in the fmi format.
 * @param filepath the path to the fmi file
 * @param limit the maximum number of lines to read from the file, -1 for no limit only use if with_edges is false
 * @param with_edges if true, the function will also load the edges of the graph, otherwise it will only load the nodes
 * @return
 */
std::tuple<std::vector<Pointc>,Graph> load_fmi(const std::string &filepath, int limit, bool with_edges) {
    std::cout << "Loading fmi format from " << filepath << std::endl;
    std::vector<Pointc> systems; // id, x, y, z, jump multiplicator

    // print current working directory
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filepath << std::endl;
        raise(SIGINT); // Raise SIGINT to terminate the program
    }

    int counter = 0;

    int number_of_nodes = 0;
    int number_of_edges = 0;
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
        if (limit > 0 && counter >= limit) {
            break;
        }
        if (technical_counter == -2) {
            // First line contains the number of nodes and edges
            number_of_nodes = std::stoi(line);
            technical_counter++;
            continue;
        }
        if (technical_counter == -1) {
            // Second line contains the number of edges
            number_of_edges = std::stoi(line);
            technical_counter++;
            cout << "Loading Nodes \n>";
            continue;
        }
        if (technical_counter < number_of_nodes) {
            if (number_of_nodes > 100 && technical_counter != 0 && counter % (number_of_nodes/100) == 0) {
                std::cout << "|";
                cout.flush();
            }

            // Read node data
            technical_counter++;
            std::string delimiter = " ";

            auto arr = split(line, delimiter);
            // this creates the point. the data is in the format: id nonesense x y ...
            systems.emplace_back(stoi(arr[2]), stod(arr[3]), stod(arr[0]));
        } else {
            if (!with_edges) {
                // If we are not loading edges, we can stop here
                std::cout << endl << "Skipping Edges\n>";
                break;
            }
            if (technical_counter == number_of_nodes) {
                // We have read all nodes
                std::cout << endl << "Loading Edges\n>";
            }
            technical_counter++;
            if (technical_counter >= number_of_nodes + number_of_edges) {
                // We have read all nodes and edges
                break;
            }
            if (number_of_nodes > 100 && number_of_edges > 100 && (counter-number_of_nodes) % (number_of_edges/100) == 0) {
                std::cout << "|";
                cout.flush();
            }
            // Read edge data
            std::string delimiter = " ";
            auto arr = split(line, delimiter);
            if (!graph_created) {
                g = Graph(systems);
                graph_created = true;
            }
            g.addEdge(stoi(arr[0]), stoi(arr[1]), stod(arr[2]), true);
        }
    }


    return make_tuple(systems, g);
}

Graph load_coastline(const std::string &filepath) {
    std::cout << "Loading coastline format from " << filepath << std::endl;
    std::vector<Pointc> points;
    std::vector<tuple<int, int>> edges;

    // print current working directory
    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filepath << std::endl;
        raise(SIGINT); // Raise SIGINT to terminate the program
    }

    int counter = 0;
    int number_of_nodes = 0;
    int number_of_edges = 0;

    bool coastline_active = false;
    int last = -1;
    int first = -1;

    for (std::string line; std::getline(file, line);){
        // defensive programming
        if (line.empty() || line[0] == '#') {
            continue; // Skip empty lines
        }

        // if the line is just a number start a new coastline
        auto line_v = split(line, " ");
        if (line_v.size() == 1) {
            coastline_active = false;
            if (last != -1 && first != -1) {
                // add edge from last point to first point
                edges.emplace_back(last, first);
                number_of_edges++;
            }
            continue;
        }
        if (line_v.size() == 2) {
            // check if the point exists in the map
            long double x = stod(line_v[0]);
            long double y = stod(line_v[1]);
            points.emplace_back(x, y, counter);
            number_of_nodes++;

            if (coastline_active) {
                // add edge from last point
                edges.emplace_back(last, counter);
                number_of_edges++;
                last = counter;
            }else{
                last = counter;
                first = counter;
                coastline_active = true;
            }
            counter++;
        }
    }
    Graph g = Graph(points);
    for (const auto& edge : edges) {
        g.addEdge(get<0>(edge), get<1>(edge), 1.0); // Assuming weight of 1.0 for coastline edges
    }
    return g;
}

Graph load_gamemap(const std::string &filepath) {

}

