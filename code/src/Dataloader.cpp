//
// Created by sebastian on 05.05.25.
//

#include "Dataloader.hpp"
#include <csignal>
#include <iostream>
#include <fstream>

#include "Quadtree.hpp"
#include "../lib/json.hpp"

using json = nlohmann::json;
using namespace std;

Star::Star(int id, double x, double y, double z, double jump_multiplicator)
    : id(id), x(x), y(y), z(z), jump_multiplicator(jump_multiplicator) {}


/**
 * Load a JSON file containing system data.
 * Each line in the file should be a valid JSON object.
 * The function will parse each line and extract the system ID, coordinates, and jump multiplicator.
 *
 * @param filepath Path to the JSON file.
 * @param limit Maximum number of lines to read from the file. Default is -1 (no limit).
 * @return A vector of tuples containing system data (id, x, y, z, jump multiplicator).
 */
std::vector<Star> load_json(const std::string& filepath, int limit) {
    std::cout << "Loading json from " << filepath << std::endl;
    std::vector<Star> systems; // id, x, y, z, jump multiplicator

    std::ifstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filepath << std::endl;
        raise(SIGINT); // Raise SIGINT to terminate the program
    }

    int counter = 0;

    try {
        for (std::string line; std::getline(file, line);) {
            // update counter
            if (counter % 1000 == 0) {
                std::cout << "Loaded " << counter << " lines" << std::endl;
            }
            counter++;

            // Skip empty lines and comments
            if (line.empty() || line[0] == '#' || line[0] == '[' || line[0] == ']') {
                std::cout << line << std::endl;
                continue; // Skip empty lines
            }
            if (limit > 0 && counter >= limit) {
                break;
            }

            // remove trailing commas
            if (line.back() == ',') {
                line.pop_back();
            }

            auto data = json::parse(line);
            int id = data["id64"];
            int x = data["coords"]["x"];
            int y = data["coords"]["y"];
            int z = data["coords"]["z"];
            int jump_multiplicator = data["mainStar"] == "Neutron Star" ? 4 : 1;

            systems.emplace_back(id, x, y, z, jump_multiplicator);
        }
    } catch (const json::parse_error& e) {
        std::cerr << "JSON parse error: " << e.what() << std::endl;
        raise(SIGINT); // Raise SIGINT to terminate the program
    }

    return systems;
}

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
 * This function loads a CH graph in the fmi format.
 * @param filepath the path to the fmi file
 * @param limit the maximum number of lines to read from the file, -1 for no limit only use if with_edges is false
 * @param with_edges if true, the function will also load the edges of the graph, otherwise it will only load the nodes
 * @return
 */
std::tuple<std::vector<Point>,Graph> load_fmi_ch(const std::string &filepath, int limit, bool with_edges) {
    std::cout << "Loading fmi format from " << filepath << std::endl;
    std::vector<Point> systems; // id, x, y, z, jump multiplicator

    // print current working directory
    std::cout << "Current working directory: " << std::filesystem::current_path() << std::endl;
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
            cout << "The graph has " << number_of_nodes << " nodes\n>";
            technical_counter++;
            continue;
        }
        if (technical_counter == -1) {
            // Second line contains the number of edges
            number_of_edges = std::stoi(line);
            cout << "The graph has " << number_of_edges << " edges\n>";
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
            // this creates the point. the data is in the format: id nonesense x y nonsense CH-level...
            systems.emplace_back(stod(arr[2]), stod(arr[3]), stoi(arr[0]), stoi(arr[5]));
            //systems.emplace_back(stod(arr[2]), stod(arr[3]), stoi(arr[0]), -1);
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
            int u_id = stoi(arr[0]);
            int v_id = stoi(arr[1]);
            Point u = g.id_point_map[u_id];
            Point v = g.id_point_map[v_id];
            double dist = stod(arr[3]); //euklidian_distance(u,v);
            g.addEdge(u_id, v_id, dist);
        }
    }

    cout << endl;
    return make_tuple(systems, g);
}

/**
 * This function loads a graph in the fmi format.
 * @param filepath the path to the fmi file
 * @param limit the maximum number of lines to read from the file, -1 for no limit only use if with_edges is false
 * @param with_edges if true, the function will also load the edges of the graph, otherwise it will only load the nodes
 * @return
 */
std::tuple<std::vector<Point>,Graph> load_fmi(const std::string &filepath, int limit, bool with_edges) {
    std::cout << "Loading fmi format from " << filepath << std::endl;
    std::vector<Point> systems; // id, x, y, z, jump multiplicator

    // print current working directory
    std::cout << "Current working directory: " << std::filesystem::current_path() << std::endl;
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
            systems.emplace_back(stoi(arr[0]), stod(arr[2]), stod(arr[3]), -1);
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
            g.addEdge(stoi(arr[0]), stoi(arr[1]), stod(arr[2]));
        }
    }


    return make_tuple(systems, g);
}

