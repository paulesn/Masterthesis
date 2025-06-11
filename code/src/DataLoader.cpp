//
// Created by sebastian on 05.05.25.
//

#include <csignal>
#include <iostream>
#include <fstream>
#include <tuple>
#include <tuple>
#include <tuple>
#include <tuple>
#include <vector>
#include "../lib/json.hpp"
#include "Graph.cpp"

using json = nlohmann::json;

/**
 * Load a JSON file containing system data.
 * Each line in the file should be a valid JSON object.
 * The function will parse each line and extract the system ID, coordinates, and jump multiplicator.
 *
 * @param filepath Path to the JSON file.
 * @param limit Maximum number of lines to read from the file. Default is -1 (no limit).
 * @return A vector of tuples containing system data (id, x, y, z, jump multiplicator).
 */

struct Star {
    int id;
    double  x, y, z, jump_multiplicator;
    Star(int id, double x, double y, double z, double jump_multiplicator)
        : id(id), x(x), y(y), z(z), jump_multiplicator(jump_multiplicator) {}
};


std::vector<Star> load_json(const std::string& filepath, int limit = -1) {
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


std::tuple<std::vector<Star>,Graph> load_fmi(const std::string &filepath, int limit = -1) {
    std::cout << "Loading fmi format from " << filepath << std::endl;
    std::vector<Star> systems; // id, x, y, z, jump multiplicator

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
            g = Graph(number_of_nodes);
            continue;
        }
        if (technical_counter == -1) {
            // Second line contains the number of edges
            number_of_edges = std::stoi(line);
            technical_counter++;
            cout << ">";
            continue;
        }
        if (technical_counter < number_of_nodes) {
            if (counter % (number_of_nodes/100) == 0) {
                std::cout << "|";
                cout.flush();
            }

            // Read node data
            technical_counter++;
            std::string delimiter = " ";

            auto arr = split(line, delimiter);
            systems.emplace_back(stoi(arr[0]), stod(arr[2]), stod(arr[3]), 0, 1);
        } else {
            if (counter == number_of_nodes) {
                // We have read all nodes
                std::cout << endl << ">";
            }
            if (technical_counter >= number_of_nodes + number_of_edges) {
                // We have read all nodes and edges
                break;
            }
            if ((counter-number_of_nodes) % (number_of_edges/100) == 0) {
                std::cout << "|";
                cout.flush();
            }
            // Read edge data
            std::string delimiter = " ";
            auto arr = split(line, delimiter);
            g.addEdge(stoi(arr[0]), stoi(arr[1]), stod(arr[2]));
        }
    }


    return make_tuple(systems, g);
}

