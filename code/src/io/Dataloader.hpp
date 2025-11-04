#ifndef DATALOADER_HPP
#define DATALOADER_HPP

#include <string>
#include <vector>
#include <tuple>
#include "../structure/Graph.hpp"

// Star data structure
struct Star {
    int id;
    double x, y, z, jump_multiplicator;

    Star(int id, double x, double y, double z, double jump_multiplicator);
};

// Load star systems from a JSON file
std::vector<Star> load_json(const std::string& filepath, int limit = -1);

// Split a string by a given delimiter
std::vector<std::string> split(const std::string& s, const std::string& delimiter);

// Load a star system and graph from a file in FMI format
std::tuple<std::vector<Pointc>, Graph> load_fmi(const std::string& filepath, int limit = -1, bool with_edges = true);

Graph load_coastline(const std::string &filepath);
#endif // DATALOADER_HPP

Graph load_gamemap(const std::string &filepath);
