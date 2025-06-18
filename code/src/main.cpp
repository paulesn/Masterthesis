
#include "Dataloader.hpp"
#include "Quadtree.hpp"
#include <ostream>
#include <iostream>
#include <set>
#include <chrono>
#include <csignal>


using namespace std;


int main() {
    auto tup = load_fmi("../../data/0200.32.fmi", -1);
    auto systems = get<0>(tup);
    Graph graph = get<1>(tup);
    std::cout << "Loaded " << systems.size() << " systems." << std::endl;


    auto tree = new Quadtree(32);
    int counter = 0;
    std::vector<Point> nodes;

    std::cout << "Inserting systems into the quadtree..." << std::endl;
    Point id_1 = Point(0,0,-1);
    for (auto& system : systems) {
        counter++;
        if (counter % 10000 == 0) {
            std::cout << "Inserted " << counter << " systems into the quadtree." << std::endl;
        }
        auto p = Point(system.x, system.y, counter);
        tree->insert(p);
        if (counter == 1) {
            id_1 = p; // Store the first point for later checks
        }
        if (id_1.id !=-1 && !tree->contains(id_1)) {
            raise(SIGINT); // Raise SIGINT to terminate the program
        }
        nodes.emplace_back(p);

    }
    for (Point node: nodes) {
        if (!tree->contains(node)) {
            std::cout << "Error: Quadtree does not contain point" << node.id <<": (" << node.x << ", " << node.y << ")." << std::endl;
            //raise(SIGINT); // Raise SIGINT to terminate the program
        }
    }
}
