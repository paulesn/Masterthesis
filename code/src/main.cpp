
//#include "WSPD.cpp"
#include "DataLoader.cpp"
#include "Quadtree.cpp"
#include <ostream>
#include <iostream>



int main() {
    auto tup = load_fmi("data/0100.32.fmi", -1);
    auto systems = get<0>(tup);
    Graph graph = get<1>(tup);
    std::cout << "Loaded " << systems.size() << " systems." << std::endl;
    auto tree = new Quadtree(32);
    int counter = 0;
    std::vector<Point> nodes;
    std::cout << "Inserting systems into the quadtree..." << std::endl;
    for (auto& system : systems) {
        counter++;
        if (counter % 100 == 0) {
            std::cout << "Inserted " << counter << " systems into the quadtree." << std::endl;
        }
        tree->insert(Point(system.x, system.y));
        nodes.emplace_back(system.x, system.y);
    }

    std::cout << "Inserted all systems into the quadtree." << std::endl;
    auto result = tree->wspd(2);
    std::cout << "Found " << result.size() << " well-separated pairs." << std::endl;
    std::cout << "For   " << systems.size() << " nodes." << std::endl;
}
