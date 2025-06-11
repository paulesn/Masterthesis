
//#include "WSPD.cpp"
#include "DataLoader.cpp"
#include "Quadtree.cpp"
#include <ostream>
#include <iostream>



int main() {
    auto tup = load_fmi("data/0200.32.fmi", -1);
    auto systems = get<0>(tup);
    Graph graph = get<1>(tup);
    std::cout << "Loaded " << systems.size() << " systems." << std::endl;
    auto tree = new Quadtree(32);
    int counter = 0;
    std::vector<Point> nodes;
    std::cout << "Inserting systems into the quadtree..." << std::endl;
    for (auto& system : systems) {
        counter++;
        if (counter % 10000 == 0) {
            std::cout << "Inserted " << counter << " systems into the quadtree." << std::endl;
        }
        tree->insert(Point(system.x, system.y));
        nodes.emplace_back(system.x, system.y);
    }

    std::cout << "Inserted all systems into the quadtree." << std::endl;
    auto start = std::chrono::high_resolution_clock::now();
    auto result = tree->wspd(2);
    auto end = std::chrono::high_resolution_clock::now();
    std::cout << "wspd(2) took "
              << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
              << " ms." << std::endl;
    std::cout << "Found " << result.size() << " well-separated pairs." << std::endl;
    std::cout << "For   " << systems.size() << " nodes." << std::endl;
    set <std::tuple<Point, Point>> uniqeue_pairs;
    for (const auto& pair : result) {
        if (get<0>(pair) == get<1>(pair)) continue; // skip self-pairs
        if (get<0>(pair) < get<1>(pair)) {
            uniqeue_pairs.insert(pair);
        } else {
            // ensure unique pairs by always inserting in a consistent order
            std::pair<Point, Point> reversed_pair(get<1>(pair), get<0>(pair));
            uniqeue_pairs.insert(reversed_pair);
        }
    }
    std::cout << "Unique pairs: " << uniqeue_pairs.size() << std::endl;
}
