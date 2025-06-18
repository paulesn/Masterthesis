#include <set>
#include <unordered_set>
#include "../src/Quadtree.cpp" // Adjust to your filename/path
#include "../src/DataLoader.cpp" // Adjust to your filename/path

void REQUIRE(bool condition) {
    if (!condition) {
        throw std::runtime_error("Test failed (This is a custom REQUIRE implementation)");
    }
    cout << "Test passed" << std::endl;
}
void REQUIRE_FALSE(bool condition) {
    if (condition) {
        throw std::runtime_error("Test failed");
    }
    cout << "Test passed" << std::endl;
}

void TEST_CASE_1() {
    std::cout << ("Point equality and hashing") << std::endl;
    Point a(1.0, 2.0);
    Point b(1.0, 2.0);
    Point c(2.0, 3.0);

    REQUIRE(a == b);
    REQUIRE_FALSE(a == c);

    std::unordered_set<Point, PointHash> set;
    set.insert(a);
    set.insert(b); // Should not insert again
    set.insert(c);

    REQUIRE(set.size() == 2);
}

void TEST_CASE_2() {
    std::cout << ("Quadtree insertion and depth") << std::endl;
    Quadtree qt(2.0);

    qt.insert(Point(0.5, 0.5));
    qt.insert(Point(-0.5, -0.5));


    // Check duplicate insertion
    qt.insert(Point(0.5, 0.5));
    REQUIRE(qt.root.points.size() == 2); // Should not double-insert
}

void TEST_CASE_3() {
    std::cout << ("Quadtree handles out-of-bound points") << std::endl;
    Quadtree qt(1.0);

    qt.insert(Point(5.0, 5.0)); // Out of initial bounds

    // Ensure point was still inserted
    REQUIRE(qt.pointSet.find(Point(5.0, 5.0)) != qt.pointSet.end());
}

void TEST_CASE_4() {
    std::cout << ("Quadtree correctly splits") << std::endl;
    Quadtree qt(2.0);
    qt.insert(Point(0.1, 0.1));
    qt.insert(Point(0.2, 0.2)); // Should trigger split

    REQUIRE_FALSE(qt.root.is_leaf); // Node should have split
}

void TEST_CASE_5() {
    std::cout << ("Area containment") << std::endl;
    Quadtree qt(1.0);
    Point inside(0.5, 0.5);
    Point outside(2.0, 2.0);

    REQUIRE(qt.root.area_contains(inside));
    REQUIRE_FALSE(qt.root.area_contains(outside));
}

void TEST_CASE_6() {
    std::cout << ("WSPD produces well-separated pairs") << std::endl;
    Quadtree qt(2.0);

    // Insert non-overlapping points
    qt.insert(Point(-1.0, -1.0));
    qt.insert(Point(1.0, 1.0));

    auto pairs = qt.wspd(1.0); // Separation factor

    REQUIRE_FALSE(pairs.empty());
    for (const auto& [a, b] : pairs) {
        REQUIRE(!(a == b));
    }
}

void TEST_CASE_7() {
    std::cout << ("WSPD produces correct results with test.fmi") << std::endl;

    auto tup = load_fmi("data/0100.32.fmi", -1, false);
    auto systems = get<0>(tup);
    Graph graph = get<1>(tup);

    auto tree = new Quadtree(32);
    std::vector<Point> nodes;
    for (auto& system : systems) {

        tree->insert(Point(system.x, system.y));
        nodes.emplace_back(system.x, system.y);
    }

    auto result = tree->wspd(2);
    cout << "Found " << result.size() << " well-separated pairs." << std::endl;

    REQUIRE(result.size() == 537705);
}

void TEST_CASE_8() {
    std::cout << ("WSPD runs on a random set of nodes.fmi") << std::endl;

    auto tree = new Quadtree(32);
    std::vector<Point> nodes;
    for (int i=0; i < 1000; i++) {
        tree->insert(Point(rand() % 1000 + 1, rand() % 1000 + 1));
    }

    auto result = tree->wspd(2);

}

int main() {
    TEST_CASE_1();
    TEST_CASE_2();
    TEST_CASE_3();
    TEST_CASE_4();
    TEST_CASE_5();
    TEST_CASE_6();
    TEST_CASE_7();
    TEST_CASE_8();

    cout << "####################################################################################" << std::endl;
    cout << "All tests passed!" << std::endl;
    cout << "####################################################################################" << std::endl;
}