
#include "Quadtree.cpp"
#include <cmath>
#import <vector>
#import <tuple>
#include <ostream>
#include <iostream>

using namespace std;

struct Rect {
    int x, y, z; // center coordinates
    int l; // length of the side
    int level; // level of the area
    int id;
    explicit Rect(int x, int y, int z, int l, int level)
        : x(x), y(y), z(z), l(l), level(level) {
        static int global_id = 0; // static variable to keep track of the number of rectangles
        id = global_id++;
    }
};


void update_log() {
    static int counter = 0;            // persists across calls, initialized once
    static int counter2 = 0;           // persists across calls, initialized once
    if (++counter >= 10000) {          // increment and check
        std::cout << ++counter2 << '\n';     // print when the count reaches 1000
        counter = 0;                  // reset for the next batch of 1000
    }
}


/**
 * This is a recursive implementation of wspd that stops going deeper into a pair if it is already well seperated
 * @param a a rectangle
 * @param b a rectangle
 * @return the amount of pairs of rectangles that are in the same area
 */
tuple<int, int> rec_wspd_3d(Rect area, Rect target, double s, int max_level, double jump_range) {
    update_log();
    if (max_level <= 0) {
        return tuple(0,0);
    }
    vector<Rect> areas =  {};
    // add sub areas to the vector
    if (area.level > target.level) {
        auto temp = area;
        area = target;
        target = area;
    }
    // target now has the larger level ==> therefore is deeper
    int level = area.level;
    areas.emplace_back(area);
    areas.emplace_back(area.x+(area.l/2), area.y+(area.l/2), area.z+(area.l/2), area.l/2, level+1);
    areas.emplace_back(area.x+(area.l/2), area.y+(area.l/2), area.z-(area.l/2), area.l/2, level+1);
    areas.emplace_back(area.x+(area.l/2), area.y-(area.l/2), area.z+(area.l/2), area.l/2, level+1);
    areas.emplace_back(area.x+(area.l/2), area.y-(area.l/2), area.z-(area.l/2), area.l/2, level+1);
    areas.emplace_back(area.x-(area.l/2), area.y+(area.l/2), area.z+(area.l/2), area.l/2, level+1);
    areas.emplace_back(area.x-(area.l/2), area.y+(area.l/2), area.z-(area.l/2), area.l/2, level+1);
    areas.emplace_back(area.x-(area.l/2), area.y-(area.l/2), area.z+(area.l/2), area.l/2, level+1);
    areas.emplace_back(area.x-(area.l/2), area.y-(area.l/2), area.z-(area.l/2), area.l/2, level+1);

    int sum = 0;
    int jump = 0;
    for (int i = 0; i < areas.size(); ++i) {
        auto a = areas[i];
        for (int o = 0; o < areas.size(); ++o) {
            auto b = areas[o];
            if (a.id < b.id) {
                continue; // skip the areas so that no pair appears double
            }
            if (a.id == area.id && b.id == target.id) {
                continue;
                // skip the areas so that no pair appears double
                // this was already checked in the previous recursion
            }
            if (b.id == area.id && a.id == target.id) {
                continue;
                // skip the areas so that no pair appears double
                // this is a duplicate of the check before
            }
            double dist_squared = pow(a.x - b.x, 2) + pow(a.y - b.y, 2) + pow(a.z - b.z, 2);
            double radius_squared = max((2*a.l*a.l), (2*b.l*b.l));
            if (dist_squared >= s*radius_squared) {
                sum++;
                if (dist_squared > jump_range) {
                    jump++;
                }
                continue; // skip the areas that are well separated
            }
            sum += get<0>(rec_wspd_3d(a, b, s, max_level-1, jump_range));
            jump += get<1>(rec_wspd_3d(a, b, s, max_level-1, jump_range));
        }
    }
    return tuple(sum, jump);
}

/**
 * This is a recursive implementation of wspd that stops going deeper into a pair if it is already well seperated
 * @param a a rectangle
 * @param b a rectangle
 * @return the amount of pairs of rectangles that are in the same area
 */
tuple<int, int> rec_wspd_2d(Rect area, Rect target, double s, int max_level, double jump_range) {
    update_log();
    if (max_level <= 0) {
        return tuple(0,0);
    }
    vector<Rect> areas =  {};
    // add sub areas to the vector
    if (area.level > target.level) {
        auto temp = area;
        area = target;
        target = area;
    }
    // target now has the larger level ==> therefore is deeper
    int level = area.level;
    areas.emplace_back(area);
    areas.emplace_back(area.x+(area.l/2), area.y+(area.l/2), 0, area.l/2, level+1);
    areas.emplace_back(area.x+(area.l/2), area.y-(area.l/2), 0, area.l/2, level+1);
    areas.emplace_back(area.x-(area.l/2), area.y+(area.l/2), 0, area.l/2, level+1);
    areas.emplace_back(area.x-(area.l/2), area.y-(area.l/2), 0, area.l/2, level+1);

    int sum = 0;
    int jump = 0;
    for (int i = 0; i < areas.size(); ++i) {
        auto a = areas[i];
        for (int o = 0; o < areas.size(); ++o) {
            auto b = areas[o];
            if (a.id < b.id) {
                continue; // skip the areas so that no pair appears double
            }
            if (a.id == area.id && b.id == target.id) {
                continue;
                // skip the areas so that no pair appears double
                // this was already checked in the previous recursion
            }
            if (b.id == area.id && a.id == target.id) {
                continue;
                // skip the areas so that no pair appears double
                // this is a duplicate of the check before
            }
            double dist_squared = pow(a.x - b.x, 2) + pow(a.y - b.y, 2) + pow(a.z - b.z, 2);
            double radius_squared = max((2*a.l*a.l), (2*b.l*b.l));
            if (dist_squared >= s*radius_squared) {
                sum++;
                if (dist_squared > jump_range) {
                    jump++;
                }
                continue; // skip the areas that are well separated
            }
            sum += get<0>(rec_wspd_2d(a, b, s, max_level-1, jump_range));
            jump += get<1>(rec_wspd_2d(a, b, s, max_level-1, jump_range));
        }
    }
    return tuple(sum, jump);
}

/**
 * This is a recursive implementation of wspd that stops going deeper into a pair if it is already well seperated
 * @param a a rectangle
 * @param b a rectangle
 * @return the amount of pairs of rectangles that are in the same area
 */
tuple<int, int> rec_wspd_graph(Rect area, Rect target, double s, int max_level, Quadtree* tree) {
    update_log();
    if (max_level <= 0) {
        return tuple(0,0);
    }
    vector<Rect> areas =  {};
    // add sub areas to the vector
    if (area.level > target.level) {
        auto temp = area;
        area = target;
        target = area;
    }
    // skip area if it is empty in octtree:
    if (tree->get_all_in_rect(area.x, area.y, area.l/2).empty()) {
        return tuple(0, 0);
    }
    if (tree->get_all_in_rect(target.x, target.y, target.l/2).empty()) {
        return tuple(0, 0);
    }
    // target now has the larger level ==> therefore is deeper
    int level = area.level;
    areas.emplace_back(area);
    areas.emplace_back(area.x+(area.l/2), area.y+(area.l/2), 0, area.l/2, level+1);
    areas.emplace_back(area.x+(area.l/2), area.y-(area.l/2), 0, area.l/2, level+1);
    areas.emplace_back(area.x-(area.l/2), area.y+(area.l/2), 0, area.l/2, level+1);
    areas.emplace_back(area.x-(area.l/2), area.y-(area.l/2), 0, area.l/2, level+1);

    int sum = 0;
    int jump = 0;
    for (int i = 0; i < areas.size(); ++i) {
        auto a = areas[i];
        for (int o = 0; o < areas.size(); ++o) {
            auto b = areas[o];
            if (a.id < b.id) {
                continue; // skip the areas so that no pair appears double
            }
            if (a.id == area.id && b.id == target.id) {
                continue;
                // skip the areas so that no pair appears double
                // this was already checked in the previous recursion
            }
            if (b.id == area.id && a.id == target.id) {
                continue;
                // skip the areas so that no pair appears double
                // this is a duplicate of the check before
            }
            double dist_squared = pow(a.x - b.x, 2) + pow(a.y - b.y, 2) + pow(a.z - b.z, 2);
            double radius_squared = max((2*a.l*a.l), (2*b.l*b.l));
            if (dist_squared >= s*radius_squared) {
                sum++;
                continue; // skip the areas that are well separated
            }
            auto temp = rec_wspd_graph(a, b, s, max_level-1, tree);
            sum += get<0>(temp);
            jump += get<1>(temp);
        }
    }
    return tuple(sum, jump);
}

tuple<int, int> run_wspd(
    double s,
    int max_depth,
    double max_j,
    int size_of_largest_rect
) {
    vector<Rect> areas =  {};  // vector of areas
    vector<tuple<Rect, Rect>> pairs = {};  // vector of pairs

    Rect all = Rect(0, 0, 0, size_of_largest_rect, 0); // all areas
    return rec_wspd_2d(all, all, s, max_depth, max_j); // start the recursion
}

tuple<int, int> run_wspd_graph(
    double s,
    vector<Star> systems,
    int max_depth
) {
    vector<Rect> areas =  {};  // vector of areas
    vector<tuple<Rect, Rect>> pairs = {};  // vector of pairs
    auto tree = new Quadtree(10, 1000000); // create a new octtree with a maximum of 10 objects per node and a maximum of 1 million nodes
    for (auto& system : systems) {
        tree->insert(&system); // insert the systems into the octtree
    }

    Rect all = Rect(0, 0, 0, 100000, 0); // all areas
    return rec_wspd_graph(all, all, s, max_depth, tree); // start the recursion
}