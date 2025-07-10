//
// Created by sebastian on 10.07.25.
//


#include "../src/ThetaSpanner.hpp"
int main() {
    std::cout << "ThetaSpanner Test" << std::endl;
    std::vector<Point> points = {
        Point(0, 0, 0, 0),
        Point(1, 1, 1, 0),
        Point(0, 1, 2, 0),
        Point(1, 0, 3, 0),
        Point(-1, 0, 4, 0),
        Point(10000, 0, 3, 0),
        Point(9999, 9999, 3, 0),
        Point(0, 0, 3, 0),
    };
    for (const auto& p : points) {
        for (const auto& p2: points)
        std::cout << "Point: " << p << " - " << p2 << ": " << angle_between(p2, p) <<std::endl;
    }
}