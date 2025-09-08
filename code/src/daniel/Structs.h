//
// Created by kochda on 31.01.24.
//

#ifndef EXACT_SHIPROUTING_STRUCTS_H
#define EXACT_SHIPROUTING_STRUCTS_H

#include <limits>
#include <map>

#include <string>
#include <iostream>
enum ImprovementAlgorithms{
    None, Greedy, Full, Funnel
};



const double maxWeight = std::numeric_limits<double>::infinity();
const unsigned int maxUInt = std::numeric_limits<unsigned int>::infinity();




#endif //EXACT_SHIPROUTING_STRUCTS_H
