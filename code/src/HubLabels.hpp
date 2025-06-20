//
// Created by sebastian on 18.06.25.
//

#ifndef HUBLABELS_H
#define HUBLABELS_H

# include <vector>
# include <tuple>

using namespace std;


//////////////////////////////////////////////////////////////////////////////////////
/// QUERIES
//////////////////////////////////////////////////////////////////////////////////////
std::tuple<int, int> optimized_query(
    std::vector<std::tuple<int,int>>& forward_labels,
    std::vector<std::tuple<int,int>>& backward_labels
);


std::tuple<int, int> query(
std::tuple<
    std::vector<std::vector<std::tuple<int,int>>>, // wenn es nur die relevanten label bekommt ist besser
    std::vector<std::vector<std::tuple<int,int>>> // und bootstrab muss nicht mehr so viel kopieren
    > hub_labels,
int start,
int target
);

// Function to bootstrap labels based on forward and backward hub labels
std::vector<std::tuple<int,int>> bootstrap(
    int node,
    std::vector<std::vector<std::tuple<int,int>>>& forward_hub_labels,
    std::vector<std::vector<std::tuple<int,int>>>& backward_hub_labels
);

//////////////////////////////////////////////////////////////////////////////////////
/// HUB LABEL GENERATOR
//////////////////////////////////////////////////////////////////////////////////////


/**
 * This function returns all nodes reachable from the start node given the upwards edges
 * @param start  a node (int)
 * @param upwards_edges a map of sources to a vector of targets
 * @return a set of node ids
 */
// Function to find reachable labels with their minimal costs
std::vector<std::tuple<int, int>> reachable(
    int start,
    const std::vector<std::vector<std::tuple<int, int>>>& edges,
    const std::vector<std::vector<std::tuple<int, int>>>& existing_labels
);

// Optimized function to generate hub labels
std::tuple<std::vector<std::vector<std::tuple<int,int>>>,std::vector<std::vector<std::tuple<int,int>>>>
generate_hub_labels(
    const std::vector<std::vector<int>>& order,                             // the oder in which the nodes are processed as CH a list of nodes for each level
    const std::vector<std::vector<std::tuple<int,int>>>& upwards_edges,     // all edges that go upwards in the order
    const std::vector<std::vector<std::tuple<int,int>>>& downwards_edges,   // all edges that go downwards in the order
    int max_threads =-1                                                  // limit the number of threads
    );

int main2();
#endif //HUBLABELS_H
