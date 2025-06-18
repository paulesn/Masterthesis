//
// Created by sebastian on 18.06.25.
//

#ifdef PYTHON
# include <pybind11/pybind11.h>
# include <pybind11/stl.h>
# include <pybind11/complex.h>
#endif

# include <vector>
# include <tuple>
# include <omp.h>
# include <queue>
# include <limits>
# include <algorithm>
# include <iostream>
#include <cstdio>
#include <fstream>

#include "HubLabels.hpp"

using namespace std;


//////////////////////////////////////////////////////////////////////////////////////
/// QUERIES
//////////////////////////////////////////////////////////////////////////////////////
std::tuple<int, int> optimized_query(
    std::vector<std::tuple<int,int>>& forward_labels,
    std::vector<std::tuple<int,int>>& backward_labels
)
{
    // both hub label lists are sorted by hub labels
    int forward_index = 0;
    int backwards_index = 0;
    int forward;
    int backward;
    int distance = std::numeric_limits<int>::max();
    int hub = -1;


    while(forward_index < forward_labels.size() & backwards_index < backward_labels.size()) {
        forward = get<0>(forward_labels[forward_index]);
        backward = get<0>(backward_labels[backwards_index]);
        if (forward == backward){
            int temp_dist = get<1>(forward_labels[forward_index]) + get<1>(backward_labels[backwards_index]);
            if(temp_dist < distance) {
                hub = forward;
                distance = temp_dist;
            }
        }
        if (forward > backward) {
            backwards_index++;
        } else {
            forward_index++;
        }
    }
    return {hub, distance};
}


std::tuple<int, int> query(
std::tuple<
    std::vector<std::vector<std::tuple<int,int>>>, // wenn es nur die relevanten label bekommt ist besser
    std::vector<std::vector<std::tuple<int,int>>> // und bootstrab muss nicht mehr so viel kopieren
    > hub_labels,
int start,
int target
)
{
    std::vector<std::vector<std::tuple<int,int>>>& forward_labels = std::get<0>(hub_labels);
    std::vector<std::vector<std::tuple<int,int>>>& backward_labels = std::get<1>(hub_labels);

    return optimized_query(forward_labels[target], backward_labels[start]);
}

// Function to bootstrap labels based on forward and backward hub labels
std::vector<std::tuple<int,int>> bootstrap(
    int node,
    std::vector<std::vector<std::tuple<int,int>>>& forward_hub_labels,
    std::vector<std::vector<std::tuple<int,int>>>& backward_hub_labels
) {
    // Optimize forward labels
    vector<tuple<int, int>> final_forwards_labels;
    for (auto& label : forward_hub_labels[node]) {
        if (std::get<0>(label) == node) {
            final_forwards_labels.push_back(label);
            continue;
        }
        auto query_result = optimized_query(forward_hub_labels[node], backward_hub_labels[std::get<0>(label)]);
        if (std::get<1>(query_result) == std::get<1>(label)) {
            final_forwards_labels.push_back(label);
        }

    }

    return {final_forwards_labels};
}

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
) {
    std::vector<tuple<int, int>> costs;  // Initialize label cost as tuple (label, cost)
    // add first cost as (label,cost) = (start,0)
    costs.push_back({start, 0}); // Start with the assumption of zero cost for the initial position

    // Process all edges leaving the 'start' node
    for (const auto& edge : edges[start]) {
        int target = std::get<0>(edge);
        int edge_cost = std::get<1>(edge);

        // Check all labels at the target node
        for (const auto& label_info : existing_labels[target]) {
            int label = std::get<0>(label_info);
            int label_cost = std::get<1>(label_info);
            int total_cost = edge_cost + label_cost;

            // Update the label cost if it's cheaper
            costs.push_back({label, total_cost});
        }
    }

    // sort costs and only take the one with the lowes value
    sort(costs.begin(), costs.end());
    vector<tuple<int, int>> results;
    int current_hub = get<0>(costs[0]);
    int current_lowest_value = get<1>(costs[0]);

    for (int index = 1; index < costs.size(); index++) {
        if (current_hub == get<0>(costs[index])) {
            // if this entry has the same hub as the last one
            current_lowest_value = min(current_lowest_value, get<1>(costs[index]));
            continue;
        }
        // this entry represents a new hub
        results.push_back({current_hub, current_lowest_value});
        current_hub = get<0>(costs[index]);
        current_lowest_value = get<1>(costs[index]);
    }
    // we have to add the element if it is the last one in the costs vector
    results.push_back({current_hub, current_lowest_value});

    return results;
}

// Optimized function to generate hub labels
std::tuple<std::vector<std::vector<std::tuple<int,int>>>,std::vector<std::vector<std::tuple<int,int>>>>
generate_hub_labels(
    const std::vector<std::vector<int>>& order,                             // the oder in which the nodes are processed as CH a list of nodes for each level
    const std::vector<std::vector<std::tuple<int,int>>>& upwards_edges,     // all edges that go upwards in the order
    const std::vector<std::vector<std::tuple<int,int>>>& downwards_edges,   // all edges that go downwards in the order
    int max_threads                                                  // limit the number of threads
    ) {
    std::vector<std::vector<std::tuple<int,int>>> forward_hub_labels;
    std::vector<std::vector<std::tuple<int,int>>> backward_hub_labels;

    //cout << "Memory at start of function" << endl;
    //printMemoryUsage();

    // calculate amount of nodes over all levels
    int total_length = 0;
    for (const auto& level : order) {
        total_length += level.size();
    }

    // identify number of threads for parallelization
    int num_threads = max(1,omp_get_num_procs()-1);  // Get the number of available processors
    if (max_threads > 0) {
        num_threads = min(num_threads, max_threads);
    }

    // Reserve memory for maps and sets to minimize reallocations
    forward_hub_labels.reserve(total_length);
    backward_hub_labels.reserve(total_length);
    for (size_t i = 0; i < total_length; ++i) {
        forward_hub_labels.emplace_back();
        backward_hub_labels.emplace_back();
    }

    //cout << "Memory after reservig space" <<  endl;
    //printMemoryUsage();

    for (const vector<int>& level : order) {
        #pragma omp parallel for num_threads(num_threads) shared(forward_hub_labels, backward_hub_labels, upwards_edges, downwards_edges)
        for (const int& node : level) {

            forward_hub_labels[node] = reachable(node, upwards_edges, forward_hub_labels);
            backward_hub_labels[node] = reachable(node, downwards_edges, backward_hub_labels);

            forward_hub_labels[node] = bootstrap(node, forward_hub_labels, backward_hub_labels);
            backward_hub_labels[node] = bootstrap(node, backward_hub_labels, forward_hub_labels);

        }
    }
    return {forward_hub_labels, backward_hub_labels};
}
