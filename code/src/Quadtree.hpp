#ifndef QUADTREE_HPP
#define QUADTREE_HPP

#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <unordered_set>
#include "Graph.hpp" // You may want to replace this with Graph.hpp if possible

// Hash function for Point
struct PointHash {
    std::size_t operator()(const Point& p) const;
};

// Euclidean distance between two points
double euklidian_distance(Point a, Point b);

// Node of the Quadtree
struct QuadtreeNode {
    std::vector<Point> points;
    double cX, cY, height;
    QuadtreeNode *NW, *NO, *SW, *SO;
    bool is_leaf;
    QuadtreeNode *parent = nullptr; // Pointer to the parent node, if needed
    double inner_sp_distance = -1; // Used for WSPD with shortest path distance

    QuadtreeNode(double cX_, double cY_, double height_, QuadtreeNode* parent);
    QuadtreeNode(const QuadtreeNode& other);
    QuadtreeNode& operator=(const QuadtreeNode& other);
    ~QuadtreeNode();

    void insert(Point s);
    [[nodiscard]] bool area_contains(Point s) const;
    [[nodiscard]] bool intercept_rect(double x, double y, double h) const;
    [[nodiscard]] std::string rec_string(int level = 0) const;

    [[nodiscard]] std::vector<Point> get_all_points() const;

    [[nodiscard]] bool contains(const Point & p) const;

private:
    void insertInternal(Point s);
};

// Quadtree structure
struct Quadtree {
    QuadtreeNode root;
    double height;
    std::unordered_set<Point, PointHash> pointSet;

    explicit Quadtree(double height_ = 1);
    ~Quadtree() = default;

    void insert(Point p);
    void print() const;
    std::vector<Point> get_all_points() const;
    bool contains(Point p) const{return root.contains(p);}
};

// Normalize the ordering of two QuadtreeNodes
std::tuple<QuadtreeNode*, QuadtreeNode*> normalize_pair(QuadtreeNode* a, QuadtreeNode* b);

// Well-separated pair decomposition (WSPD) using Euclidean distance
std::vector<std::tuple<Point, Point>> wspd(const Quadtree* tree, double s);

// WSPD using shortest-path distances from the Graph
std::vector<std::tuple<Point, Point>> wspd_spd(const Quadtree* tree, double s, Graph g);


#endif // QUADTREE_HPP
