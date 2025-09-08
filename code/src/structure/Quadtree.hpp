#ifndef QUADTREE_HPP
#define QUADTREE_HPP

#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <unordered_set>
#include "Graph.hpp" // You may want to replace this with Graph.hpp if possible

/**
 * Hash function for Point
 * for use in unordered_set
 */
struct PointHash {
    std::size_t operator()(const Pointc& p) const;
};

/**
 * Euclidean distance between two points
 */
double euklidian_distance(Pointc a, Pointc b);

/**
 * Node of the Quadtree with defind 4 or 0 children
 */
struct QuadtreeNode {
    /**
     * A list of points in this cell. This list is never cleared so that each cell knows all points that were inserted into it.
     * TODO evaluate if this is necessary, maybe we can just use the points in the leaf nodes.
     */
    std::vector<Pointc> points;
    /**
     * The position and size information of the node.
     * cX and cY are the center coordinates of the node, height is half the width/height of the node;
     * Effectivly the distance from the center to any edge.
     */
    double cX, cY, height;
    /**
    * The four children of the node, if it is not a leaf.
    * NW = North-West, NO = North-East, SW = South-West, SO = South-East (derived from German)
    * The children are nullptr if the node is a leaf.
    *
    */
    QuadtreeNode *NW, *NO, *SW, *SO;

    /**
     * The information if the node is a leaf or not. True == no children, false == has children.
     */
    bool is_leaf;
    /**
     * Pointer to the parent node, if needed. This is useful for traversing back up the tree.
     * It is used, for example, in nearest neighbor search or when calculating the well-separated pair decomposition (WSPD).
     */
    QuadtreeNode *parent = nullptr; // Pointer to the parent node, if needed

    /**
     * The longest shortest path between any two points in this node.
     * this is calculated on demand and is used for the WSPD with shortest path distance.
     * It is stored here to avoid recalculating it multiple times.
     */
    double inner_sp_distance = -1; // Used for WSPD with shortest path distance

    /**
     * "Default" constructor for a QuadtreeNode.
     * @param cX_ center of the new node in x direction
     * @param cY_ center of the new node in y direction
     * @param height_ height of the root node, which is half the width/height of the node
     * @param parent optional pointer to the parent node, default is nullptr
     */
    QuadtreeNode(double cX_, double cY_, double height_, QuadtreeNode* parent=nullptr);

    /**
     * Copy constructor for a QuadtreeNode.
     * @param other the QuadtreeNode to copy from. As the children are pointers, they are copied as well but reference still the old children.
     */
    QuadtreeNode(const QuadtreeNode& other);
    
    QuadtreeNode& operator=(const QuadtreeNode& other);
    ~QuadtreeNode();

    void insert(Pointc s);
    [[nodiscard]] bool area_contains(Pointc s) const;
    [[nodiscard]] bool intercept_rect(double x, double y, double h) const;
    [[nodiscard]] std::string rec_string(int level = 0) const;
    [[nodiscard]] std::vector<Pointc> get_all_points() const;

    [[nodiscard]] bool contains(const Pointc & p) const;

    std::vector<int> circle_intersect(double x, double y, double r) const;

private:

    /**
     * this function returns all points that are in the area between two rays taking the source point as origin.
     * and the rays are defined by the angles angle_a and angle_b. relative to the x-axis.
     * between the two rays it is always assumed that angle_a > angle_b.
     * @param source_x the x coordinate of the source point
     * @param source_y the y coordinate of the source point
     * @param angle_a the angle of the first ray
     * @param angle_b the angle of the second ray
     * @return a vector of node ids that are in the area between the two rays.
     */
    bool internal_angles_intersect(Pointc source, double angle_a, double angle_b) const ;

    void insertInternal(Pointc s);
};

// Quadtree structure
struct Quadtree {
    QuadtreeNode root;
    double height;
    std::unordered_set<Pointc, PointHash> pointSet;

    explicit Quadtree(double height_ = 1);
    ~Quadtree() = default;



    bool insert(Pointc p);
    void print() const;
    std::vector<Pointc> get_all_points() const;
    bool contains(Pointc p) const{return root.contains(p);}


    std::vector<int> circle_intersect(double x, double y, double r) const {return root.circle_intersect(x,y,r);}
};

// Normalize the ordering of two QuadtreeNodes
std::tuple<QuadtreeNode*, QuadtreeNode*> normalize_pair(QuadtreeNode* a, QuadtreeNode* b);

// Well-separated pair decomposition (WSPD) using Euclidean distance
Graph wspd(const Quadtree* tree, double s, Graph* g);

// WSPD using shortest-path distances from the Graph
Graph wspd_spd(const Quadtree* tree, double s, Graph g);



#endif // QUADTREE_HPP
