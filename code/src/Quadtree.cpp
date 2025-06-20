#include "Quadtree.hpp"

#include <cassert>
#include <chrono>
#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>
#include <string>
#include <iosfwd>
#include <set>
#include <unordered_set>
#include "Graph.hpp"


using namespace std;

std::size_t PointHash::operator()(const Point &p) const {
    return std::hash<int>()(p.x) ^ (std::hash<int>()(p.y) << 1);
}

double euklidian_distance(Point a, Point b) {
    auto dx = abs(a.x - b.x);
    auto dy = abs(a.y - b.y);
    return sqrt(dx*dx + dy*dy);
}

::QuadtreeNode::QuadtreeNode(double cX_, double cY_, double height_, QuadtreeNode* parent_): cX(cX_), cY(cY_), height(height_), parent(parent_) {
    is_leaf = true;
    NW = NO = SW = SO = nullptr;
}

::QuadtreeNode::QuadtreeNode(const QuadtreeNode &other): points(other.points),
                                                         cX(other.cX), cY(other.cY), height(other.height),
                                                         is_leaf(other.is_leaf) {
    if (other.NW) NW = new QuadtreeNode(*other.NW); else NW = nullptr;
    if (other.NO) NO = new QuadtreeNode(*other.NO); else NO = nullptr;
    if (other.SW) SW = new QuadtreeNode(*other.SW); else SW = nullptr;
    if (other.SO) SO = new QuadtreeNode(*other.SO); else SO = nullptr;
    parent = other.parent; // Copy the parent pointer
    //cout << "Copy constructor called for QuadtreeNode at (" << cX << ", " << cY << ") with height " << height << endl;
}

::QuadtreeNode & ::QuadtreeNode::operator=(const QuadtreeNode &other) {
    if (this == &other) return *this;

    // 1) Clean up existing children
    delete NW; delete NO; delete SW; delete SO;

    // 2) Copy primitive/data members
    points    = other.points;
    cX        = other.cX;
    cY        = other.cY;
    height    = other.height;
    is_leaf   = other.is_leaf;

    // 3) Clone children if they exist
    if (other.NW) NW = new QuadtreeNode(*other.NW); else NW = nullptr;
    if (other.NO) NO = new QuadtreeNode(*other.NO); else NO = nullptr;
    if (other.SW) SW = new QuadtreeNode(*other.SW); else SW = nullptr;
    if (other.SO) SO = new QuadtreeNode(*other.SO); else SO = nullptr;

    return *this;
}

::QuadtreeNode::~QuadtreeNode() {
    delete NW;
    delete NO;
    delete SW;
    delete SO;
}

void ::QuadtreeNode::insert(Point s) {
    points.emplace_back(s);
    insertInternal(s);
}

bool ::QuadtreeNode::area_contains(Point s) const {
    const double x = s.x;
    const double y = s.y;
    return x >= (cX - height) && x <= (cX + height)
           && y >= (cY - height) && y <= (cY + height);
}

bool ::QuadtreeNode::intercept_rect(double x, double y, double h) const {
    // return true if the cell contains a corner of the other rect
    if (area_contains(Point(x+h, y+h))) return true;
    if (area_contains(Point(x+h, y-h))) return true;
    if (area_contains(Point(x-h, y+h))) return true;
    if (area_contains(Point(x-h, y-h))) return true;

    // check if the other rect contains one corner of the cell
    if (x - h > cX + height || x + h < cX - height) return false;
    if (y - h > cY + height || y + h < cY - height) return false;

    return true;
}

std::string QuadtreeNode::rec_string(const int level) const {
    std::string text;
    for (int i = 0; i < level; ++i) text += "|";
    if (is_leaf) {
        if (points.empty()) text += "0\n";
        else text += "1\n";
    }
    else {
        text += "-[" + std::to_string(points.size()) + "]\n";

        text += NW->rec_string(level + 1);
        text += NO->rec_string(level + 1);
        text += SW->rec_string(level + 1);
        text += SO->rec_string(level + 1);
    }
    return text;
}

vector<Point> QuadtreeNode::get_all_points() const {
    if (is_leaf) {
        return points;
    }
    vector<Point> all_points;
    if (NW) {
        auto nw_points = NW->get_all_points();
        all_points.insert(all_points.end(), nw_points.begin(), nw_points.end());
    }
    if (NO) {
        auto no_points = NO->get_all_points();
        all_points.insert(all_points.end(), no_points.begin(), no_points.end());
    }
    if (SW) {
        auto sw_points = SW->get_all_points();
        all_points.insert(all_points.end(), sw_points.begin(), sw_points.end());
    }
    if (SO) {
        auto so_points = SO->get_all_points();
        all_points.insert(all_points.end(), so_points.begin(), so_points.end());
    }
    return all_points;
}

bool QuadtreeNode::contains(const Point &p) const {
    if (!is_leaf) {
        if (NW->area_contains(p)) return NW->contains(p);
        if (NO->area_contains(p)) return NO->contains(p);
        if (SW->area_contains(p)) return SW->contains(p);
        if (SO->area_contains(p)) return SO->contains(p);
        return false;
    }
    if (points.empty()) {
        return false; // Leaf is empty, no points to contain
    } else {
        // Check if the point is in the current leaf's points
        for (const auto& pt : points) {
            if (pt == p) {
                return true; // Point found in this leaf
            }
        }
        return false; // Point not found in this leaf
    }
}

void QuadtreeNode::insertInternal(Point s) {
    // 1) If already split, forward to exactly one child
    if (!is_leaf) {
        if (NW->area_contains(s)) {NW->insert(s); return;}
        if (NO->area_contains(s)) {NO->insert(s); return;}
        if (SW->area_contains(s)) {SW->insert(s); return;}
        if (SO->area_contains(s)) {SO->insert(s); return;}

        cout << "No fitting subarea found THIS IS AN ERROR" << endl;
        return;
    }

    //We are in a leaf, but we do not have to store anything as this already happened in insert()
    if (points.size() > 1) {
        //split the node
        const double h2 = height * 0.5;
        NW = new QuadtreeNode(cX + h2, cY + h2, h2, this);
        NO = new QuadtreeNode(cX + h2, cY - h2, h2, this);
        SW = new QuadtreeNode(cX - h2, cY + h2, h2, this);
        SO = new QuadtreeNode(cX - h2, cY - h2, h2, this);
        is_leaf = false;
        // reinsert the points to let them flow down a layer
        for (const auto pt: points) {
            insertInternal(pt);
        }
        // points.clear(); //TODO DO I NEED THE CLEAR? clear the points as they are now stored in the children
    }
}


::Quadtree::Quadtree(double height_): root(0,0, height_, nullptr), height(height_), pointSet() {}

bool ::Quadtree::insert(Point p) {
    static int counter = 0;
    //std::cout << counter << std::endl;
    counter++;
    auto it = pointSet.find(p);
    if (it != pointSet.end()) {
        auto found = *it;
        //std::cout << "Point (" << p.x << ", " << p.y << ") already exists in the quadtree as point (" << found.x << ", " << found.y << ")." << std::endl;
        return false; // point already in tree, do not insert again
    }
    pointSet.insert(p);
    if (root.QuadtreeNode::area_contains(p)) {
        root.QuadtreeNode::insert(p);
    } else {
        std::cout << "Point (" << p.x << ", " << p.y << ") is outside the quadtree bounds." << std::endl;
        const auto points = (root.QuadtreeNode::get_all_points());
        // create new quadtree large enough to contain the point
        const double newHeight = std::max(height, std::sqrt(p.x * p.x + p.y * p.y));
        auto newRoot = QuadtreeNode(0.0, 0.0, newHeight, nullptr);
        height = newHeight;
        newRoot.insert(p);
        root = newRoot;
        for (const auto old_point: points) {root.insert(old_point);}
    }
    return true; // point successfully inserted
}

void Quadtree::print() const { std::cout << root.QuadtreeNode::rec_string(); }

std::vector<Point> Quadtree::get_all_points() const {
    return root.get_all_points();
}

/**
 *
 * @param a a quadtree node
 * @param b another quadtree node
 * @return a tuple of the two nodes, normalized such that a is always the one with smaller cX or cY
 */
std::tuple<QuadtreeNode*, QuadtreeNode*> normalize_pair(QuadtreeNode* a, QuadtreeNode* b) {
    if (a->cX < b->cX || (a->cX == b->cX && a->cY < b->cY)) {
        return std::make_tuple(a, b);
    } else {
        return std::make_tuple(b, a);
    }
}

/**
 * a well separated pair decomposition (WSPD) algorithm that uses the euklidian distance to check if two sets of points are well-separated.
 */
Graph wspd(const Quadtree* tree, const double s, Graph* g) {

    auto root = tree->root;
    Graph spanner = Graph(g->id_point_map.size());

    std::vector<std::tuple<QuadtreeNode*, QuadtreeNode*>> parings_todo;
    std::set<std::pair<QuadtreeNode*, QuadtreeNode*>> seen;
    parings_todo.emplace_back(&root, &root);


    int skipped_because_nlptr = 0;
    int skipped_because_duplicate = 0;
    int skipped_add_because_nlptr = 0;

    int counter = 0;
    while (true){
        counter++;
        if (counter %100000 == 0) {
            cout << "Processed " << counter << " pairs. Current size of parings_todo: " << parings_todo.size() << "                                      \r";
        }

        // loop end condition
        if (parings_todo.empty()) break;

        // get next pair
        auto self = get<0>(parings_todo.back());
        auto other = get<1>(parings_todo.back());

        parings_todo.pop_back();


        // check the self case
        if (self == other) {
            if (self->is_leaf) {
                // if the two nodes are the same return nothing
                continue;
            }
            parings_todo.emplace_back(self->NW, self->NW);
            parings_todo.emplace_back(self->NW, self->NO);
            parings_todo.emplace_back(self->NW, self->SW);
            parings_todo.emplace_back(self->NW, self->SO);

            parings_todo.emplace_back(self->NO, self->NO);
            parings_todo.emplace_back(self->NO, self->SW);
            parings_todo.emplace_back(self->NO, self->SO);

            parings_todo.emplace_back(self->SW, self->SW);
            parings_todo.emplace_back(self->SW, self->SO);

            parings_todo.emplace_back(self->SO, self->SO);
            continue;
        }

        // check for well separation
        // Case 1 - both nodes are leaf nodes
        if (self->is_leaf && other->is_leaf) {
            // if the two nodes are the same return nothing
            if (self != other && !self->points.empty() && !other->points.empty()) {
                spanner.addEdge(self->points[0].id, other->points[0].id, euklidian_distance(self->points[0], other->points[0]));
            }
            continue;
        }

        // Case 2 - nodes are well separated
        // use the default euklidian distance check
        double radius = 0.0;
        double other_radius = 0.0;
        if (self->is_leaf) radius = sqrt(self->height*self->height)*sqrt(2);
        if (other->is_leaf) other_radius = sqrt(other->height*other->height)*sqrt(2);
        auto dist = euklidian_distance(Point(self->cX, self->cY), Point(other->cX, other->cY));
        if (max(radius, other_radius)*s < dist) {
            if (!self->points.empty() && !other->points.empty()) {
                // if the two nodes are well separated, add the first point of each set as a pair
                spanner.addEdge(self->points[0].id, other->points[0].id, euklidian_distance(self->points[0], other->points[0]));
            }
            continue;
        }



        // Choose larger for further splitting
        QuadtreeNode* larger = (self->height > other->height) ? self : other;
        QuadtreeNode* smaller = (self->height > other->height) ? other : self;
        if (larger->is_leaf) std::swap(larger, smaller);

        for (auto* child : {larger->NW, larger->NO, larger->SW, larger->SO}) {
            if (child) parings_todo.emplace_back(smaller, child);
        }
    }
    cout << parings_todo.size() << endl;
    cout << "Skipped because null pointer: " << skipped_because_nlptr << endl;
    cout << "Skipped because duplicate: " << skipped_because_duplicate << endl;
    cout << "Skipped because null pointer in add: " << skipped_add_because_nlptr << endl;
    return spanner;
}

/**
 * This function is a variant of wspd that uses the shortest path distance to check if two sets of points are well-separated.
 * It is used for the WSPD algorithm.
 * @param s
 * @return
 */
Graph wspd_spd(const Quadtree* tree, const double s, Graph g) {

    auto root = tree->root;
    Graph spanner = Graph(g.id_point_map.size());// new graph with same nodes

    //std::vector<std::tuple<Point, Point>> pairs; NO longer needed
    std::vector<std::tuple<QuadtreeNode*, QuadtreeNode*>> parings_todo;
    std::set<std::pair<QuadtreeNode*, QuadtreeNode*>> seen;
    parings_todo.emplace_back(&root, &root);


    int skipped_because_nlptr = 0;
    int skipped_because_duplicate = 0;
    int skipped_add_because_nlptr = 0;

    int counter = 0;
    while (true){
        counter++;
        if (counter %100000 == 0) {
            cout << "Processed " << counter << " pairs. Current size of parings_todo: " << parings_todo.size() << "\r";
        }

        // loop end condition
        if (parings_todo.empty()) break;

        // get next pair
        auto self = get<0>(parings_todo.back());
        auto other = get<1>(parings_todo.back());

        parings_todo.pop_back();


        // check the self case
        if (self == other) {
            if (self->is_leaf) {
                // if the two nodes are the same return nothing
                continue;
            }
            parings_todo.emplace_back(self->NW, self->NW);
            parings_todo.emplace_back(self->NW, self->NO);
            parings_todo.emplace_back(self->NW, self->SW);
            parings_todo.emplace_back(self->NW, self->SO);

            parings_todo.emplace_back(self->NO, self->NO);
            parings_todo.emplace_back(self->NO, self->SW);
            parings_todo.emplace_back(self->NO, self->SO);

            parings_todo.emplace_back(self->SW, self->SW);
            parings_todo.emplace_back(self->SW, self->SO);

            parings_todo.emplace_back(self->SO, self->SO);
            continue;
        }

        // check for well separation
        // Case 1 - both nodes are leaf nodes
        auto temp = g.dijkstra(698, 686, -1);

        if (self->is_leaf && other->is_leaf) {
            // if the two nodes are the same return nothing
            if (self != other && !self->points.empty() && !other->points.empty()) {
                // TODO replace with hublabel query
                spanner.addEdge(self->points[0].id,other->points[0].id, get<1>(g.dijkstra(self->points[0].id,other->points[0].id)));
            }
            continue;
        }

        // Case 2 - nodes are well separated
        auto set1 (self->points);
        auto set2 (other->points);

        auto result = g.wspdCheck(set1, set2,self->inner_sp_distance, other->inner_sp_distance, s );
        auto path = std::get<0>(result);
        double dist = std::get<1>(result);

        if (!path.empty()){
            // TODO CHOOSE CORRECT SNIPPET
            /* this code snippet would have added the shortest path between the two sets of points
            for (int i = 0; i< path.size()-1; i++) {
                auto p1 = Point(0,0, path[i]); // Attention: the points loose their coordiantes in this process
                auto p2 = Point(0,0, path[i]);
                pairs.emplace_back(p1, p2);
            }
            */
            // instead this code snipped adds the start and end point of the path as direct edge
            auto p1 = g.id_point_map[0]; // Attention: the points loose their coordiantes in this process
            auto p2 = g.id_point_map[path.size()-1];
            spanner.addEdge(p1.id,p2.id, dist, true);

            continue;
        }

        //----------------------------



        // Choose larger for further splitting
        QuadtreeNode* larger = (self->height > other->height) ? self : other;
        QuadtreeNode* smaller = (self->height > other->height) ? other : self;
        if (larger->is_leaf) std::swap(larger, smaller);

        for (auto* child : {larger->NW, larger->NO, larger->SW, larger->SO}) {
            if (child) parings_todo.emplace_back(smaller, child);
        }
    }
    cout << parings_todo.size() << endl;
    cout << "Skipped because null pointer: " << skipped_because_nlptr << endl;
    cout << "Skipped because duplicate: " << skipped_because_duplicate << endl;
    cout << "Skipped because null pointer in add: " << skipped_add_because_nlptr << endl;
    return spanner;

}
