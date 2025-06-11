#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>
#include <string>
#include <iosfwd>
#include <set>
#include <unordered_set>

using namespace std;

struct Point {
    double x, y;
    Point(const double x_, const double y_) : x(x_), y(y_) {}

    // Overload the equality operator for Point
    bool operator==(const Point& other) const {
        return (x == other.x && y == other.y);
    }
    // Overload the less than operator for Point
    bool operator<(const Point& other) const {
        return (x < other.x || (x == other.x && y < other.y));
    }

    // Overload the output operator for Point
    friend std::ostream& operator<<(std::ostream& os, const Point& p) {
        os << "(" << p.x << ", " << p.y << ")";
        return os;
    }
};

struct PointHash {
    std::size_t operator()(const Point& p) const {
        return std::hash<int>()(p.x) ^ (std::hash<int>()(p.y) << 1);
    }
};

double euklidian_distance(Point a, Point b) {
    auto dx = abs(a.x - b.x);
    auto dy = abs(a.y - b.y);
    return sqrt(dx*dx + dy*dy);
}

struct QuadtreeNode {
    std::vector<Point> points;                                                     // A reference to all points in this subtree
    double cX, cY, height;                                                          // the center and height of the regular rectangle this node represents
    QuadtreeNode* NW;
    QuadtreeNode* NO;
    QuadtreeNode* SW;
    QuadtreeNode* SO;
    // the four children of this node
    bool is_leaf;                                                                   // for fast identification of leafs

    QuadtreeNode(const double cX_, const double cY_, const double height_)
        : cX(cX_), cY(cY_), height(height_) {
        is_leaf = true;
        NW = NO = SW = SO = nullptr;
    }

    // Deep-copy copy constructor
    QuadtreeNode(const QuadtreeNode& other)
      : points(other.points),
        cX(other.cX), cY(other.cY), height(other.height),
        is_leaf(other.is_leaf)
    {
        if (other.NW) NW = new QuadtreeNode(*other.NW); else NW = nullptr;
        if (other.NO) NO = new QuadtreeNode(*other.NO); else NO = nullptr;
        if (other.SW) SW = new QuadtreeNode(*other.SW); else SW = nullptr;
        if (other.SO) SO = new QuadtreeNode(*other.SO); else SO = nullptr;
    }

    // Deep-copy copy-assignment operator
    QuadtreeNode& operator=(const QuadtreeNode& other) {
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



    ~QuadtreeNode() {
        delete NW;
        delete NO;
        delete SW;
        delete SO;
    }

    // Public insert: count this new point exactly once
    void insert(Point s) {
        points.emplace_back(s);
        insertInternal(s);
    }

    // Check if rectangle contains (x,y)
    [[nodiscard]] bool area_contains(const Point s) const {
        const double x = s.x;
        const double y = s.y;
        return x >= (cX - height) && x <= (cX + height)
            && y >= (cY - height) && y <= (cY + height);
    }

    [[nodiscard]] bool intercept_rect(const double x, const double y, const double h) const {
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

    [[nodiscard]] std::string rec_string(const int level = 0) const {
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


    [[nodiscard]] int getDepth(int call_level) const {
        if (is_leaf) return call_level;
        return max({
            NW->getDepth(call_level + 1),
            NO->getDepth(call_level + 1),
            SW->getDepth(call_level + 1),
            SO->getDepth(call_level + 1)
            });
    }

private:

    // Internal helper — does NOT increment totalObjects
    void insertInternal(const Point s) {
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
            NW = new QuadtreeNode(cX + h2, cY + h2, h2);
            NO = new QuadtreeNode(cX + h2, cY - h2, h2);
            SW = new QuadtreeNode(cX - h2, cY + h2, h2);
            SO = new QuadtreeNode(cX - h2, cY - h2, h2);
            is_leaf = false;
            // reinsert the points to let them flow down a layer
            for (const auto pt: points) {
                insertInternal(pt);
            }
        }
    }
};

struct Quadtree {
    QuadtreeNode root;
    double height;
    std::unordered_set<Point, PointHash> pointSet;

    explicit Quadtree(const double height_=1)
        : root(0,0, height_), height(height_), pointSet() {}
    ~Quadtree() = default;

    void insert(const Point p) {
        if (pointSet.find(p) != pointSet.end()) {
            return; // point already in tree, do not insert again
        }
        pointSet.insert(p);
        if (root.area_contains(p)) {
            root.insert(p);
        } else {
            std::cout << "Point (" << p.x << ", " << p.y << ") is outside the quadtree bounds." << std::endl;
            const auto points = (root.points);
            // create new quadtree large enough to contain the point
            const double newHeight = std::max(height, std::sqrt(p.x * p.x + p.y * p.y));
            auto newRoot = QuadtreeNode(0.0, 0.0, newHeight);
            height = newHeight;
            newRoot.insert(p);
            root = newRoot;
            for (const auto old_point: points) {insert(old_point);}
        }
    }

    void print() { std::cout << root.rec_string(); }

    int getDepth() {
        return root.getDepth(1);
    }

    auto normalize_pair(QuadtreeNode* a, QuadtreeNode* b) {
        if (a->cX < b->cX || (a->cX == b->cX && a->cY < b->cY)) {
            return std::make_tuple(a, b);
        } else {
            return std::make_tuple(b, a);
        }
    }

    /**
     *
     * @param s
     * @param wsp_check
     * @return
     */
    auto wspd(
            const double s
            //bool (*wsp_check)(vector<Point>*, vector<Point>*, double) = euklidian
        ) {
        const int depth = getDepth(); // the depth of the tree
        const int cells = pow(4, depth)*depth;
        const int cell_pairs = cells*cells;

        std::vector<std::tuple<Point, Point>> pairs;
        std::vector<std::tuple<QuadtreeNode*, QuadtreeNode*>> parings_todo;
        std::set<std::pair<QuadtreeNode*, QuadtreeNode*>> seen;

        parings_todo.reserve(cell_pairs);
        parings_todo.emplace_back(&root, &root);


        int skipped_because_nlptr = 0;
        int skipped_because_duplicate = 0;
        int skipped_add_because_nlptr = 0;

        //for (int i=0; i < cell_pairs; i++) {
        while (true){

            // loop end condition
            if (parings_todo.empty()) break;

            // get next pair
            auto self = get<0>(parings_todo.back());
            auto other = get<1>(parings_todo.back());

            parings_todo.pop_back();

            // defensive null check
            if (other == nullptr) {skipped_because_nlptr++; continue;}
            if (self == nullptr) {skipped_because_nlptr++; continue;}

            //normalize pair
            auto pair = normalize_pair(self, other);

            // check if we already processed this pair
            if (seen.find(std::make_pair(self, other)) != seen.end()) {
                skipped_because_duplicate++;
                continue;
            }
            seen.insert(std::make_pair(self, other));


            // check if one of the subtrees is empty, if yes return nothing
            if (self->points.empty() || other->points.empty()) {
                continue;
            }

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
                if (self != other) {
                    pairs.emplace_back(self->points[0], other->points[0]);
                }
                continue;
            }

            // Case 2 - nodes are well separated
            // TODO easier wspd_check remove later
            auto radius = sqrt(self->height*self->height + self->height/2*self->height/2);
            auto other_radius = sqrt(other->height*other->height + other->height/2*other->height/2);
            auto dist = euklidian_distance(Point(self->cX, self->cY), Point(other->cX, other->cY));
            if (max(radius, other_radius)*s < dist) {
                pairs.emplace_back(self->points[0], other->points[0]);
                continue;
            }

            //if (wsp_check(&(self->points), &(other->points), s)) {
            //    pairs.emplace_back(self->points[0], other->points[0]);
            //    continue;
            //}

            // Choose larger for further splitting
            QuadtreeNode* larger = (self->height > other->height) ? self : other;
            QuadtreeNode* smaller = (self->height > other->height) ? other : self;
            if (larger->is_leaf) std::swap(larger, smaller);

            auto lNW = larger->NW;
            if (lNW == nullptr) {
                skipped_add_because_nlptr++;
                continue;
            }
            parings_todo.emplace_back(smaller, lNW);

            auto lNO = larger->NO;
            if (lNO == nullptr) {
                skipped_add_because_nlptr++;
                continue;
            }
            parings_todo.emplace_back(smaller, lNO);

            auto lSW = larger->SW;
            if (lSW == nullptr) {
                skipped_add_because_nlptr++;
                continue;
            }
            parings_todo.emplace_back(smaller, lSW);

            auto lSO = larger->SO;
            if (lSO == nullptr) {
                skipped_add_because_nlptr++;
                continue;
            }
            parings_todo.emplace_back(smaller, lSO);
        }
        cout << pairs.size() << endl;
        cout << parings_todo.size() << endl;
        cout << "Skipped because null pointer: " << skipped_because_nlptr << endl;
        cout << "Skipped because duplicate: " << skipped_because_duplicate << endl;
        cout << "Skipped because null pointer in add: " << skipped_add_because_nlptr << endl;
        return pairs;
    }

};
