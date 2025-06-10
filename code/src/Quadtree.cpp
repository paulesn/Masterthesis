#include <iostream>
#include <vector>
#include <cmath>
#include <csignal>
#include <tuple>
#include <string>
#include <functional>
#include <iosfwd>
#include <iosfwd>
#include <vector>
#include <vector>

using namespace std;

struct Point {
    double x, y;
    Point(const double x_, const double y_) : x(x_), y(y_) {}
};

double euklidian_distance(Point a, Point b) {
    auto dx = abs(a.x - b.x);
    auto dy = abs(a.y - b.y);
    return sqrt(dx*dx + dy*dy);
}

/**
 * TODO maybe not calculating the correct circle
 * @param points a vector of points
 * @return the radius and center of a circle enclosing the points
 */
tuple<Point, double> euklidian_radius(vector<Point> points) {
    double dim_max = 0;
    double centerX = points[0].x;
    double centerY = points[0].y;
    for (int i = 0; i < points.size(); i++) {
        for (int j = i + 1; j < points.size(); j++) {
            auto dim = euklidian_distance(points[i], points[j]);
            if (dim > dim_max) {
                dim_max = dim;
                centerX = abs(points[i].x-points[j].x);
                centerY = abs(points[i].y-points[j].y);
            }
        }
    }
    auto rad = dim_max/2;
    return make_tuple(Point(centerX, centerY), rad);
}

bool euklidian(vector<Point> pA, vector<Point> pB, double s) {
    auto A = euklidian_radius(pA);
    auto radA = get<1>(A);
    auto centerA = get<0>(A);
    auto B = euklidian_radius(pB);
    auto radB = get<1>(B);
    auto centerB = get<0>(B);

    auto dist = euklidian_distance(centerA, centerB);
    auto max_rad = max(abs(radA), abs(radB));
    return dist-(radA+radB) < max_rad*s;

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
        for (auto p : points) {
            if (p.x == s.x && p.y == s.y) {
                cout<<"Already added" << endl;
                return;
            }
        } // do not add any point twice
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


    std::vector<std::tuple<Point, Point>> rec_wspd(
        QuadtreeNode *other,
        const double s,
        double (*dist_function)(Point a, Point b) = euklidian_distance,
        tuple<Point, double> (*radius_function)(vector<Point>) = euklidian_radius
    ){
        // defensive null check
        if (other == nullptr) return {};

        // check if one of the subtrees is empty, if yes return nothing
        if (points.empty() || other->points.empty()) {
            cout << "empty subttree" << endl;
            return {};
        }

        // if both of the nodes are leaf nodes, it is always well separated
        if (is_leaf && other->is_leaf) {
            // if the two nodes are the same return nothing
            if (points[0].x == other->points[0].x && points[0].y == other->points[0].y) {return {};}

            //TODO remove security check if any is empty
            if (points.empty() || other->points.empty()) {return {};}

            cout << "Both leaf nodes reached" << endl;

            return { std::make_tuple(points[0], other->points[0]) };
        }

        //self-case
        if (this ==  other) {
            // return all pairs of children
            std::vector<std::tuple<Point, Point>> pairs;
            for (const auto child1: {NW,NO,SW,SO}) {
                if (!child1->points.empty()) { // defensive code to prevent to many checks
                    for (const auto child2: {NW,NO,SW,SO}) {
                        if (!child2->points.empty()) {
                            auto temp_pairs = child1->rec_wspd(child2, s, dist_function, radius_function);
                            if (!temp_pairs.empty()) pairs.insert(pairs.end(), temp_pairs.begin(),temp_pairs.end());
                        }
                    }
                }
            }
        };

        const double dist = dist_function(Point(cX,cY), Point(other->cX,other->cY));
        const double rad = std::max(get<1>(radius_function(points)), get<1>(radius_function(other->points)));
        if (dist >= s * rad) {
            cout << "WSP found with dist = " << dist << endl;
            std::vector<std::tuple<Point, Point>> pairs;
            pairs.emplace_back(points[0], other->points[0]);
        }

        // Choose larger for further splitting
        QuadtreeNode* larger = (height > other->height) ? this : other;
        QuadtreeNode* smaller = (height > other->height) ? other : this;
        if (larger->is_leaf) std::swap(larger, smaller);

        std::vector<std::tuple<Point, Point>> pairs;

        auto lNW = larger->NW;
        auto lnw_pairs = smaller->rec_wspd(lNW, s, dist_function, radius_function);
        pairs.insert(pairs.end(), lnw_pairs.begin(),lnw_pairs.end());

        auto lNO = larger->NO;
        auto lno_pairs = smaller->rec_wspd(lNO, s, dist_function, radius_function);
        pairs.insert(pairs.end(), lno_pairs.begin(),lno_pairs.end());

        auto lSW = larger->SW;
        auto lsw_pairs = smaller->rec_wspd(lSW, s, dist_function, radius_function);
        pairs.insert(pairs.end(), lsw_pairs.begin(),lsw_pairs.end());

        auto lSO = larger->SO;
        auto lso_pairs = smaller->rec_wspd(lSO, s, dist_function, radius_function);
        pairs.insert(pairs.end(), lso_pairs.begin(),lso_pairs.end());

        return pairs;
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

    explicit Quadtree(const double height_=1)
        : root(0,0, height_), height(height_) {}
    ~Quadtree() = default;

    void insert(const Point p) {
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

    auto rec_wspd(const double s) {
        return root.rec_wspd(&root, s);
    }

    /**
     *
     * @param s
     * @param wsp_check
     * @return
     */
    auto wspd(
            const double s,
            bool (*wsp_check)(vector<Point>, vector<Point>, double) = euklidian
        ) {
        const int depth = getDepth(); // the depth of the tree
        const int cells = 4^(depth)+1;
        const int cell_pairs = cells*cells;

        std::vector<std::tuple<Point, Point>> pairs;
        std::vector<std::tuple<QuadtreeNode*, QuadtreeNode*>> parings_todo;
        parings_todo.reserve(cells);
        parings_todo.emplace_back(&root, &root);

        //for (int i=0; i < cell_pairs; i++) {
        while (true){
            cout << "----------------------------------------------------------------"  << endl;
            if (parings_todo.empty()) break;

            auto self = get<0>(parings_todo.back());
            auto other = get<1>(parings_todo.back());
            parings_todo.pop_back();

            // defensive null check
            if (other == nullptr) continue;
            if (self == nullptr) continue;

            // check if one of the subtrees is empty, if yes return nothing
            if (self->points.empty() || other->points.empty()) {
                cout << "empty subttree" << endl;
                continue;
            }

            // if both of the nodes are leaf nodes, it is always well separated
            if (self->is_leaf && other->is_leaf) {
                // if the two nodes are the same return nothing
                if (self->points[0].x == other->points[0].x && self->points[0].y == other->points[0].y) {continue;}

                //TODO remove security check if any is empty
                if (self->points.empty() || other->points.empty()) {continue;}

                cout << "Both leaf nodes reached" << endl;
                pairs.emplace_back(self->points[0], other->points[0]);
            }

            //self-case
            if (self == other) {
                cout << "Self Case" << endl;
                // add all pairs of children
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
            };

            if (wsp_check(self->points, other->points, s)) {
                cout << "WSP found with dist" << endl;
                pairs.emplace_back(self->points[0], other->points[0]);
                continue;
            }

            // Choose larger for further splitting
            QuadtreeNode* larger = (self->height > other->height) ? self : other;
            QuadtreeNode* smaller = (self->height > other->height) ? other : self;
            if (larger->is_leaf) std::swap(larger, smaller);


            cout << "no wspd found -> keep separating" << endl;
            auto lNW = larger->NW;
            parings_todo.emplace_back(smaller, lNW);

            auto lNO = larger->NO;
            parings_todo.emplace_back(smaller, lNO);

            auto lSW = larger->SW;
            parings_todo.emplace_back(smaller, lSW);

            auto lSO = larger->SO;
            parings_todo.emplace_back(smaller, lSO);
        }
        cout << pairs.size() << endl;
        cout << parings_todo.size() << endl;
        return pairs;
    }

};
