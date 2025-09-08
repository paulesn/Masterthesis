//
// Created by kochda on 01.02.24.
//

#ifndef EXACT_SHIPROUTING_GEOMETRY_H
#define EXACT_SHIPROUTING_GEOMETRY_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Point_3 Point_3;
typedef K::Segment_2 Segment_2;
typedef K::Line_2 Line_2;

typedef u_int ID;
typedef std::pair<double, double> WGS84Coordinates;

namespace boost {
    namespace serialization {

        template<class Archive>
        void serialize(Archive & ar, Point_2 & p, const unsigned int version)
        {
            double x = p.x(), y = p.y();
            ar & boost::serialization::make_nvp("x", x);
            ar & boost::serialization::make_nvp("y", y);

            if (Archive::is_loading::value) {
                p = Point_2(x, y);
            }
        }

    } // namespace serialization
} // namespace boost


class Point {
private:
    Point_2 point;
public:
    friend class boost::serialization::access;
    Point(double x, double y);
    Point(Point_2& point);
    Point();
    bool operator==(const Point& p2) const;
    Point operator+(const Point& p2) const;
    friend std::ostream& operator<<(std::ostream& os, const Point& pt);
    double x() const;
    double y() const;
    double distance(const Point& p2) const;
    static double distance(const Point& p1, const Point& p2);
    double squaredDistance(const Point& p2) const;
    static double squaredDistance(const Point& p1, const Point& p2);
    static CGAL::Orientation orientation(const Point& p1, const Point& p2, const Point& checkPoint);
    // static Point intersectTwoLines(const Point& p1, const Point& p2, const Point& q1, const Point& q2);
    WGS84Coordinates pointInWGS84() const;
    Point_2 getCGALPoint() const;
    template<class Archive>
    void serialize(Archive &ar, const unsigned int version) {
        ar & this->point;
    }
};

typedef std::vector<Point> Points;

class ClosestPointSearch {
private:
    std::shared_ptr<Points> points;
    std::vector<std::vector<std::vector<ID>>> nodesPerCell;
    double xMin, yMin, xCellWidth, yCellWidth;
public:
    ClosestPointSearch(std::shared_ptr<Points> pointsPointer, int numCellsPerDimension, const std::vector<bool>& activePoints);
    ClosestPointSearch(std::shared_ptr<Points> pointsPointer, int numCellsPerDimension);
    void setPoints(std::shared_ptr<Points> pointsPointer, int numCellsPerDimension, const std::vector<bool>& activePoints);
    void setPoints(std::shared_ptr<Points> pointsPointer, int numCellsPerDimension);
    ID getClosestPoint(const Point& point) const;
    void printStats() const;
};

double distanceBetweenSegmentAndPoint(const Point& segmentStart, const Point& segmentEnd, const Point& point);
bool twoLinesIntersectingEachOther(const Point& l1Start, const Point& l1End, const Point& l2Start, const Point& l2End);
std::shared_ptr<Points> readPointsFromGraph(const std::string& filename);
std::shared_ptr<Points> readPoints(const std::string& filename);
void writePointsToFile(const std::shared_ptr<Points>& points, const std::string& filename);
std::vector<int> readPointsInTriangles(const std::string& filename);
std::vector<std::pair<ID, double>> readNearestNeighbors(const std::string& filename);
std::vector<std::vector<int>> readTriangleIntersections(const std::string& filename);
double haversineDistance(const WGS84Coordinates& p1, const WGS84Coordinates& p2);
#endif //EXACT_SHIPROUTING_GEOMETRY_H
