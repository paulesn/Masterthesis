//
// Created by kochda on 31.01.24.
//

#ifndef EXACT_SHIPROUTING_TRIANGULATION_H
#define EXACT_SHIPROUTING_TRIANGULATION_H

#include <memory>
#include <vector>
#include "Geometry.h"
#include <queue>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/array.hpp>
#include <boost/archive/binary_oarchive.hpp>

using namespace std;

class Comparator {
public:
    bool operator() (std::tuple<std::pair<ID, ID>, int, int, double>& a, std::tuple<std::pair<ID, ID>, int, int, double>& b) const {
        return get<3>(a) > get<3>(b);
    }
};

typedef ID InternalID;
typedef ID GlobalID;
typedef int TriangleIndex;
typedef std::array<TriangleIndex, 3> Neighbors;

class Triangulation{
private:

public:
    friend class boost::serialization::access;
    class Triangle {
    private:
        std::array<ID, 3> triangleIndices;
        Triangulation* parentTriangulation;
    public:
        friend class boost::serialization::access;
        friend class Triangulation;
        Triangle();
        Triangle(const std::array<ID, 3>& indices, Triangulation* parent): triangleIndices(indices), parentTriangulation(parent) {};
        Triangle(ID x, ID y, ID z, Triangulation* parent): parentTriangulation(parent) {this->triangleIndices = {x,y,z};};
        bool isCCW();
        ID getID(int i) const;
        std::array<ID, 3> getTriangleIndices() const;
        int getIndex(ID pointID) const;
        int indexEdge(ID p1, ID p2) const;
        static int CCW(int i);
        static int CW(int i);
        bool pointInsideTriangle(const Point& point) const;
        double maxSideLength() const;
        double minDistanceToTriangle(const Point& point) const;
        double maxDistanceToTriangle(const Point& point) const;
        bool lineIntersectsTriangle(const Point& startLine, const Point& endLine) const;
        bool intersectsAnotherTriangle(const Triangulation& otherTriang, TriangleIndex triangleID) const;
        double area() const;
        template <class Archive> void serialize(Archive &ar, const unsigned int version);
    };
    // Constructors
    Triangulation();
    Triangulation(const std::string triangulationFile, const std::string& pointsFile);
    Triangulation(const std::string triangulationFile, std::shared_ptr<Points> points);
    Triangulation(const std::string triangulationFile);

    // Functionality
    bool pointInTriangulation(ID pointID) const;
    bool directVisibility(InternalID p1, InternalID p2) const;
    bool directVisibility(const Point& p1, const Point& p2, const ClosestPointSearch& cps) const;
    std::vector<GlobalID> oneToAllVisibility(const InternalID startPoint, const bool steinerPointsAsTargets,
        const std::vector<bool>& markedTriangles = std::vector<bool>()) const;
    std::vector<GlobalID> oneToAllVisibility(const Point& startPoint, const ClosestPointSearch& closestPointSearch,
        const bool steinerPointsAsTargets) const;
    ID computeSpannerNode(ID startPoint, const Point& leftEndCone, const Point& rightEndCone);
    ID computeSpannerNode(const Point& startPoint, const Point& leftEndCone, const Point& rightEndCone, const ClosestPointSearch& cps);
    ID computeSpannerNode(ID startPoint, const Point& leftEndCone, const Point& rightEndCone, const std::vector<bool>& shrunkNodes);
    void createVisibilityGraph(std::vector<std::vector<GlobalID>> &edges, std::vector<std::vector<double>> &costs,
        const bool steinerPointsAsTargets, const std::vector<uint>& startPoints, const std::vector<bool>& markedTriangles = std::vector<bool>()) const;
    void createVisibilityGraph(std::vector<std::vector<GlobalID>> &edges, std::vector<std::vector<double>> &costs,
        const bool steinerPointsAsTargets, const std::vector<bool>& markedTriangles = std::vector<bool>()) const;
    TriangleIndex findTriangleContainingPoint(const Point& point, const ClosestPointSearch& closestPointSearch) const;
    bool isConstrainedPoint(const int &point) const;
    std::vector<uint> getConstrainedPoints() const;
    std::vector<uint> getMarkedConstrainedPoints(const std::vector<bool>& markedTriangles) const;
    template <class Archive> void serialize(Archive &ar, const unsigned int version);
    void calculateConstrainedAndAmbigiousPoints();
    void createThetaGraph(std::vector<std::vector<GlobalID>> &edges, std::vector<std::vector<double>> &costs, uint numIntervals);
    void createThetaGraph(std::vector<std::vector<GlobalID>> &edges, std::vector<std::vector<double>> &costs, uint numIntervals, const std::vector<bool>& shrunkNodes);
    void createThetaGraph(std::vector<std::vector<GlobalID>> &edges, std::vector<std::vector<double>> &costs, uint numIntervals, const std::vector<Point>& sourcePoints, const ClosestPointSearch& cps);

    // Read stuff
    void readFileTopologyFromGraph(const std::string& input);
    void readFromGraph(const std::string& input);
    void readTriangulationFile(const std::string& input);
    void readFromBinary(const std::string& file);

    // write stuff
    void outputAsGL(const std::string& filename, const std::vector<bool>& highlightedTriangles = std::vector<bool>(),
        bool highlightCoastlines = false);
    void writeTrianglesToFile(const std::string& filename, const std::vector<ID>& IDMapping = std::vector<ID>());
    void saveVisibilityGraph(const std::vector<std::vector<GlobalID>> &edges,
                                    const std::vector<std::vector<double>> &costs,
                                    const std::string& filename) const;
    void saveTriangulationAsGraph(const std::string& filename) const;
    void saveTriangulationAsGraph(const std::string& filename, const std::vector<bool>& markedTriangles) const;
    void saveTriangulationAsFMI(const std::string& filename, int edgeWeightScaling) const;
    void writeToBinary(const std::string& file) const;

    // Getter
    std::vector<bool> getPointsInTriangulation() const;
    Point getPoint(InternalID id) const;
    std::shared_ptr<Points> getPoints() const;
    int getNumTriangles() const;
    Triangle getTriangle(TriangleIndex index) const;
    std::vector<TriangleIndex> getAdjacentTriangles(ID point) const;
    Neighbors getNeighbors(TriangleIndex triangle) const;
    const std::vector<Triangle>& getTriangles() const;
    uint getAmountPoints() const;

    // Setter
    void setPoints(std::shared_ptr<Points> pts);
private:
    std::shared_ptr<Points> points;
    std::vector<bool> pointsInTriangulation;
    std::vector<Triangle> triangles;
    std::vector<Neighbors> neighbors;
    std::vector<std::vector<TriangleIndex>> trianglesAtNode;
    std::vector<bool> constrainedPoints;
    std::vector<bool> ambigiousPoints; // This name is taken from the polyanya code. It means that a node has more than one obstacle adjacent
    u_int amountPoints;
    void stepOneTriangle(int &oldTriangle, int &newTriangle, const Point &p1cart, const Point &p2cart) const;
    void oneToAllVisibilityInnerWorkings(const Point& startPoint, std::vector<GlobalID>& seenPoints,
        std::queue<std::tuple<std::pair<int,int>, int, int>>& segmentQueue, const bool steinerPointsAsTargets,
        std::vector<bool> markedTriangles = std::vector<bool>()) const;
    void computeSpannerNodeInnerWorkings(const Point& startPoint, std::priority_queue<std::tuple<std::pair<ID, ID>, int, int, double>, std::vector<std::tuple<std::pair<ID, ID>, int, int, double>>, Comparator>& segmentQueue, ID& spannerNode, double& distanceToClosestPoint);
    void computeSpannerNodeInnerWorkings(ID startPoint, std::priority_queue<std::tuple<std::pair<ID, ID>, int, int, double>, std::vector<std::tuple<std::pair<ID, ID>, int, int, double>>, Comparator>& segmentQueue, ID& spannerNode, double& distanceToClosestPoint, const std::vector<bool>& shrunkNodes);
};

#endif //EXACT_SHIPROUTING_TRIANGULATION_H
