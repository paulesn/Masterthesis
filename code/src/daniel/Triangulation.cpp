//
// Created by kochda on 31.01.24.
//

#include <queue>
#include "Triangulation.h"
#include "Structs.h"
#include "Timer.h"
#include <boost/archive/binary_iarchive.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "Progressbar.h"
#include <CGAL/squared_distance_2.h>

using namespace std;


Triangulation::Triangulation(){}

Triangulation::Triangulation(const std::string triangulationFile, const std::string& pointsFile) {
    Triangulation::points = readPoints(pointsFile);
    this->readTriangulationFile(triangulationFile);
    this->calculateConstrainedAndAmbigiousPoints();
}

Triangulation::Triangulation(const std::string triangulationFile, std::shared_ptr<Points> points) {
    Triangulation::points = points;
    this->readTriangulationFile(triangulationFile);
    this->calculateConstrainedAndAmbigiousPoints();
}

Triangulation::Triangulation(const std::string triangulationFile) {
    this->readTriangulationFile(triangulationFile);
}


bool Triangulation::Triangle::isCCW() {
    const Point& p0 = parentTriangulation->getPoint(this->getID(0));
    const Point& p1 = parentTriangulation->getPoint(this->getID(1));
    const Point& p2 = parentTriangulation->getPoint(this->getID(2));
    CGAL::Sign orientation = Point::orientation(p0, p1, p2);
    if (orientation == CGAL::POSITIVE)
    {
        return true;
    }
    else if (orientation == CGAL::NEGATIVE)
    {
        return false;
    }
    else if (orientation == CGAL::COPLANAR)
    {
        cout << "COPLANAR TRIANGLE" << endl;
        cout << setprecision(20);
        cout << this->getID(0) << " " << this->getID(1) << " " << this->getID(2) << " " << endl;
        cout << p0 << " " << p1 << " " << p2 << endl;
        return false;
    }
    else
    {
        cout << "This should not happen ... " << endl;
        return false;
    }
}

ID Triangulation::Triangle::getID(int i) const {
    return triangleIndices[(i%3)];
}

bool Triangulation::Triangle::lineIntersectsTriangle(const Point &startLine, const Point &endLine) const {
    // For this function we assume that the line starts AND ends outside of the triangle. Therefore only two sides
    // have to be checked for intersection

    const Point& p0 = parentTriangulation->getPoint(this->getID(0));
    const Point& p1 = parentTriangulation->getPoint(this->getID(1));
    const Point& p2 = parentTriangulation->getPoint(this->getID(2));

    if (twoLinesIntersectingEachOther(startLine, endLine, p0, p1)) {
        return true;
    }
    if (twoLinesIntersectingEachOther(startLine, endLine, p1, p2)) {
        return true;
    }
    return false;
}


void Triangulation::readFileTopologyFromGraph(const std::string& input)
{
    this->pointsInTriangulation.resize(points->size());
    std::cout << "Reading Topology of Triangulation" << std::endl;
    std::ifstream inFile(input);
    std::string line;
    int numberVertices = 0;
    int numberTriangles = 0;
    int counter = 0;
    int cwCounter = 0;
    while (getline(inFile, line))
    {
        istringstream iss(line);
        if (counter == 0)
        {
            iss >> numberVertices;
            trianglesAtNode.resize(numberVertices);
        }
        else if (counter == 1)
        {
            iss >> numberTriangles;
        }
        else if (counter > 1 && counter <= numberVertices + 1)
        {
            double lat, lon;
            iss >> lat >> lon;
        }
        else
        {
            int p1, p2, p3;
            iss >> p1 >> p2 >> p3;
            this->pointsInTriangulation[p1] = true;
            this->pointsInTriangulation[p2] = true;
            this->pointsInTriangulation[p3] = true;
            Triangle triangle(p1,p2,p3, this);
            if (!triangle.isCCW()) {
                cwCounter++;
                triangle.triangleIndices[1] = p3;
                triangle.triangleIndices[2] = p2;
            }
            this->triangles.push_back(triangle);
            trianglesAtNode[p1].push_back(this->triangles.size() - 1);
            trianglesAtNode[p2].push_back(this->triangles.size() - 1);
            trianglesAtNode[p3].push_back(this->triangles.size() - 1);
        }
        counter++;
    }
    cout << cwCounter << " triangles of " << this->triangles.size() << " were cw"<< endl;
    cout << "Read Topology from file. Now looking for neighbor triangles." << endl;
    neighbors.resize(triangles.size());
    for (int i = 0; i < this->triangles.size(); ++i)
    {
        // Indizierung ist so gemacht, dass Kante Dreieck so indiziert ist, dass es gegenüber vom Punkt mit gleichem Index liegt
        // Paar Optimierungen im Algorithmus sind noch möglich
        neighbors[i] = {-1, -1, -1};
        int p0 = this->triangles[i].getID(0), p1 = this->triangles[i].getID(1), p2 = this->triangles[i].getID(2);
        for (int j = 0; j < trianglesAtNode[p0].size(); ++j)
        { // Iterate over triangles adjacent to p0
            // In each triangle check whether there is edge p0-p1
            if (this->triangles[trianglesAtNode[p0][j]].getID(0) == p1 || this->triangles[trianglesAtNode[p0][j]].getID(1) == p1 || this->triangles[trianglesAtNode[p0][j]].getID(2) == p1)
            {
                if (i != this->trianglesAtNode[p0][j])
                {
                    neighbors[i][2] = trianglesAtNode[p0][j];
                    break;
                }
            }
        }
        for (int j = 0; j < trianglesAtNode[p1].size(); ++j)
        { // Iterate over triangles adjacent to p1
            // In each triangle check whether there is edge p1-p2
            if (this->triangles[trianglesAtNode[p1][j]].getID(0) == p2 || this->triangles[trianglesAtNode[p1][j]].getID(1) == p2 || this->triangles[trianglesAtNode[p1][j]].getID(2) == p2)
            {
                if (i != this->trianglesAtNode[p1][j])
                {
                    neighbors[i][0] = trianglesAtNode[p1][j];
                    break;
                }
            }
        }
        for (int j = 0; j < trianglesAtNode[p2].size(); ++j)
        { // Iterate over triangles adjacent to p2
            // In each triangle check whether there is edge p2-p0
            if (this->triangles[trianglesAtNode[p2][j]].getID(0) == p0 || this->triangles[trianglesAtNode[p2][j]].getID(1) == p0 || this->triangles[trianglesAtNode[p2][j]].getID(2) == p0)
            {
                if (i != this->trianglesAtNode[p2][j])
                {
                    neighbors[i][1] = trianglesAtNode[p2][j];
                    break;
                }
            }
        }
    }
    //    cout << "Sorting triangles at each point counter clockwise" << endl;
    //    this->sortTrianglesCounterClockWise();
    this->amountPoints = 0;
    for (bool point: this->pointsInTriangulation) {
        if (point)
            ++this->amountPoints;
    }
    cout << "File reading done" << endl;
}

void Triangulation::readTriangulationFile(const std::string& input)
{
    std::cout << "Reading Topology of Triangulation" << std::endl;
    this->pointsInTriangulation.resize(points->size());
    std::ifstream inFile(input);
    std::string line;
    int numberTriangles = 0;
    int counter = 0;
    int cwCounter = 0;
    trianglesAtNode.resize(points->size());
    while (getline(inFile, line))
    {
        istringstream iss(line);
        if (counter == 0)
        {
            iss >> numberTriangles;
        }
        else
        {
            int p1, p2, p3;
            iss >> p1 >> p2 >> p3;
            this->pointsInTriangulation[p1] = true;
            this->pointsInTriangulation[p2] = true;
            this->pointsInTriangulation[p3] = true;
            Triangle triangle(p1,p2,p3, this);
            if (!triangle.isCCW()) {
                cwCounter++;
                triangle.triangleIndices[1] = p3;
                triangle.triangleIndices[2] = p2;
                ;
            }
            this->triangles.push_back(triangle);
            trianglesAtNode[p1].push_back(this->triangles.size() - 1);
            trianglesAtNode[p2].push_back(this->triangles.size() - 1);
            trianglesAtNode[p3].push_back(this->triangles.size() - 1);
        }
        counter++;
    }
    cout << cwCounter << " triangles of " << this->triangles.size() << " were cw"<< endl;
    cout << "Read Topology from file. Now looking for neighbor triangles." << endl;
    neighbors.resize(triangles.size());
    for (int i = 0; i < this->triangles.size(); ++i)
    {
        // Indizierung ist so gemacht, dass Kante Dreieck so indiziert ist, dass es gegenüber vom Punkt mit gleichem Index liegt
        // Paar Optimierungen im Algorithmus sind noch möglich
        neighbors[i] = {-1, -1, -1};
        int p0 = this->triangles[i].getID(0), p1 = this->triangles[i].getID(1), p2 = this->triangles[i].getID(2);
        for (int j = 0; j < trianglesAtNode[p0].size(); ++j)
        { // Iterate over triangles adjacent to p0
            // In each triangle check whether there is edge p0-p1
            if (this->triangles[trianglesAtNode[p0][j]].getID(0) == p1 || this->triangles[trianglesAtNode[p0][j]].getID(1) == p1 || this->triangles[trianglesAtNode[p0][j]].getID(2) == p1)
            {
                if (i != this->trianglesAtNode[p0][j])
                {
                    neighbors[i][2] = trianglesAtNode[p0][j];
                    break;
                }
            }
        }
        for (int j = 0; j < trianglesAtNode[p1].size(); ++j)
        { // Iterate over triangles adjacent to p1
            // In each triangle check whether there is edge p1-p2
            if (this->triangles[trianglesAtNode[p1][j]].getID(0) == p2 || this->triangles[trianglesAtNode[p1][j]].getID(1) == p2 || this->triangles[trianglesAtNode[p1][j]].getID(2) == p2)
            {
                if (i != this->trianglesAtNode[p1][j])
                {
                    neighbors[i][0] = trianglesAtNode[p1][j];
                    break;
                }
            }
        }
        for (int j = 0; j < trianglesAtNode[p2].size(); ++j)
        { // Iterate over triangles adjacent to p2
            // In each triangle check whether there is edge p2-p0
            if (this->triangles[trianglesAtNode[p2][j]].getID(0) == p0 || this->triangles[trianglesAtNode[p2][j]].getID(1) == p0 || this->triangles[trianglesAtNode[p2][j]].getID(2) == p0)
            {
                if (i != this->trianglesAtNode[p2][j])
                {
                    neighbors[i][1] = trianglesAtNode[p2][j];
                    break;
                }
            }
        }
    }
    //    cout << "Sorting triangles at each point counter clockwise" << endl;
    //    this->sortTrianglesCounterClockWise();
    this->amountPoints = 0;
    for (bool point: this->pointsInTriangulation) {
        if (point)
            ++this->amountPoints;
    }
    cout << "File reading done" << endl;
}

int Triangulation::Triangle::getIndex(ID pointID) const {
    if (this->triangleIndices[0] == pointID) {
        return 0;
    } else if (this->triangleIndices[1] == pointID) {
        return 1;
    } else if (this->triangleIndices[2] == pointID) {
        return 2;
    } else {
        return -1;
    }
}

int Triangulation::Triangle::CCW(int i) {
    return (i+1) % 3;
}

int Triangulation::Triangle::CW(int i) {
    return (i+2) % 3;
}

bool Triangulation::directVisibility(InternalID p1, InternalID p2) const {
    // DEBUG
    //    cout << "Points: " << p1 << " " << p2 << endl;
    //    this->visitedTriangles.clear();
    //    this->visitedTriangles.resize(triangles.size(), false);
    // DEBUG END
    if (p1 == p2)
    {
        // Convention that a point does not see itself to avoid self loops in resulting graphs
        return false;
    }
    Point p1Point = (*this->points)[p1];
    Point p2Point = (*this->points)[p2];
    // Find start triangle.
    int forwardTriangle = -1;
    int forwardTriangleNext = -1;
    // Coplanar edge case is not yet implemented
    for (TriangleIndex triangleIndex : this->trianglesAtNode[p1])
    {
        // find indexPoint of p1 in the triangle
        int p1Index = this->triangles[triangleIndex].getIndex(p1);
        // First check whether p1 and p2 are direct neighbors
        if (triangles.at(triangleIndex).getID(Triangulation::Triangle::CCW(p1Index)) == p2 ||
            triangles.at(triangleIndex).getID(Triangulation::Triangle::CW(p1Index)) == p2)
        {
            return true;
        }
        // Convert the needed triangle points into cartesian coordinates

        Point ccwp1cart = (*this->points)[triangles[triangleIndex].getID(Triangle::CCW(p1Index))];
        Point cwp1cart = (*this->points)[triangles[triangleIndex].getID(Triangle::CW(p1Index))];

        // orientation tests whether we are in the start triangle
        //        cout << orientation(p1Point, p2Point, ccwp1cart) << " " << orientation(p1Point, p2Point, cwp1cart) << endl;
        if (Point::orientation(p1Point, p2Point, ccwp1cart) == CGAL::NEGATIVE && // NEGATIVE means its on the side of the plain where the points are counterclockwise
            Point::orientation(p1Point, p2Point, cwp1cart) == CGAL::POSITIVE)
        {
            forwardTriangle = triangleIndex;
            //            this->visitedTriangles.at(forwardTriangle) = true;
            // As for the moment we ignore the coplanar case, we get the second triangle for free
            forwardTriangleNext = this->neighbors[forwardTriangle][p1Index];
        }
    }
    // Now traverse through the triangles until we hit land or the target
    while (forwardTriangleNext != -1)
    {
        // Erst überprüfe, ob Zielpunkt im Dreieck liegt
        // TODO: Überlegen, ob ich nur gegenüberliegenden checken muss?
        int targetIndex = triangles[forwardTriangleNext].getIndex(p2);
        if (targetIndex != -1) {
            return true;
        }

        stepOneTriangle(forwardTriangle, forwardTriangleNext, p1Point, p2Point);
        //        this->visitedTriangles.at(forwardTriangle) = true;
    }
    //    cout << endl;
    return false;
}

bool Triangulation::directVisibility(const Point& p1, const Point& p2, const ClosestPointSearch& cps) const {
    if (p1 == p2) {
        return false;
    }
    Point p1Copy = p1, p2Copy = p2;
    ID closest1 = cps.getClosestPoint(p1Copy);
    ID closest2 = cps.getClosestPoint(p2Copy);
    bool startInNode = false;
    // If source and target are on the triangulation we use the old directVisibility function for that case
    if (this->getPoint(closest1) == p1Copy && this->getPoint(closest2) == p2Copy) {
        return this->directVisibility(closest1, closest2);
    }
    if (this->getPoint(closest1) == p1Copy) {
        startInNode = true;
    }
    if (!(this->getPoint(closest1) == p1Copy) && this->getPoint(closest2) == p2Copy) {
        // We only want to implement the special case where p1 is on the triangulation and p2 is not. We can achieve that by switching p1 and p2, should it be the other way around.
        p1Copy = p2;
        p2Copy = p1;
        ID tempId = closest1;
        closest1 = closest2;
        closest2 = tempId;
        startInNode = true;
    }
    // First find the triangles containing p1 and p2
    TriangleIndex t1Index = this->findTriangleContainingPoint(p1Copy, cps);
    TriangleIndex t2Index = this->findTriangleContainingPoint(p2Copy, cps);
    if (t1Index == -1 || t2Index == -1) {
        return false;
    }
    if (t1Index == t2Index) {
        return true;
    }
    // We traverse the triangulation until we hit an obstacle or the triangle we step into is t2
    // We need to initialize
    TriangleIndex forwardTriangle = t1Index;
    TriangleIndex forwardTriangleNext = -1;
    // We need to distinguish the case where the start is on a node and the case, where the start is in a triangle
    if (startInNode) {
        for (TriangleIndex triangleIndex : this->trianglesAtNode[closest1])
        {
            if (triangleIndex == t2Index) {
                return true;
            }
            // find indexPoint of p1 in the triangle
            int p1Index = this->triangles[triangleIndex].getIndex(closest1);

            // Convert the needed triangle points into cartesian coordinates

            Point ccwp1cart = (*this->points)[triangles[triangleIndex].getID(Triangle::CCW(p1Index))];
            Point cwp1cart = (*this->points)[triangles[triangleIndex].getID(Triangle::CW(p1Index))];

            // orientation tests whether we are in the start triangle
            //        cout << orientation(p1Point, p2Point, ccwp1cart) << " " << orientation(p1Point, p2Point, cwp1cart) << endl;
            if (Point::orientation(p1Copy, p2Copy, ccwp1cart) == CGAL::NEGATIVE && // NEGATIVE means its on the side of the plain where the points are counterclockwise
                Point::orientation(p1Copy, p2Copy, cwp1cart) == CGAL::POSITIVE)
            {
                forwardTriangle = triangleIndex;
                //            this->visitedTriangles.at(forwardTriangle) = true;
                // As for the moment we ignore the coplanar case, we get the second triangle for free
                forwardTriangleNext = this->neighbors[forwardTriangle][p1Index];
            }
        }
    } else {
        Triangulation::Triangle t1 = this->getTriangle(t1Index);
        for (int i = 0; i < 3; ++i) {
            const Point& trianglePoint1 = this->getPoint(t1.getID(i));
            const Point& trianglePoint2 = this->getPoint(t1.getID(i + 1));
            CGAL::Orientation orientation1 = Point::orientation(p1Copy, p2Copy, trianglePoint1);
            CGAL::Orientation orientation2 = Point::orientation(p1Copy, p2Copy, trianglePoint2);
            if (orientation1 == CGAL::NEGATIVE && orientation2 == CGAL::POSITIVE) {
                forwardTriangleNext = this->neighbors[t1Index][(i+2) % 3];
                break;
            }
            // I hope this is correct.
            // if (orientation1 == CGAL::ZERO || orientation2 == CGAL::ZERO) {
            //     return true;
            // }
        }
    }

    while (forwardTriangleNext != -1) {
        if (forwardTriangleNext == t2Index) {
            return true;
        }
        stepOneTriangle(forwardTriangle, forwardTriangleNext, p1Copy, p2Copy);
    }
    return false;
}


void Triangulation::stepOneTriangle(int &oldTriangle, int &newTriangle, const Point &p1cart,
                                    const Point &p2cart) const {
    // Find the index of the old triangle in the new one, as it is the same indexPoint as the point we have to check
    auto it = find(this->neighbors[newTriangle].begin(), this->neighbors[newTriangle].end(), oldTriangle);
    int pointIndex = it - this->neighbors[newTriangle].begin();
    Point checkPoint = (*this->points)[(this->triangles[newTriangle]).getID(pointIndex)];
    int outgoingTriangleIndex;
    // Again, Coplanar is not looked at.
    CGAL::Sign cpOrientation = Point::orientation(p1cart, p2cart, checkPoint);
    if (cpOrientation == CGAL::POSITIVE)
    { // positive is "left"
        outgoingTriangleIndex = Triangulation::Triangle::CCW(pointIndex);
    }
    else if (cpOrientation == CGAL::NEGATIVE)
    { // Negative is "right"
        outgoingTriangleIndex = Triangulation::Triangle::CW(pointIndex);
    }
    else if (cpOrientation == CGAL::ZERO)
    {
        //        cout << "ZERO CASE" << endl;
        int point = this->triangles[newTriangle].getID(pointIndex);
        if (this->isConstrainedPoint(point))
        {
            oldTriangle = newTriangle;
            newTriangle = -1;
            return;
        }
        for (int triangleIndex : this->trianglesAtNode[point])
        {
            // find indexPoint of p1 in the triangle
            int p1Index = this->triangles[triangleIndex].getIndex(point);
            // Convert the needed triangle points into cartesian coordinates
            Point ccwp1cart = (*points)[triangles[triangleIndex].getID(Triangle::CCW(p1Index))];
            Point cwp1cart = (*points)[triangles[triangleIndex].getID(Triangle::CW(p1Index))];

            // orientation tests whether we are in the start triangle
            //        cout << orientation(p1Point, p2Point, ccwp1cart) << " " << orientation(p1Point, p2Point, cwp1cart) << endl;
            if (Point::orientation(p1cart, p2cart, ccwp1cart) == CGAL::NEGATIVE && // NEGATIVE means its on the side of the plain where the points are counterclockwise
                Point::orientation(p1cart, p2cart, cwp1cart) == CGAL::POSITIVE)
            {
                oldTriangle = triangleIndex;
                // As for the moment we ignore the coplanar case, we get the second triangle for free
                newTriangle = this->neighbors[oldTriangle][p1Index];
                return;
            }
        }
        // If no outgoing triangle is found, we return newTriangle -1, such that the algorithm terminates with false
        oldTriangle = newTriangle;
        newTriangle = -1;
        return;
    }
    oldTriangle = newTriangle;
    newTriangle = this->neighbors[newTriangle][outgoingTriangleIndex];
}

bool Triangulation::isConstrainedPoint(const int& point) const {

    // for (const int &adjacentTriangle : this->trianglesAtNode[point])
    // {
    //     for (int i = 0; i < this->neighbors[adjacentTriangle].size(); ++i)
    //     {
    //         if (this->neighbors[adjacentTriangle][i] == -1)
    //         {
    //             if (this->triangles[adjacentTriangle].getID(Triangulation::Triangle::CCW(i)) == point ||
    //                     this->triangles[adjacentTriangle].getID(Triangulation::Triangle::CW(i)) == point) {
    //                 return true;
    //             }
    //         }
    //     }
    // }
    // return false;
    return this->constrainedPoints[point];
}

void Triangulation::calculateConstrainedAndAmbigiousPoints() {
    this->constrainedPoints.clear();
    this->constrainedPoints.resize(points->size());
    this->ambigiousPoints.clear();
    this->ambigiousPoints.resize(points->size(), false);
    cout << "Calculating ambigious and constrained points" << endl;
    for (int point = 0; point < points->size(); ++point) {
        int obstacleCount = 0;
        for (const int &adjacentTriangle : this->trianglesAtNode[point])
        {
            for (int i = 0; i < this->neighbors[adjacentTriangle].size(); ++i)
            {
                if (this->neighbors[adjacentTriangle][i] == -1)
                {
                    if (this->triangles[adjacentTriangle].getID(Triangulation::Triangle::CCW(i)) == point ||
                            this->triangles[adjacentTriangle].getID(Triangulation::Triangle::CW(i)) == point) {
                        obstacleCount++;
                            }
                }
            }
        }
        if (obstacleCount != 0) {
            this->constrainedPoints[point] = true;
            if (obstacleCount > 2) {
                this->ambigiousPoints[point] = true;
            }
        }
    }
}


void Triangulation::oneToAllVisibilityInnerWorkings(const Point& startPoint, std::vector<GlobalID>& seenPoints,
        std::queue<std::tuple<std::pair<int,int>, int, int>>& segmentQueue, const bool steinerPointsAsTargets,
        std::vector<bool> markedTriangles) const {
    typedef std::pair<int, int> Segment;
    while (!segmentQueue.empty())
    {
//                cout << segmentQueue.size() << endl;
        tuple<Segment, int, int> currentStuff = segmentQueue.front();
        segmentQueue.pop();
        Segment currentSegment = get<0>(currentStuff); //.first is the ccw end, .second the cw one
        int originTriangle = get<1>(currentStuff);
        int edgeID = get<2>(currentStuff);
//        this->visitedTriangles[originTriangle] = true;
//        cout << originTriangle << endl;
        if (!markedTriangles.empty() && !markedTriangles[originTriangle]) {
            continue;
        }
        Point segment2Cart = (*points)[currentSegment.second];
        Point segment1Cart = (*points)[currentSegment.first];

        if (this->neighbors[originTriangle][edgeID] == -1)
        {
            auto orientation1 = Point::orientation(startPoint, segment2Cart, (*points)[triangles[originTriangle].getID((edgeID + 2) % 3)]);
            auto orientation2 = Point::orientation(startPoint, segment1Cart, (*points)[triangles[originTriangle].getID((edgeID + 1) % 3)]);
            if (orientation1 == CGAL::ZERO || orientation1 == CGAL::COPLANAR)
            {
                seenPoints.push_back(currentSegment.second);
            }
            if (orientation2 == CGAL::ZERO || orientation2 == CGAL::COPLANAR)
            {
                seenPoints.push_back(currentSegment.first);
            }
            continue;
        }

        int neighborTriangle = this->neighbors[originTriangle][edgeID];
        int entryEdge = this->triangles[neighborTriangle].indexEdge(this->triangles[originTriangle].getID((edgeID + 1) % 3),
                                                                    this->triangles[originTriangle].getID((edgeID + 2) % 3));

        Point checkPointCart = (*points)[triangles[this->neighbors[originTriangle][edgeID]].getID(entryEdge)];

        auto orientationS1 = Point::orientation(startPoint, segment1Cart, checkPointCart);
        auto orientationS2 = Point::orientation(startPoint, segment2Cart, checkPointCart);

        if (orientationS1 == CGAL::POSITIVE && orientationS2 == CGAL::NEGATIVE)
        {
            Segment segment1 = make_pair(currentSegment.first, this->triangles[neighborTriangle].getID(entryEdge));
            segmentQueue.push(make_tuple(segment1, neighborTriangle, (entryEdge + 1) % 3));
            Segment segment2 = make_pair(this->triangles[neighborTriangle].getID(entryEdge), currentSegment.second);
            segmentQueue.push(make_tuple(segment2, neighborTriangle, (entryEdge + 2) % 3));
            if (steinerPointsAsTargets){
                seenPoints.push_back(triangles[this->neighbors[originTriangle][edgeID]].getID(entryEdge));
            }
        }
        else if (orientationS1 == CGAL::POSITIVE && (orientationS2 == CGAL::POSITIVE || orientationS2 == CGAL::COPLANAR))
        {
            segmentQueue.push(make_tuple(currentSegment, neighborTriangle, (entryEdge + 1) % 3));
        }
        else if ((orientationS1 == CGAL::NEGATIVE || orientationS1 == CGAL::COPLANAR) && orientationS2 == CGAL::NEGATIVE)
        {
            segmentQueue.push(make_tuple(currentSegment, neighborTriangle, (entryEdge + 2) % 3));
        }
        else
        {
            cout << "This should not happen ...." << endl;
            cout << orientationS1 << " " << orientationS2 << endl;
        }
    }
    sort(seenPoints.begin(), seenPoints.end());
    seenPoints.erase(unique(seenPoints.begin(), seenPoints.end()), seenPoints.end());
}

std::vector<GlobalID> Triangulation::oneToAllVisibility(const InternalID startPoint,
    const bool steinerPointsAsTargets, const std::vector<bool>& markedTriangles) const {
    typedef pair<int, int> Segment;
    vector<GlobalID> seenPoints;
    Point startPointCart = (*this->points)[startPoint];
    std::queue<tuple<Segment, int, int>> segmentQueue;
    for (const int &neighborTriangle : this->trianglesAtNode[startPoint])
    {
        if (!markedTriangles.empty() && !markedTriangles[neighborTriangle]) {
            continue;
        }
        int startIndex = this->triangles[neighborTriangle].getIndex(startPoint);
        Segment segment = make_pair(this->triangles[neighborTriangle].getID(Triangle::CCW(startIndex)),
                                    this->triangles[neighborTriangle].getID(Triangle::CW(startIndex)));
        segmentQueue.push(make_tuple(segment, neighborTriangle, startIndex));
        // Now check if there is a constrained area adjacent to the triangle and if yes add the corresponding point
        Neighbors adjacentTriangles = this->neighbors[neighborTriangle];
        if (adjacentTriangles[(startIndex + 1) % 3] == -1)
        {
            seenPoints.push_back(triangles[neighborTriangle].getID((startIndex + 2) % 3));
        }
        if (adjacentTriangles[(startIndex + 2) % 3] == -1)
        {
            seenPoints.push_back(triangles[neighborTriangle].getID((startIndex + 1) % 3));
        }
        if (steinerPointsAsTargets) {
            seenPoints.push_back(triangles[neighborTriangle].getID((startIndex + 2) % 3));
            seenPoints.push_back(triangles[neighborTriangle].getID((startIndex + 1) % 3));
        }
    }
    this->oneToAllVisibilityInnerWorkings(startPointCart, seenPoints, segmentQueue, steinerPointsAsTargets, markedTriangles);
    return seenPoints;
}

std::vector<GlobalID> Triangulation::oneToAllVisibility(const Point &startPoint,
    const ClosestPointSearch& closestPointSearch, const bool steinerPointsAsTargets) const {
    ID closestPoint = closestPointSearch.getClosestPoint(startPoint);
    if ((*this->points)[closestPoint] == startPoint) {
        return this->oneToAllVisibility(closestPoint, steinerPointsAsTargets);
    }
    //TODO Hier für nicht doppelte ClosestPoint Suche sorgen (Doppelt in findTriangleContainingPoint)
    TriangleIndex triangleContainingStartPointIndex = this->findTriangleContainingPoint(startPoint, closestPointSearch);
    const Triangle& triangleContainingStartPoint = this->triangles[triangleContainingStartPointIndex];
    vector<GlobalID> seenPoints;
    if (triangleContainingStartPointIndex == -1) {
        return seenPoints;
    }
    typedef pair<int, int> Segment;
    std::queue<tuple<Segment, int, int>> segmentQueue;
    for (int i = 0; i < 3; ++i) {
        Segment segment = std::make_pair(triangleContainingStartPoint.getID(i+1), triangleContainingStartPoint.getID(i+2));
        segmentQueue.push(make_tuple(segment, triangleContainingStartPointIndex, i));
    }
    this->oneToAllVisibilityInnerWorkings(startPoint, seenPoints, segmentQueue, steinerPointsAsTargets);
    return seenPoints;
}

ID Triangulation::computeSpannerNode(const Point& startPoint, const Point &leftEndCone, const Point &rightEndCone, const ClosestPointSearch& cps) {
    ID closestPointToStart = cps.getClosestPoint(startPoint);
    if ((*this->points)[closestPointToStart] == startPoint) {
        return this->computeSpannerNode(closestPointToStart, leftEndCone, rightEndCone);
    }
    points->push_back(leftEndCone); // points->size() - 2 Remember to delete
    points->push_back(rightEndCone); // points->size() - 1
    ID leftEndConeID = points->size() - 2, rightEndConeID = points->size() - 1;
    typedef pair<ID, ID> Segment; // A tuple only determines the direction of the cone. For distance calculations one needs the actual segment
    ID closestPoint = -1;
    double distanceToClosestPoint = maxWeight;
    std::priority_queue<tuple<Segment, int, int, double>, std::vector<tuple<Segment, int, int, double>>, Comparator> segmentQueue;
    TriangleIndex triangleContainingStartPointIndex = this->findTriangleContainingPoint(startPoint, cps);
    if (triangleContainingStartPointIndex == -1) {
        return -1;
    }
    const Triangle& triangleContainingStartPoint = this->triangles[triangleContainingStartPointIndex];
    for (int i = 0; i < 3; ++i) {
        std::pair<ID, ID> currentEdge = std::make_pair(triangleContainingStartPoint.getID(i+1), triangleContainingStartPoint.getID(i+2));
        Segment_2 cgalEdge(this->getPoint(currentEdge.first).getCGALPoint(), this->getPoint(currentEdge.second).getCGALPoint());
        double distance = squared_distance(startPoint.getCGALPoint(), cgalEdge);
        CGAL::Orientation leftConeWithLeftEnd = Point::orientation(startPoint, this->getPoint(currentEdge.second), leftEndCone);
        CGAL::Orientation leftConeWithRightEnd = Point::orientation(startPoint, this->getPoint(currentEdge.first), leftEndCone);
        CGAL::Orientation rightConeWithLeftEnd = Point::orientation(startPoint, this->getPoint(currentEdge.second), rightEndCone);
        CGAL::Orientation rightConeWithRightEnd = Point::orientation(startPoint, this->getPoint(currentEdge.first), rightEndCone);
        bool leftConeIsInsideEdge = (leftConeWithLeftEnd != CGAL::LEFT_TURN && leftConeWithRightEnd != CGAL::RIGHT_TURN);
        bool rightConeIsInsideEdge = (rightConeWithLeftEnd != CGAL::LEFT_TURN && rightConeWithRightEnd != CGAL::RIGHT_TURN);
        if (leftConeIsInsideEdge && rightConeIsInsideEdge) {
            Segment segment = std::make_pair(rightEndConeID, leftEndConeID);
            segmentQueue.push(make_tuple(segment, triangleContainingStartPointIndex, i, distance));
            break;
        } else if (leftConeIsInsideEdge && !rightConeIsInsideEdge) {
            Segment segment = std::make_pair(currentEdge.first, leftEndConeID);
            segmentQueue.push(make_tuple(segment, triangleContainingStartPointIndex, i, distance));
            double dist = Point::squaredDistance(startPoint, this->getPoint(currentEdge.first));
            if (closestPoint == -1) {
                closestPoint = currentEdge.first;
                distanceToClosestPoint = dist;
            } else if (dist < distanceToClosestPoint) {
                closestPoint = currentEdge.first;
                distanceToClosestPoint = dist;
            }
        } else if (rightConeIsInsideEdge && !leftConeIsInsideEdge) {
            Segment segment = std::make_pair(rightEndConeID, currentEdge.second);
            segmentQueue.push(make_tuple(segment, triangleContainingStartPointIndex, i, distance));
            double dist = Point::squaredDistance(startPoint, this->getPoint(currentEdge.first));
            if (closestPoint == -1) {
                closestPoint = currentEdge.second;
                distanceToClosestPoint = dist;
            } else if (dist < distanceToClosestPoint) {
                closestPoint = currentEdge.second;
                distanceToClosestPoint = dist;
            }
        }
    }
    this->computeSpannerNodeInnerWorkings(startPoint, segmentQueue, closestPoint, distanceToClosestPoint);
    points->pop_back();
    points->pop_back();
    return closestPoint;
}

ID Triangulation::computeSpannerNode(ID startPoint, const Point &leftEndCone, const Point &rightEndCone) {
    points->push_back(leftEndCone); // points->size() - 2 Remember to delete
    points->push_back(rightEndCone); // points->size() - 1
    ID leftEndConeID = points->size() - 2, rightEndConeID = points->size() - 1;
    typedef pair<ID, ID> Segment; // A tuple only determines the direction of the cone. For distance calculations one needs the actual segment
    ID closestPoint = -1;
    double distanceToClosestPoint = maxWeight;
    Point startPointCart = (*this->points)[startPoint];

    std::priority_queue<tuple<Segment, int, int, double>, std::vector<tuple<Segment, int, int, double>>, Comparator> segmentQueue;
    for (const int &neighborTriangle : this->trianglesAtNode[startPoint])
    {
        int startIndex = this->triangles[neighborTriangle].getIndex(startPoint);
        CGAL::Orientation leftConeWithLeftOrientation = Point::orientation((*points)[startPoint], (*points)[this->triangles[neighborTriangle].getID(Triangle::CW(startIndex))], leftEndCone);
        CGAL::Orientation leftConeWithRightOrientation = Point::orientation((*points)[startPoint], (*points)[this->triangles[neighborTriangle].getID(Triangle::CCW(startIndex))], leftEndCone);
        CGAL::Orientation rightConeWithRightOrientation = Point::orientation((*points)[startPoint], (*points)[this->triangles[neighborTriangle].getID(Triangle::CCW(startIndex))], rightEndCone);
        CGAL::Orientation rightConeWithLeftOrientation = Point::orientation((*points)[startPoint], (*points)[this->triangles[neighborTriangle].getID(Triangle::CW(startIndex))], rightEndCone);
        bool leftEndConeInsideTriangle = leftConeWithLeftOrientation == CGAL::NEGATIVE && leftConeWithRightOrientation == CGAL::POSITIVE;
        bool rightEndConeInsideTriangle = rightConeWithLeftOrientation == CGAL::NEGATIVE && rightConeWithRightOrientation == CGAL::POSITIVE;
        if (leftEndConeInsideTriangle && !rightEndConeInsideTriangle) {
            Segment segment = make_pair(this->triangles[neighborTriangle].getID(Triangle::CCW(startIndex)), leftEndConeID);
            // We do not compute the real intersection, but just use the edge of the triangle to calculate the distance. This gives a distance, which is too small, but that's okay
            Segment_2 cgalSegment((*this->points)[this->triangles[neighborTriangle].getID(Triangle::CW(startIndex))].getCGALPoint(),
                (*points)[ this->triangles[neighborTriangle].getID(Triangle::CCW(startIndex))].getCGALPoint());
            double distance = squared_distance(cgalSegment, startPointCart.getCGALPoint());
            segmentQueue.emplace(make_tuple(segment, neighborTriangle, startIndex, distance));
            double dist = Point::squaredDistance(startPointCart, (*points)[this->triangles[neighborTriangle].getID(Triangle::CCW(startIndex))]);
            if (closestPoint == -1) {
                closestPoint = this->triangles[neighborTriangle].getID(Triangle::CCW(startIndex));
                distanceToClosestPoint = dist;
            } else {
                if (dist < distanceToClosestPoint) {
                    closestPoint = this->triangles[neighborTriangle].getID(Triangle::CCW(startIndex));
                    distanceToClosestPoint = dist;
                }
            }

        } else if (!leftEndConeInsideTriangle && rightEndConeInsideTriangle) {
            Segment segment = make_pair(rightEndConeID, this->triangles[neighborTriangle].getID(Triangle::CW(startIndex)));

            Segment_2 cgalSegment((*this->points)[this->triangles[neighborTriangle].getID(Triangle::CW(startIndex))].getCGALPoint(),
                (*points)[ this->triangles[neighborTriangle].getID(Triangle::CCW(startIndex))].getCGALPoint());
            double distance = squared_distance(cgalSegment, startPointCart.getCGALPoint());
            segmentQueue.emplace(make_tuple(segment, neighborTriangle, startIndex, distance));
            double dist = Point::squaredDistance(startPointCart, (*points)[this->triangles[neighborTriangle].getID(Triangle::CW(startIndex))]);
            if (closestPoint == -1) {
                closestPoint = this->triangles[neighborTriangle].getID(Triangle::CW(startIndex));
                distanceToClosestPoint = dist;
            } else {
                if (dist < distanceToClosestPoint) {
                    distanceToClosestPoint = dist;
                    closestPoint = this->triangles[neighborTriangle].getID(Triangle::CW(startIndex));
                }
            }
        } else if (leftEndConeInsideTriangle && rightEndConeInsideTriangle) {
            Segment segment = make_pair(rightEndConeID, leftEndConeID);
            Segment_2 cgalSegment((*this->points)[this->triangles[neighborTriangle].getID(Triangle::CW(startIndex))].getCGALPoint(),
                (*points)[ this->triangles[neighborTriangle].getID(Triangle::CCW(startIndex))].getCGALPoint());
            double distance = squared_distance(cgalSegment, startPointCart.getCGALPoint());
            segmentQueue.emplace(make_tuple(segment, neighborTriangle, startIndex, distance));
        } else if (leftConeWithLeftOrientation == CGAL::LEFT_TURN && rightConeWithRightOrientation == CGAL::RIGHT_TURN) {
            // This check is for the case where the triangle is completely contained in the cone.
            // So far the test is also positive when the triangle is on the opposite side and really wide.
            // So check that the triangle is on the right side.
            if (rightConeWithLeftOrientation == CGAL::RIGHT_TURN && leftConeWithRightOrientation == CGAL::LEFT_TURN) {
                Segment segment = make_pair(this->triangles[neighborTriangle].getID(Triangle::CCW(startIndex)), this->triangles[neighborTriangle].getID(Triangle::CW(startIndex)));
                Segment_2 cgalSegment((*this->points)[this->triangles[neighborTriangle].getID(Triangle::CW(startIndex))].getCGALPoint(),
                    (*points)[ this->triangles[neighborTriangle].getID(Triangle::CCW(startIndex))].getCGALPoint());
                double distance = squared_distance(cgalSegment, startPointCart.getCGALPoint());
                segmentQueue.emplace(make_tuple(segment, neighborTriangle, startIndex, distance));
                double dist = Point::squaredDistance(startPointCart, (*points)[this->triangles[neighborTriangle].getID(Triangle::CCW(startIndex))]);
                if (closestPoint == -1) {
                    closestPoint = this->triangles[neighborTriangle].getID(Triangle::CCW(startIndex));
                    distanceToClosestPoint = dist;
                } else {
                    if (dist < distanceToClosestPoint) {
                        distanceToClosestPoint = dist;
                        closestPoint = this->triangles[neighborTriangle].getID(Triangle::CCW(startIndex));
                    }
                }
                dist = Point::squaredDistance(startPointCart, (*points)[this->triangles[neighborTriangle].getID(Triangle::CW(startIndex))]);
                if (dist < distanceToClosestPoint) {
                    distanceToClosestPoint = dist;
                    closestPoint = this->triangles[neighborTriangle].getID(Triangle::CW(startIndex));
                }
            }

        }
    }
    this->computeSpannerNodeInnerWorkings(this->getPoint(startPoint), segmentQueue, closestPoint, distanceToClosestPoint);
    points->pop_back();
    points->pop_back();
    return closestPoint;
}

void Triangulation::computeSpannerNodeInnerWorkings(const Point& startPoint,
    std::priority_queue<std::tuple<std::pair<ID, ID>, int, int, double>, std::vector<std::tuple<std::pair<ID, ID>, int, int, double> >, Comparator> &segmentQueue,
    ID &spannerNode, double &distanceToClosestPoint) {
    typedef std::pair<ID, ID> Segment;
    while (!segmentQueue.empty())
    {
//                cout << segmentQueue.size() << endl;
        tuple<Segment, int, int, double> currentStuff = segmentQueue.top();
        segmentQueue.pop();
        Segment currentSegment = get<0>(currentStuff); //.first is the ccw end, .second the cw one
        int originTriangle = get<1>(currentStuff);
        int edgeID = get<2>(currentStuff);
        double currentDistance = get<3>(currentStuff);
        if (currentDistance > distanceToClosestPoint) {
            return;
        }
//        this->visitedTriangles[originTriangle] = true;
//        cout << originTriangle << endl;
        Point segment2Cart = (*points)[currentSegment.second];
        Point segment1Cart = (*points)[currentSegment.first];

        int neighborTriangle = this->neighbors[originTriangle][edgeID];
        if (neighborTriangle == -1)
        {
            auto orientation1 = Point::orientation(startPoint, segment2Cart, (*points)[triangles[originTriangle].getID((edgeID + 2) % 3)]);
            auto orientation2 = Point::orientation(startPoint, segment1Cart, (*points)[triangles[originTriangle].getID((edgeID + 1) % 3)]);
            if (orientation1 == CGAL::ZERO || orientation1 == CGAL::COPLANAR)
            {
                double distance = Point::squaredDistance(startPoint, (*points)[triangles[originTriangle].getID((edgeID + 2) % 3)]);
                if (distance < distanceToClosestPoint) {
                    distanceToClosestPoint = distance;
                    spannerNode = triangles[originTriangle].getID((edgeID + 2) % 3);
                }
            }
            if (orientation2 == CGAL::ZERO || orientation2 == CGAL::COPLANAR)
            {
                double distance = Point::squaredDistance(startPoint, (*points)[triangles[originTriangle].getID((edgeID + 1) % 3)]);
                if (distance < distanceToClosestPoint) {
                    distanceToClosestPoint = distance;
                    spannerNode = triangles[originTriangle].getID((edgeID + 1) % 3);
                }
            }
            continue;
        }

        int entryEdge = this->triangles[neighborTriangle].indexEdge(this->triangles[originTriangle].getID((edgeID + 1) % 3),
                                                                    this->triangles[originTriangle].getID((edgeID + 2) % 3));

        ID checkPointID = triangles[neighborTriangle].getID(entryEdge);
        Point checkPointCart = (*points)[triangles[neighborTriangle].getID(entryEdge)];

        auto orientationS1 = Point::orientation(startPoint, segment1Cart, checkPointCart);
        auto orientationS2 = Point::orientation(startPoint, segment2Cart, checkPointCart);


        if (orientationS1 == CGAL::POSITIVE && orientationS2 == CGAL::NEGATIVE)
        {
            Segment segment1 = make_pair(currentSegment.first, checkPointID);
            Segment_2 cgalSegment1((*points)[triangles[neighborTriangle].getID(Triangle::CW(entryEdge))].getCGALPoint(), (*points)[triangles[neighborTriangle].getID(entryEdge)].getCGALPoint());
            double distance = squared_distance(cgalSegment1, startPoint.getCGALPoint());
            segmentQueue.emplace(make_tuple(segment1, neighborTriangle, (entryEdge + 1) % 3, distance));

            Segment segment2 = make_pair(checkPointID, currentSegment.second);
            Segment_2 cgalSegment2((*points)[triangles[neighborTriangle].getID(entryEdge)].getCGALPoint(), (*points)[triangles[neighborTriangle].getID(Triangle::CCW(entryEdge))].getCGALPoint());
            double distance2 = squared_distance(cgalSegment2, startPoint.getCGALPoint());
            segmentQueue.push(make_tuple(segment2, neighborTriangle, (entryEdge + 2) % 3, distance2));
        }
        else if (orientationS1 == CGAL::POSITIVE && (orientationS2 == CGAL::POSITIVE || orientationS2 == CGAL::COPLANAR))
        {
            Segment_2 cgalSegment((*points)[triangles[neighborTriangle].getID(Triangle::CW(entryEdge))].getCGALPoint(), checkPointCart.getCGALPoint());
            double distance = squared_distance(cgalSegment, startPoint.getCGALPoint());
            segmentQueue.push(make_tuple(currentSegment, neighborTriangle, (entryEdge + 1) % 3, distance));
        }
        else if ((orientationS1 == CGAL::NEGATIVE || orientationS1 == CGAL::COPLANAR) && orientationS2 == CGAL::NEGATIVE)
        {
            Segment_2 cgalSegment((*points)[triangles[neighborTriangle].getID(Triangle::CCW(entryEdge))].getCGALPoint(), checkPointCart.getCGALPoint());
            double distance = squared_distance(cgalSegment, startPoint.getCGALPoint());
            segmentQueue.push(make_tuple(currentSegment, neighborTriangle, (entryEdge + 2) % 3, distance));
        }
        else
        {
            cout << "This should not happen ...." << endl;
            cout << orientationS1 << " " << orientationS2 << endl;
        }
    }
}

ID Triangulation::computeSpannerNode(ID startPoint, const Point &leftEndCone, const Point &rightEndCone, const vector<bool>& shrunkNodes) {
    points->push_back(leftEndCone); // points->size() - 2 Remember to delete
    points->push_back(rightEndCone); // points->size() - 1
    ID leftEndConeID = points->size() - 2, rightEndConeID = points->size() - 1;
    typedef pair<ID, ID> Segment; // A tuple only determines the direction of the cone. For distance calculations one needs the actual segment
    ID closestPoint = -1;
    double distanceToClosestPoint = maxWeight;
    Point startPointCart = (*this->points)[startPoint];

    std::priority_queue<tuple<Segment, int, int, double>, std::vector<tuple<Segment, int, int, double>>, Comparator> segmentQueue;
    for (const int &neighborTriangle : this->trianglesAtNode[startPoint])
    {
        int startIndex = this->triangles[neighborTriangle].getIndex(startPoint);
        CGAL::Orientation leftConeWithLeftOrientation = Point::orientation((*points)[startPoint], (*points)[this->triangles[neighborTriangle].getID(Triangle::CW(startIndex))], leftEndCone);
        CGAL::Orientation leftConeWithRightOrientation = Point::orientation((*points)[startPoint], (*points)[this->triangles[neighborTriangle].getID(Triangle::CCW(startIndex))], leftEndCone);
        CGAL::Orientation rightConeWithRightOrientation = Point::orientation((*points)[startPoint], (*points)[this->triangles[neighborTriangle].getID(Triangle::CCW(startIndex))], rightEndCone);
        CGAL::Orientation rightConeWithLeftOrientation = Point::orientation((*points)[startPoint], (*points)[this->triangles[neighborTriangle].getID(Triangle::CW(startIndex))], rightEndCone);
        bool leftEndConeInsideTriangle = leftConeWithLeftOrientation == CGAL::NEGATIVE && leftConeWithRightOrientation == CGAL::POSITIVE;
        bool rightEndConeInsideTriangle = rightConeWithLeftOrientation == CGAL::NEGATIVE && rightConeWithRightOrientation == CGAL::POSITIVE;
        if (leftEndConeInsideTriangle && !rightEndConeInsideTriangle) {
            Segment segment = make_pair(this->triangles[neighborTriangle].getID(Triangle::CCW(startIndex)), leftEndConeID);
            // We do not compute the real intersection, but just use the edge of the triangle to calculate the distance. This gives a distance, which is too small, but that's okay
            Segment_2 cgalSegment((*this->points)[this->triangles[neighborTriangle].getID(Triangle::CW(startIndex))].getCGALPoint(),
                (*points)[ this->triangles[neighborTriangle].getID(Triangle::CCW(startIndex))].getCGALPoint());
            double distance = squared_distance(cgalSegment, startPointCart.getCGALPoint());
            segmentQueue.emplace(make_tuple(segment, neighborTriangle, startIndex, distance));
            double dist = Point::squaredDistance(startPointCart, (*points)[this->triangles[neighborTriangle].getID(Triangle::CCW(startIndex))]);
            ID closestPointCandidate = this->triangles[neighborTriangle].getID(Triangle::CCW(startIndex));
            if (closestPoint == -1 && shrunkNodes[closestPointCandidate]) {
                closestPoint = closestPointCandidate;
                distanceToClosestPoint = dist;
            } else if (shrunkNodes[closestPointCandidate]){
                if (dist < distanceToClosestPoint) {
                    closestPoint = closestPointCandidate;
                    distanceToClosestPoint = dist;
                }
            }

        } else if (!leftEndConeInsideTriangle && rightEndConeInsideTriangle) {
            Segment segment = make_pair(rightEndConeID, this->triangles[neighborTriangle].getID(Triangle::CW(startIndex)));

            Segment_2 cgalSegment((*this->points)[this->triangles[neighborTriangle].getID(Triangle::CW(startIndex))].getCGALPoint(),
                (*points)[ this->triangles[neighborTriangle].getID(Triangle::CCW(startIndex))].getCGALPoint());
            double distance = squared_distance(cgalSegment, startPointCart.getCGALPoint());
            segmentQueue.emplace(make_tuple(segment, neighborTriangle, startIndex, distance));
            double dist = Point::squaredDistance(startPointCart, (*points)[this->triangles[neighborTriangle].getID(Triangle::CW(startIndex))]);
            ID closestPointCandidate = this->triangles[neighborTriangle].getID(Triangle::CW(startIndex));
            if (closestPoint == -1 && shrunkNodes[closestPointCandidate]) {
                closestPoint = closestPointCandidate;
                distanceToClosestPoint = dist;
            } else if (shrunkNodes[closestPointCandidate]) {
                if (dist < distanceToClosestPoint) {
                    distanceToClosestPoint = dist;
                    closestPoint = closestPointCandidate;
                }
            }
        } else if (leftEndConeInsideTriangle && rightEndConeInsideTriangle) {
            Segment segment = make_pair(rightEndConeID, leftEndConeID);
            Segment_2 cgalSegment((*this->points)[this->triangles[neighborTriangle].getID(Triangle::CW(startIndex))].getCGALPoint(),
                (*points)[ this->triangles[neighborTriangle].getID(Triangle::CCW(startIndex))].getCGALPoint());
            double distance = squared_distance(cgalSegment, startPointCart.getCGALPoint());
            segmentQueue.emplace(make_tuple(segment, neighborTriangle, startIndex, distance));
        } else if (leftConeWithLeftOrientation == CGAL::LEFT_TURN && rightConeWithRightOrientation == CGAL::RIGHT_TURN) {
            // This check is for the case where the triangle is completely contained in the cone.
            // So far the test is also positive when the triangle is on the opposite side and really wide.
            // So check that the triangle is on the right side.
            if (rightConeWithLeftOrientation == CGAL::RIGHT_TURN && leftConeWithRightOrientation == CGAL::LEFT_TURN) {
                Segment segment = make_pair(this->triangles[neighborTriangle].getID(Triangle::CCW(startIndex)), this->triangles[neighborTriangle].getID(Triangle::CW(startIndex)));
                Segment_2 cgalSegment((*this->points)[this->triangles[neighborTriangle].getID(Triangle::CW(startIndex))].getCGALPoint(),
                    (*points)[ this->triangles[neighborTriangle].getID(Triangle::CCW(startIndex))].getCGALPoint());
                double distance = squared_distance(cgalSegment, startPointCart.getCGALPoint());
                segmentQueue.emplace(make_tuple(segment, neighborTriangle, startIndex, distance));
                double dist = Point::squaredDistance(startPointCart, (*points)[this->triangles[neighborTriangle].getID(Triangle::CCW(startIndex))]);
                ID closestPointCandidate = this->triangles[neighborTriangle].getID(Triangle::CCW(startIndex));
                if (closestPoint == -1 && shrunkNodes[closestPointCandidate]) {
                    closestPoint = closestPointCandidate;
                    distanceToClosestPoint = dist;
                } else if (shrunkNodes[closestPointCandidate]) {
                    if (dist < distanceToClosestPoint) {
                        distanceToClosestPoint = dist;
                        closestPoint = closestPointCandidate;
                    }
                }
                dist = Point::squaredDistance(startPointCart, (*points)[this->triangles[neighborTriangle].getID(Triangle::CW(startIndex))]);
                if (dist < distanceToClosestPoint && shrunkNodes[closestPointCandidate]) {
                    distanceToClosestPoint = dist;
                    closestPoint = this->triangles[neighborTriangle].getID(Triangle::CW(startIndex));
                }
            }

        }
    }
    this->computeSpannerNodeInnerWorkings(startPoint, segmentQueue, closestPoint, distanceToClosestPoint, shrunkNodes);
    points->pop_back();
    points->pop_back();
    return closestPoint;
}

void Triangulation::computeSpannerNodeInnerWorkings(ID startPointID,
    std::priority_queue<std::tuple<std::pair<ID, ID>, int, int, double>, std::vector<std::tuple<std::pair<ID, ID>, int, int, double> >, Comparator> &segmentQueue,
    ID &spannerNode, double &distanceToClosestPoint, const vector<bool>& shrunkNodes) {
    typedef std::pair<ID, ID> Segment;
    Point startPoint = (*points)[startPointID];
    while (!segmentQueue.empty())
    {
//                cout << segmentQueue.size() << endl;
        tuple<Segment, int, int, double> currentStuff = segmentQueue.top();
        segmentQueue.pop();
        Segment currentSegment = get<0>(currentStuff); //.first is the ccw end, .second the cw one
        int originTriangle = get<1>(currentStuff);
        int edgeID = get<2>(currentStuff);
        double currentDistance = get<3>(currentStuff);
        if (currentDistance > distanceToClosestPoint) {
            return;
        }
//        this->visitedTriangles[originTriangle] = true;
//        cout << originTriangle << endl;
        Point segment2Cart = (*points)[currentSegment.second];
        Point segment1Cart = (*points)[currentSegment.first];

        int neighborTriangle = this->neighbors[originTriangle][edgeID];
        if (neighborTriangle == -1)
        {
            auto orientation1 = Point::orientation(startPoint, segment2Cart, (*points)[triangles[originTriangle].getID((edgeID + 2) % 3)]);
            auto orientation2 = Point::orientation(startPoint, segment1Cart, (*points)[triangles[originTriangle].getID((edgeID + 1) % 3)]);
            if (orientation1 == CGAL::ZERO || orientation1 == CGAL::COPLANAR)
            {
                double distance = Point::squaredDistance(startPoint, (*points)[triangles[originTriangle].getID((edgeID + 2) % 3)]);
                ID spannerNodeCandidate = triangles[originTriangle].getID((edgeID + 2) % 3);
                if (distance < distanceToClosestPoint && shrunkNodes[spannerNodeCandidate]) {
                    distanceToClosestPoint = distance;
                    spannerNode = spannerNodeCandidate;
                }
            }
            if (orientation2 == CGAL::ZERO || orientation2 == CGAL::COPLANAR)
            {
                ID spannerNodeCandidate = triangles[originTriangle].getID((edgeID + 1) % 3);
                double distance = Point::squaredDistance(startPoint, (*points)[triangles[originTriangle].getID((edgeID + 1) % 3)]);
                if (distance < distanceToClosestPoint && shrunkNodes[spannerNodeCandidate]) {
                    distanceToClosestPoint = distance;
                    spannerNode = spannerNodeCandidate;
                }
            }
            continue;
        }

        int entryEdge = this->triangles[neighborTriangle].indexEdge(this->triangles[originTriangle].getID((edgeID + 1) % 3),
                                                                    this->triangles[originTriangle].getID((edgeID + 2) % 3));

        ID checkPointID = triangles[neighborTriangle].getID(entryEdge);
        Point checkPointCart = (*points)[triangles[neighborTriangle].getID(entryEdge)];

        auto orientationS1 = Point::orientation(startPoint, segment1Cart, checkPointCart);
        auto orientationS2 = Point::orientation(startPoint, segment2Cart, checkPointCart);


        if (orientationS1 == CGAL::POSITIVE && orientationS2 == CGAL::NEGATIVE)
        {
            Segment segment1 = make_pair(currentSegment.first, checkPointID);
            Segment_2 cgalSegment1((*points)[triangles[neighborTriangle].getID(Triangle::CW(entryEdge))].getCGALPoint(), (*points)[triangles[neighborTriangle].getID(entryEdge)].getCGALPoint());
            double distance = squared_distance(cgalSegment1, startPoint.getCGALPoint());
            segmentQueue.emplace(make_tuple(segment1, neighborTriangle, (entryEdge + 1) % 3, distance));

            Segment segment2 = make_pair(checkPointID, currentSegment.second);
            Segment_2 cgalSegment2((*points)[triangles[neighborTriangle].getID(entryEdge)].getCGALPoint(), (*points)[triangles[neighborTriangle].getID(Triangle::CCW(entryEdge))].getCGALPoint());
            double distance2 = squared_distance(cgalSegment2, startPoint.getCGALPoint());
            segmentQueue.push(make_tuple(segment2, neighborTriangle, (entryEdge + 2) % 3, distance2));
        }
        else if (orientationS1 == CGAL::POSITIVE && (orientationS2 == CGAL::POSITIVE || orientationS2 == CGAL::COPLANAR))
        {
            Segment_2 cgalSegment((*points)[triangles[neighborTriangle].getID(Triangle::CW(entryEdge))].getCGALPoint(), checkPointCart.getCGALPoint());
            double distance = squared_distance(cgalSegment, startPoint.getCGALPoint());
            segmentQueue.push(make_tuple(currentSegment, neighborTriangle, (entryEdge + 1) % 3, distance));
        }
        else if ((orientationS1 == CGAL::NEGATIVE || orientationS1 == CGAL::COPLANAR) && orientationS2 == CGAL::NEGATIVE)
        {
            Segment_2 cgalSegment((*points)[triangles[neighborTriangle].getID(Triangle::CCW(entryEdge))].getCGALPoint(), checkPointCart.getCGALPoint());
            double distance = squared_distance(cgalSegment, startPoint.getCGALPoint());
            segmentQueue.push(make_tuple(currentSegment, neighborTriangle, (entryEdge + 2) % 3, distance));
        }
        else
        {
            cout << "This should not happen ...." << endl;
            cout << orientationS1 << " " << orientationS2 << endl;
        }
    }
}


int Triangulation::Triangle::indexEdge(ID p1, ID p2) const {
    for (int i = 0; i < 3; ++i)
    {
        if ((this->triangleIndices[i] == p1 && this->triangleIndices[(i + 1) % 3] == p2) ||
            (this->triangleIndices[i] == p2 && this->triangleIndices[(i + 1) % 3] == p1))
        {
            return (i + 2) % 3;
        }
    }
    return -1;
}

void Triangulation::createVisibilityGraph(std::vector<std::vector<GlobalID>> &edges, std::vector<std::vector<double>> &costs,
    const bool steinerPointsAsTargets, const std::vector<bool>& markedTriangles) const
{
    // We're only interested in the points on constrained edges. But to keep indices the same, we're going to keep every point, but won't
    // have any incoming or outgoing edges to those on high sea
    edges.resize((*points).size());
    cout << "Creating all the visibility" << endl;
    Timer timer;
    timer.start();
    std::vector<uint> pointsToStartFrom(this->points->size(), 0);
    if (markedTriangles.empty()) {
        // pointsToStartFrom = this->getConstrainedPoints();
        iota(pointsToStartFrom.begin(), pointsToStartFrom.end(), 0);
    } else {
        pointsToStartFrom = this->getMarkedConstrainedPoints(markedTriangles);
    }

    timer.mark();
    cout << "Time to find the start points of the visibility creation: " << timer.msecs() << std::endl;
    ProgressBar pb;
    pb.start(points->size());
#pragma omp parallel for shared(pb)
    for (int i = 0; i < pointsToStartFrom.size(); ++i) {
        pb.update(i);
        uint point = pointsToStartFrom[i];
        edges[point] = this->oneToAllVisibility(point, steinerPointsAsTargets, markedTriangles);
    }
    cout << "\n Now calculating the edge costs" << endl;
    costs.resize(edges.size());
    ProgressBar pb2;
    pb2.start(edges.size());
    for (int i = 0; i < edges.size(); ++i)
    {
        pb2.update(i);
//        if (i % 100 == 0)
//        {
//            timer.stop();
//            cout << float(i) / float((*points).size()) * 100. << "% Time so far: " << timer.secs() << " Remaining: "
//                 << timer.secs() * float((*points).size()) / float(i) - timer.secs() << endl;
//            timer.cont();
//        }
        costs[i].resize(edges[i].size());
        for (int j = 0; j < edges[i].size(); ++j)
        {
            costs[i][j] = Point::distance((*points)[i], (*points)[edges[i][j]]);
        }
    }
    timer.stop();
    cout << "\n Creating Visibility Graph took " << timer.secs() << " seconds." << endl;
}

void Triangulation::createVisibilityGraph(std::vector<std::vector<GlobalID>> &edges, std::vector<std::vector<double>> &costs,
    const bool steinerPointsAsTargets, const vector<uint>& startPoints, const std::vector<bool>& markedTriangles) const
{
    // We're only interested in the points on constrained edges. But to keep indices the same, we're going to keep every point, but won't
    // have any incoming or outgoing edges to those on high sea
    edges.resize((*points).size());
    Timer timer;
    timer.start();
    ProgressBar pb;
    pb.start(points->size());
#pragma omp parallel for shared(pb)
    for (int i = 0; i < startPoints.size(); ++i) {
        pb.update(i);
        uint point = startPoints[i];
        edges[point] = this->oneToAllVisibility(point, steinerPointsAsTargets, markedTriangles);
    }
    cout << "\n Now calculating the edge costs" << endl;
    costs.resize(edges.size());
    ProgressBar pb2;
    pb2.start(edges.size());
    for (int i = 0; i < edges.size(); ++i)
    {
        pb2.update(i);
//        if (i % 100 == 0)
//        {
//            timer.stop();
//            cout << float(i) / float((*points).size()) * 100. << "% Time so far: " << timer.secs() << " Remaining: "
//                 << timer.secs() * float((*points).size()) / float(i) - timer.secs() << endl;
//            timer.cont();
//        }
        costs[i].resize(edges[i].size());
        for (int j = 0; j < edges[i].size(); ++j)
        {
            costs[i][j] = Point::distance((*points)[i], (*points)[edges[i][j]]);
        }
    }
    timer.stop();
    cout << "\n Creating Visibility Graph took " << timer.secs() << " seconds." << endl;
}

void Triangulation::saveTriangulationAsGraph(const std::string &filename) const {
    std::ofstream writer(filename);
    writer << this->points->size() << "\n" << this->triangles.size() << "\n";
    writer << std::setprecision(20);
    for (const Point& point: *this->points) {
        writer << point.x() << " " << point.y() << "\n";
    }
    for (const Triangle& triangle: this->triangles) {
        writer << triangle.getID(0) << " " << triangle.getID(1) << " " << triangle.getID(2) << "\n";
    }
    writer.close();
}

void Triangulation::saveTriangulationAsGraph(const std::string &filename, const std::vector<bool> &markedTriangles) const {
    // First count the marked triangles
    uint counter = 0;
    for (bool markedTriangle : markedTriangles) {
        if (markedTriangle) {
            ++counter;
        }
    }
    std::ofstream writer(filename);
    writer << this->points->size() << "\n" << counter << "\n";
    writer << std::setprecision(20);
    for (const Point& point: *this->points) {
        writer << point.x() << " " << point.y() << "\n";
    }
    for (int i = 0; i < this->triangles.size(); ++i) {
        if (markedTriangles[i]) {
            const Triangle& triangle = this->triangles[i];
            writer << triangle.getID(0) << " " << triangle.getID(1) << " " << triangle.getID(2) << "\n";
        }
    }
}

void Triangulation::saveVisibilityGraph(const std::vector<std::vector<GlobalID>> &edges,
                                        const std::vector<std::vector<double>> &costs, const std::string &filename) const {
    ofstream outfile;
    outfile.open(filename);
    outfile << "#\n#\n#\n#\n\n";
    outfile << (this->points)->size() << "\n";
    unsigned long edgesCount = 0;
    for (int i = 0; i < edges.size(); ++i)
    {
        edgesCount += edges[i].size();
    }
    outfile << edgesCount << "\n";
    std::cout << "Wrote points to file." << std::endl;
    outfile << setprecision(20);
    for (int i = 0; i < this->points->size(); ++i)
    {
        outfile << i << " " << i << " " << (*this->points)[i].x() << " " << (*this->points)[i].y() << " 0\n";
    }
    for (int i = 0; i < edges.size(); ++i)
    {
        for (int j = 0; j < edges[i].size(); ++j)
        {
            outfile << i << " " << edges[i][j] << " " << costs[i][j] << " 0 0\n";
        }
    }
    outfile.close();
    cout << "Wrote to file" << endl;
}

bool Triangulation::Triangle::pointInsideTriangle(const Point &point) const {
    // Returns whether point lays inside this triangle

    CGAL::Orientation orientation1 = Point::orientation(parentTriangulation->getPoint(this->getID(0)),
        parentTriangulation->getPoint(this->getID(1)), point);
    CGAL::Orientation orientation2 = Point::orientation(parentTriangulation->getPoint(this->getID(1)),
        parentTriangulation->getPoint(this->getID(2)), point);
    CGAL::Orientation orientation3 = Point::orientation(parentTriangulation->getPoint(this->getID(2)),
        parentTriangulation->getPoint(this->getID(0)), point);
    if ((orientation1 == CGAL::Sign::LEFT_TURN || orientation1 == CGAL::Sign::COLLINEAR) &&
        (orientation2 == CGAL::Sign::LEFT_TURN || orientation2 == CGAL::Sign::COLLINEAR) &&
        (orientation3 == CGAL::Sign::LEFT_TURN || orientation3 == CGAL::Sign::COLLINEAR)) {
        return true;
    }
    return false;
}


TriangleIndex Triangulation::findTriangleContainingPoint(const Point &point, const ClosestPointSearch &closestPointSearch) const {
    // first find point closest to the search point. Most likely the triangle we're looking for is adjacent
    // If not we use these triangles as starting points for a directed search, greedily picking the triangle closest to
    // the point we're looking for
    ID closestPoint = closestPointSearch.getClosestPoint(point);
    for (TriangleIndex neighborTriang: this->trianglesAtNode[closestPoint]) {
        if (this->triangles[neighborTriang].pointInsideTriangle(point)) {
            return neighborTriang;
        }
    }
    // If we reach this point in the code, the point was not in an adjacent triangle of the nearest neighbor
    // Then we start a search. We use a PQ, checking the closest triangle first
    vector<bool> lookedAtTriangles(this->triangles.size(), false);
    std::priority_queue<std::pair<double, TriangleIndex>, std::vector<std::pair<double, TriangleIndex>>,
        std::greater<std::pair<double, TriangleIndex>>> pq;
    for (TriangleIndex neighborTriang: this->trianglesAtNode[closestPoint]) {
        const Point& triangleCorner = (*this->points)[this->triangles[neighborTriang].getID(0)];
        pq.push(std::make_pair(triangleCorner.distance(point), neighborTriang));
        if (neighborTriang != -1) {
            lookedAtTriangles.at(neighborTriang) = true;
        }
    }
    while (!pq.empty()) {
        TriangleIndex currentTriangleIndex = pq.top().second;
        pq.pop();
        if (this->triangles[currentTriangleIndex].pointInsideTriangle(point)) {
            return currentTriangleIndex;
        }
        for (const TriangleIndex& neighbor: this->neighbors[currentTriangleIndex]) {
            if (neighbor == -1) {
                continue;
            }
            if (!lookedAtTriangles.at(neighbor)) {
                const Point& triangleCorner = (*this->points)[this->triangles[neighbor].getID(0)];
                pq.push(std::make_pair(triangleCorner.distance(point), neighbor));
                if (neighbor != -1) {
                    lookedAtTriangles.at(neighbor) = true;
                }
            }
        }
    }
    return -1;
}

std::shared_ptr<Points> Triangulation::getPoints() const {
    return Triangulation::points;
}

int Triangulation::getNumTriangles() const {
    return this->triangles.size();
}

void Triangulation::outputAsGL(const std::string& filename, const std::vector<bool> &highlightedTriangles,
    bool highlightCoastlines) {
    cout << "Outputting file" << endl;
    ofstream myfile;
    myfile.open(filename);

    int amountPoints = this->points->size();
    int amountEdges = this->triangles.size() * 3;

    myfile << amountPoints << "\n"
           << amountEdges << "\n";
    myfile << setprecision(20);
    for (int i = 0; i < this->points->size(); ++i)
    {
        WGS84Coordinates coords = (*this->points)[i].pointInWGS84();
        myfile << coords.first << " " << coords.second << "\n";

    }
    for (int i = 0; i < this->triangles.size(); ++i) {
        int thickness;
        int color;
        if (!highlightedTriangles.empty() && highlightedTriangles[i])
        {
            thickness = 1;
            color = 2;
        }
        else
        {
            thickness = 1;
            color = 1;
        }
        if (!highlightCoastlines) {
            myfile << triangles[i].getID(0) << " " << triangles[i].getID(1) << " " << thickness << " " << color << "\n";
            myfile << triangles[i].getID(1) << " " << triangles[i].getID(2) << " " << thickness << " " << color << "\n";
            myfile << triangles[i].getID(2) << " " << triangles[i].getID(0) << " " << thickness << " " << color << "\n";
        } else {
            int highlightColor = 5;
            int thisColor;
            for (int j = 0; j < 3; ++j) {
                if (neighbors[i][(j + 2)%3] == -1) {
                    thisColor = highlightColor;
                } else {
                    thisColor = color;
                }
                myfile << triangles[i].getID(j) << " " << triangles[i].getID(j + 1) << " " << thickness << " " << thisColor << "\n";
            }
        }
    }
    myfile.close();
    cout << "Wrote to file" << endl;
}

Triangulation::Triangle Triangulation::getTriangle(TriangleIndex index) const {
    return this->triangles[index];
}

void Triangulation::writeTrianglesToFile(const std::string& filename, const std::vector<ID>& IDMapping) {
    std::ofstream writer(filename);
    writer << this->triangles.size() << "\n";
    for (const Triangle& triangle: this->triangles) {
        if (IDMapping.empty()) {
            writer << triangle.getID(0) << " " << triangle.getID(1) << " " << triangle.getID(2) << "\n";
        } else {
            writer << IDMapping[triangle.getID(0)] << " " << IDMapping[triangle.getID(1)] << " "
                << IDMapping[triangle.getID(2)] << "\n";
        }
    }
    writer.close();
}

void Triangulation::setPoints(std::shared_ptr<Points> pts) {
    points = pts;
}

bool Triangulation::pointInTriangulation(ID pointID) const {
    return this->pointsInTriangulation[pointID];
}

std::vector<bool> Triangulation::getPointsInTriangulation() const {
    return this->pointsInTriangulation;
}

void Triangulation::saveTriangulationAsFMI(const std::string &filename, int edgeWeightScaling) const {
    long num_nodes = points->size();

    struct Edge
    {
        ID target;
        double dist;
        Edge() {}
        Edge(ID target_, double dist_) : target(target_), dist(dist_) {}
    };
    long num_edges = 0;
    vector<vector<Edge>> edges_per_node(num_nodes);
    {
        struct FullEdge
        {
            ID source, target;
            FullEdge() {}
            FullEdge(ID source_, ID target_) : source(source_), target(target_) {}
        };
        vector<FullEdge> all_edges;
        for (const Triangle &triang : triangles)
        {
            for (int i = 0; i < 3; i++)
            {
                for (int j = i + 1; j < 3; j++)
                {
                    all_edges.push_back(FullEdge(triang.getID(i), triang.getID(j)));
                    all_edges.push_back(FullEdge(triang.getID(j), triang.getID(i)));
                }
            }
        }
        sort(all_edges.begin(), all_edges.end(), [](const FullEdge &e1, const FullEdge &e2)
             { return e1.source < e2.source || (e1.source == e2.source && e1.target < e2.target); });
        for (long i = 0; i < all_edges.size(); i++)
        {
            if (i == 0 || all_edges.at(i).source != all_edges.at(i - 1).source || all_edges.at(i).target != all_edges.at(i - 1).target)
            {
                double dist = Point::distance((*points)[all_edges.at(i).source], (*points)[all_edges.at(i).target]);
                edges_per_node.at(all_edges.at(i).source).push_back(Edge(all_edges.at(i).target, dist));
                num_edges++;
            }
        }
    }
    vector<bool> is_reached(num_nodes, false);
    // DEBUG
    if (num_nodes > 0)
    { // check if graph is connected
        vector<ID> queue;
        long queue_index = 0;
        queue.push_back(0);
        is_reached.at(0) = true;
        while(queue_index < queue.size())
        {
            ID current_node = queue.at(queue_index++);
            for(const Edge& edge : edges_per_node.at(current_node))
            {
                if(!is_reached.at(edge.target))
                {
                    is_reached.at(edge.target) = true;
                    queue.push_back(edge.target);
                }
            }
        }
        if(queue.size() != num_nodes)
        {
            cout << "WARNING: only " << queue.size() << " of " << num_nodes << " nodes are reachable from node 0." << endl;
            /*vector<NodeID> map_node_ids(num_nodes, CC::NO_ENTRY);

            long counter = 0;
            for(long i = 0; i < queue.size(); i++)
            {
                map_node_ids.at(queue.at(i)) = i;
            }
            num_edges = 0;
            vector<vector<Edge>> new_edges_per_node(queue.size());
            for(long i = 0; i < queue.size(); i++)
            {
                for(const Edge& edge : edges_per_node.at(queue.at(i)))
                {
                    if(is_reached.at(edge.target))
                    {
                        //DEBUG
                        if(map_node_ids.at(edge.target) == CC::NO_ENTRY)
                        {
                            cout << "ERROR: inconsistent data in search" << endl;
                        }
                        //DEBUG END
                        new_edges_per_node.at(i).push_back(Edge(map_node_ids.at(edge.target), edge.dist));
                        num_edges++;
                    }
                }
                sort(new_edges_per_node.at(i).begin(), new_edges_per_node.at(i).end(), [](const Edge& e1, const Edge& e2){return e1.target < e2.target;});
            }
            edges_per_node = new_edges_per_node;
            num_nodes = queue.size();*/
        }
    }
    // DEBUG END
    long counter_zero_edges = 0;
    ofstream writer(filename);
    writer << setprecision(20);
    writer << "#\n#\n#\n#\n\n";
    writer << num_nodes << "\n"
           << num_edges << "\n";
    for (long i = 0; i < num_nodes; i++)
    {
        writer << i << " " << i << " " << this->points->at(i).x() << " " << this->points->at(i).y() << " 0\n";
    }
    if (edgeWeightScaling != -1)
    {
        for (long i = 0; i < num_nodes; i++)
        {
            for (const Edge &edge : edges_per_node.at(i))
            {
                long edge_weight = round(edge.dist * edgeWeightScaling);
                if(edge_weight > INT_MAX)
                {
                    cout << "ERROR: edge weight is too large, exiting " << edge_weight << endl;
                    exit(1);
                }
                if(edge_weight == 0)
                {
                    counter_zero_edges++;
                }
                if(edge_weight < 0) {
                    cout << "ERROR: edge weight is far too large. " << edge_weight << endl;
                }
                writer << i << " " << edge.target << " " << edge_weight << " 0 0\n";
            }
        }
    }
    else
    {
        for (long i = 0; i < num_nodes; i++)
        {
            for (const Edge &edge : edges_per_node.at(i))
            {
                writer << i << " " << edge.target << " " << edge.dist << " 0 0\n";
            }
        }
    }
    writer.close();
    if(counter_zero_edges > 0)
    {
        cout << "WARNING: There are " << counter_zero_edges << " edges with weight 0 after scaling, this is " << (100 * double(counter_zero_edges) / num_edges) << "\% of all edges" << endl;
    }
}

std::vector<TriangleIndex> Triangulation::getAdjacentTriangles(ID point) const {
    return this->trianglesAtNode.at(point);
}

Neighbors Triangulation::getNeighbors(TriangleIndex triangle) const {
    return this->neighbors.at(triangle);
}

Point Triangulation::getPoint(InternalID id) const {
    return (*Triangulation::points)[id];
}

double Triangulation::Triangle::maxSideLength() const {
    double maxLength = 0;
    double currentLength;
    Point p1,p2;
    for (int i = 0; i < 3; ++i) {
        p1 = parentTriangulation->getPoint(this->getID(i));
        p2 = parentTriangulation->getPoint(this->getID(i+1));
        currentLength = Point::distance(p1,p2);
        maxLength = std::max(maxLength, currentLength);
    }
    return maxLength;
}

double Triangulation::Triangle::minDistanceToTriangle(const Point &point) const {
    double minDistance = maxWeight;
    double dist;
    for (int i = 0; i < 3; ++i) {
        const Point& p1 = parentTriangulation->getPoint(this->getID(i));
        const Point& p2 = parentTriangulation->getPoint(this->getID(i+1));
        Segment_2 segment(p1.getCGALPoint(), p2.getCGALPoint());
        dist = sqrt(CGAL::squared_distance(segment, point.getCGALPoint()));
        minDistance = min(minDistance, dist);
    }
    return minDistance;
}

double Triangulation::Triangle::maxDistanceToTriangle(const Point &point) const {
    double maxDistance = - maxWeight;
    double dist;
    for (int i = 0; i < 3; ++i) {
        const Point& p1 = parentTriangulation->getPoint(this->getID(i));
        dist = point.distance(p1);
        maxDistance = max(maxDistance, dist);
    }
    return maxDistance;
}

const std::vector<Triangulation::Triangle>& Triangulation::getTriangles() const {
    return this->triangles;
}


std::array<ID, 3> Triangulation::Triangle::getTriangleIndices() const {
    return this->triangleIndices;
}

bool Triangulation::Triangle::intersectsAnotherTriangle(const Triangulation &otherTriang, TriangleIndex triangleID) const {
    Triangulation::Triangle otherTriangle = otherTriang.getTriangle(triangleID);
    auto otherPoints = otherTriang.getPoints();
    // First check if a corner lays inside the shrunk triangle
    for (int i = 0; i < 3; ++i) {
        Point corner = (*otherPoints)[otherTriangle.getID(i)];
        if (this->pointInsideTriangle(corner)) {
            return true;
        }
    }

    // Now it is still possible that there are only intersections of edges
    // Iterate over the shrunk triangle edges and check whether they intersect the ref triangle
    for (int i = 0; i < 3; ++i) {
        Point sourceEdge = parentTriangulation->getPoint(this->getID(i)), targetEdge = parentTriangulation->getPoint(this->getID(i+1));
        if (otherTriangle.lineIntersectsTriangle(sourceEdge, targetEdge)) {
            return true;
        }
    }
    // Should no check return true, there is no intersection
    return false;
}

double Triangulation::Triangle::area() const {
    const Point& p0 = parentTriangulation->getPoint(this->getID(0));
    const Point& p1 = parentTriangulation->getPoint(this->getID(1));
    const Point& p2 = parentTriangulation->getPoint(this->getID(2));
    CGAL::Triangle_2<K> triangle(p0.getCGALPoint(), p1.getCGALPoint(), p2.getCGALPoint());
    return triangle.area();
}


uint Triangulation::getAmountPoints() const {
    return this->amountPoints;
}

std::vector<uint> Triangulation::getMarkedConstrainedPoints(const std::vector <bool> &markedTriangles) const {
    std::vector<uint> markedConstrainedPoints;
    markedConstrainedPoints.reserve(this->points->size());
    std::vector<bool> lookedAtPoints(this->points->size(), false);
    for (TriangleIndex triangleIndex = 0; triangleIndex < this->getNumTriangles(); ++triangleIndex) {
        if (markedTriangles[triangleIndex]) {
            for (int i = 0; i < 3; ++i) {
                uint pointIndex = this->getTriangle(triangleIndex).getID(i);
                if (lookedAtPoints[pointIndex]) {
                    continue;
                }
                lookedAtPoints[pointIndex] = true;
                if (this->isConstrainedPoint(pointIndex)) {
                    markedConstrainedPoints.push_back(pointIndex);
                }
            }
        }
    }
    markedConstrainedPoints.shrink_to_fit();
    return markedConstrainedPoints;
}

std::vector<uint> Triangulation::getConstrainedPoints() const {
    std::vector<uint> constrainedPoints;
    constrainedPoints.reserve(this->points->size());
    for (int i = 0; i < this->points->size(); ++i) {
        if (this->isConstrainedPoint(i)) {
            constrainedPoints.push_back(i);
        }
    }
    constrainedPoints.shrink_to_fit();
    return constrainedPoints;
}

template<class Archive>
void Triangulation::serialize(Archive &ar, const unsigned int version) {
    ar & this->points;
    ar & this->pointsInTriangulation;
    ar & this->triangles;
    ar & this->neighbors;
    ar & this->trianglesAtNode;
    ar & this->amountPoints;
    ar & this->constrainedPoints;
    ar & this->ambigiousPoints;
}


void Triangulation::writeToBinary(const std::string &file) const {
    std::cout << "Start binary write" << std::endl;
    std::ofstream myHLfile(file, std::ios_base::out | std::ios_base::binary);
    boost::archive::binary_oarchive oa(myHLfile);
    oa &(*this);
}

void Triangulation::readFromBinary(const std::string &file) {
    std::cout << "Start binary read" << std::endl;
    std::ifstream myHLfile(file, std::ios_base::in | std::ios_base::binary);
    boost::archive::binary_iarchive ia(myHLfile);
    ia &(*this);
}

template<class Archive>
void Triangulation::Triangle::serialize(Archive &ar, const unsigned int version) {
    ar & this->triangleIndices;
    ar & this->parentTriangulation;
}

void Triangulation::readFromGraph(const std::string &input) {
    this->points = readPointsFromGraph(input);
    this->readFileTopologyFromGraph(input);
    this->calculateConstrainedAndAmbigiousPoints();
}

Triangulation::Triangle::Triangle() {

}

void Triangulation::createThetaGraph(std::vector<std::vector<GlobalID> > &edges, std::vector<std::vector<double> > &costs, uint numIntervals) {
    double maxX = -maxWeight, maxY = - maxWeight, minX = maxWeight, minY = maxWeight;
    for (Point& point : *this->points) {
        maxX = max(maxX, point.x());
        maxY = max(maxY, point.y());
        minX = min(minX, point.x());
        minY = min(minY, point.y());
    }
    // double radius = sqrt((maxX - minX)/2 * (maxX - minX)/2 + (maxY - minY)/2 * (maxY - minY)/2); // Think about this again when not tired
    double radius = 1;
    vector<Point> pointsDefiningAnglesGlobal;
    for (uint i = 0; i < numIntervals; i++) {
        double angle = 2*M_PI*i/numIntervals;
        double offset =  0.0001;
        pointsDefiningAnglesGlobal.emplace_back(Point(radius*cos(angle + offset), radius*sin(angle + offset)));
    }
    Timer timer;
    timer.start();
    ProgressBar pb;
    pb.start(points->size());
// #pragma omp parallel for shared(pb)
    for (int node = 0; node < points->size(); ++node) {
        pb.update(node);
        vector<Point> pointsDefiningAnglesLocal = pointsDefiningAnglesGlobal;
        for (Point& angle: pointsDefiningAnglesLocal) {
            angle = angle + (*this->points)[node];
        }
        vector<GlobalID> outgoingEdges;
        vector<double> outgoingCosts;
        outgoingEdges.reserve(numIntervals);
        for (int i = 0; i < numIntervals; ++i) {
            Point rightCone =  pointsDefiningAnglesLocal[i];
            Point leftCone =  pointsDefiningAnglesLocal[(i + 1)%numIntervals];
            ID spannerNode = computeSpannerNode(node, leftCone, rightCone);
            if (spannerNode == -1) {
                continue;
            }
            outgoingEdges.push_back(spannerNode);
            double distance = (*points)[node].distance((*points)[spannerNode]);
            outgoingCosts.push_back(distance);
        }
        edges.push_back(outgoingEdges);
        costs.push_back(outgoingCosts);
    }
    // Create reverse edges
    std::vector<std::vector<GlobalID> > reverseEdges(edges.size());
    std::vector<std::vector<double> > reverseCosts(edges.size());
    for (int source = 0; source < edges.size(); ++source) {
        for (int i = 0; i < edges[source].size(); ++i) {
            ID target = edges[source][i];
            double cost = costs[source][i];
            reverseEdges[target].push_back(source);
            reverseCosts[target].push_back(cost);
        }
    }
    for (int source = 0; source < reverseEdges.size(); ++source) {
        for (int i = 0; i < reverseEdges[source].size(); ++i) {
            ID target = reverseEdges[source][i];
            double cost = reverseCosts[source][i];
            edges[source].push_back(target);
            costs[source].push_back(cost);
        }
    }
    timer.stop();
    cout << "\n Creating Theta Graph took " << timer.secs() << " seconds." << endl;
}

void Triangulation::createThetaGraph(std::vector<std::vector<GlobalID> > &edges, std::vector<std::vector<double> > &costs, uint numIntervals, const std::vector<Point> &sourcePoints, const ClosestPointSearch& cps) {
    double maxX = -maxWeight, maxY = - maxWeight, minX = maxWeight, minY = maxWeight;
    for (Point& point : *this->points) {
        maxX = max(maxX, point.x());
        maxY = max(maxY, point.y());
        minX = min(minX, point.x());
        minY = min(minY, point.y());
    }
    // double radius = sqrt((maxX - minX)/2 * (maxX - minX)/2 + (maxY - minY)/2 * (maxY - minY)/2); // Think about this again when not tired
    double radius = 1;
    vector<Point> pointsDefiningAnglesGlobal;
    for (uint i = 0; i < numIntervals; i++) {
        double angle = 2*M_PI*i/numIntervals;
        double offset =  0.0001;
        pointsDefiningAnglesGlobal.emplace_back(Point(radius*cos(angle + offset), radius*sin(angle + offset)));
    }
    Timer timer;
    timer.start();
    ProgressBar pb;
    pb.start(points->size());
// #pragma omp parallel for shared(pb)
    for (int node = 0; node < sourcePoints.size(); ++node) {
        pb.update(node);
        vector<Point> pointsDefiningAnglesLocal = pointsDefiningAnglesGlobal;
        for (Point& angle: pointsDefiningAnglesLocal) {
            angle = angle + sourcePoints[node];
        }
        vector<GlobalID> outgoingEdges;
        vector<double> outgoingCosts;
        outgoingEdges.reserve(numIntervals);
        for (int i = 0; i < numIntervals; ++i) {
            Point rightCone =  pointsDefiningAnglesLocal[i];
            Point leftCone =  pointsDefiningAnglesLocal[(i + 1)%numIntervals];
            ID spannerNode = computeSpannerNode(sourcePoints[node], leftCone, rightCone, cps);
            if (spannerNode == -1) {
                continue;
            }
            outgoingEdges.push_back(spannerNode);
            double distance = sourcePoints[node].distance((*points)[spannerNode]);
            outgoingCosts.push_back(distance);
        }
        edges.push_back(outgoingEdges);
        costs.push_back(outgoingCosts);
        // DEBUG
        ofstream writer("theta_graph.txt");
        writer << setprecision(20);
        writer << 1 + numIntervals + outgoingEdges.size() << "\n";
        writer << numIntervals + outgoingEdges.size() << "\n";
        auto sourceWGS84 = sourcePoints[node].pointInWGS84();
        writer << sourceWGS84.first << " " << sourceWGS84.second << "\n";
        for (int i = 0; i < numIntervals; ++i) {
            auto intervalWGS84 = pointsDefiningAnglesLocal[i].pointInWGS84();
            writer << intervalWGS84.first << " " << intervalWGS84.second << "\n";
        }
        for (int i = 0; i < outgoingEdges.size(); ++i) {
            auto targetWGS84 = (*points)[outgoingEdges[i]].pointInWGS84();
            writer << targetWGS84.first << " " << targetWGS84.second << "\n";
        }
        for (int i = 0; i < numIntervals; ++i) {
            writer << "0 " << i + 1<< " 4 3\n";
        }
        for (int i = 0; i < outgoingEdges.size(); ++i) {
            writer << "0 " << i + 1 + numIntervals << " 4 4\n";
        }
        writer.close();
        string stuff;
        std::cin >> stuff;
        // DEBUG END
    }
    // Create reverse edges
    std::vector<std::vector<GlobalID> > reverseEdges(edges.size());
    std::vector<std::vector<double> > reverseCosts(edges.size());
    for (int source = 0; source < edges.size(); ++source) {
        for (int i = 0; i < edges[source].size(); ++i) {
            ID target = edges[source][i];
            double cost = costs[source][i];
            reverseEdges[target].push_back(source);
            reverseCosts[target].push_back(cost);
        }
    }
    for (int source = 0; source < reverseEdges.size(); ++source) {
        for (int i = 0; i < reverseEdges[source].size(); ++i) {
            ID target = reverseEdges[source][i];
            double cost = reverseCosts[source][i];
            edges[source].push_back(target);
            costs[source].push_back(cost);
        }
    }
    timer.stop();
    cout << "\n Creating Theta Graph took " << timer.secs() << " seconds." << endl;
}

void Triangulation::createThetaGraph(std::vector<std::vector<GlobalID> > &edges, std::vector<std::vector<double> > &costs, uint numIntervals, const vector<bool>& shrunkNodes) {
    double maxX = -maxWeight, maxY = - maxWeight, minX = maxWeight, minY = maxWeight;
    for (Point& point : *this->points) {
        maxX = max(maxX, point.x());
        maxY = max(maxY, point.y());
        minX = min(minX, point.x());
        minY = min(minY, point.y());
    }
    // double radius = sqrt((maxX - minX)/2 * (maxX - minX)/2 + (maxY - minY)/2 * (maxY - minY)/2); // Think about this again when not tired
    double radius = 1;
    vector<Point> pointsDefiningAnglesGlobal;
    for (uint i = 0; i < numIntervals; i++) {
        double angle = 2*M_PI*i/numIntervals;
        double offset =  0.0001;
        pointsDefiningAnglesGlobal.emplace_back(Point(radius*cos(angle + offset), radius*sin(angle + offset)));
    }
    Timer timer;
    timer.start();
    ProgressBar pb;
    pb.start(points->size());
// #pragma omp parallel for shared(pb)
    for (int node = 0; node < points->size(); ++node) {
        pb.update(node);
        vector<Point> pointsDefiningAnglesLocal = pointsDefiningAnglesGlobal;
        for (Point& angle: pointsDefiningAnglesLocal) {
            angle = angle + (*this->points)[node];
        }
        vector<GlobalID> outgoingEdges;
        vector<double> outgoingCosts;
        outgoingEdges.reserve(numIntervals);
        for (int i = 0; i < numIntervals; ++i) {
            Point rightCone =  pointsDefiningAnglesLocal[i];
            Point leftCone =  pointsDefiningAnglesLocal[(i + 1)%numIntervals];
            ID spannerNode = computeSpannerNode(node, leftCone, rightCone, shrunkNodes);
            if (spannerNode == -1) {
                continue;
            }
            outgoingEdges.push_back(spannerNode);
            double distance = (*points)[node].distance((*points)[spannerNode]);
            outgoingCosts.push_back(distance);
        }
        edges.push_back(outgoingEdges);
        costs.push_back(outgoingCosts);
    }
    // Create reverse edges
    std::vector<std::vector<GlobalID> > reverseEdges(edges.size());
    std::vector<std::vector<double> > reverseCosts(edges.size());
    for (int source = 0; source < edges.size(); ++source) {
        for (int i = 0; i < edges[source].size(); ++i) {
            ID target = edges[source][i];
            double cost = costs[source][i];
            reverseEdges[target].push_back(source);
            reverseCosts[target].push_back(cost);
        }
    }
    for (int source = 0; source < reverseEdges.size(); ++source) {
        for (int i = 0; i < reverseEdges[source].size(); ++i) {
            ID target = reverseEdges[source][i];
            double cost = reverseCosts[source][i];
            edges[source].push_back(target);
            costs[source].push_back(cost);
        }
    }
    timer.stop();
    cout << "\n Creating Theta Graph took " << timer.secs() << " seconds." << endl;
}

