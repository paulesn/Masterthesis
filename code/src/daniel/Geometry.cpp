//
// Created by kochda on 01.02.24.
//

#include <memory>
#include "Geometry.h"
#include <algorithm>
#include "Structs.h"

#include "Progressbar.h"
#include "Timer.h"

Point::Point() {
    this->point = Point_2();
}

Point::Point(double x, double y) {
    this->point = Point_2(x,y);
}

Point::Point(Point_2 &point) {
    this->point = point;
}


bool Point::operator==(const Point &p2) const {
    return (this->x() == p2.x() && this->y() == p2.y());
}

Point Point::operator+(const Point &p2) const {
    return Point(this->x() + p2.x(), this->y() + p2.y());
}

std::ostream &operator<<(std::ostream &os, const Point& pt) {
    os << pt.x() << " " << pt.y();
    return os;
}

double Point::x() const {
    return this->point.x();
}

double Point::y() const {
    return this->point.y();
}

double Point::distance(const Point &p2) const {
    double xdiff = this->point.x() - p2.x();
    double ydiff = this->point.y() - p2.y();
    double dist = sqrt(xdiff * xdiff + ydiff * ydiff);
    return dist;
}

double Point::distance(const Point &p1, const Point &p2) {
    return p1.distance(p2);
}

double Point::squaredDistance(const Point &p2) const {
    double xdiff = this->point.x() - p2.x();
    double ydiff = this->point.y() - p2.y();
    double dist = xdiff * xdiff + ydiff * ydiff;
    return dist;
}

double Point::squaredDistance(const Point &p1, const Point &p2) {
    return p1.squaredDistance(p2);
}

CGAL::Orientation Point::orientation(const Point &p1, const Point &p2, const Point &checkPoint) {
    return CGAL::orientation(p1.point, p2.point, checkPoint.point);
}

WGS84Coordinates Point::pointInWGS84() const {
    double lon = (this->x() / 20037508.34) * 180;
    double lat = (this->y() / 20037508.34) * 180;

    lat = 180 / M_PI * (2 * atan(exp(lat * M_PI / 180)) - M_PI / 2);
    return std::make_pair(lat, lon);
}

Point_2 Point::getCGALPoint() const {
    return this->point;
}

std::shared_ptr<Points> readPointsFromGraph(const std::string& filename) {
    // Reads points from a .graph file
    std::cout << "Reading Points." << std::endl;
    std::string temp;
    std::ifstream inputFile(filename);
    std::shared_ptr<Points> points = std::make_shared<Points>();
    std::string line;
    inputFile >> line;
    int numNodes = stoi(line);
    points->resize(numNodes);
    inputFile >> temp;

    for (int i = 0; i < numNodes; ++i) {
        double x,y;
        inputFile >> x >> y;
        Point point(x,y);
        (*points)[i] = point;
    }
    return points;
}

std::shared_ptr<Points> readPoints(const std::string& filename) {
    // Reads points from a .points file
    std::cout << "Reading Points." << std::endl;
    std::string temp;
    std::ifstream inputFile(filename);
    std::shared_ptr<Points> points = std::make_shared<Points>();
    std::string line;
    inputFile >> line;
    int numNodes = stoi(line);
    points->resize(numNodes);
    for (int i = 0; i < numNodes; ++i) {
        double x,y;
        inputFile >> x >> y;
        Point point(x,y);
        (*points)[i] = point;
    }
    return points;
}

ClosestPointSearch::ClosestPointSearch(std::shared_ptr<Points> pointsPointer, int numCellsPerDimension, const std::vector<bool>& activePoints) {
    this->setPoints(pointsPointer, numCellsPerDimension, activePoints);
}

ClosestPointSearch::ClosestPointSearch(std::shared_ptr<Points> pointsPointer, int numCellsPerDimension) {
    this->setPoints(pointsPointer, numCellsPerDimension);
}

void ClosestPointSearch::setPoints(std::shared_ptr<Points> pointsPointer, int numCellsPerDimension, const std::vector<bool> &activePoints) {
    Timer timer;
    timer.start();
    this->points = pointsPointer;
    double xMax = -maxWeight, yMax = -maxWeight;
    this->yMin = maxWeight;
    this->xMin = maxWeight;
    for (int i = 0; i < (this->points)->size(); ++i) {
        if (activePoints[i]) {
            Point& point = (*this->points)[i];
            this->xMin = std::min(this->xMin, point.x());
            xMax = std::max(xMax, point.x());
            this->yMin = std::min(this->yMin, point.y());
            yMax = std::max(yMax, point.y());
        }
    }
    this->nodesPerCell = std::vector<std::vector<std::vector<ID>>>(numCellsPerDimension, std::vector<std::vector<ID>>(numCellsPerDimension));
    this->xCellWidth = (xMax - this->xMin)/(numCellsPerDimension - 1);
    this->yCellWidth = (yMax - this->yMin)/(numCellsPerDimension - 1);
    ProgressBar progressBar;
    progressBar.start(this->points->size());
    for (int i = 0; i < this->points->size(); ++i) {
        progressBar.update(i);
        if (activePoints[i]) {
            const Point& point = (*this->points)[i];
            int xIndex = (point.x() - xMin)/xCellWidth;
            int yIndex = (point.y() - yMin)/yCellWidth;
            this->nodesPerCell[xIndex][yIndex].push_back(i);
        }
    }
    progressBar.stop();
    timer.stop();
    std::cout << "Created Closest Point Search. It took " << timer.secs() << "seconds." << std::endl;
}

void ClosestPointSearch::setPoints(std::shared_ptr<Points> pointsPointer, int numCellsPerDimension) {
    Timer timer;
    timer.start();
    this->points = pointsPointer;
    double xMax = -maxWeight, yMax = -maxWeight;
    this->yMin = maxWeight;
    this->xMin = maxWeight;
    for (int i = 0; i < (this->points)->size(); ++i) {
        Point& point = (*this->points)[i];
        this->xMin = std::min(this->xMin, point.x());
        xMax = std::max(xMax, point.x());
        this->yMin = std::min(this->yMin, point.y());
        yMax = std::max(yMax, point.y());
    }
    this->nodesPerCell = std::vector<std::vector<std::vector<ID>>>(numCellsPerDimension, std::vector<std::vector<ID>>(numCellsPerDimension));
    this->xCellWidth = (xMax - this->xMin)/(numCellsPerDimension - 1);
    this->yCellWidth = (yMax - this->yMin)/(numCellsPerDimension - 1);
    ProgressBar progressBar;
    progressBar.start(this->points->size());
    for (int i = 0; i < this->points->size(); ++i) {
        progressBar.update(i);
        const Point& point = (*this->points)[i];
        int xIndex = (point.x() - xMin)/xCellWidth;
        int yIndex = (point.y() - yMin)/yCellWidth;
        this->nodesPerCell[xIndex][yIndex].push_back(i);
    }
    progressBar.stop();
    timer.stop();
    std::cout << "Created Closest Point Search. It took " << timer.secs() << "seconds." << std::endl;
}


ID ClosestPointSearch::getClosestPoint(const Point &point) const {
    int xIndex = (point.x() - this->xMin)/this->xCellWidth;
    int yIndex = (point.y() - this->yMin)/this->yCellWidth;
    double minDist = maxWeight;
    int minIndex = -1;
    int layer = 0;
    int nodeCounter = 0; // Counts how many points we have seen. As long as it is zero, we should not terminate. As soon as we have seen a single point, we should terminate.
    while (nodeCounter == 0) {
        for (int xCell = xIndex - layer; xCell < xIndex + layer; ++xCell) {
            for (int yCell = yIndex - layer; yCell < yIndex + layer; ++yCell) {
                if (xCell < 0 || yCell < 0 || xCell >= nodesPerCell.size() || yCell >= nodesPerCell.at(0).size()) {
                    continue;
                }
                for (ID pointIndex: this->nodesPerCell[xCell][yCell]) {
                    double dist = std::max(std::abs((*this->points)[pointIndex].x() - point.x()), std::abs((*this->points)[pointIndex].y() - point.y()));
                    if (dist < minDist) {
                        minDist = dist;
                        minIndex = pointIndex;
                    }
                }
                nodeCounter += nodesPerCell[xCell][yCell].size();
            }
        }
        layer++;
    }
    return minIndex;
}

void ClosestPointSearch::printStats() const {
    int maxCellSize = 0;
    int minCellSize = std::numeric_limits<int>::max();
    int avgCellSize = 0;
    int avgCellSizeWithoutZeros = 0;
    int withoutZerosCounter = 0;
    for (int xIndex = 0; xIndex < this->nodesPerCell.size(); ++xIndex) {
        for (int yIndex = 0; yIndex < this->nodesPerCell[xIndex].size(); ++yIndex) {
            maxCellSize = std::max(maxCellSize, int(this->nodesPerCell[xIndex][yIndex].size()));
            minCellSize = std::min(minCellSize, int(this->nodesPerCell[xIndex][yIndex].size()));
            avgCellSize += this->nodesPerCell[xIndex][yIndex].size();
            if (this->nodesPerCell[xIndex][yIndex].size() != 0) {
                avgCellSizeWithoutZeros += this->nodesPerCell[xIndex][yIndex].size();
                withoutZerosCounter++;
            }
        }
    }
    std::cout << "Max Cell Size: " << maxCellSize << std::endl;
    std::cout << "Min Cell Size: " << minCellSize << std::endl;
    std::cout << "Avg Cell Size: " << avgCellSize/(this->nodesPerCell.size() * this->nodesPerCell.size())<< std::endl;
    std::cout << "Average Cell Size without empty cells: " << avgCellSizeWithoutZeros/withoutZerosCounter << std::endl;
}


void writePointsToFile(const std::shared_ptr<Points> &points, const std::string& filename) {
    std::ofstream writer(filename);
    writer << std::setprecision(20);
    writer << points->size() << "\n";
    for (const Point& point: *points) {
        writer << point.x() << " " << point.y() << "\n";
    }
    writer.close();
}

std::vector<int> readPointsInTriangles(const std::string &filename) {
    std::ifstream reader(filename);
    int numPoints;
    reader >> numPoints;
    std::vector<int> pointsInTriangles(numPoints, -1);
    int temp;
    for (int i = 0; i < numPoints; ++i) {
        reader >> temp;
        pointsInTriangles[i] = temp;
    }
    return pointsInTriangles;
}

bool twoLinesIntersectingEachOther(const Point& l1Start, const Point& l1End, const Point& l2Start, const Point& l2End) {
    if (Point::orientation(l1Start, l1End, l2Start) != Point::orientation(l1Start, l1End, l2End)) {
        if (Point::orientation(l2Start, l2End, l1Start) != Point::orientation(l2Start, l2End, l1End)) {
            return true;
        }
    }
    return false;
}

double distanceBetweenSegmentAndPoint(const Point& segmentStart, const Point& segmentEnd, const Point& point) {
    K::Segment_2 segment(segmentStart.getCGALPoint(), segmentEnd.getCGALPoint());
    return sqrt(CGAL::squared_distance(segment, point.getCGALPoint()));
}

std::vector<std::pair<ID, double>> readNearestNeighbors(const std::string& filename) {
    std::vector<std::pair<ID, double>> nearestNeighbors;
    std::ifstream reader(filename);
    std::string line;
    ID closestPointID;
    double distance;
    while(getline(reader, line)) {
        std::istringstream iss(line);
        iss >> closestPointID >> distance;
        nearestNeighbors.push_back(std::make_pair(closestPointID, distance));
    }
    reader.close();
    return nearestNeighbors;
}

std::vector<std::vector<int> > readTriangleIntersections(const std::string &filename) {
    std::ifstream reader(filename);
    std::string line;
    u_int shrunkTriangle, intersectedTriangle;
    std::vector<std::vector<int>> refTrianglesIntersectedByShrunkTriangles;
    while (std::getline(reader, line)) {
        refTrianglesIntersectedByShrunkTriangles.push_back(std::vector<int>());
        std::istringstream iss(line);
        iss >> shrunkTriangle;
        while (iss >> intersectedTriangle) {
            refTrianglesIntersectedByShrunkTriangles.back().push_back(intersectedTriangle);
        }
    }
    reader.close();
    return refTrianglesIntersectedByShrunkTriangles;
}

double haversineDistance(const WGS84Coordinates &p1, const WGS84Coordinates &p2) {
    // distance between latitudes
    // and longitudes
    double lon1 = p1.second;
    double lon2 = p2.second;
    double lat1 = p1.first;
    double lat2 = p2.first;

    double dLat = (lat2 - lat1) *
                  M_PI / 180.0;
    double dLon = (lon2 - lon1) *
                  M_PI / 180.0;

    // convert to radians
    lat1 = (lat1)*M_PI / 180.0;
    lat2 = (lat2)*M_PI / 180.0;

    double sinusdlat = sin(dLat / 2);
    double sinusdlon = sin(dLon / 2);
    // apply formulae
    double a = sinusdlat * sinusdlat +
               sinusdlon * sinusdlon *
                   cos(lat1) * cos(lat2);
    double rad = 6371;
    double c = 2 * asin(sqrt(a));
    return rad * c;
}

// Point Point::intersectTwoLines(const Point &p1, const Point &p2, const Point &q1, const Point &q2) {
//     Line_2 pSegment(p1.getCGALPoint(), p2.getCGALPoint());
//     Line_2 qSegment(q1.getCGALPoint(), q2.getCGALPoint());
//     auto resIntersection = CGAL::intersection(pSegment, qSegment);
//     Point_2* intersection = std::get_if<Point_2 >(&*resIntersection);
//     return Point(*intersection);
// }
