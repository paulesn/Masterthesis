

#include <math.h>
#include <vector>

#include "DataLoader.cpp"

struct OcttreeNode {
    std::vector<OcttreeNode> children;
    int objectCount;
    std::vector<Star*> objects; // list of object pointer
    double cX, cY, cZ, height;// bounding box it is always an edge-aligned cube the coordinates provide the center

    explicit OcttreeNode(double cX, double cY, double cZ, double height)
        : cX(cX), cY(cY), cZ(cZ), height(height), objectCount(0) {}

    bool intercept(
        double x, double y, double z, double r
    ) {
        // if the center of the cube and the center of the sphere are closer than the radius the sphere intersects the cube
        // else:
        //     compute the point on the sphere closest to the center of the cube
        //     if that point is inside the cube, the sphere intersects the cube
        // else:
        //     the sphere does not intersect the cube

        double centerDist = x*cX + y*cY + z*cZ;
        if (centerDist < r*r) {return true;}

        // calculate the closest point on the sphere to the center of the cube
        /**
         *c is the center of the cube
         *s is the center of the sphere (in variable no prefix)
         *s + d = c ==> d = c-s = x
         * to normalize x in its length, it is divided by its own length and then
         * multiplied with the radius of the sphere to provide a point on the sphere
         */
        double distX = cX-x;
        double distY = cY-y;
        double distZ = cZ-z;
        const double lDist = sqrt(distX*distX + distY*distY + distZ*distZ);
        distX /= lDist*r;
        distY /= lDist*r;
        distZ /= lDist*r;

        // check if the point is in the cube:
        if (distX < x-height || distX > x+height) {return false;}
        if (distY < y-height || distY > y+height) {return false;}
        if (distZ < z-height || distZ > z+height) {return false;}
        return true;
    }

    std::vector<Star*> find(
        double x, double y, double z, double r
    ) {
    /**
     * identify all objects in the octtree that are within the radius of the sphere
     */
        std::vector<Star*> result;
        if (intercept(x, y, z, r)) {
            if (std::size(children) > 0) {
                for (int i = 0; i < children.size(); ++i) {
                    auto childResult = children[i].find(x, y, z, r);
                    result.insert(result.end(), childResult.begin(), childResult.end());
                }
            } else {
                for (int i = 0; i < objects.size(); ++i) {
                    // check if the object is in the sphere
                    if (objects[i]->x < x-r || objects[i]->x > x+r) {continue;}
                    if (objects[i]->y < y-r || objects[i]->y > y+r) {continue;}
                    if (objects[i]->z < z-r || objects[i]->z > z+r) {continue;}
                    result.push_back(objects[i]);
                }
            }
        }
        return result;
    }

    void insert(
        const Star* s, int maxCount
    ) {
        int x = s->x;
        int y = s->y;
        int z = s->z;

        if (std::size(children) > 0) {
            for (auto child : children) {
                const double h = child.height;
                const double cX = child.cX;
                double cY = child.cY;
                double cZ = child.cZ;
                if (x < cX-h || x > cX+h) {continue;}
                if (y < cY-h || y > cY+h) {continue;}
                if (z < cZ-h || z > cZ+h) {continue;}
                child.insert(s, maxCount);
                child.objectCount++;
                return;
            }
        } else {
            // this is a leaf node we have to add the object here
            objects.push_back(s);
            // now we have to check if the leaf has gotten to big
            if (std::size(children) > maxCount) {
                // create 8 children and insert all objects again
                children.push_back(OcttreeNode(cX+height/2, cY+height/2, cZ+height/2, height/2));
                children.push_back(OcttreeNode(cX+height/2, cY+height/2, cZ-height/2, height/2));
                children.push_back(OcttreeNode(cX+height/2, cY-height/2, cZ+height/2, height/2));
                children.push_back(OcttreeNode(cX+height/2, cY-height/2, cZ-height/2, height/2));
                children.push_back(OcttreeNode(cX-height/2, cY+height/2, cZ+height/2, height/2));
                children.push_back(OcttreeNode(cX-height/2, cY+height/2, cZ-height/2, height/2));
                children.push_back(OcttreeNode(cX-height/2, cY-height/2, cZ+height/2, height/2));
                children.push_back(OcttreeNode(cX-height/2, cY-height/2, cZ-height/2, height/2));
                for (auto obj: objects) {
                    insert(obj, maxCount); // TODO get coords of the objects
                }
                objects.clear();
                objectCount = 0;
            }

        }
    }
};

struct Octtree
{
    OcttreeNode* root;
    int maxObjects;
    int height;

    Octtree(int maxObjects, int height)
        : root(nullptr), maxObjects(maxObjects), height(height) {
        root = new OcttreeNode(0, 0, 0, height);
    }
    ~Octtree() {
        delete root;
    }
    void insert(const Star* s) {
        root->insert(s, maxObjects);
    }
    std::vector<Star*> find(double x, double y, double z, double r) {
        return root->find(x, y, z, r);
    }
    std::vector<Star*> find(Star* s, double r) {
        return root->find(s->x, s->y, s->z, r);
    }
};