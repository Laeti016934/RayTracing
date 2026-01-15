#ifndef KDTREE_H
#define KDTREE_H

#include <vector>
#include "AABB.h"
#include "PrimitiveRef.h"
#include "Hit.h"

class Scene;   // forward declaration
class Ray;     // forward declaration

struct Node
{
    // Spatial data : AABB
    AABB boundingBox;

    // If internal node
    int splittingAxis;    // Splitting axis (X, Y or Z)
    float cuttingPos;       // Cutting position
    Node* leftChild;        // Left child 
    Node* rightChild;       // Right child

    // If leaf node
    std::vector<PrimitiveRef> primitives;    // List of primitives 

    bool isLeaf; // Node State 

    Node()
        : splittingAxis(-1),
          cuttingPos(0.f),
          leftChild(nullptr),
          rightChild(nullptr),
          isLeaf(false)
    {}
};

class KDTree
{
private:
    Node* root;                 // Root node
    int maxDepth;               // Maximum depth
    int maxPrimitivesPerLeaf;   // Max number of primitives per leaf
    float epsilon;              // Epsilon of minimal size

    // Construction récursive
    Node* buildNode(AABB &aabb, std::vector<PrimitiveRef> &prim, const Scene& scene, int depth);

public:
    KDTree(/* args */);
    ~KDTree();

    // Construction globale
    void buildKDTree(const Scene& scene);

    bool intersectPrimitive(const PrimitiveRef& prim, const Scene& scene, const Ray& ray, Hit& hit) const;

    bool intersectNode(Node* node, const Scene& scene, const Ray& ray, float t_entry, float t_exit, Hit& hit) const;

};

KDTree::KDTree(/* args */)
{
}

KDTree::~KDTree()
{
}





#endif