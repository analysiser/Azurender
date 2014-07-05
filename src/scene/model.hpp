/**
 * @file model.hpp
 * @brief Model class
 *
 * @author Eric Butler (edbutler)
 */

#ifndef _462_SCENE_MODEL_HPP_
#define _462_SCENE_MODEL_HPP_

#include "scene/scene.hpp"
#include "scene/mesh.hpp"
#include <iostream>
#include <set>
#include <vector>
#include <map>
#include <list>

namespace _462 {

    
    /**
     * @brief   BoxNode is for storing the bounding volumn of primitives
     *          It also acts as temporary linked list node that would be converted
     *          into bvh tree after initializing model meshs.
     */
    struct BoxNode
    {
        Box boundingBox;        // bounding box of the node
        size_t triangleIndex;   // mesh triangle index it refers
        Vector3 midPoint;       // bounding box mid point, used for tree subdivision

        
        BoxNode *prev;
        BoxNode *next;

        
        BoxNode(){
            boundingBox = Box(Vector3::Zero(), Vector3::Zero());
            midPoint = Vector3::Zero();
            prev = NULL;
            next = NULL;
        }
    };

    
    /**
     * @brief   LinkedList that stores BoxNode. It was used for intersection test, 
     *          then I reorganized it into a BVH tree
     */
    class LinkedList
    {
    public:
        BoxNode *head;
        BoxNode *tail;

        
        LinkedList();
        LinkedList(const LinkedList &other);
        ~LinkedList();

        
        // Create a BoxNode with box and triangle mesh index
        BoxNode *create(Box box, size_t idx);

        
        // Append a BoxNode, append it to the tail of the linkedlist
        void append(BoxNode *target);

        
        // Detach a BoxNode, but not release the memory, memory releasing would
        // be handled by BVHNode
        void detach(BoxNode *target);

        
        // size of the linked list
        size_t count;
    };

    
    /**
     * @brief   BVH Node on BVH Tree (Bounding Volumn Hierarchical node)
     *          It reads in a linked list of BoxNodes, make the list organized
     *          as BVH tree. Only the leaf node of BVH tree stores the actual 
     *          data, while other nodes has the bounding box only for intersection
     *          tests.
     */
    class BVHNode
    {
    public:
        BVHNode();
        BVHNode(LinkedList *linkedList, int axis);
        ~BVHNode();

        
        LinkedList *list;           // Linked list for BoxNode
        Box bbox;                   // Bounding box of the BVH node
        BoxNode *data;              // BoxNode data that actually stored in BVHNode

        
        BVHNode *left;              // Pointer for left sub tree
        BVHNode *right;             // Pointer for right sub tree

        
        // Get the bounding box for current BVHNode, called AFTER initializing left and right subtree
        void combineBoundingBox();

        
        // Test ray - BVH node intersection, put intersect's triangle indices to a vector
        void nodeIntersect(Ray r, real_t t0, real_t t1, std::vector<size_t> &idxList);

        
        // Subdivide BVH tree node according to given axis, applied 3D kd-tree subdivision algorithm
        // No SAH, used median
        void subdivideList(LinkedList *leftList, LinkedList *rightList, int axis);

    };

    

    
    /**
     * A mesh of triangles.
     */
    class Model : public Geometry
    {
    public:

        
        const Mesh* mesh;
        const Material* material;

        
        LinkedList *modelBoxNodes;
        Box bBox;

        
        BVHNode *root;          // this should be used as const...

        
        Model();

        
        virtual ~Model();

        
        virtual void render() const;

        
        //
        virtual void createBoundingBox() const;

        
        // Override of virtual function from Geometry
        virtual bool hit(Ray ray, real_t t0, real_t t1, HitRecord &rec) const;

        
        // Create a BVH Tree for optimization
        void createBVHTree();

                
    };
} /* _462 */

#endif /* _462_SCENE_MODEL_HPP_ */
