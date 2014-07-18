/**
 * @file model.cpp
 * @brief Model class
 *
 * @author Eric Butler (edbutler)
 * @author Zeyang Li (zeyangl)
 */

#include "scene/model.hpp"
#include "scene/material.hpp"
#include "application/opengl.hpp"
#include "scene/triangle.hpp"
#include <iostream>
#include <cstring>
#include <string>
#include <fstream>
#include <sstream>


namespace _462 {
    
    LinkedList::LinkedList() {
        
        head = NULL;
        tail = NULL;
        
        count = 0;
    }
    
    LinkedList::LinkedList(const LinkedList &other)
    {
        head = NULL;
        tail = NULL;
        count = 0;
        
        BoxNode *tmp = other.head;
        while (tmp != NULL) {
            BoxNode *newNode = create(tmp->boundingBox, tmp->triangleIndex);
            this->append(newNode);
            tmp = tmp->next;
        }
    }
    
    LinkedList::~LinkedList() {
        while (head != NULL) {
            BoxNode *tmp = head;
            head = head->next;
            assert(tmp != NULL);
            tmp->next = NULL;
            delete tmp;
        }
        tail = NULL;
    }
    
    BoxNode* LinkedList::create(Box box, size_t idx)
    {
        BoxNode *tmp = new BoxNode();
        tmp->boundingBox = box;
        tmp->midPoint = (box.bounds[0] + box.bounds[1]) / 2.0;
        tmp->triangleIndex = idx;
        tmp->prev = NULL;
        tmp->next = NULL;
        
        return tmp;
    }
    
    void LinkedList::append(BoxNode *target)
    {
        assert(target != NULL);
        
        if (target == NULL) {
            return;
        }
        
        if (head == NULL) {
            
            head = target;
            tail = target;
        }
        else {
            
            assert(tail != NULL);
            
            tail->next = target;
            target->prev = tail;
            tail = target;
            
        }
        count++;
    }
    
    void LinkedList::detach(BoxNode *target)
    {
        assert(target != NULL);
        
        if (target == NULL) {
            return;
        }
        
        BoxNode *prevNode = target->prev;
        BoxNode *nextNode = target->next;
        
        if (prevNode != NULL) {
            prevNode->next = nextNode;
        }
        if (nextNode != NULL) {
            nextNode->prev = prevNode;
        }
        
        if (head == target) {
            head = nextNode;
        }
        if (tail == target) {
            tail = prevNode;
        }
        
        target->prev = NULL;
        target->next = NULL;
        
        count--;
    }

    BVHNode::BVHNode()
    {
        this->list = NULL;
        this->left = NULL;
        this->right = NULL;
        this->data = NULL;
    }
    
    BVHNode::BVHNode(LinkedList *linkedList, int axis)
    {
//        this->list = new LinkedList(*linkedList);
        this->list = linkedList;
        this->left = NULL;
        this->right = NULL;
        this->data = NULL;
        
        assert(this->list != NULL);
        
        // Set as leaf node
        if (this->list->count == 1) {
            data = this->list->head;
            this->list->detach(data);
            bbox = data->boundingBox;
            
            assert(this->list->count == 0);
            this->list = NULL;
        }
        else if (this->list->count == 2) {
            LinkedList *ll = new LinkedList();
            LinkedList *rl = new LinkedList();
            
            BoxNode *head = this->list->head;
            BoxNode *tail = this->list->tail;
            
            this->list->detach(head);
            this->list->detach(tail);
            
            ll->append(head);
            rl->append(tail);
            
            this->left = new BVHNode(ll, (axis+1)%3);
            this->right = new BVHNode(rl, (axis+1)%3);
            
            delete ll;
            delete rl;
            
            assert(this->list->count == 0);
            this->list = NULL;
            
            combineBoundingBox();
            
        }
        else {
            // Naive subdivision
            LinkedList *ll = new LinkedList();
            LinkedList *rl = new LinkedList();
            
            // Subdivision along a given axis
            subdivideList(ll, rl, axis);
            
            this->left = new BVHNode(ll, (axis+1)%3);
            this->right = new BVHNode(rl, (axis+1)%3);
            
            delete ll;
            delete rl;
            
            assert(this->list->count == 0);
            this->list = NULL;
            
            combineBoundingBox();
        }
    }
    
    BVHNode::~BVHNode()
    {
        delete this->data;
        delete this->list;
        delete this->left;
        delete this->right;
    }
    
    /**
     * @brief   Combine the bouning box of left and right subtree, 
     *          set it to 'This' node's bounding box
     */
    void BVHNode::combineBoundingBox()
    {
        Box box1 = left->bbox;
        Box box2 = right->bbox;
        
        real_t minx = box1.bounds[0].x < box2.bounds[0].x ? box1.bounds[0].x :
box2.bounds[0].x;
        real_t miny = box1.bounds[0].y < box2.bounds[0].y ? box1.bounds[0].y :
box2.bounds[0].y;
        real_t minz = box1.bounds[0].z < box2.bounds[0].z ? box1.bounds[0].z :
box2.bounds[0].z;
        bbox.bounds[0] = Vector3(minx, miny, minz);
        
        real_t maxx = box1.bounds[1].x > box2.bounds[1].x ? box1.bounds[1].x :
box2.bounds[1].x;
        real_t maxy = box1.bounds[1].y > box2.bounds[1].y ? box1.bounds[1].y :
box2.bounds[1].y;
        real_t maxz = box1.bounds[1].z > box2.bounds[1].z ? box1.bounds[1].z :
box2.bounds[1].z;
        bbox.bounds[1] = Vector3(maxx, maxy, maxz);
        
    }
    
    /**
     * @brief   Test if the ray intersects with BVH Tree node, if there is, push
     *          the model triangle meshes' index into a vector for futher testing
     * @param   r           Incoming ray
     * @param   t0          lower limit of t
     * @param   t1          upeer limit of t
     * @param   &idxList    reference of vector, stores the indices that might hit
     * @return  BVHNode     the node that intersects
     */
    void BVHNode::nodeIntersect(Ray r, real_t t0, real_t t1, std::vector<size_t>
&idxList)
    {
        if (this->bbox.intersect(r, t0, t1)) {
            
            // is leaf node
            if (this->data != NULL) {
                idxList.push_back(this->data->triangleIndex);
            }
            // is subtree node
            else {
                if (this->left != NULL) {
                    this->left->nodeIntersect(r, t0, t1, idxList);
                }
                
                if (this->right != NULL) {
                    this->right->nodeIntersect(r, t0, t1, idxList);
                }
            }
        }
    }
    
    // Magic code
    void BVHNode::subdivideList(LinkedList *leftList, LinkedList *rightList, int axis)
    {
        std::set<real_t> distribution;
        BoxNode *tmp = this->list->head;
        while (tmp != NULL) {
            real_t var = 0;
            switch (axis) {
                case 0:
                    var = tmp->midPoint.x;
                    break;
                case 1:
                    var = tmp->midPoint.y;
                    break;
                case 2:
                    var = tmp->midPoint.z;
                    break;
                default:
                    break;
            }
            
            distribution.insert(var);
            tmp = tmp->next;
        }
        
        // find median number
        real_t median = 0;
        size_t mid = distribution.size()/2;
        size_t idx = 0;
        std::set<real_t>::iterator i;
        for (i = distribution.begin(); i != distribution.end(); i++) {
            if (idx == mid) {
                median = *i;
                break;
            }
            ++idx;
        }
        
        // subdivision
        tmp = this->list->head;
        while (tmp != NULL) {
            real_t var = 0;
            switch (axis) {
                case 0:
                    var = tmp->midPoint.x;
                    break;
                case 1:
                    var = tmp->midPoint.y;
                    break;
                case 2:
                    var = tmp->midPoint.z;
                    break;
                default:
                    break;
            }
            if (var <= median) {
                BoxNode *target = tmp;
                tmp = tmp->next;
                
                this->list->detach(target);
                leftList->append(target);
            }
            else {
                BoxNode *target = tmp;
                tmp = tmp->next;
                
                this->list->detach(target);
                rightList->append(target);
            }
        }
        
        // special condition
        // when all nodes are either on left or on right
        if ((leftList->count == 0) && (rightList->count > 0)) {
            BoxNode *target = rightList->tail;
            rightList->detach(target);
            leftList->append(target);
        }
        
        if ((rightList->count == 0) && (leftList->count > 0)) {
            BoxNode *target = leftList->tail;
            leftList->detach(target);
            rightList->append(target);
        }
    }

    Model::Model() : mesh( 0 ), material( 0 ) {
        modelBoxNodes = new LinkedList();
        bBox = Box(Vector3::Zero(), Vector3::Zero());
        type = eModel;
        root = NULL;
        bvhTree = nullptr;
    }
    Model::~Model() {
        delete modelBoxNodes;
        delete root;
    }
    
    void Model::render() const
    {
        if ( !mesh )
            return;
        if ( material )
            material->set_gl_state();
        mesh->render();
        if ( material )
            material->reset_gl_state();
    }
    
    void Model::createBoundingBox() const
    {
        // if already have bounding boxes, then no need to recalculate again
        // reduced about 10 MB for each time a new ray tracing happens
//        if (this->bBox.bounds[0] == Vector3::Zero() && this->bBox.bounds[1] == Vector3::Zero()) {
//            MeshTriangle const *triangles = mesh->get_triangles();
//            
//            for (size_t i = 0; i < mesh->num_triangles(); i++) {
//                
//                Vector3 A = mesh->vertices[triangles[i].vertices[0]].position;
//                Vector3 B = mesh->vertices[triangles[i].vertices[1]].position;
//                Vector3 C = mesh->vertices[triangles[i].vertices[2]].position;
//                
//                Box boundingBox = getBoundingBoxForTriangle(A, B, C);
//                
//                BoxNode *newNode = modelBoxNodes->create(boundingBox, i);
//                modelBoxNodes->append(newNode);
//            }
//        }
        
        if (!bvhTree) {
            
            MeshTriangle const *triangles = mesh->get_triangles();
            bvhTree = new azBVHTree(mesh->num_triangles());

            for (size_t i = 0; i < mesh->num_triangles(); i++) {
                
                Vector3 A = mesh->vertices[triangles[i].vertices[0]].position;
                Vector3 B = mesh->vertices[triangles[i].vertices[1]].position;
                Vector3 C = mesh->vertices[triangles[i].vertices[2]].position;
                
                BndBox bbox(A);
                bbox.include(B);
                bbox.include(C);
                
//                bvhTree->leafNodes_[i] = azBVHTree::azBVNode(bbox, i);
                bvhTree->setLeaf(azBVHTree::azBVNode(bbox, i), i);
            }
            
            bvhTree->buildBVHTree();
        }
    }
    
    bool Model::hit(Ray ray, real_t t0, real_t t1, HitRecord &rec) const
    {
        Ray r = Ray(invMat.transform_point(ray.e), invMat.transform_vector(ray.d));

//        if (!bBox.intersect(r, t0, t1)) {
//            return false;
//        }
        
        // Indices list that stores all the indices of triangles that might hit
//        std::vector<size_t> indexList;
//        root->nodeIntersect(r, t0, t1, indexList);
        
        auto rayTriangleIntersectionTest = [this](Ray rr, real_t tt0, real_t tt1, real_t &tt, INT64 triIndex){
            MeshTriangle const *triangles = mesh->get_triangles();
            
            MeshVertex A = mesh->vertices[triangles[triIndex].vertices[0]];
            MeshVertex B = mesh->vertices[triangles[triIndex].vertices[1]];
            MeshVertex C = mesh->vertices[triangles[triIndex].vertices[2]];
            
            // result.x = beta, result.y = gamma, result.z = t
            Vector3 result = getResultTriangleIntersection(rr, A.position, B.position, C.position);
            
            if (result.z < tt0 || result.z > tt1) {
                return false;
            }
            
            if (result.y < 0 || result.y > 1) {
                return false;
            }
            
            if (result.x < 0 || result.x > 1 - result.y) {
                return false;
            }
            
            tt = result.z;
            return true;
        };

        INT64 idx = -1;
        
        if (bvhTree->getFirstIntersectIndex(r, t0, t1, idx, rayTriangleIntersectionTest)) {
            
            MeshTriangle const *triangles = mesh->get_triangles();
            MeshVertex A = mesh->vertices[triangles[idx].vertices[0]];
            MeshVertex B = mesh->vertices[triangles[idx].vertices[1]];
            MeshVertex C = mesh->vertices[triangles[idx].vertices[2]];
            
            // result.x = beta, result.y = gamma, result.z = t
            Vector3 result = getResultTriangleIntersection(r, A.position, B.position, C.position);
            
            rec.type = eTriangle;
            
            rec.position = ray.e + result.z * ray.d;
            
            real_t beta = result.x;
            real_t gamma = result.y;
            real_t alpha = 1 - beta - gamma;
            
            rec.normal = normalize(alpha * (normMat * A.normal) + beta * (normMat * B.normal) + gamma * (normMat * C.normal));
            
            // For texture mapping adjustment
            A.tex_coord = getAdjustTexCoord(A.tex_coord);
            B.tex_coord = getAdjustTexCoord(B.tex_coord);
            C.tex_coord = getAdjustTexCoord(C.tex_coord);
            
            Vector2 tex_cood_interpolated =
            alpha * A.tex_coord + beta * B.tex_coord + gamma * C.tex_coord;
            
            int width = 0, height = 0;
            if (material) {
                material->get_texture_size(&width, &height);
            }
            
            rec.diffuse = material->diffuse;
            rec.ambient = material->ambient;
            rec.specular = material->specular;
            rec.phong = material->phong;
            
            rec.texture = material->get_texture_pixel(tex_cood_interpolated.x * width, tex_cood_interpolated.y * height);
            
            rec.t = result.z;
            
            rec.refractive_index = material->refractive_index;
            
//            rec.isLight = this->isLight;

            
            return true;
        }
        
        return false;
        
//        if (indexList.size() == 0) {
//            return false;
//        }
        
//        MeshTriangle const *triangles = mesh->get_triangles();
//        real_t tt = INFINITY;
//        
//        for (size_t i = 0; i < indexList.size(); i++) {
//            
//            size_t idx = indexList[i];
//            
//            MeshVertex A = mesh->vertices[triangles[idx].vertices[0]];
//            MeshVertex B = mesh->vertices[triangles[idx].vertices[1]];
//            MeshVertex C = mesh->vertices[triangles[idx].vertices[2]];
//            
//            // result.x = beta, result.y = gamma, result.z = t
//            Vector3 result = getResultTriangleIntersection(r, A.position,
//B.position, C.position);
//            
//            if (result.z < t0 || result.z > t1) {
//                continue;
//            }
//            
//            if (result.y < 0 || result.y > 1) {
//                continue;
//            }
//            
//            if (result.x < 0 || result.x > 1 - result.y) {
//                continue;
//            }
//            
//            if (result.z > tt) {
//                continue;
//            }
//            
//            rec.type = eTriangle;
//            
//            rec.position = ray.e + result.z * ray.d;
//            
//            real_t beta = result.x;
//            real_t gamma = result.y;
//            real_t alpha = 1 - beta - gamma;
//            
//            rec.normal = normalize(alpha * (normMat * A.normal) + beta * (normMat * B.normal) + gamma * (normMat * C.normal));
//
//            // For texture mapping adjustment
//            A.tex_coord = getAdjustTexCoord(A.tex_coord);
//            B.tex_coord = getAdjustTexCoord(B.tex_coord);
//            C.tex_coord = getAdjustTexCoord(C.tex_coord);
//            
//            Vector2 tex_cood_interpolated =
//            alpha * A.tex_coord + beta * B.tex_coord + gamma * C.tex_coord;
//            
//            int width = 0, height = 0;
//            if (material) {
//                material->get_texture_size(&width, &height);
//            }
//            
//            rec.diffuse = material->diffuse;
//            rec.ambient = material->ambient;
//            rec.specular = material->specular;
//            rec.phong = material->phong;
//            
//            rec.texture = material->get_texture_pixel(tex_cood_interpolated.x * width, tex_cood_interpolated.y * height);
//            
//            rec.t = result.z;
//            
//            rec.refractive_index = material->refractive_index;
//            tt = result.z;
//            
////            rec.isLight = this->isLight;
//        }
//        
//        if (tt < INFINITY) {
//            return true;
//        }
//        
//        return false;
    }
    
    /*
    // @brief To create a BVH Tree for model
    void Model::createBVHTree()
    {
        // Attention here, this function is called multiple times but root is only 
        // intiated once
        // this is to avoid memory leak. modelBoxNodes are only initiated once and it 
        // would be empty then
//        if (root == NULL) {
//            root = new BVHNode(modelBoxNodes, 0);
//            this->bBox = root->bbox;
//        }
        
    }*/
} /* _462 */
