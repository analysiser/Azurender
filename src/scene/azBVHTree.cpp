//
//  azBVHTree.cpp
//  Azurender
//
//  Created by Xiao Li on 7/11/14.
//
//

#include "azBVHTree.hpp"

namespace _462 {
    
//    template <class FN>
//    bool azBVHTree::insectFirstRay(const Ray& ray, real_t& t0, real_t &t1, FN& rayIntersectionFunction) {
//        
//    }
    
    void azBVHTree::azBVNode::buildBoundingBox(std::vector<azBVNode>::iterator leavesBegin,
                                               std::vector<azBVNode>::iterator leavesEnd) {
        
        std::vector<azBVNode>::iterator it;
        for (it = leavesBegin; it != leavesEnd; it++) {
//            const azBVNode node = *it;
            this->include(*it);
        }
        
//        UINT size = leavesEnd - leavesBegin;
//        for (UINT i = 0; i < size; i++) {
//            const azBVNode node = *(leavesBegin + i);
//            this->include(node);
//        }
    }
    
    void azBVHTree::azBVNode::buildDown(std::vector<azBVNode>::iterator leavesBegin,
                                        std::vector<azBVNode>::iterator leavesEnd) {
        
        assert(leavesBegin != leavesEnd);
        
        UINT size = leavesEnd - leavesBegin;
        std::cout<<size<<std::endl;
        
        assert(size > 1);
        
        // this node is a leaf node
//        if (size > 1) {
//            this->isLeaf_ = false;
//            this->buildBoundingBox(leavesBegin, leavesEnd);
//            this->d_ = this->getLongestEdge();
//            
//            // get the middle position along the longest edge of this bounding box
//            std::vector<azBVNode>::iterator leavesMiddle =
//            std::partition(leavesBegin, leavesEnd, [this](const azBVNode &aNode) {
//                return ( (aNode.pMin[this->d_] + aNode.pMax[this->d_]) <= (this->pMin[this->d_] + this->pMax[this->d_]) );
//            });
//            
//            
//            
//        }
        
        
        // this is the node that is directly above two leaves
        if (size == 2) {
            this->isLeaf_ = false;
            this->buildBoundingBox(leavesBegin, leavesEnd);
            this->leftChild_ = &(*leavesBegin);
            this->rightChild_ = &(*(leavesBegin + 1));
        }
        else if (size > 2) {
            
            this->isLeaf_ = false;
            
            this->buildBoundingBox(leavesBegin, leavesEnd);
            d_ = this->getLongestEdge();
            
            // get the middle position along the longest edge of this bounding box
            std::vector<azBVNode>::iterator leavesMiddle =
            std::partition(leavesBegin, leavesEnd, [this](const azBVNode &aNode) {
                return ( (aNode.pMin[this->d_] + aNode.pMax[this->d_]) < (this->pMin[this->d_] + this->pMax[this->d_]) );
            });
            
            UINT leftSize = leavesMiddle - leavesBegin;
            UINT rightSize = leavesEnd - leavesMiddle;
            
//            std::cout<<leftSize<<" "<<rightSize<<std::endl;
            
            if ((leftSize == 0) or (rightSize == 0)) {
                
                UINT median = (leftSize + rightSize - 1)/2 + 1;
                std::nth_element(leavesBegin, leavesBegin + median, leavesEnd,
                                 [this](const azBVNode &nodeA, const azBVNode &nodeB) {
                                     
                            return ( (nodeA.pMin[this->d_] + nodeA.pMax[this->d_]) <=
                                     (nodeB.pMin[this->d_] + nodeB.pMax[this->d_]) );
                });
                
                leavesMiddle = leavesBegin + median;
                leftSize = leavesMiddle - leavesBegin;
                rightSize = leavesEnd - leavesMiddle;
                
//                std::cout<<"adjusted: "<<leftSize<<" "<<rightSize<<std::endl;
                assert((leftSize != 0) and (rightSize != 0));
            }
            
            if ((leftSize > 0) and (rightSize > 0)) {
                
                // always make sure mode nodes on left side
                if (leftSize >= rightSize) {
                    
                    // if there is only one node on right
                    if (rightSize == 1) {
                        this->leftChild_ = this + 1;
                        this->rightChild_ = &(*leavesMiddle);
                        
                        this->leftChild_->buildDown(leavesBegin, leavesMiddle);
                    }
                    else {
                        this->leftChild_ = this + 1;
                        this->rightChild_ = this + leftSize;
                        
                        this->leftChild_->buildDown(leavesBegin, leavesMiddle);
                        this->rightChild_->buildDown(leavesMiddle, leavesEnd);
                    }
                }
                else {
                    
                    // if there is only one node the left
                    if (leftSize == 1) {
                        this->leftChild_ = this + 1;
                        this->rightChild_ = &(*leavesMiddle);
                        
                        this->leftChild_->buildDown(leavesMiddle, leavesEnd);
                    }
                    else {
                        this->leftChild_ = this + 1;
                        this->rightChild_ = this + rightSize;
                        
                        this->rightChild_->buildDown(leavesBegin, leavesMiddle);
                        this->leftChild_->buildDown(leavesMiddle, leavesEnd);
                    }
                }
            }
            
        }
        
        
        
    }
    
//    template <class FN>
    void azBVHTree::azBVNode::intersectRay(const Ray& r,
                                           real_t& t0,
                                           real_t& t1,
                                           UINT &index,
                                           std::vector<size_t> &indices)
    {
        
        if (this->intersect(r, t0, t1)) {
            
            if (this->isLeaf()) {
                indices.push_back(this->idx2_);
            }
            else {
                if (this->leftChild_ != nullptr) {
                    this->leftChild_->intersectRay(r, t0, t1, index, indices);
                }
                
                if (this->rightChild_ != nullptr) {
                    this->rightChild_->intersectRay(r, t0, t1, index, indices);
                }
            }
        }
    }
    
    
}