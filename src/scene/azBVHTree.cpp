//
//  azBVHTree.cpp
//  Azurender
//
//  Created by Xiao Li on 7/11/14.
//
//

#include "azBVHTree.hpp"

namespace _462 {
    
    void azBVHTree::azBVNode::buildBoundingBox(std::vector<azBVNode>::iterator leavesBegin,
                                               std::vector<azBVNode>::iterator leavesEnd) {
        
        for (auto it = leavesBegin; it != leavesEnd; it++) {
            this->include(it->pMin);
            this->include(it->pMax);
        }
    }
    
    void azBVHTree::azBVNode::buildDown(std::vector<azBVNode>::iterator leavesBegin,
                                        std::vector<azBVNode>::iterator leavesEnd) {
        
        assert(leavesBegin != leavesEnd);
        
        UINT size = leavesEnd - leavesBegin;
        assert(size > 1);
        
        // this node is a leaf node
        if (size > 1) {
            this->isLeaf_ = false;
            this->buildBoundingBox(leavesBegin, leavesEnd);
            this->d_ = this->getLongestEdge();
            
            // get the middle position along the longest edge of this bounding box
            std::vector<azBVNode>::iterator leavesMiddle =
            std::partition(leavesBegin, leavesEnd, [this](const azBVNode &aNode) {
                return ( (aNode.pMin[this->d_] + aNode.pMax[this->d_]) < (this->pMin[this->d_] + this->pMax[this->d_]) );
            });
            
            
            UINT firstSize, secondSize;
            firstSize = leavesMiddle - leavesBegin;
            secondSize = leavesEnd - leavesMiddle;
            
            assert( size == (firstSize + secondSize) );
            
            if (firstSize == 0 or secondSize == 0) {
                
                UINT median = (firstSize + secondSize - 1)/2 + 1;
                std::nth_element(leavesBegin, leavesBegin + median, leavesEnd,
                                 [this](const azBVNode &nodeA, const azBVNode &nodeB) {
                                     
                                     return ( (nodeA.pMin[this->d_] + nodeA.pMax[this->d_]) <
                                              (nodeB.pMin[this->d_] + nodeB.pMax[this->d_]) );
                                 });
                
                leavesMiddle = leavesBegin + median;
                firstSize = leavesMiddle - leavesBegin;
                secondSize = leavesEnd - leavesMiddle;
            }
            
            assert( size == (firstSize + secondSize) );
            
            if (firstSize > 1) {
                
                if (secondSize > 1) {
                    
                    if (firstSize >= secondSize) {
                        this->leftChild_ = this + 1;
                        this->rightChild_ = this + firstSize;
                        
                        this->leftChild_->buildDown(leavesBegin, leavesMiddle);
                        this->rightChild_->buildDown(leavesMiddle, leavesEnd);
                    }
                    else {
                        this->leftChild_ = this + 1;
                        this->rightChild_ = this + secondSize;
                        
                        this->leftChild_->buildDown(leavesMiddle, leavesEnd);
                        this->rightChild_->buildDown(leavesBegin, leavesMiddle);
                    }
                }
                else if (secondSize == 1) {
                    this->leftChild_ = this + 1;
                    this->rightChild_ = &*leavesMiddle;
                    this->leftChild_->buildDown(leavesBegin, leavesMiddle);
                    
                }
                else {
                    assert(0);
                }
            }
            else {
                
                if (firstSize == 1) {
                    if (secondSize > firstSize) {
                        this->leftChild_ = this + 1;
                        this->rightChild_ = &*leavesBegin;
                        this->leftChild_->buildDown(leavesMiddle, leavesEnd);
                    }
                    else if (secondSize == firstSize){
                        this->leftChild_ = &*leavesBegin;
                        this->rightChild_ = &*leavesMiddle;
                    }
                    else {
                        assert(0);
                    }
                }
                else {
                    assert(0);
                }
            }
            
            
        }
        
    }
    
    
}