//
//  azBVHTree.hpp
//  Azurender
//
//  Created by Xiao Li on 7/11/14.
//
//

#ifndef __Azurender__azBVHTree__
#define __Azurender__azBVHTree__

#include <iostream>

#include "math/matrix.hpp"
#include "math/quaternion.hpp"
#include "math/vector.hpp"

#include "scene.hpp"
#include "scene/model.hpp"
#include "scene/BndBox.hpp"
#include "scene/ray.hpp"

namespace _462 {
    
    typedef std::int64_t INT64;
    typedef std::uint32_t UINT32;
    typedef std::uint8_t  UINT8;
    
    class azBVHTree
    {
    public:
        
        class azBVNode;
        typedef std::vector<azBVNode> azBVNodesArray;
        
        azBVHTree () {}
        
        azBVHTree (const UINT32 &leafSize) {
            
            leafsize_ = leafSize;
            branchsize_ = leafSize - 1;
            size_ = leafsize_ + branchsize_;
            
            leafNodes_ = azBVNodesArray(leafsize_);
            branchNodes_ = azBVNodesArray(branchsize_);
            root_ = &branchNodes_[0];
        }
        
        azBVNodesArray getLeafNodes() { return leafNodes_; }
        UINT32 getLeafSize() { return leafsize_; }
        
        azBVNodesArray getBranchNodes() { return branchNodes_; }
        UINT32 getBranchSize() { return branchsize_; }
        
        void setLeaf(BndBox bbox, UINT32 index) {
            assert(index < leafsize_);
            leafNodes_[index] = azBVNode(bbox, index);
        }
        
        azBVNode *root() {
            
            assert(branchsize_ > 0);
            return &branchNodes_[0];
        }
        
        void buildBVHTree() {
            
            assert(root_ != nullptr);
            assert(leafsize_ > 0);
            assert(branchsize_ > 0);
            
            root_ = &*branchNodes_.begin();
            root_->buildDown(leafNodes_.begin(), leafNodes_.end());
        }
        
        template <typename FN>
        bool getFirstIntersectIndex(const Ray& r,
                                    real_t& t0,
                                    real_t& t1,
                                    INT64& index,
                                    FN &func) {
            INT64 idx = -1;
            this->root()->intersectRayTest(r, t0, t1, idx, func);
            if (idx != -1) {
                index = idx;
                return true;
            }
            return false;
        }
        
        
        // Nested class azBVNode, supporting structure for BVHTree
        class azBVNode : public BndBox {

        public:
            /**
             * The number of dimensions.
             */
            static const size_t DIM = 3;
            
            azBVNode() : idx1_(0),
            idx2_(std::numeric_limits<unsigned int>::max()),
            edgeValue_(0), isLeaf_(0), leftChild_(0), rightChild_(0)
            { }
            
            azBVNode(const BndBox &bbox, UINT32 index) {
                
                this->d_ = 0;
                this->pMin = bbox.pMin;
                this->pMax = bbox.pMax;
                this->idx1_ = 0;
                this->idx2_ = index;
                this->isLeaf_ = true;
                this->leftChild_ = nullptr;
                this->rightChild_ = nullptr;
            }
            
            azBVNode(const azBVNode *other)
            {
                this->d_ = other->d_;
                this->pMin = other->pMin;
                this->pMax = other->pMax;
                this->idx1_ = other->idx1_;
                this->idx2_ = other->idx2_;
                this->isLeaf_ = other->isLeaf_;
                this->leftChild_ = other->leftChild_;
                this->rightChild_ = other->rightChild_;
                
            }
            
            UINT8 getLongestEdge() {
                UINT8 d = 0;
                real_t longest = getEdgeValue(d);
                for (UINT8 i = 1; i < DIM; i++) {
                    real_t edgeValue_i = getEdgeValue(i);
                    if (edgeValue_i > longest) {
                        longest = edgeValue_i;
                        d = i;
                    }
                }
                return d;
            }
            
            real_t getEdgeValue(UINT8 d) {
                return (pMax[d] - pMin[d]);
            }
            
            UINT32 getLeafIndex() {
                assert(idx1_ == 0);
                return idx2_;
            }
            
            bool isLeaf(){ return isLeaf_; }
            
            // Build the bounding box for all nodes from leavesBegin to leavesEnd
            void buildBoundingBox(std::vector<azBVNode>::iterator leavesBegin,
                                  std::vector<azBVNode>::iterator leavesEnd);
            
            // Build down
            void buildDown(std::vector<azBVNode>::iterator leavesBegin,
                           std::vector<azBVNode>::iterator leavesEnd);
            
            // Ray intersect with box
            template <typename FN>
            void intersectRayTest(const Ray& r,
                                  real_t& t0,
                                  real_t& t1,
                                  int64_t& index,
                                  FN &func) {
                
                if (this->intersect(r, t0, t1)) {
                    if (this->isLeaf()) {
                        real_t tt;
                        if (func(r, t0, t1, tt, this->idx2_)) {
                            t1 = tt;
                            index = this->idx2_;
                        }

                    }
                    else {
                        if (this->leftChild_ != nullptr) {
                            this->leftChild_->intersectRayTest(r, t0, t1, index, func);
                        }
                        if (this->rightChild_ != nullptr) {
                            this->rightChild_->intersectRayTest(r, t0, t1, index, func);
                        }
                        
                        
                    }
                    
                }
                
                
                
            }
            
//        private:
//            UINT8 edge_;
            UINT8 d_;
            UINT32 idx1_, idx2_;
            real_t edgeValue_;
            bool isLeaf_;
            
            azBVNode *leftChild_, *rightChild_;
        
        };
        
    private:
        UINT32 size_, leafsize_, branchsize_;
        azBVNode *root_;
        azBVNodesArray leafNodes_;
        azBVNodesArray branchNodes_;
        
        
    };
    
    
    
}


#endif /* defined(__Azurender__azBVHTree__) */
