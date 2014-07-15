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
    
    typedef std::uint32_t UINT;
    typedef std::uint8_t  UINT8;
    typedef bool (*rayIntersectionFunction)(const Ray& r, real_t t0, real_t t1, real_t& tt, UINT index);
    
    class azBVHTree
    {
    public:
        
        class azBVNode;
        typedef std::vector<azBVNode> azBVNodesArray;
//        typedef std::vector<azBVNode> azBVNodesList;
        
        azBVHTree () : size_(0),
        leafsize_(0),
        branchsize_(0),
        branchNodes_(0),
        leafNodes_(0)
        { }
        
        azBVHTree (const UINT &leafSize) {
            
            leafsize_ = leafSize;
            branchsize_ = leafSize - 1;
            size_ = leafsize_ + branchsize_;
            
            leafNodes_ = azBVNodesArray(leafsize_);
            branchNodes_ = azBVNodesArray(branchsize_);
        }
        
//        template <class FN>
//        bool insectFirstRay(const Ray& ray, real_t& t0, real_t &t1, FN& rayIntersectionFunction);
        
        
        
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
            
            azBVNode(const BndBox &bbox, UINT index) {
                
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
            
            UINT getLeafIndex() {
                assert(idx1_ == 0);
                return idx2_;
            }
            
            bool isLeaf(){ return isLeaf_; }
            
            void buildBoundingBox(std::vector<azBVNode>::iterator leavesBegin,
                                  std::vector<azBVNode>::iterator leavesEnd);
            
            // Build down
            void buildDown(std::vector<azBVNode>::iterator leavesBegin,
                           std::vector<azBVNode>::iterator leavesEnd);
            
            // Ray intersect with box
//            template <class FN>
            void intersectRay(const Ray& r,
                              real_t& t0,
                              real_t& t1,
                              UINT &index,
                              std::vector<size_t> &indices);
            
//        private:
//            UINT8 edge_;
            UINT8 d_;
            UINT idx1_, idx2_;
            real_t edgeValue_;
            bool isLeaf_;
            
            azBVNode *leftChild_, *rightChild_;
        };
        
//    private:
        UINT size_, leafsize_, branchsize_;
        azBVNode *root_;
        azBVNodesArray leafNodes_;
        azBVNodesArray branchNodes_;
        
    private:
        
        
        
    };
    
    
    
}


#endif /* defined(__Azurender__azBVHTree__) */
