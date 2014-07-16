//
//  BndBox.cpp
//  Azurender
//
//  Created by Xiao Li on 7/10/14.
//
//

#include "BndBox.hpp"

namespace _462 {
    
    void BndBox::include(const _462::Vector3 &p)
    {
        pMin.x = std::min(pMin.x, p.x);
        pMin.y = std::min(pMin.y, p.y);
        pMin.z = std::min(pMin.z, p.z);
        
        pMax.x = std::max(pMax.x, p.x);
        pMax.y = std::max(pMax.y, p.y);
        pMax.z = std::max(pMax.z, p.z);
    }
    
    void BndBox::include(const BndBox &bbox)
    {
        include(bbox.pMin);
        include(bbox.pMax);
//        pMin.x = std::min(pMin.x, bbox.pMin.x);
//        pMin.y = std::min(pMin.y, bbox.pMin.y);
//        pMin.z = std::min(pMin.z, bbox.pMin.z);
//        
//        pMax.x = std::max(pMax.x, bbox.pMax.x);
//        pMax.y = std::max(pMax.y, bbox.pMax.y);
//        pMax.z = std::max(pMax.z, bbox.pMax.z);
    }
    
    BndBox BndBox::expand(const _462::BndBox &b, const _462::Vector3 &p)
    {
        BndBox ret = b;
        
        ret.pMin.x = std::min(b.pMin.x, p.x);
        ret.pMin.y = std::min(b.pMin.y, p.y);
        ret.pMin.z = std::min(b.pMin.z, p.z);
        
        ret.pMax.x = std::max(b.pMax.x, p.x);
        ret.pMax.y = std::max(b.pMax.y, p.y);
        ret.pMax.z = std::max(b.pMax.z, p.z);
        
        return ret;
    }
    
    BndBox BndBox::expand(const _462::BndBox &b, const _462::BndBox &b2)
    {
        BndBox ret = b;
        
        ret.pMin.x = std::min(b.pMin.x, b2.pMin.x);
        ret.pMin.y = std::min(b.pMin.y, b2.pMin.y);
        ret.pMin.z = std::min(b.pMin.z, b2.pMin.z);
        
        ret.pMax.x = std::max(b.pMax.x, b2.pMax.x);
        ret.pMax.y = std::max(b.pMax.y, b2.pMax.y);
        ret.pMax.z = std::max(b.pMax.z, b2.pMax.z);
        
        return ret;
    }
    
    bool BndBox::intersect(const Ray &r, real_t t0, real_t t1) const
    {
        real_t mint = t0, maxt = t1;
        
        for (int i = 0; i < 3; i++) {
            float invRayDir = 1.f/r.d[i];
            float tNear = (pMin[i] - r.e[i]) * invRayDir;
            float tFar  = (pMax[i] - r.e[i]) * invRayDir;
            
            if (tNear > tFar) std::swap(tNear, tFar);
            mint = tNear > mint ? tNear : mint;
            maxt = tFar  < maxt ? tFar  : maxt;
            
            if (mint > maxt) return false;
        }
        
        return true;
    }
    
}
