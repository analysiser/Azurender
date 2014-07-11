//
//  BndBox.h
//  Azurender
//
//  Created by Xiao Li on 7/10/14.
//
//

#ifndef __Azurender__BndBox__
#define __Azurender__BndBox__

#include <iostream>

#include "math/matrix.hpp"
#include "math/quaternion.hpp"
#include "math/vector.hpp"

#include "scene/ray.hpp"

namespace _462 {
    
    class BndBox
    {
    public:
        BndBox () {
            pMin = Vector3(INFINITY, INFINITY, INFINITY);
            pMax = Vector3(-INFINITY, -INFINITY, -INFINITY);
        }
        
        BndBox (const Vector3 &p) : pMin(p), pMax(p) { }
        BndBox (const Vector3 &p1, const Vector3 &p2) {
            pMin = Vector3(std::min(p1.x, p2.x), std::min(p1.y, p2.y), std::min(p1.z, p2.z));
            pMax = Vector3(std::max(p1.x, p2.x), std::max(p1.y, p2.y), std::max(p1.z, p2.z));
        }
        
        bool isInside(const Vector3 &pt) const {
            return (pt.x >= pMin.x && pt.x <= pMax.x &&
                    pt.y >= pMin.y && pt.y <= pMax.y &&
                    pt.z >= pMin.z && pt.z <= pMax.z);
        }
        
        float surfaceArea() const {
            Vector3 d = pMax - pMin;
            return 2.f * (d.x * d.y + d.x * d.z + d.y * d.z);
        }
        
        float volume() const {
            Vector3 d = pMax - pMin;
            return d.x * d.y * d.z;
        }
        
        void include(const Vector3 &p);
        BndBox expand(const BndBox &b, const Vector3 &p);
        BndBox expand(const BndBox &b, const BndBox &b2);
        
        // Test if a given ray hits a bounding box
        bool intersect(const Ray &r, real_t t0, real_t t1) const;
        
        Vector3 pMin, pMax;
    };
}



#endif /* defined(__Azurender__BndBox__) */
