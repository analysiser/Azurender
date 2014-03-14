//
//  Photon.h
//  P3
//
//  Created by Xiao Li on 3/8/14.
//
//

#ifndef __P3__Photon__
#define __P3__Photon__

#include <iostream>
#include "math/color.hpp"
#include "math/vector.hpp"

using namespace _462;

enum ePhotonType {
    ePhotonTypeNone = 0,
    ePhotonTypeDirect = 1,
    ePhotonTypeIndirect = 2,
    ePhotonTypeCaustics = 3,
    };

class Photon {
    
public:
    
    // Naive implementation
    // Will be optimized later
    Vector3 position;       // Hit position
    Vector3 direction;      // Incoming direction
    Vector3 normal;         // Normal of the intersection point, for not sampling a sphere but a disk
    
    Color3 color;           // Photon color, also flux
    
    real_t squaredDistance;
    
    ePhotonType photonType; // Photon type, update when photon tracing
    
    char mask;              // Mask to indicate where the photon belongs to
    // 0x1 : Hit on diffusive surface, 0x1 : Hit before, 0x0 : Not hit before
    // 0x2 : Hit on specular surface, 0x2 : Hit before, 0x0 : Not hit before
    
    Photon(){}
    
    Photon(Color3 lcolor) {
        position = Vector3::Zero();
        direction = Vector3::Zero();
        normal = Vector3::Zero();
        
        color = lcolor;
        
        photonType = ePhotonTypeNone;
        mask = 0x0;
        squaredDistance = 0;
    }
    
    // Photon squared distance greater
    bool operator()(const Photon& a,const Photon& b) const{
        return (a.squaredDistance < b.squaredDistance);
    }
    
    bool operator<(const Photon &p2) const {
        return squaredDistance < p2.squaredDistance;
    }
    
//    bool operator==(const Photon &p2) const {
//        return (photon == p2.photon) && (squaredDistance == p2.squaredDistance);
//    }
};

#endif /* defined(__P3__Photon__) */

