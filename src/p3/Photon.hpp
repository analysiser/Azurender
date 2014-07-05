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

#define COLOR_COEEF (0.004)

using namespace _462;

class Photon {
    
public:
    
    // Naive implementation
    // Will be optimized later
    Vector3 position;       // Hit position
//    Vector3 direction;      // Incoming direction
//    Vector3 normal;         // Normal of the intersection point, for not sampling a sphere but a disk
    
//    Color3 color;           // Photon color, also flux
    
//    real_t squaredDistance;
    
    // compressed color
    unsigned char ccolor[4];
    
    unsigned char phi, theta; // compressed incident direction

    char mask;              // Mask to indicate where the photon belongs to
    // 0x1 : Hit on diffusive surface, 0x1 : Hit before, 0x0 : Not hit before
    // 0x2 : Hit on specular surface, 0x2 : Hit before, 0x0 : Not hit before
    char splitAxis;
    
    Photon(){}
    
    Photon(Color3 lcolor) {
        position = Vector3::Zero();
//        direction = Vector3::Zero();
//        normal = Vector3::Zero();
        
        setColor(lcolor);
//        color = lcolor;
        
        mask = 0x0;
//        squaredDistance = 0;
        splitAxis = -1;
    }
    
    void setColor(Color3 color)
    {
        color.to_array(ccolor);
    }
    
    Color3 getColor()
    {
        Color3 color = Color3::Black();
        color.r = (ccolor[0]) * COLOR_COEEF;
        color.g = (ccolor[1]) * COLOR_COEEF;
        color.b = (ccolor[2]) * COLOR_COEEF;
        return color;
    }
    
    // Photon squared distance greater
//    bool operator()(const Photon& a,const Photon& b) const{
//        return (a.squaredDistance < b.squaredDistance);
//    }
    
//    bool operator<(const Photon &p2) const {
//        return squaredDistance < p2.squaredDistance;
//    }
};

struct ClosePhoton
{
    size_t kdtreeIndex;
    real_t x, y, z;
    real_t squaredDistance;
    short splitAxis;
    
    ClosePhoton()
    {
        kdtreeIndex = -1;
        x = y = z =0;
        splitAxis = -1;
    }
    
    ClosePhoton(Photon photon, size_t idx)
    {
        kdtreeIndex = idx;
        x = photon.position.x;
        y = photon.position.y;
        z = photon.position.z;
        squaredDistance = -1;
        splitAxis = photon.splitAxis;
    }
    
    bool operator<(const ClosePhoton &p2) const {
        return squaredDistance < p2.squaredDistance;
    }
};

#endif /* defined(__P3__Photon__) */

