//
//  montecarlo.h
//  Azurender
//
//  Created by Xiao Li on 9/7/14.
//
//

#include "montecarlo.hpp"

namespace _462 {
    
    void ConcentricSampleDisk(float u1, float u2, float *dx, float *dy) {
        float r, theta;
        // Map uniform random numbers to $[-1,1]^2$
        float sx = 2 * u1 - 1;
        float sy = 2 * u2 - 1;
        
        // Map square to $(r,\theta)$
        
        // Handle degeneracy at the origin
        if (sx == 0.0 && sy == 0.0) {
            *dx = 0.0;
            *dy = 0.0;
            return;
        }
        if (sx >= -sy) {
            if (sx > sy) {
                // Handle first region of disk
                r = sx;
                if (sy > 0.0) theta = sy/r;
                else          theta = 8.0f + sy/r;
            }
            else {
                // Handle second region of disk
                r = sy;
                theta = 2.0f - sx/r;
            }
        }
        else {
            if (sx <= sy) {
                // Handle third region of disk
                r = -sx;
                theta = 4.0f - sy/r;
            }
            else {
                // Handle fourth region of disk
                r = -sy;
                theta = 6.0f + sx/r;
            }
        }
        theta *= M_PI / 4.f;
        *dx = r * cosf(theta);
        *dy = r * sinf(theta);
    }
    
    Vector3 UniformSampleHemisphere(float u1, float u2)
    {
        float z = u1;
        float r = sqrtf(std::max(0.f, 1.f - z * z));
        float phi = 2 * M_PI * u2;
        float x = r * cosf(phi);
        float y = r * sinf(phi);
        return Vector3(x, y, z);
    }
    
}


