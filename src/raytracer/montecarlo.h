//
//  montecarlo.h
//  Azurender
//
//  Created by Xiao Li on 9/7/14.
//
//

#include "math/math.hpp"
#include "math/vector.hpp"

namespace _462 {
    
    void ConcentricSampleDisk(float u1, float u2, float *dx, float *dy);
    
    inline Vector3 CosineSampleHemisphere(float u1, float u2) {
        Vector3 ret;
        ConcentricSampleDisk(u1, u2, &ret.x, &ret.y);
        ret.z = sqrtf(max(0.f, 1.f - ret.x*ret.x - ret.y*ret.y));
        return ret;
    }
}