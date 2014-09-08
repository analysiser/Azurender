//
//  montecarlo.h
//  Azurender
//
//  Created by Xiao Li on 9/7/14.
//
//

#include <iostream>

#include "math/math.hpp"
#include "math/vector.hpp"



namespace _462 {
    
    void ConcentricSampleDisk(float u1, float u2, float *dx, float *dy);
    Vector3 UniformSampleHemisphere(float u1, float u2);
    
    inline Vector3 CosineSampleHemisphere(float u1, float u2) {
        
        float x, y, z;
        ConcentricSampleDisk(u1, u2, &x, &y);
        z = std::sqrtf(std::max(0.f, 1.f - x * x - y * y));
        return Vector3(x, y, z);
    }
}