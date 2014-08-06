//
//  azReflection.cpp
//  Azurender
//
//  Created by Xiao Li on 8/4/14.
//
//

#ifndef __Azurender__azReflection__
#define __Azurender__azReflection__



#include "math/color.hpp"
#include "math/math.hpp"
#include "math/vector.hpp"

namespace _462 {
    
    class azReflection {
        
    public:

        azReflection() {}
        ~azReflection() {}        
        
        static Vector3 reflect(Vector3 i, Vector3 n)
        {
            return i - 2.0 * n * dot(n, i);
        }
        
        // tempoaryt
        static Vector3 refract(Vector3 d, Vector3 n, real_t ratio)
        {
            real_t dot_dn = dot(d, n);
            real_t square = real_t(1) - (pow(ratio, 2) * (real_t(1) - pow(dot_dn, 2)));
            
            if (square < 0) {
                return Vector3().Zero();
            }
            
            return (ratio * (d - (n * dot_dn))) - (n * sqrt(square));
        }
        
        // Incoming direction d
        static Vector3 refract(Vector3 d, Vector3 n, const float &etai, const float &etat)
        {
            float cosi = dot(d, n);
            cosi = clamp(cosi, -1.f, 1.f);
            
            // Compute indices of refraction for dielectric
            // From air to dieletric
            bool entering = cosi < 0.;
            float ei = etai, et = etat;
            
            // From dieletric to air
            if (!entering)
            {
                std::swap(ei, et);
                n = -n;
                cosi = -cosi;
            }
            
            float square = 1.f - (pow(ei/et, 2.f) * (1.f - pow(cosi, 2.f)));
            
            if (square < 0) {
                return Vector3().Zero();
            }
            
            return (ei/et * (d - (n * cosi))) - (n * sqrt(square));
        }
    };
    
    class azFresnel {
        
    public:
        
        azFresnel() {}
        ~azFresnel() {}
        
        // PBRT version
        static float FrDiel(float cosi, float cost, const float &etai, const float &etat)
        {
//            Color3 cEtai = Color3(etai, etai, etai);
//            Color3 cEtat = Color3(etat, etat, etat);
            
            float Rparl = ((etat * cosi) - (etai * cost)) /
                          ((etat * cosi) + (etai * cost));
            float Rperp = ((etai * cosi) - (etat * cost)) /
                          ((etai * cosi) + (etat * cost));
            
            return (Rparl * Rparl + Rperp * Rperp) / 2.f;
        }
        
        // PBRT version
        static float FresnelDielectricEvaluate(float cosi, const float &etai, const float &etat)
        {
            // Compute Fresnel reflectance for dielectric
            cosi = clamp(cosi, -1.f, 1.f);
            
            // Compute indices of refraction for dielectric
            bool entering = cosi < 0.;
            float ei = etai, et = etat;
            
            if (!entering)  std::swap(ei, et);
            
            // Compute _sint_ using Snell's law
            float sint = ei/et * sqrtf(std::max(0.f, 1.f - cosi * cosi));
            if (sint >= 1.) {
                // Handle total internal reflection
                return 1.;
            }
            else {
                float cost = sqrtf(std::max(0.f, 1.f - sint*sint));
                return FrDiel(fabsf(cosi), cost, ei, et);
            }
        }
        
        // source: https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
        static float FresnelDielectricDielectric(float Eta, float CosTheta)
        {
            float SinTheta2 = 1 - CosTheta * CosTheta;
            
            float t0 = sqrt(1 - (SinTheta2 / (Eta * Eta)));
            float t1 = Eta * t0;
            float t2 = Eta * CosTheta;
            
            float rs = (CosTheta - t1) / (CosTheta + t1);
            float rp = (t0 - t2) / (t0 + t2);
            
            return 0.5 * (rs * rs + rp * rp);
        }
        
    };
}


#endif /* defined(__Azurender__azReflection__) */