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
    
    // source: https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
    
    class azReflection {
        
    public:

        azReflection() {}
        ~azReflection() {}        
        
        static Vector3 reflect(Vector3 i, Vector3 n)
        {
            return i - 2.0 * n * dot(n, i);
        }
        
//        static Vector3 refract( Vector3 i, Vector3 n, float inv_eta )
//        {
//            float cosi = dot(-i, n);
//            float cost2 = 1.0f - inv_eta * inv_eta * (1.0f - cosi * cosi);
//            Vector3 t = inv_eta * i + ((inv_eta * cosi - sqrt(abs(cost2))) * n);
//            //            return t * (Vector3)(cost2 > 0);
//            return t * (cost2 > 0 ? 1.0f : 0.0f);
//        }
        
        static Vector3 refract(Vector3 d, Vector3 n, real_t ratio)
        {
            real_t dot_dn = dot(d, n);
            real_t square = real_t(1) - (pow(ratio, 2) * (real_t(1) - pow(dot_dn, 2)));
            
            if (square < 0) {
                return Vector3().Zero();
            }
            
            return (ratio * (d - (n * dot_dn))) - (n * sqrt(square));
        }
    };
    
    class azFresnel {
        
    public:
        
        azFresnel() {}
        ~azFresnel() {}
        
        
        // Xiao 8/4/2014
        static float FresnelDielectricDielectric(float Eta, float CosTheta)
        {
            float c = CosTheta;
            float temp = Eta * Eta + c * c - 1;
            
            if (temp < 0)
                return 1;
            
            float g = sqrt(temp);
            return 0.5 * Square((g - c) / (g + c)) *
            (1 + Square(( (g + c)  * c - 1) / ((g - c) * c + 1)));
        }
        
//        static float FresnelDielectricDielectric(float Eta, float CosTheta)
//        {
//            float SinTheta2 = 1 - CosTheta * CosTheta;
//            
//            float t0 = sqrt(1 - (SinTheta2 / (Eta * Eta)));
//            float t1 = Eta * t0;
//            float t2 = Eta * CosTheta;
//            
//            float rs = (CosTheta - t1) / (CosTheta + t1);
//            float rp = (t0 - t2) / (t0 + t2);
//            
//            return 0.5 * (rs * rs + rp * rp);
//        }
        
    };
}


#endif /* defined(__Azurender__azReflection__) */