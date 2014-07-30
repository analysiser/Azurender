//
//  azLights.cpp
//  Azurender
//
//  Created by Xiao Li on 7/17/14.
//
//

#include "azLights.hpp"

#include "math/random462.hpp"



namespace _462 {
    
    Light::Light() : nSamples(1)
    {
        
    }
    
    Light::~Light() { }
    
    bool Light::initialize()
    {
        return true;
    }
    
    PointLight::PointLight()
    {
        position = Vector3::Zero();
        color = Color3::Black();
        intensity = 1.0;
        radius = 0.0;
    }
    
    
    // p and light in world position, need xform in future
    Color3 PointLight::SampleLight(const Vector3 &p, const Vector3 &normal, float t0, float t1) const
    {
        Vector3 dir;
        Vector3 surfP = SamplePointOnLight(dir);//position + ran * radius;
        Vector3 diff  = surfP - p;
        Vector3 diffn = normalize(diff);
        float tl = diff.x / diffn.x;
        
        // could sample
        if (tl > t0 && tl < t1) {
            
            Vector3 l = normalize(surfP - p);
            Vector3 n = normal;
            
            real_t dot_nl = dot(n, l);
            
//            return (color * intensity) * (1.0 / squared_distance(surfP, p)) * (dot_nl > 0 ? dot_nl : 0);
            // debug
            return (color * intensity) * (dot_nl > 0 ? dot_nl : 0);
        }
        
        return Color3::Black();
    }
    
    Vector3 PointLight::SamplePointOnLight(Vector3 &dir) const
    {
        // TODO: do not count in matrix at this point
        real_t x = _462::random_gaussian();
        real_t y = _462::random_gaussian();
        real_t z = _462::random_gaussian();
        
        Vector3 ran = Vector3(x, y, z);
        normalize(ran);
        dir = ran;
        
        return position + ran * radius;
    }
    
    real_t PointLight::Power() const
    {
        return (4.f * M_PI * intensity);
    }
    
}