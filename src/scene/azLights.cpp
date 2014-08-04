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
    
    
#pragma mark - PointLight
    bool PointLight::initialize() const
    {
        return true;
    }
    
    PointLight::PointLight()
    {
        position = Vector3::Zero();
        color = Color3::Black();
        intensity = 1.0;
        radius = 0.0;
        attenuation_constant = 0.0;
    }
    
    
    // p and light in world position, need xform in future
    Color3 PointLight::SampleLight(const Vector3 &p, const Vector3 &normal, float t0, float t1) const
    {
        Vector3 surfP = SamplePointOnLight();//position + ran * radius;
        Vector3 diff  = surfP - p;
        Vector3 diffn = normalize(diff);
        float tl = diff.x / diffn.x;
        
        // could sample
        if (tl > t0 && tl < t1) {
            
            Vector3 l = diffn;
            Vector3 n = normal;
            
            real_t dot_nl = dot(n, l);
            
            if (attenuation_constant != 0.0) {
                return ( color * intensity * std::max(dot_nl, 0.0) * (real_t(1) / (attenuation_constant)) );
            }
            else {
                Color3 resultColor = ( color * (intensity * (1.0 / squared_distance(surfP, p))) * std::max(dot_nl, 0.0) );
                resultColor = clamp(resultColor, 0, 1.0);
                
                // TODO: is this right?
                return resultColor;//( color * (intensity * (1.0 / squared_distance(surfP, p))) * (dot_nl > 0 ? dot_nl : 0) );
            }
            
        }
        
        return Color3::Black();
    }
    
    Vector3 PointLight::SamplePointOnLight() const
    {
        real_t x = _462::random_gaussian();
        real_t y = _462::random_gaussian();
        real_t z = _462::random_gaussian();
        
        Vector3 ran = Vector3(x, y, z);
        normalize(ran);
        
        return position + ran * radius;
    }
    
    Ray PointLight::getRandomRayFromLight() const
    {
        real_t x = _462::random_gaussian();
        real_t y = _462::random_gaussian();
        real_t z = _462::random_gaussian();
        
        Vector3 ran = Vector3(x, y, z);
        normalize(ran);
        
        return Ray(position, ran);
    }
    
    Vector3 PointLight::getLightEmissionDirection(const Vector3 &sampleOnLight) const
    {
        return normalize(sampleOnLight - position);
    }
    
    Vector3 PointLight::getPointToLightDirection(const Vector3 &incidentPos, const Vector3 &lightPos) const
    {
        return normalize(lightPos - incidentPos);
    }
    
    real_t PointLight::Power() const
    {
        return (4.f * M_PI * intensity);
    }
    
#pragma mark - DistantLight
    DistantLight::DistantLight()
    {
        position = Vector3::Zero();
        color = Color3::Black();
        intensity = 1.0;
        direction = Vector3::Zero();
        inversed_direction = Vector3::Zero();
    }
    
    bool DistantLight::initialize() const
    {
        direction = normalize(direction);
        inversed_direction = -direction;
        inversed_direction = normalize(inversed_direction);
        return true;
    }
    
    Color3 DistantLight::SampleLight(const Vector3 &/* p */, const Vector3 &normal, float t0, float t1) const
    {
        float tl = 2;
        // visible
        if (tl > t0 && tl <= t1) {
            Vector3 l = inversed_direction;
            Vector3 n = normal;
            real_t dot_nl = dot(n, l);
            
            // no fall off
            Color3 resultColor = ( color * intensity ) * (dot_nl > 0 ? dot_nl : 0 );
//            Color3 resultColor = ( color * (intensity * (1.0 / (tl * squared_distance(direction, Vector3::Zero())))) * (dot_nl > 0 ? dot_nl : 0) );
            resultColor = clamp(resultColor, 0, 0.99);
//            float clampV = 0.1;
//            for (int i = 0; i < 3; i++) {
//                if (resultColor[i] > clampV) {
//                    return resultColor;
//                }
//            }
//            resultColor = clamp(resultColor, clampV, 1.0);
            
            return resultColor;
        }
        return Color3::Black();
    }
    
    // sample a point on light surface
    Vector3 DistantLight::SamplePointOnLight() const
    {
        return Vector3::Zero();
    }
    
    Ray DistantLight::getRandomRayFromLight() const
    {
        // This should not be called
        std::cout<<"No photon mapping for directional light!"<<std::endl;
        return Ray(Vector3::Zero(), direction);
    }
    
    // given a sample point on light, get the emission direction from light source to the point
    // return normalized direction
    Vector3 DistantLight::getLightEmissionDirection(const Vector3 &/* sampleOnLight */) const
    {
        return direction;
    }
    
    // given a incident point on light, get the direction from incident point to light source
    // return normalized direction
    Vector3 DistantLight::getPointToLightDirection(const Vector3 &/* incidentPos */, const Vector3 &/* lightPos */) const
    {
        return inversed_direction;
    }
    
    real_t DistantLight::Power() const
    {
        // TODO: estimate according to the whole scene
        return (4.f * M_PI * intensity);
    }
    
    
}