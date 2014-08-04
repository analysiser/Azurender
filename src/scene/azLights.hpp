//
//  azLights.hpp
//  Azurender
//
//  Created by Xiao Li on 7/17/14.
//
//

#ifndef __Azurender__azLights__
#define __Azurender__azLights__

#include <iostream>

#include "math/camera.hpp"
#include "math/color.hpp"
#include "math/matrix.hpp"
#include "math/quaternion.hpp"
#include "math/vector.hpp"

#include "ray.hpp"

namespace _462 {
    
    // Light class based on PBRT implementation
    // Base class for all lights
    class Light {
    public:
        
        Light();
        virtual ~Light();
        
        // SampleLight, returns the radiance arriving at point p due to this light
        virtual Color3 SampleLight(const Vector3 &p, const Vector3 &normal, float t0, float t1) const = 0;
        
        // return a sample point on light surface
        virtual Vector3 SamplePointOnLight() const = 0;
        
        virtual Ray getRandomRayFromLight() const = 0;
        
        virtual Vector3 getLightEmissionDirection(const Vector3 &sampleOnLight) const = 0;
        
        virtual Vector3 getPointToLightDirection(const Vector3 &incidentPos, const Vector3 &lightPos) const = 0;
        
        virtual real_t Power() const = 0;
//
//        virtual bool isDeltaLight() const = 0;
        
        virtual bool initialize() const = 0;
        
        const int nSamples;
        
        
//    protected:
        
        // Transformation matrix
        const Matrix4 matLightToWorld;
        
        // Inverse transformation matrix
        const Matrix4 matWorldToLight;
        
        // Normal transformation matrix
        const Matrix3 matNormal;
        
        Quaternion orientation;
        
        Vector3 position;
        
        Color3 color;
        
        float intensity;
        
//        float radius;
        
//    private:

    };
    
    
    // Basic point light source
    class PointLight : public Light {
    public:
        
        PointLight();
//        ~PointLight();
        
//        PointLight(const Matrix4 &light2world, const Color3 &intensity) {
//            
//        }
        bool initialize() const;
        
        // TODO: add falloff operator
        Color3 SampleLight(const Vector3 &p, const Vector3 &normal, float t0, float t1) const;
        
        // sample a point on light surface
        Vector3 SamplePointOnLight() const;
        
        Ray getRandomRayFromLight() const;
        
        // given a sample point on light, get the emission direction from light source to the point
        // return normalized direction
        Vector3 getLightEmissionDirection(const Vector3 &sampleOnLight) const;
        
        // given a incident point on light, get the direction from incident point to light source
        // return normalized direction
        Vector3 getPointToLightDirection(const Vector3 &incidentPos, const Vector3 &lightPos) const;
        
        real_t Power() const;
        
//        Color3 Power() const {};
//        
//        bool isDeltaLight() const { return false; }
        
        
//    private:
        float radius;
        
        int attenuation_constant;
        
    };
    
    class DistantLight : public Light {
    public:
        DistantLight();
        
        bool initialize() const;
        
        Color3 SampleLight(const Vector3 &p, const Vector3 &normal, float t0, float t1) const;
        
        // sample a point on light surface
        Vector3 SamplePointOnLight() const;
        
        Ray getRandomRayFromLight() const;
        
        // given a sample point on light, get the emission direction from light source to the point
        // return normalized direction
        Vector3 getLightEmissionDirection(const Vector3 &sampleOnLight) const;
        
        // given a incident point on light, get the direction from incident point to light source
        // return normalized direction
        Vector3 getPointToLightDirection(const Vector3 &incidentPos, const Vector3 &lightPos) const;
        
        real_t Power() const;
        
        mutable Vector3 direction;
        mutable Vector3 inversed_direction;
        
    };
    
}

#endif /* defined(__Azurender__azLights__) */
