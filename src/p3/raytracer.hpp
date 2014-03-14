/**
 * @file raytacer.hpp
 * @brief Raytracer class
 *
 * Implement these functions for project 3.
 *
 * @author H. Q. Bovik (hqbovik)
 * @bug Unimplemented
 */

#ifndef _462_RAYTRACER_HPP_
#define _462_RAYTRACER_HPP_

#define MAX_DEPTH 5

#include "math/color.hpp"
#include "math/random462.hpp"
#include "math/vector.hpp"
#include "scene/scene.hpp"
#include "p3/Photon.hpp"
#include <stack>

namespace _462 {
    
    class Scene;
    class Ray;
    struct Intersection;
    class Raytracer
    {
    public:
        
        Raytracer();
        
        ~Raytracer();
        
        bool initialize(Scene* scene, size_t num_samples,
                        size_t width, size_t height);
        
        bool raytrace(unsigned char* buffer, real_t* max_time);
        
        // indirect and caustics list for photons to trace
        std::vector<Photon> photon_indirect_list;
        std::vector<Photon> photon_caustic_list;
        
        
    private:
        
        Color3 trace_pixel(const Scene* scene,
                           size_t x,
                           size_t y,
                           size_t width,
                           size_t height);
        
        // Photon Mapping
        void photonScatter(const Scene* scene);
        
        // Photon Tracing for global illumination
        void photonTrace(Ray ray, real_t t0, real_t t1, int depth);
        
        // Raytracing helper function, to decide if there is a hit on a surface to shade
        Color3 trace(Ray ray, real_t t0, real_t t1, int depth);
        
        // Shading function, shades the hit record from a surface
        Color3 shade(Ray ray, HitRecord record, real_t t0, real_t t1, int depth);
        
        // Shading of direct illumination
        Color3 shade_direct_illumination(HitRecord &record, real_t t0, real_t t1);
        
        // TODO: seperate this to independent functions
        // Shading of photon radiance gathering
        Color3 shade_photon_illumination(HitRecord &record);
        
        // Shading of caustics
        Color3 shade_caustics(HitRecord &record);
        
        // Shading of Indirect lightning
        Color3 shade_indirect_illumination(HitRecord &record);
        
        // retrieve the closest hit record
        HitRecord getClosestHit(Ray r, real_t t0, real_t t1, bool *isHit, bool *isLight);
        
        // sample a Light source on the volumn of sphere light
        Vector3 sampleLightSource(SphereLight light);
        
        Ray getPhotonEmissionRayFromLight(SphereLight light);
        
        // helper function for sampling a point on a given unit sphere
        Vector3 samplePointOnUnitSphere();
        
        // sample a random direction along the hemisphere defined at the normal
        Vector3 uniformSampleHemisphere(const Vector3& normal);
        
        // the scene to trace
        Scene* scene;
        
        // the dimensions of the image to trace
        size_t width, height;
        
        // the next row to raytrace
        size_t current_row;
        
        //
        unsigned int num_samples;
        
//        // Current refraction index
//        real_t current_refractive_index;
        
        
        
    };
    
    inline real_t getGaussianFilterWeight(real_t dist_sqr, real_t radius_sqr);
    /*!
     @brief get the gaussian filter weight, jensen P67, this value is the only weight
     @param dist_x_p    distance between x and photon p
     @param dist_max    r is the maximum distance
    */
    inline real_t getConeFilterWeight(real_t dist_x_p, real_t dist_max);
    
} /* _462 */

#endif /* _462_RAYTRACER_HPP_ */
