/**
 * @file raytacer.cpp
 * @brief Raytracer class
 *
 * Implement these functions for project 4.
 *
 * @author H. Q. Bovik (hqbovik)
 * @bug Unimplemented
 */

#include "raytracer.hpp"
#include "scene/scene.hpp"
#include "scene/model.hpp"

#include <SDL_timer.h>
#include <iostream>

#ifdef OPENMP // just a defense in case OpenMP is not installed.

#ifndef __APPLE__
#include <omp.h>
#endif

#endif

#define RAYTRACE_DEPTH          5
#define PHOTON_TRACE_DEPTH      5

#define EPSILON                     0.00001
#define TMAX                        300

#define PROB_DREFLECTION            0.5F
#define PROB_DABSORB                0.7F
#define INDIRECT_PHOTON_NEEDED      200000      // 200000
#define CAUSTICS_PHOTON_NEEDED      50000       // 50000

// Gaussian filter constants
#define E                               (2.718)
#define ALPHA                           (0.918)
#define BETA                            (1.953)
#define ONE_MINUS_E_TO_MINUS_BETA       (0.858)

// Cone filter constants    jensen P67, k >= 1 is a filter constant characterizing the filter
#define CONE_K                      (1)

namespace _462 {
    
    // max number of threads OpenMP can use. Change this if you like.
#define MAX_THREADS 8
    
    static const unsigned STEP_SIZE = 8;
    
    Raytracer::Raytracer()
    : scene(0), width(0), height(0) { }
    
    // random real_t in [0, 1)
    static inline real_t random()
    {
        return real_t(rand())/RAND_MAX;
    }
    
    Raytracer::~Raytracer() { }
    
    /**
     * Initializes the raytracer for the given scene. Overrides any previous
     * initializations. May be invoked before a previous raytrace completes.
     * @param scene The scene to raytrace.
     * @param width The width of the image being raytraced.
     * @param height The height of the image being raytraced.
     * @return true on success, false on error. The raytrace will abort if
     *  false is returned.
     */
    bool Raytracer::initialize(Scene* scene, size_t num_samples,
                               size_t width, size_t height)
    {
        /*
         * omp_set_num_threads sets the maximum number of threads OpenMP will
         * use at once.
         */
#ifndef __APPLE__
#ifdef OPENMP
        omp_set_num_threads(MAX_THREADS);
#endif
#endif
        this->scene = scene;
        this->num_samples = num_samples;
        this->width = width;
        this->height = height;
        
        current_row = 0;
        
        Ray::init(scene->camera);
        scene->initialize();
        
//        const SphereLight* lights = scene->get_lights();
        
        // Initialization or precompuation before the trace
        for (size_t i = 0; i < scene->num_geometries(); i++) {
            scene->get_geometries()[i]->createBoundingBox();
            if (scene->get_geometries()[i]->type == eModel) {
                Model *model = (Model *)scene->get_geometries()[i];
                // For models, create BVH tree for optimization
                model->createBVHTree();
            }
        }
        
        // TODO: Photon Mapping: build kd tree
        // scatter photons
        photonScatter(scene);
        
        
        return true;
    }
    
    // scatter photons
    void Raytracer::photonScatter(const Scene* scene)
    {
        // assuming there is only one light source
        SphereLight l = scene->get_lights()[0];
        // scatter indirect
        while ((photon_indirect_list.size() < INDIRECT_PHOTON_NEEDED))
        {
//            Vector3 p = samplePointOnUnitSphere();
//            Ray photonRay = Ray(p + l.position, p);
            Ray photonRay = getPhotonEmissionRayFromLight(scene->get_lights()[0]);
            // Make color only the portion of light
            photonRay.photon.color = l.color;// * (real_t(1)/(real_t)(INDIRECT_PHOTON_NEEDED + CAUSTICS_PHOTON_NEEDED));
            photonTrace(photonRay, EPSILON, TMAX, PHOTON_TRACE_DEPTH);
        }
        
        // scatter indirect
        while ((photon_caustic_list.size() < CAUSTICS_PHOTON_NEEDED))
        {
//            Vector3 p = samplePointOnUnitSphere();
//            Ray photonRay = Ray(p + l.position, p);
            Ray photonRay = getPhotonEmissionRayFromLight(scene->get_lights()[0]);
            // Make color only the portion of light
            photonRay.photon.color = l.color;// * (real_t(1)/(real_t)(CAUSTICS_PHOTON_NEEDED + INDIRECT_PHOTON_NEEDED));
            photonTrace(photonRay, EPSILON, TMAX, PHOTON_TRACE_DEPTH);
        }
        
        printf("Finished Scattering, indirect num = %ld, caustic num = %ld\n", photon_indirect_list.size(), photon_caustic_list.size());
        
        // balance
        std::vector<Photon> tmp_indirect(photon_indirect_list.size() + 1);
        std::vector<Photon> tmp_caustic(photon_caustic_list.size() + 1);

        balance(1, tmp_indirect, photon_indirect_list);
        balance(1, tmp_caustic, photon_caustic_list);
        
        kdtree_photon_indirect_list = tmp_indirect;
        kdtree_photon_caustic_list = tmp_caustic;
        
        pClosePhotonIndirect = new ClosePhoton[kdtree_photon_indirect_list.size()];
        pClosePhotonCaustics = new ClosePhoton[kdtree_photon_caustic_list.size()];
        
        // initialize close photon kdtree
        for (size_t i = 0; i < kdtree_photon_indirect_list.size(); i++)
        {
            Photon p = kdtree_photon_indirect_list[i];
            pClosePhotonIndirect[i] = ClosePhoton(p, i);
        }
        
        for (size_t i = 0; i < kdtree_photon_caustic_list.size(); i++)
        {
            Photon p = kdtree_photon_caustic_list[i];
            pClosePhotonCaustics[i] = ClosePhoton(p, i);
        }
        
        printf("Finished Balancing!\n");
    }
    
    /**
     * Performs a raytrace on the given pixel on the current scene.
     * The pixel is relative to the bottom-left corner of the image.
     * @param scene The scene to trace.
     * @param x The x-coordinate of the pixel to trace.
     * @param y The y-coordinate of the pixel to trace.
     * @param width The width of the screen in pixels.
     * @param height The height of the screen in pixels.
     * @return The color of that pixel in the final image.
     */
    Color3 Raytracer::trace_pixel(const Scene* scene,
                                  size_t x,
                                  size_t y,
                                  size_t width,
                                  size_t height)
    {
        assert(x < width);
        assert(y < height);
        
        real_t dx = real_t(1)/width;
        real_t dy = real_t(1)/height;
        
        Color3 res = Color3::Black();
        unsigned int iter;
        for (iter = 0; iter < num_samples; iter++)
        {
            // pick a point within the pixel boundaries to fire our
            // ray through.
            real_t i = real_t(2)*(real_t(x)+random())*dx - real_t(1);
            real_t j = real_t(2)*(real_t(y)+random())*dy - real_t(1);
            
//            current_refractive_index = scene->refractive_index;
            
            Ray r = Ray(scene->camera.get_position(), Ray::get_pixel_dir(i, j));
            res += trace(r, EPSILON, TMAX, RAYTRACE_DEPTH);
            
        }
        return res*(real_t(1)/num_samples);
    }
    
    /**
     * Raytraces some portion of the scene. Should raytrace for about
     * max_time duration and then return, even if the raytrace is not copmlete.
     * The results should be placed in the given buffer.
     * @param buffer The buffer into which to place the color data. It is
     *  32-bit RGBA (4 bytes per pixel), in row-major order.
     * @param max_time, If non-null, the maximum suggested time this
     *  function raytrace before returning, in seconds. If null, the raytrace
     *  should run to completion.
     * @return true if the raytrace is complete, false if there is more
     *  work to be done.
     */
    bool Raytracer::raytrace(unsigned char* buffer, real_t* max_time)
    {
        // TODO Add any modifications to this algorithm, if needed.
        
        static const size_t PRINT_INTERVAL = 64;
        
        // the time in milliseconds that we should stop
        unsigned int end_time = 0;
        bool is_done;
        
        if (max_time)
        {
            // convert duration to milliseconds
            unsigned int duration = (unsigned int) (*max_time * 1000);
            end_time = SDL_GetTicks() + duration;
        }
        
        // until time is up, run the raytrace. we render an entire group of
        // rows at once for simplicity and efficiency.
        for (; !max_time || end_time > SDL_GetTicks(); current_row += STEP_SIZE)
        {
            // we're done if we finish the last row
            is_done = current_row >= height;
            // break if we finish
            if (is_done) break;
            
            int loop_upper = std::min(current_row + STEP_SIZE, height);
            
            // This tells OpenMP that this loop can be parallelized.
#ifndef __APPLE__
#pragma omp parallel for
#endif
            for (int c_row = current_row; c_row < loop_upper; c_row++)
            {
                /*
                 * This defines a critical region of code that should be
                 * executed sequentially.
                 */
#ifndef __APPLE__
#pragma omp critical
#endif
                {
                    if (c_row % PRINT_INTERVAL == 0)
                        printf("Raytracing (Row %d)\n", c_row);
                }
                
                for (size_t x = 0; x < width; x++)
                {
                    // trace a pixel
                    Color3 color = trace_pixel(scene, x, c_row, width, height);
                    // write the result to the buffer, always use 1.0 as the alpha
                    color.to_array(&buffer[4 * (c_row * width + x)]);
                }
            }
        }
        
        if (is_done) {
            printf("Done raytracing!");
            // TODO: shading
            
        }
        
        return is_done;
    }
    
    /**
     * Do Photon tracing, store photons according to their type into indirect map or caustics map
     * @param   ray         incoming ray
     * @param   t0          lower limit of t
     * @param   t1          upper limit of t
     * @param   depth       maximun depth
     */
    void Raytracer::photonTrace(Ray ray, real_t t0, real_t t1, int depth)
    {
        if (depth == 0) {
            return;
        }
        
        bool isHit = false;
        bool isLight = false;
        HitRecord record = getClosestHit(ray, t0, t1, &isHit, &isLight);
        
        if (isHit && !isLight) {
            // refractive
            if (record.refractive_index > 0) {
//                real_t R = 1;
                
                // Trace for specularly transmissive
                real_t cos_theta = dot((ray.d), record.normal);
                Vector3 rd;     //refract ray direction
                real_t idxRatio = real_t(1)/record.refractive_index;
                real_t c = 0;
                
                bool isRefraction = false;
                if (cos_theta <= 0) {
                    // From air to dielectric material
                    isRefraction = isRefract(ray.d, record.normal, idxRatio, rd);
                    if(isRefraction)
                        c = dot(-ray.d, record.normal);
                }
                else {
                    // From dielectric material to air
                    isRefraction = isRefract(ray.d, -record.normal, real_t(1) / idxRatio, rd);
                    if (isRefraction)
                        c = dot(rd, record.normal);
                }
                
                if (isRefraction) {
                    // create refractive reflection for photon
                    Ray refractRay = Ray(record.position + EPSILON * rd , normalize(rd));
                    refractRay.photon = ray.photon;
                    refractRay.photon.color *= record.specular;
                    refractRay.photon.mask |= 0x2;
                    photonTrace(refractRay, t0, t1, depth-1);
                }
            }
            // specular reflective
            else if (record.specular != Color3::Black())
            {
                // Pure reflective surface
                if (record.diffuse == Color3::Black()) {

                    Ray reflectRay = Ray(record.position, normalize(ray.d - 2 * dot(ray.d, record.normal) * record.normal));
                    reflectRay.photon = ray.photon;
                    reflectRay.photon.color *= record.specular;
                    reflectRay.photon.mask |= 0x2;
                    photonTrace(reflectRay, t0, t1, depth - 1);
                    
                }
                // Hit on a surface that is both reflective and diffusive
                else {
                    real_t prob = random();
                    // Then there is a possibility of whether reflecting or absorbing
                    if (prob < 0.5) {
                        Ray reflectRay = Ray(record.position, normalize(ray.d - 2 * dot(ray.d, record.normal) * record.normal));
                        reflectRay.photon = ray.photon;
                        reflectRay.photon.color *= record.specular;
                        reflectRay.photon.mask |= 0x2;
                        photonTrace(reflectRay, t0, t1, depth - 1);
                    }
                    else {
                        // absorb
                        if (photon_indirect_list.size() < INDIRECT_PHOTON_NEEDED) {
                            ray.photon.position = record.position;
                            ray.photon.direction = -ray.d;
//                            ray.photon.color = ray.photon.color;// * record.diffuse;
                            ray.photon.normal = record.normal;
                            photon_indirect_list.push_back(ray.photon);
                        }
                    }
                }
            }
            // diffusive
            else {
                // direct illumination, do not store
                if (ray.photon.mask == 0x0) {
                    real_t prob = random();
                    if (prob > PROB_DABSORB) {
                        Ray photonRay = Ray(record.position, uniformSampleHemisphere(record.normal));
                        photonRay.photon = ray.photon;
                        photonRay.photon.mask |= 0x1;
                        photonRay.photon.color = ray.photon.color * record.diffuse;
                        photonTrace(photonRay, t0, t1, depth-1);
                    }
                }
                // caustics
                else if (ray.photon.mask == 0x2) {
                    if (photon_caustic_list.size() < CAUSTICS_PHOTON_NEEDED) {
                        ray.photon.position = record.position;
                        ray.photon.direction = -ray.d;
                        ray.photon.normal = (record.normal);
                        photon_caustic_list.push_back(ray.photon);
                        
                    }
                }
                // indirect illumination
                else {
                    real_t prob = random();
                    if (prob < PROB_DABSORB) {
                        // Store photon in indirect illumination map
                        if (photon_indirect_list.size() < INDIRECT_PHOTON_NEEDED) {
                            ray.photon.position = (record.position);
                            ray.photon.direction = (-ray.d);
                            ray.photon.normal = (record.normal);
//                            ray.photon.color = ray.photon.color * record.diffuse;
                            photon_indirect_list.push_back(ray.photon);
                        }
                    }
                    else {
                        // Generate a diffusive reflect
                        Ray photonRay = Ray(record.position, uniformSampleHemisphere(record.normal));
                        photonRay.photon = ray.photon;
                        photonRay.photon.mask |= 0x1;
                        photonRay.photon.color = (ray.photon.color * record.diffuse);
                        photonTrace(photonRay, t0, t1, depth-1);
                    }
                }
            }
        }
    }
    
    /**
     * @brief   Ray tracing a given ray, return the color of the hiting point on surface
     * @param   ray     Ray to trace
     * @param   t0      Lower limit of t
     * @param   t1      Upper limit of t
     * @param   depth   Recursion depth
     * @return  Color3  Shading color of hit point, if there is no hit, 
     *                  then background color
     */
    Color3 Raytracer::trace(Ray ray, real_t t0, real_t t1, int depth)
    {
        bool isHit;
        bool isLight;
        HitRecord record = getClosestHit(ray, t0, t1, &isHit, &isLight);
        
        if (isHit && !isLight) {
            return shade(ray, record, t0, t1, depth);
        }
        else if (isLight) {
            // TODO: lixiao debug force color to be white...
//            return Color3::White();//scene->get_lights()[0].color;
            return scene->get_lights()[0].color;
        }
        else
            return scene->background_color;
        
    }
    
    /**
     * @brief   Shading a hit point on surface
     * @param   ray     Ray to trace
     * @param   record  Hit record
     * @param   t0      Lower limit of t
     * @param   t1      Upper limit of t
     * @param   depth   Recursion depth
     * @return  Color3  Shading color of hit point, after calculating diffusive, ambient, reflection etc.
     */
    Color3 Raytracer::shade(Ray ray, HitRecord record, real_t t0, real_t t1, int depth)
    {
        Color3 radiance = Color3::Black();
        
        if (depth == 0) {
            return radiance;
        }
        
        // refraction
        if (record.refractive_index > 0) {
            real_t R = 1;
            
            Color3 reflColor = Color3::Black();
            Color3 refrColor = Color3::Black();
            
            // Trace reflection color
            Vector3 reflectDirection = ray.d - 2 * dot(ray.d, record.normal) * record.normal;
            Ray reflectRay = Ray(record.position + EPSILON * reflectDirection, normalize(reflectDirection));
            reflColor = record.specular * trace(reflectRay, t0, t1, depth - 1);
            
            // Trace for specularly transmissive
            real_t cos_theta = dot((ray.d), record.normal);
            Vector3 rd;     //refract ray direction
            real_t idxRatio = real_t(1)/record.refractive_index;
            real_t c = 0;
            
            bool isRefraction = false;
            
            if (cos_theta <= 0) {
                // From air to dielectric material
                isRefraction = isRefract(ray.d, record.normal, idxRatio, rd);
                if(isRefraction)
                    c = dot(-ray.d, record.normal);
            }
            else {
                // From dielectric material to air
                isRefraction = isRefract(ray.d, -record.normal, real_t(1) / idxRatio, rd);
                if (isRefraction) {
                    c = dot(rd, record.normal);
                }
            }
            
            // If there is refraction, recursively calculate color
            if (isRefraction) {
                
                Ray refractRay = Ray(record.position + EPSILON * rd , normalize(rd));
                real_t R0 = ((idxRatio - real_t(1)) * (idxRatio - real_t(1))) / ((idxRatio + real_t(1)) * (idxRatio + real_t(1)));
                R = R0 + (real_t(1) - R0) * pow((real_t(1) - c), 5);
                refrColor = trace(refractRay, t0, t1, depth - 1);
            }
            
            radiance += R * reflColor + (1.0 - R) * refrColor;
            
            
            
            
        }
        // reflection
        else if (record.specular != Color3::Black())
        {
            // Trace for specularly reflective
            Ray reflectRay = Ray(record.position, normalize(ray.d - 2 * dot(ray.d, record.normal) * record.normal));
            radiance += record.specular * trace(reflectRay, t0, t1, depth - 1);
            
//            if (record.diffuse != Color3::Black())
//                radiance *= shade_direct_illumination(record, t0, t1);
            
        }
        else {
            
            // Trace each light source for direct illumination
            radiance += shade_direct_illumination(record, t0, t1);
            
//            radiance += shade_photon_illumination(record);
            // Shading caustics
            radiance += shade_caustics(record);
            
            // Shading global illumination
//            radiance += shade_indirect_illumination(record);
        }
        
        // for raytracing ambient
//        if (record.diffuse != Color3::Black() && record.refractive_index == 0) {
//            radiance += scene->ambient_light * record.ambient;
//        }
        
        // TODO:
        // try to add high light with blinn-phong model
        if (record.phong > 0) {
            SphereLight light = scene->get_lights()[0];
            Vector3 lightpos = sampleLightSource(light);
            Vector3 l = normalize(lightpos - record.position);
            Ray rayToLight = Ray(record.position, l);
            bool isHit, isHitLight;
            getClosestHit(rayToLight, EPSILON, TMAX, &isHit, &isHitLight);
            Color3 highLightColor = light.color;
            if (isHit && isHitLight) {
                // should be high light
                Vector3 r = -l + 2 * dot(l, record.normal) * record.normal;
                Vector3 e = -ray.d;
                real_t coef = dot(e, r);
                coef = coef > 0 ? coef : 0;
                coef = pow(coef, record.phong);
                highLightColor *= coef;
            }
            else
                highLightColor = Color3::Black();
            
            radiance += highLightColor;
        }
        
        return record.texture * radiance;
        // consider diffusive color
//        return record.texture * record.diffuse * radiance;
        
    }
    
    /**
     * @brief Do direct illumination for ray tracing
     */
    Color3 Raytracer::shade_direct_illumination(HitRecord &record, real_t t0, real_t t1)
    {
        Color3 res = Color3::Black();
        
        bool isSoftShadow = true;
        
        for (size_t i = 0; i < scene->num_lights(); i++) {
            
            if (isSoftShadow)
            {
                // make random points on a light, this would make egdes of shadows smoother
                size_t sample_num_per_light = 16;    // todo: add samples
                for (size_t j = 0; j < sample_num_per_light; j++) {
                    
                    // Get a random point on light sphere
                    Vector3 sampleLightPos = sampleLightSource(scene->get_lights()[i]);
                    
                    // scene->get_lights()[i].position; // For hard shadows
                    
                    Vector3 d_shadowRay = sampleLightPos - record.position;
                    Vector3 d_shadowRay_normolized = normalize(d_shadowRay);
                    Ray shadowRay = Ray(record.position, d_shadowRay_normolized);
                    real_t tl = d_shadowRay.x / d_shadowRay_normolized.x;
                    
                    bool isInShadow = false;
                    bool isHitLight = false;
                    tl = tl < t1 ? tl : t1;
                    getClosestHit(shadowRay, t0, tl, &isInShadow, &isHitLight);
                    
                    // if the obejct is not in shadow and is opaque
                    if (!isInShadow && (record.refractive_index == 0) && !isHitLight) {
                        res += (record.getPixelColorFromLight(scene->get_lights()[i])) * (real_t(1)/real_t(sample_num_per_light));
                    }
                    else if (isInShadow && isHitLight) {
                        res += (record.getPixelColorFromLight(scene->get_lights()[i])) * (real_t(1)/real_t(sample_num_per_light));
                    }
                    else
                    {
                        // in shadow
                    }
                }
            }
            else
            {
                SphereLight light = scene->get_lights()[i];
                Vector3 sampleLightPos;
                if (light.type == 1) {
                    // Get a random point on light sphere
                    sampleLightPos = light.position;
                    
                   
                }
                else if (light.type == 2)
                {
                    sampleLightPos = (light.vertex1 + light.vertex2)/2;
                }
                else
                    printf("error light type\n");
                
                Vector3 d_shadowRay = sampleLightPos - record.position;
                Vector3 d_shadowRay_normolized = normalize(d_shadowRay);
                Ray shadowRay = Ray(record.position, d_shadowRay_normolized);
                real_t tl = d_shadowRay.x / d_shadowRay_normolized.x;
                
                bool isInShadow = false;
                bool isHitLight = false;
                tl = tl < t1 ? tl : t1;
                getClosestHit(shadowRay, t0, tl, &isInShadow, &isHitLight);
                
                // if the obejct is not in shadow and is opaque
                if (!isInShadow && (record.refractive_index == 0) && !isHitLight) {
                    res += (record.getPixelColorFromLight(scene->get_lights()[i]));
                }
                else if (isInShadow && isHitLight) {
                    res += (record.getPixelColorFromLight(scene->get_lights()[i]));
                }
                else
                {
                    // in shadow
                }
            }
            
            
        }
        
        return res;
    }
    
    // Shading of caustics
    Color3 Raytracer::shade_caustics(HitRecord &record)
    {
        Color3 causticsColor = Color3::Black();
        std::vector<Photon> nearestCausticPhotons;
        
        if (record.refractive_index == 0)
        {
            // sample radius
            real_t sampleSquaredRadiusCaustics = 0.001;   // 0.1 // 0.001 // 0.0225
            real_t maxSquaredDistCaustics = 0.001;      // 0.001 // 0.001
            size_t maxPhotonsEstimate = 60;
            
            // gather samples
            // TODO: doing nearest neighbor search with KD Tree
            locatePhotons(1, record.position, kdtree_photon_caustic_list, nearestCausticPhotons, sampleSquaredRadiusCaustics, maxPhotonsEstimate);
            
//            for (std::vector<Photon>::iterator it = kdtree_photon_caustic_list.begin(); it != kdtree_photon_caustic_list.end(); it++)
//            {
//                real_t squaredDist = squared_distance(record.position, it->position);
//                if (squaredDist < sampleSquaredRadiusCaustics)
//                {
////                    if (dot(record.normal, it->normal) == 1)
//                    {
//                        Photon tmp = *it;
//                        tmp.squaredDistance = squaredDist;
//                        
//                        // if photons less than needed to estimate
//                        if (nearestCausticPhotons.size() < maxPhotonsEstimate) {
//                            nearestCausticPhotons.push_back(tmp);
//                            std::push_heap(nearestCausticPhotons.begin(), nearestCausticPhotons.end());
//                        }
//                        else
//                        {
////                            std::make_heap(nearestCausticPhotons.begin(), nearestCausticPhotons.end(), Photon());
//                            if (tmp < nearestCausticPhotons.front())
//                            {
//                                std::pop_heap(nearestCausticPhotons.begin(), nearestCausticPhotons.end());
//                                nearestCausticPhotons.pop_back();
//                                nearestCausticPhotons.push_back(tmp);
//                                std::push_heap(nearestCausticPhotons.begin(), nearestCausticPhotons.end());
//                            }
//                        }
//                        
//                    }
//                }
//            }
            
            // calculate radiance
            for (size_t i = 0; i < nearestCausticPhotons.size(); i++)
            {
                causticsColor +=
                record.getPhotonLambertianColor(nearestCausticPhotons[i].direction, nearestCausticPhotons[i].color)
                * getConeFilterWeight(sqrt(nearestCausticPhotons[i].squaredDistance), sqrt(sampleSquaredRadiusCaustics));
                //* getGaussianFilterWeight(nearestCausticPhotons[i].squaredDistance, sampleSquaredRadiusCaustics);
            }
            
            // color/= PI*r^2
            causticsColor *= (real_t(1)/((real_t(1) - real_t(2)/(real_t(3) * CONE_K)) * (PI * maxSquaredDistCaustics))) * 0.001;
            //causticsColor *= (real_t(1)/(PI * maxSquaredDistCaustics));
        }
        
        nearestCausticPhotons.clear();
        
        return causticsColor;
        
    }
    
    // Shading of Indirect lightning
    Color3 Raytracer::shade_indirect_illumination(HitRecord &record)
    {
        Color3 indirectColor = Color3::Black();
        std::vector<Photon> nearestIndirectPhotons;
        
        if (record.refractive_index == 0)
        {
            real_t sampleSquaredRadiusIndirect = 0.001;       // 0.5
            real_t maxSquaredDistIndirect = 0.001;          // 0.001
            size_t maxPhotonsEstimate = 100;
            
//            locatePhotons(1, record.position, kdtree_photon_indirect_list, nearestIndirectPhotons, sampleSquaredRadiusIndirect, maxPhotonsEstimate);
            
            // gather samples
            // TODO: doing nearest neighbor search with KD Tree
            for (std::vector<Photon>::iterator it = photon_indirect_list.begin(); it != photon_indirect_list.end(); it++)
            {
                real_t squaredDist = squared_distance(record.position, it->position);
                if (squaredDist < sampleSquaredRadiusIndirect)
                {
//                    if (dot(record.normal, it->normal) == 1)
                    {
                        Photon tmp = *it;
                        tmp.squaredDistance = squaredDist;
                        
                        // if photons less than needed to estimate
                        if (nearestIndirectPhotons.size() < maxPhotonsEstimate) {
                            nearestIndirectPhotons.push_back(tmp);
                            std::push_heap(nearestIndirectPhotons.begin(), nearestIndirectPhotons.end());
                        }
                        else
                        {
                            if (tmp < nearestIndirectPhotons.front())
                            {
                                std::pop_heap(nearestIndirectPhotons.begin(), nearestIndirectPhotons.end());
                                nearestIndirectPhotons.pop_back();
                                nearestIndirectPhotons.push_back(tmp);
                                std::push_heap(nearestIndirectPhotons.begin(), nearestIndirectPhotons.end());
                            }
                        }
                        
                    }
                }
            }
            
            // cone filter to gather radiance of indirect illumination photons
            for (size_t i = 0; i < nearestIndirectPhotons.size(); i++)
            {
                indirectColor += record.getPhotonLambertianColor(nearestIndirectPhotons[i].direction, nearestIndirectPhotons[i].color)
                //* getConeFilterWeight(sqrt(nearestIndirectPhotons[i].squaredDistance), sqrt(sampleSquaredRadiusIndirect));
                * getGaussianFilterWeight(nearestIndirectPhotons[i].squaredDistance, sampleSquaredRadiusIndirect);
            }
            //indirectColor *= (real_t(1)/((real_t(1) - real_t(2)/(real_t(3) * CONE_K)) * (PI * maxSquaredDistIndirect)));
            indirectColor *= (real_t(1)/(PI * maxSquaredDistIndirect)) * 0.001;
            
        }
        
        return indirectColor;
    }
    
    Color3 Raytracer::shade_photon_illumination(HitRecord &record)
    {
//        Color3 radiance = Color3::Black();
        std::vector<Photon> nearestCausticPhotons;
        std::vector<Photon> nearestIndirectPhotons;
        
        // Calculate photon colors
        Color3 indirectColor = Color3::Black();
        Color3 causticsColor = Color3::Black();
        
        
//        if ((record.refractive_index == 0) && (record.diffuse != Color3::Black()))
        if (record.refractive_index == 0)
        {
            real_t sampleSquaredRadiusCaustics = 0.1;   // 0.1
            real_t sampleSquaredRadiusIndirect = 0.5;
            
            real_t maxSquaredDistCaustics = 0.001;      // 0.01
            real_t maxSquaredDistIndirect = 0.001;
            
//            real_t guassianMaxCaustics = -INFINITY;
//            real_t gaussianMaxIndirect = -INFINITY;
            // caustics
            for (size_t i = 0; i < photon_caustic_list.size(); i++) {
                real_t squaredDist = squared_distance(record.position, photon_caustic_list[i].position);
                if (squaredDist < sampleSquaredRadiusCaustics)
                {
                    if (dot(record.normal, photon_indirect_list[i].normal) == 1)
                    {
                        photon_caustic_list[i].squaredDistance = squaredDist;
                        nearestCausticPhotons.push_back(photon_caustic_list[i]);
                    }
                }
            }
            
            // indirect lighting
            for (size_t i = 0; i < photon_indirect_list.size(); i++) {
                real_t squaredDist = squared_distance(record.position, photon_indirect_list[i].position);
                if (squaredDist < sampleSquaredRadiusIndirect)
                {
                    if (dot(record.normal, photon_indirect_list[i].normal) == 1)
                    {
                        photon_indirect_list[i].squaredDistance = squaredDist;
                        nearestIndirectPhotons.push_back(photon_indirect_list[i]);
                    }
                    
                }
            }
            
            
            // calculate radiance
            for (size_t i = 0; i < nearestCausticPhotons.size(); i++)
            {
                causticsColor +=
                record.getPhotonLambertianColor(nearestCausticPhotons[i].direction, nearestCausticPhotons[i].color)
                * getGaussianFilterWeight(nearestCausticPhotons[i].squaredDistance, sampleSquaredRadiusCaustics);
            }
            causticsColor *= (real_t(1)/(PI * maxSquaredDistCaustics));
            
            // cone filter to gather radiance of indirect illumination photons
            for (size_t i = 0; i < nearestIndirectPhotons.size(); i++)
            {
                indirectColor += record.getPhotonLambertianColor(nearestIndirectPhotons[i].direction, nearestIndirectPhotons[i].color)
                //* getConeFilterWeight(sqrt(nearestIndirectPhotons[i].squaredDistance), sqrt(sampleSquaredRadiusIndirect));
                * getGaussianFilterWeight(nearestIndirectPhotons[i].squaredDistance, sampleSquaredRadiusIndirect);
            }
            //indirectColor *= (real_t(1)/((real_t(1) - real_t(2)/(real_t(3) * CONE_K)) * (PI * maxSquaredDistIndirect)));
            indirectColor *= (real_t(1)/(PI * maxSquaredDistIndirect));
        }
        
        nearestCausticPhotons.clear();
        nearestIndirectPhotons.clear();
        
        return indirectColor + causticsColor;
    }
    
    /**
     * @brief   Sampling light volumn by sampling on light 
     *          sphere with gaussian distribution
     * @param   light   Instance of SphereLight
     */
    Vector3 Raytracer::sampleLightSource(SphereLight light)
    {
        Vector3 ran = Vector3::Zero();
        
        // sphere light
        if (light.type == 1) {
            
            real_t x = _462::random_gaussian();
            real_t y = _462::random_gaussian();
            real_t z = _462::random_gaussian();
            
            ran = Vector3(x, y, z);
            
            ran = normalize(ran);
            ran *= light.radius;
            ran += light.position;
            
        }
        else if (light.type == 2)
        {
            real_t b = random() * real_t(2) - real_t(1);
            real_t c = random() * real_t(2) - real_t(1);//_462::random_uniform();
//            real_t d = _462::random_uniform();
            
            Vector3 vec1 = light.vertex1 - light.position;
            Vector3 vec2 = light.vertex2 - light.position;
//            Vector3 vec3 = normalize(cross(vec2, vec1));
            
            ran = light.position + (vec1) * b + (vec2) * c;// + vec3 * d;
        }
        else
        {
            printf("Wrong Light Type, must be 1 (sphere light) or 2(parallelogram light)\n");
        }
        
        return ran;
    }
    
    /**
     * @brief   Sampling light volumn by sampling on light
     *
     * @param   light   Instance of SphereLight
     * @return  Randomized photon emmision ray
     */
    Ray Raytracer::getPhotonEmissionRayFromLight(SphereLight light)
    {
        Ray photonRay;
        // sphere light
        if (light.type == 1)
        {
            Vector3 p = samplePointOnUnitSphere();
            photonRay = Ray(p + light.position, p);
        }
        // square light
        else if (light.type == 2)
        {
            Vector3 normal = light.normal;
//            Vector3 p = samplePointOnUnitSphere();//sampleLightSource(light);
            
//            real_t u = random();
//            real_t v = 2 * PI * random();
//            Vector3 d = -Vector3(cos(v) * sqrt(u), sin(v) * sqrt(u), sqrt(1 - u));
//            std::cout<<light.position<<std::endl;
            photonRay = Ray(light.position, uniformSampleHemisphere(normal));
            
        }
        else
        {
            printf("error in light type!\n");
        }
        
        return photonRay;
    }
    
    
    /**
     * Retrieve the closest hit record
     * @param r             incoming ray
     * @param t0            lower limit of t
     * @param t1            upper limit of t
     * @param *isHit        indicate if there is a hit incident
     * @return HitRecord    the closest hit record
     */
    HitRecord Raytracer::getClosestHit(Ray r, real_t t0, real_t t1, bool *isHit, bool *isLight)
    {
        HitRecord closestHitRecord;
        HitRecord tmp;
        real_t t = t1;
        
        Geometry* const* geometries = scene->get_geometries();
        *isHit = false;
        for (size_t i = 0; i < scene->num_geometries(); i++) {
            
            if (geometries[i]->hit(r, t0, t1, tmp))
            {
                bool hitLight = (geometries[i]->isLight == 1.0);//(tmp.specular == Color3(-1.0, -1.0, -1.0));
                
                // if not hitting light
                if (!hitLight)
                {
                    if (!*isHit) {
                        *isHit = true;
                        t = tmp.t;
                        closestHitRecord = tmp;
                        
                        *isLight = false;
                    }
                    else {
                        if (tmp.t < t) {
                            t = tmp.t;
                            closestHitRecord = tmp;
                            
                            *isLight = false;
                        }
                    }
                }
                // if hit the light
                else
                {
                    // if not hit anything before
                    if (!*isHit) {
                        
                        *isHit = true;
                        t = tmp.t;
                        closestHitRecord = tmp;
                        
                        *isLight = true;
                    }
                    else {
                        if (tmp.t < t) {
                            
                            t = tmp.t;
                            closestHitRecord = tmp;
                            
//                            *isHit = false;
                            *isLight = true;
                        }
                    }
                }
            }
        }
        
        closestHitRecord.t = t;
        return closestHitRecord;
    }
    
    Vector3 Raytracer::samplePointOnUnitSphere()
    {
        real_t x = _462::random_gaussian();
        real_t y = _462::random_gaussian();
        real_t z = _462::random_gaussian();
        
        Vector3 ran = Vector3(x, y, z);
        return normalize(ran);
    }
    
    Vector3 Raytracer::uniformSampleHemisphere(const Vector3& normal)
    {
        Vector3 newDir = samplePointOnUnitSphere();
        if (dot(newDir, normal) < 0.0) {
            newDir = -newDir;
        }
        return normalize(newDir);
    }
    
    inline real_t getGaussianFilterWeight(real_t dist_sqr, real_t radius_sqr)
    {
        real_t power = -(BETA * 0.5F * (dist_sqr/radius_sqr));
        real_t numerator = 1.0F - pow(E, power);
        real_t denominator = ONE_MINUS_E_TO_MINUS_BETA;
        return ALPHA * (1.0F - (numerator/denominator));
    }
    
    inline real_t getConeFilterWeight(real_t dist_x_p, real_t dist_max)
    {
        return (real_t(1) - dist_x_p/(CONE_K * dist_max));
    }
    
    inline void balance(size_t index, std::vector<Photon> &balancedKDTree, std::vector<Photon> list)
    {
        if (index == 1) {
            assert((balancedKDTree.size() == (list.size() + 1)));
        }
        
        if (list.size() == 1) {
            balancedKDTree[index] = list[0];
            return;
        }
        
        // If there is no data in photon list, do nothing and return
        if (list.size() == 0) {
            // actually program should not run to here, this should be checked on subdivision
            return;
        }
        
        // get the surrounding cube
        Vector3 max = Vector3(-INFINITY, -INFINITY, -INFINITY);
        Vector3 min = Vector3(INFINITY, INFINITY, INFINITY);
        for (std::vector<Photon>::iterator it = list.begin(); it != list.end(); it++)
        {
            // calculate box
            max.x = it->position.x > max.x ? it->position.x : max.x;
            max.y = it->position.y > max.y ? it->position.y : max.y;
            max.z = it->position.z > max.z ? it->position.z : max.z;
            
            min.x = it->position.x < min.x ? it->position.x : min.x;
            min.y = it->position.y < min.y ? it->position.y : min.y;
            min.z = it->position.z < min.z ? it->position.z : min.z;
        }
        
        char splitAxis = -1;
        Vector3 diff = Vector3(max.x - min.x, max.y - min.y, max.z - min.z);
        if (diff.x >= diff.y && diff.x >= diff.z)
            splitAxis = 0;
        else if (diff.y >= diff.x && diff.y >= diff.z)
            splitAxis = 1;
        else if (diff.z >= diff.x && diff.z >= diff.y)
            splitAxis = 2;
        
        // Sorting the vector
        bool (*comparator)(const Photon &a, const Photon &b) = NULL;
        
        switch (splitAxis) {
            case 0:
                comparator = photonComparatorX;
                break;
            case 1:
                comparator = photonComparatorY;
                break;
            case 2:
                comparator = photonComparatorZ;
                break;
                
            default:
                break;
        }
        
        std::sort(list.begin(), list.end(), comparator);
        
        // On Left-balancing Binary Trees, J. Andreas Brentzen (jab@imm.dtu.dk)
        size_t N = list.size();
        size_t exp = (size_t)log2(N);
        size_t M = pow(2, exp);
        size_t R = N - (M - 1);
        size_t LT, RT;
        if (R <= M/2)
        {
            LT = (M - 2)/2 + R;
            RT = (M - 2)/2;
        }
        else
        {
            LT = (M - 2)/2 + M/2;
            RT = (M - 2)/2 + R - M/2;
        }
        size_t const median = LT;
        
        // if more than One data
        Photon p = list[median];
        p.splitAxis = splitAxis;
        if (index >= balancedKDTree.size()) {
            printf("index = %ld, size = %ld\n", index, balancedKDTree.size());
        }
        assert(index < balancedKDTree.size());
        balancedKDTree[index] = p;
        
//        printf("LT = %d, RT = %d\n", LT, RT);
        
        if (LT > 0) {
            std::vector<Photon> leftList(list.begin(), list.begin() + median);
            balance(2 * index, balancedKDTree, leftList);
        }
        
        if (RT > 0) {
            std::vector<Photon> rightList(list.begin() + median + 1, list.end());
            balance(2 * index + 1, balancedKDTree, rightList);
        }
        
        list.clear();
        
        
        
    }
    
    inline void locatePhotons(size_t p,
                              Vector3 position,
                              std::vector<Photon> balancedKDTree,
                              std::vector<Photon> &nearestPhotons,
                              real_t &sqrDist,
                              size_t maxNum)
    {
        // examine child nodes
        Photon photon = balancedKDTree[p];
        
        if (2 * p + 1 < balancedKDTree.size())
        {
            assert(photon.splitAxis != -1);
            Vector3 diff = position - photon.position;
            real_t diffToPlane;
            switch (photon.splitAxis) {
                case 0:
                    diffToPlane = diff.x;
                    break;
                case 1:
                    diffToPlane = diff.y;
                    break;
                case 2:
                    diffToPlane = diff.z;
                    break;
                default:
                    assert(0);
                    break;
            }
            real_t sqrDiffToPlane = diffToPlane * diffToPlane;
            
            if (diffToPlane < 0)
            {
                // search left subtree
                locatePhotons(2 * p, position, balancedKDTree, nearestPhotons, sqrDist, maxNum);
                if (sqrDiffToPlane < sqrDist)
                {
                    // check right subtree
                    locatePhotons(2 * p + 1, position, balancedKDTree, nearestPhotons, sqrDist, maxNum);
                }
            }
            else
            {
                // search right subtree
                locatePhotons(2 * p + 1, position, balancedKDTree, nearestPhotons, sqrDist, maxNum);
                if (sqrDiffToPlane < sqrDist)
                {
                    // check left subtree
                    locatePhotons(2 * p, position, balancedKDTree, nearestPhotons, sqrDist, maxNum);
                }
            }
            
        }
        
        // compute true squared distance to photon
        real_t sqrDistPhoton = squared_distance(position, photon.position);
        if (sqrDistPhoton <= sqrDist)
        {
            photon.squaredDistance = sqrDistPhoton;
            if (nearestPhotons.size() < maxNum)
            {
                nearestPhotons.push_back(photon);
                std::push_heap(nearestPhotons.begin(), nearestPhotons.end());
            }
            else
            {
                if (sqrDistPhoton < nearestPhotons.front().squaredDistance)
                {
                    std::pop_heap(nearestPhotons.begin(), nearestPhotons.end());
                    nearestPhotons.pop_back();
                    nearestPhotons.push_back(photon);
                    std::push_heap(nearestPhotons.begin(), nearestPhotons.end());
                    
                    sqrDist = photon.squaredDistance;
                }
            }
        }
    }
    
//    inline void locatePhotonsCompact(size_t p,
//                                     Vector3 position,
//                                     ClosePhoton *closePhotons,
//                                     std::vector<ClosePhoton> &nearestClosePhotons,
//                                     real_t &sqrDist,
//                                     size_t maxNum)
//    {
//        // examine child nodes
//        ClosePhoton photon = closePhotons[p];
//        
//        if (2 * p + 1 < balancedKDTree.size())
//        {
//            assert(photon.splitAxis != -1);
//            Vector3 diff = position - photon.position;
//            real_t diffToPlane;
//            switch (photon.splitAxis) {
//                case 0:
//                    diffToPlane = diff.x;
//                    break;
//                case 1:
//                    diffToPlane = diff.y;
//                    break;
//                case 2:
//                    diffToPlane = diff.z;
//                    break;
//                default:
//                    assert(0);
//                    break;
//            }
//            
//            if (diffToPlane < 0)
//            {
//                // search left subtree
//                locatePhotons(2 * p, position, balancedKDTree, nearestPhotons, sqrDist, maxNum);
//                if (diffToPlane * diffToPlane < sqrDist)
//                {
//                    // check right subtree
//                    locatePhotons(2 * p + 1, position, balancedKDTree, nearestPhotons, sqrDist, maxNum);
//                }
//            }
//            else
//            {
//                // search right subtree
//                locatePhotons(2 * p + 1, position, balancedKDTree, nearestPhotons, sqrDist, maxNum);
//                if (diffToPlane * diffToPlane < sqrDist)
//                {
//                    // check left subtree
//                    locatePhotons(2 * p, position, balancedKDTree, nearestPhotons, sqrDist, maxNum);
//                }
//            }
//            
//        }
//        
//        // compute true squared distance to photon
//        real_t sqrDistPhoton = squared_distance(position, photon.position);
//        if (sqrDistPhoton < sqrDist)
//        {
//            photon.squaredDistance = sqrDistPhoton;
//            if (nearestPhotons.size() < maxNum)
//            {
//                nearestPhotons.push_back(photon);
//                std::push_heap(nearestPhotons.begin(), nearestPhotons.end());
//            }
//            else
//            {
//                if (sqrDistPhoton < nearestPhotons.front().squaredDistance)
//                {
//                    std::pop_heap(nearestPhotons.begin(), nearestPhotons.end());
//                    nearestPhotons.pop_back();
//                    nearestPhotons.push_back(photon);
//                    std::push_heap(nearestPhotons.begin(), nearestPhotons.end());
//                    
//                    sqrDist = photon.squaredDistance;
//                }
//            }
//            
//        }
//    }
    
} /* _462 */
