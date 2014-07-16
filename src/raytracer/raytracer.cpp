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
#include "scene/triangle.hpp"
#include "scene/sphere.hpp"

#include <SDL_timer.h>
#include <iostream>
//#include <pthread.h>


#ifdef OPENMP // just a defense in case OpenMP is not installed.

#ifndef __APPLE__
#include <omp.h>
#endif

#endif

#define RAYTRACE_DEPTH          8
#define PHOTON_TRACE_DEPTH      5


#define EPSILON                     1e-12
#define TMAX                        400

#define PROB_DABSORB                0.5F
#define INDIRECT_PHOTON_NEEDED      500000      // 200000   // 500000
#define CAUSTICS_PHOTON_NEEDED      200000      // 50000    // 200000

#define NUM_SAMPLE_PER_LIGHT        1           // if I do so many times of raytracing, i dont need high number of samples

// Gaussian filter constants
#define E                               (2.718)
#define ALPHA                           (0.918)
#define BETA                            (1.953)
#define ONE_MINUS_E_TO_MINUS_BETA       (0.858)

#define TOTAL_ITERATION                 200  // 300
#define SMALL_NODE_GRANULARITY          128

#define PHOTON_QUERY_RADIUS             (0.000375)     // 0.000272
#define TWO_MUL_RADIUS                  (0.00075)
#define FOUR_BY_THREE                   (1.333333)
#define SPHERE_VOLUME                   (2.21e-10)

// Cone filter constants    jensen P67, k >= 1 is a filter constant characterizing the filter
#define CONE_K                          (1)

#define ENABLE_PATH_TRACING_GI          false

#define ENABLE_PHOTON_MAPPING           false
#define C_PHOTON_MODE                   1
#define SIMPLE_SMALL_NODE               1

#define ENABLE_DOF                      false
#define DOF_T                           9
#define DOF_R                           (0.5)
#define DOF_SAMPLE                      2


namespace _462 {
    
    // max number of threads OpenMP can use. Change this if you like.
#define MAX_THREADS 8
    
    static const unsigned STEP_SIZE = 4;
    
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
        num_iteration = 1;  //
        
        Ray::init(scene->camera);
        scene->initialize();
        
//        const SphereLight* lights = scene->get_lights();
        
        //----------------------------------------
        // initialize direction conversion tables, from Jensen's implementation
        //----------------------------------------
        for (int i = 0; i < 256; i++)
        {
            double angle = double(i) * (1.0/256.0) * M_PI;
            costheta[i] = cos( angle );
            sintheta[i] = sin( angle );
            cosphi[i]   = cos( 2.0 * angle );
            sinphi[i]   = sin( 2.0 * angle );
        }
        
        // generate bounding boxes that main contain objects that may cause caustics
//        generateCausticsBoxes();
        
        // Construction of BVH tree and bounding volumes
        int start_time = SDL_GetTicks();
        // Initialization or precompuation before the trace
        for (size_t i = 0; i < scene->num_geometries(); i++) {
            scene->get_geometries()[i]->createBoundingBox();
//            if (scene->get_geometries()[i]->type == eModel) {
//                Model *model = (Model *)scene->get_geometries()[i];
//                model->createBVHTree();
//            }
        }
        
        // initialize debug variables
        radius_clear = 0;
        radius_shadow = 0;
        clear_count = 0;
        shadow_count = 0;
        acc_pass_spent = 0;
        acc_kdtree_cons = 0;
        
        
        int end_time = SDL_GetTicks();
        int time_diff = end_time - start_time;
        printf("construct BVH time: %d ms\n", time_diff);
        
        
        raytraceColorBuffer = new Color3[width * height];
        for (size_t i = 0; i < width * height; i++) {
            raytraceColorBuffer[i] = Color3::Black();
        }
        
        // test:
        printf("c photon size = %ld, photon size = %ld, Vector3 size = %ld, Color size = %ld, TP size = %ld\n", sizeof(cPhoton),sizeof(Photon), sizeof(Vector3), sizeof(Color3), sizeof(unsigned char));
        // test end
        
#if ENABLE_PHOTON_MAPPING
            parallelPhotonScatter(scene);
    #if C_PHOTON_MODE
            cPhotonKDTreeConstruction();
    #else
            kdtreeConstruction();
    #endif
        
#endif
        
        
        
        // lixiao debug
//        std::thread trace_pixel_threads[4];
//        TracePixelData tracePixelData[4];
//        for (int i = 0; i < 4; i++)
//        {
//            tracePixelData[i].row_num = 0;
//            tracePixelData[i].buffer = new unsigned char(4 * width);
//            trace_pixel_threads[i] = std::thread(&Raytracer::tracePixelWorker, *this, &tracePixelData[i]);
//        }
        
        // initialize trace pixel data structure
        for (int i = 0; i < MAX_THREADS_TRACE; i++)
        {
            tracePixelData[i].row_num = -1;
            tracePixelData[i].buffer = new Color3[width];
//            tracePixelData[i].scene = scene;
        }
        
        return true;
    }
    
    // scatter photons
    void Raytracer::photonScatter(const Scene* scene)
    {
        photon_indirect_list.clear();
        photon_caustic_list.clear();
        
        // assuming there is only one light source
        SphereLight l = scene->get_lights()[0];
        
        unsigned int start = SDL_GetTicks();
        
        // scatter indirect
        while ((photon_indirect_list.size() < INDIRECT_PHOTON_NEEDED))
        {
            //            Vector3 p = samplePointOnUnitSphere();
            //            Ray photonRay = Ray(p + l.position, p);
            Ray photonRay = getPhotonEmissionRayFromLight(scene->get_lights()[0]);
            // Make color only the portion of light
            photonRay.photon.mask = 0;
            photonRay.photon.setColor(l.color);// * (real_t(1)/(real_t)(INDIRECT_PHOTON_NEEDED + CAUSTICS_PHOTON_NEEDED));
            photonTrace(photonRay, EPSILON, TMAX, PHOTON_TRACE_DEPTH);
            
        }

        
        unsigned int end = SDL_GetTicks();
        printf("Finished Scattering Indirect : %d ms\n", end - start);
        
        start = SDL_GetTicks();
        
        // scatter cautics
        while ((photon_caustic_list.size() < CAUSTICS_PHOTON_NEEDED))
        {
            //            Vector3 p = samplePointOnUnitSphere();
            //            Ray photonRay = Ray(p + l.position, p);
            Ray photonRay = getPhotonEmissionRayFromLight(scene->get_lights()[0]);
            //            Ray photonRay = getPhotonEmissionRayForCaustics(scene->get_lights()[0]);
            // Make color only the portion of light
            photonRay.photon.mask = 0;
            photonRay.photon.setColor(l.color);// * (real_t(1)/(real_t)(CAUSTICS_PHOTON_NEEDED + INDIRECT_PHOTON_NEEDED));
            photonTrace(photonRay, EPSILON, TMAX, PHOTON_TRACE_DEPTH);
        }

        
        end = SDL_GetTicks();
        printf("Finished Scattering Caustics : %d ms\n", end - start);
        
        printf("Finished Scattering, indirect num = %ld, caustic num = %ld\n", photon_indirect_list.size(), photon_caustic_list.size());
        
        std::vector<Photon> tmp_indirect(photon_indirect_list.size() + 1);
        std::vector<Photon> tmp_caustic(photon_caustic_list.size() + 1);
        
        start = SDL_GetTicks();
        balance(1, tmp_indirect, photon_indirect_list);
        end = SDL_GetTicks();
        printf("Finished Balancing Indirect : %d ms\n", end - start);
        
        start = SDL_GetTicks();
        balance(1, tmp_caustic, photon_caustic_list);
        end = SDL_GetTicks();
        printf("Finished Balancing Caustics : %d ms\n", end - start);

        
        kdtree_photon_indirect_list = tmp_indirect;
        kdtree_photon_caustic_list = tmp_caustic;
        
        printf("Finished Balancing!\n");
    }
    
    // worker node for parallel photon scatter
    void Raytracer::photonScatterWorker(PhotonScatterData *data)
    {
//        printf("running thread!\n");
//        PhotonScatterData localData = *data;
        
        unsigned int start = SDL_GetTicks();
        
        // scatter indirect
        while ((data->worker_photon_indirect.size() < data->indirect_needed))
        {
            Ray photonRay = getPhotonEmissionRayFromLight(data->light);
            // Make color only the portion of light
            photonRay.photon = Photon(data->light.color);
//            photonRay.photon.setColor();// * (real_t(1)/(real_t)(INDIRECT_PHOTON_NEEDED + CAUSTICS_PHOTON_NEEDED));
//            photonTrace(photonRay, EPSILON, TMAX, PHOTON_TRACE_DEPTH);
            localPhotonTrace(photonRay, EPSILON, TMAX, PHOTON_TRACE_DEPTH,
                             data->worker_photon_indirect,
                             data->worker_photon_caustics,
                             data->indirect_needed,
                             data->caustics_needed);
        }
        
        // scatter caustics
        while ((data->worker_photon_caustics.size() < data->caustics_needed))
        {
            Ray photonRay = getPhotonEmissionRayFromLight(data->light);
            // Make color only the portion of light
            photonRay.photon = Photon(data->light.color);
//            photonRay.photon.setColor(data->light.color);// * (real_t(1)/(real_t)(INDIRECT_PHOTON_NEEDED + CAUSTICS_PHOTON_NEEDED));
            //            photonTrace(photonRay, EPSILON, TMAX, PHOTON_TRACE_DEPTH);
            localPhotonTrace(photonRay, EPSILON, TMAX, PHOTON_TRACE_DEPTH,
                             data->worker_photon_indirect,
                             data->worker_photon_caustics,
                             data->indirect_needed,
                             data->caustics_needed);
        }
        
        unsigned int end = SDL_GetTicks();
        
        printf("thread finish: %dms\n", end - start);
    }
    
    void Raytracer::tracePixelWorker(TracePixelData *data)
    {
        // Actual trace
        for (size_t x = 0; x < width; x++)
        {
            Color3 color = trace_pixel(scene, x, data->row_num, width, height);
            data->buffer[x] = color;
        }
    }
    
    void Raytracer::parallelPhotonScatter(const Scene* scene)
    {
        // debug
        unsigned int start = SDL_GetTicks();
        unsigned int end;
        
        photon_indirect_list.clear();
        photon_caustic_list.clear();
        
        // assuming there is only one light source
        SphereLight l = scene->get_lights()[0];
        
        // create MAX_THREADS POSIX threads
        std::thread tid[MAX_THREADS_SCATTER];
        PhotonScatterData data[MAX_THREADS_SCATTER];
        
        for (int i = 0; i < MAX_THREADS_SCATTER; i++)
        {
            data[i].indirect_needed = INDIRECT_PHOTON_NEEDED/MAX_THREADS_SCATTER;
            data[i].caustics_needed = CAUSTICS_PHOTON_NEEDED/MAX_THREADS_SCATTER;
            data[i].light = l;
            
            tid[i] = std::thread(&Raytracer::photonScatterWorker, *this, &data[i]);
        }
        
        // wait until four threads finishes
        for (int i = 0; i < MAX_THREADS_SCATTER; i++)
        {
            tid[i].join();
        }
        
        // combine thread data to global raytracer data
        photon_indirect_list.reserve(INDIRECT_PHOTON_NEEDED);
        photon_caustic_list.reserve(CAUSTICS_PHOTON_NEEDED);
        
        for (int i = 0; i < MAX_THREADS_SCATTER; i++)
        {
//            printf("worker thread has indirect: %ld, caustics: %ld\n",data[i].worker_photon_indirect.size(), data[i].worker_photon_caustics.size());
            photon_indirect_list.insert( photon_indirect_list.end(), data[i].worker_photon_indirect.begin(), data[i].worker_photon_indirect.end() );
            photon_caustic_list.insert( photon_caustic_list.end(), data[i].worker_photon_caustics.begin(), data[i].worker_photon_caustics.end() );
        }
        
        end = SDL_GetTicks();
        printf("Finished Parallel Scattering, time = %d ms, indirect num = %ld, caustic num = %ld\n", end - start, photon_indirect_list.size(), photon_caustic_list.size());
    }
    
    void Raytracer::kdtreeConstruction()
    {
        unsigned int start = 0, end = 0, acc = 0;
        
        // balancing kdtree and copy back
        std::vector<Photon> tmp_indirect(photon_indirect_list.size() + 1);
        std::vector<Photon> tmp_caustic(photon_caustic_list.size() + 1);
        
        start = SDL_GetTicks();
        balance(1, tmp_indirect, photon_indirect_list);
        end = SDL_GetTicks();
        acc += end - start;
        printf("Finished Construct Indirect KD-Tree : %d ms\n", end - start);
        
        start = SDL_GetTicks();
        balance(1, tmp_caustic, photon_caustic_list);
        end = SDL_GetTicks();
        acc += end - start;
        printf("Finished Construct Caustics KD-Tree : %d ms\n", end - start);
        
        kdtree_photon_indirect_list = tmp_indirect;
        kdtree_photon_caustic_list = tmp_caustic;
        
        printf("Finished KD-Tree Construction! Total time = %d ms\n", acc);
    }
    
    // c style photon kdtree construction
    void Raytracer::cPhotonKDTreeConstruction()
    {
        unsigned int start = 0, end = 0;
        
        start = SDL_GetTicks();
        // clear old c style kdtree
//        kdtree_cphoton_indirect.clear();
//        kdtree_cphoton_caustics.clear();
        
        cphoton_indirect_data.clear();
        cphoton_caustics_data.clear();
        
        std::vector<cPhoton> c_indirect(photon_indirect_list.size());
        std::vector<cPhoton> c_caustics(photon_caustic_list.size());
        
        // copy from normal photon to c photon
        for (size_t i = 0; i < photon_indirect_list.size(); i++)
        {
            photon_indirect_list[i].position.to_array(c_indirect[i].position);
            c_indirect[i].index = i;
            c_indirect[i].splitAxis = -1;
        }
        
        for (size_t i = 0; i < photon_caustic_list.size(); i++)
        {
            photon_caustic_list[i].position.to_array(c_caustics[i].position);
            c_caustics[i].index = i;
            c_caustics[i].splitAxis = -1;
        }
        
        
        cphoton_indirect_data = c_indirect;
        cphoton_caustics_data = c_caustics;
        
        vvh_indirect_root = new KDNode;
        vvh_caustics_root = new KDNode;
        
        // construct vvh based kd-tree
        vvhKDTreeConstruction(cphoton_indirect_data, vvh_indirect_root);
        vvhKDTreeConstruction(cphoton_caustics_data, vvh_caustics_root);
        
        // construct ispc friendly cphoton data
        ispc_cphoton_indirect_data.size = cphoton_indirect_data.size();
        ispc_cphoton_indirect_data.posx = new float[cphoton_indirect_data.size()];
        ispc_cphoton_indirect_data.posy = new float[cphoton_indirect_data.size()];
        ispc_cphoton_indirect_data.posz = new float[cphoton_indirect_data.size()];
        
        ispc_cphoton_indirect_data.bitmap = new int8_t[cphoton_indirect_data.size()];
        
        ispc_cphoton_caustics_data.size = cphoton_caustics_data.size();
        ispc_cphoton_caustics_data.posx = new float[cphoton_caustics_data.size()];
        ispc_cphoton_caustics_data.posy = new float[cphoton_caustics_data.size()];
        ispc_cphoton_caustics_data.posz = new float[cphoton_caustics_data.size()];
        ispc_cphoton_caustics_data.bitmap = new int8_t[cphoton_caustics_data.size()];
        
        
        // TODO: ispc parallel
        for (size_t i = 0; i < cphoton_indirect_data.size(); i++)
        {
            ispc_cphoton_indirect_data.posx[i] = cphoton_indirect_data[i].position[0];
            ispc_cphoton_indirect_data.posy[i] = cphoton_indirect_data[i].position[1];
            ispc_cphoton_indirect_data.posz[i] = cphoton_indirect_data[i].position[2];
            ispc_cphoton_indirect_data.bitmap[i] = 0;
        }
        
        for (size_t i = 0; i < cphoton_caustics_data.size(); i++)
        {
            ispc_cphoton_caustics_data.posx[i] = cphoton_caustics_data[i].position[0];
            ispc_cphoton_caustics_data.posy[i] = cphoton_caustics_data[i].position[1];
            ispc_cphoton_caustics_data.posz[i] = cphoton_caustics_data[i].position[2];
            ispc_cphoton_caustics_data.bitmap[i] = 0;
        }

        end = SDL_GetTicks();
        acc_kdtree_cons += end - start;
        printf("c photon KDTree Construction Time = %d ms\n", end - start);
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
            
            
            
#if ENABLE_DOF
            Ray r = Ray(scene->camera.get_position(), Ray::get_pixel_dir(i, j));
            res += trace(r, EPSILON, TMAX, RAYTRACE_DEPTH);
            
            for (int i = 0; i < DOF_SAMPLE - 1; i++) {
                double random_r = DOF_R * random_uniform();
                double random_theta = 2 * PI * random_uniform();
                double random_x = random_r * cos(random_theta);
                double random_y = random_r * sin(random_theta);
                Vector3 focus   = r.e + r.d * DOF_T;
                Vector3 right   = cross(scene->camera.get_direction(), scene->camera.get_up());
                Vector3 cam     = scene->camera.get_position() + DOF_R * random_x * right + DOF_R * random_y * scene->camera.get_up();
                Ray sample_ray = Ray(cam, normalize(focus - cam));
                res += trace(sample_ray, EPSILON, TMAX, RAYTRACE_DEPTH);
            }
            
            res *= 1.0/(float)DOF_SAMPLE;
            
            
            
#else
            Ray r = Ray(scene->camera.get_position(), Ray::get_pixel_dir(i, j));
            res += trace(r, EPSILON, TMAX, RAYTRACE_DEPTH);
#endif
            
        }
        return res*(real_t(1)/real_t(num_samples));
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
        
//        static const size_t PRINT_INTERVAL = 64;
        
        // the time in milliseconds that we should stop
        unsigned int end_time = 0;
        bool is_done;
        
        if (max_time)
        {
            // convert duration to milliseconds
            unsigned int duration = (unsigned int) (*max_time * 1000);
            end_time = SDL_GetTicks() + duration;
//            printf("entime = %ld\n",end_time);
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
//#ifndef __APPLE__
//#pragma omp parallel for
//#endif
            // launch 4 threads for actually tracing pixel, one for each row
            // NOTE!!!! MAX_THREADS_TRACE must be same with step size
//            std::thread threads[MAX_THREADS_TRACE];
//            
//            // launch threads and start tracing
//            for (int c_row = current_row; c_row < loop_upper; c_row++)
//            {
//                int index = c_row % MAX_THREADS_TRACE;
//                tracePixelData[index].row_num = c_row;
//                threads[index] = std::thread(&Raytracer::tracePixelWorker, *this, &tracePixelData[index]);
//            }
//            
//            // when all threads hit barrier
//            for (int c_row = current_row; c_row < loop_upper; c_row++)
//            {
//                int index = c_row % MAX_THREADS_TRACE;
//                threads[index].join();
//            }
//            
//            // calculate color and update the buffer
//            for (int c_row = current_row; c_row < loop_upper; c_row++)
//            {
//                int index = c_row % MAX_THREADS_TRACE;
//                for (size_t x = 0; x < width; x++)
//                {
//                    Color3 color = tracePixelData[index].buffer[x];
//                    raytraceColorBuffer[(c_row * width + x)] += color;
//                    Color3 progressiveColor = raytraceColorBuffer[(c_row * width + x)] * ((1.0)/(num_iteration));
//                    progressiveColor.to_array(&buffer[4 * (c_row * width + x)]);
//                }
//            }
//            #pragma omp parallel for
            for (int c_row = current_row; c_row < loop_upper; c_row++)
            {
                
                /*
                 * This defines a critical region of code that should be
                 * executed sequentially.
                 */
//#ifndef __APPLE__
//#pragma omp critical
//#endif
//                {
////                    if (c_row % PRINT_INTERVAL == 0)
////                        printf("Raytracing (Row %d)\n", c_row);
//                }
                
                for (size_t x = 0; x < width; x++)
                {
                    // Xiao: debugging
                    // Measuring time
                    if (c_row == 0 && x == 0)
                    {
                        // master start
                        if (num_iteration == 1)
                        {
                            master_start = SDL_GetTicks();
                            
                        }
                        
                        pass_start = SDL_GetTicks();
                        
                        acc_cphoton_search_time = 0;
                        acc_iphoton_search_time = 0;
                    }
                    
                    // trace a pixel
                    Color3 color = trace_pixel(scene, x, c_row, width, height);
//                    applyGammaHDR(color);
                    raytraceColorBuffer[(c_row * width + x)] += color;
                    Color3 progressiveColor = raytraceColorBuffer[(c_row * width + x)] * ((1.0)/(num_iteration));
                    progressiveColor = clamp(progressiveColor, 0.0, 1.0);
                    
                    
                    // write the result to the buffer, always use 1.0 as the alpha
//                    color.to_array(&buffer[4 * (c_row * width + x)]);
                    progressiveColor.to_array(&buffer[4 * (c_row * width + x)]);
                }
            }
        }
        
        if (is_done)
        {
            if (num_iteration < TOTAL_ITERATION)
            {
                pass_end = SDL_GetTicks();
                acc_pass_spent += pass_end - pass_start;
                printf("Done One Pass! Iteration = %d, Pass spent = %dms\n", num_iteration, (pass_end - pass_start));
                printf("Pass photon search spent: indirect = %dms, caustics = %dms\n", acc_iphoton_search_time, acc_cphoton_search_time);
                // add postprocessing kernal to raytraceColorBuffer
                
                current_row = 0;
#if ENABLE_PHOTON_MAPPING
    #if C_PHOTON_MODE
                    parallelPhotonScatter(scene);

                    cPhotonKDTreeConstruction();
    #else
                    kdtreeConstruction();
    #endif
#endif
                
                is_done = false;
                num_iteration++;
                
                
            }
            else
            {
                pass_end = SDL_GetTicks();
                master_end = SDL_GetTicks();
                
                printf("Done Progressive Photon Mapping! Iteration = %d, Total Spent = %dms\n", num_iteration, master_end - master_start);
//                perPixelRender(buffer);
                
                // debug varibale update
                printf("average clear radius = %f, average shadow radius = %f, Average KDTree Construction %dms , Average Pass Spent = %dms\n", radius_clear/(float)clear_count, radius_shadow/(float)shadow_count, acc_kdtree_cons/num_iteration, acc_pass_spent/num_iteration);
            }
            
        }
        
        
        
//        if (is_done) {
//            printf("Done raytracing!");
//            // TODO: shading
//            
//            for (size_t i = 0; i < height; i++)
//            {
//                for (size_t j = 0; j < width; j++)
//                {
//                    Color3 color = raytraceColorBuffer[(i * width + j)] * (1.0/totalNum);
//                    color.to_array(&buffer[4 * (i * width + j)]);
//                }
//            }
//            
//        }
        
        return is_done;
    }
    
    void Raytracer::perPixelRender(unsigned char* buffer)
    {
        printf("Final Rendering\n");
        for (size_t y = 0; y < height; y++)
        {
            for (size_t x = 0; x < width; x++)
            {
                Color3 color = raytraceColorBuffer[(y * width + x)] * (1.0/(real_t)(TOTAL_ITERATION));
                color.to_array(&buffer[4 * (y * width + x)]);
            }
        }
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
        HitRecord record = getClosestHit(ray, t0, t1, &isHit);
        
        if (isHit && !record.isLight)
        {
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
                    refractRay.photon.setColor(refractRay.photon.getColor() * record.specular);
//                    refractRay.photon.color *= record.specular;
                    refractRay.photon.mask |= 0x2;
                    photonTrace(refractRay, t0, t1, depth-1);
                }
                
                // TODO: this generates a ring over the ceiling
//                if (isRefraction) {
//                    real_t R0 = ((idxRatio - real_t(1)) * (idxRatio - real_t(1))) / ((idxRatio + real_t(1)) * (idxRatio + real_t(1)));
//                    R = R0 + (real_t(1) - R0) * pow((real_t(1) - c), 5);
//                }
//                
//                real_t prob = random();
//                if (prob < R) {
//                    // create specular reflection for photon
//                    Vector3 reflectDirection = ray.d - 2 * dot(ray.d, record.normal) * record.normal;
//                    Ray reflectRay = Ray(record.position + EPSILON * reflectDirection, normalize(reflectDirection));
//                    reflectRay.photon = ray.photon;
//                    reflectRay.photon.mask |= 0x2;
////                    reflectRay.photon.color *= record.specular;
//                    reflectRay.photon.setColor(reflectRay.photon.getColor() * record.specular);
//                    photonTrace(reflectRay, t0, t1, depth-1);
//                    
//                }
//                else {
//                    // create refractive reflection for photon
//                    Ray refractRay = Ray(record.position + EPSILON * rd , normalize(rd));
//                    refractRay.photon = ray.photon;
//                    refractRay.photon.mask |= 0x2;
////                    refractRay.photon.color *= record.specular;
//                    refractRay.photon.setColor(refractRay.photon.getColor() * record.specular);
//                    photonTrace(refractRay, t0, t1, depth-1);
//                }

            }
            // specular reflective
            else if (record.specular != Color3::Black())
            {
                // Pure reflective surface
                if (record.diffuse == Color3::Black()) {

                    Ray reflectRay = Ray(record.position, normalize(ray.d - 2 * dot(ray.d, record.normal) * record.normal));
                    reflectRay.photon = ray.photon;
                    reflectRay.photon.setColor(reflectRay.photon.getColor() * record.specular);
//                    reflectRay.photon.color *= record.specular;
                    reflectRay.photon.mask |= 0x2;
                    photonTrace(reflectRay, t0, t1, depth - 1);
                    
                }
                // Hit on a surface that is both reflective and diffusive
                else
                {
                    real_t prob = random();
                    // Then there is a possibility of whether reflecting or absorbing
                    if (prob < 0.5) {
                        Ray reflectRay = Ray(record.position, normalize(ray.d - 2 * dot(ray.d, record.normal) * record.normal));
                        reflectRay.photon = ray.photon;
                        reflectRay.photon.setColor(reflectRay.photon.getColor() * record.specular);
//                        reflectRay.photon.color *= record.specular;
                        reflectRay.photon.mask |= 0x2;
                        photonTrace(reflectRay, t0, t1, depth - 1);
                    }
                    else {
                        // absorb
                        if (photon_indirect_list.size() < INDIRECT_PHOTON_NEEDED)
                        {
                            ray.photon.position = record.position;
//                            ray.photon.direction = -ray.d;
                            setPhotonDirection(ray.photon, -ray.d);
//                            ray.photon.setDirection(-ray.d);
//                            ray.photon.color = ray.photon.color;// * record.diffuse;
//                            ray.photon.color = ray.photon.color/
//                            ray.photon.normal = record.normal;
                            photon_indirect_list.push_back(ray.photon);
                        }

                    }
                }
            }
            // diffusive
            else {
                // direct illumination, do not store
                if (ray.photon.mask == 0x0) {
                    // consider don't do direct illumination
                    // if remove this, global photons could be faster but caustics are getting far slower
                    real_t prob = random();
                    if (prob > PROB_DABSORB) {
                        Ray photonRay = Ray(record.position, uniformSampleHemisphere(record.normal));
                        photonRay.photon = ray.photon;
                        photonRay.photon.mask |= 0x1;
//                        photonRay.photon.color = ray.photon.color * record.diffuse;
                        photonRay.photon.setColor(ray.photon.getColor() * record.diffuse);
//                        photonRay.photon.setColor(ray.photon.getColor()  * (real_t(1)/real_t(1.0 - PROB_DABSORB)));
                        
//                        ray.photon.color *= real_t(1)/real_t(1.0 - PROB_DABSORB);
//                        ray.photon.setColor(ray.photon.getColor() * (real_t(1)/real_t(1.0 - PROB_DABSORB)));
                        
                        photonTrace(photonRay, t0, t1, depth-1);
                    }
//                    else
//                    {
//                        // absorb
//                        if (photon_indirect_list.size() < INDIRECT_PHOTON_NEEDED)
//                        {
//                            ray.photon.position = record.position;
//                            ray.photon.direction = -ray.d;
//                            ray.photon.normal = record.normal;
//                            photon_indirect_list.push_back(ray.photon);
//                        }
//                    }
                }
                // caustics
                else if (ray.photon.mask == 0x2) {
//                    printf("nice mask!\n");
                    if (photon_caustic_list.size() < CAUSTICS_PHOTON_NEEDED) {
                        ray.photon.position = record.position;
//                        ray.photon.direction = -ray.d;
                        setPhotonDirection(ray.photon, -ray.d);
//                        ray.photon.normal = (record.normal);
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
//                            ray.photon.direction = (-ray.d);
                            setPhotonDirection(ray.photon, -ray.d);
//                            ray.photon.normal = (record.normal);
//                            ray.photon.color *= real_t(1)/real_t(PROB_DABSORB);
                            ray.photon.setColor(ray.photon.getColor() * (real_t(1)/real_t(PROB_DABSORB)));
                            photon_indirect_list.push_back(ray.photon);
                        }
                    }
                    else {
                        // Generate a diffusive reflect
                        Ray photonRay = Ray(record.position, uniformSampleHemisphere(record.normal));
                        photonRay.photon = ray.photon;
                        photonRay.photon.mask |= 0x1;
//                        photonRay.photon.color = (ray.photon.color * record.diffuse);
                        photonRay.photon.setColor(ray.photon.getColor() * record.diffuse);
//                        ray.photon.color *= real_t(1)/real_t(1.0 - PROB_DABSORB);
                        ray.photon.setColor(ray.photon.getColor() * (real_t(1)/real_t(1.0 - PROB_DABSORB)));
                        photonTrace(photonRay, t0, t1, depth-1);
                    }
                }
            }
        }
    }
    
    // Photon Tracing function for a single thread purposes
    void Raytracer::localPhotonTrace(Ray ray, real_t t0, real_t t1, int depth,
                                     std::vector<Photon> &indirect_list, std::vector<Photon> &caustics_list,
                                     size_t indirect_needed, size_t caustics_needed)
    {
        if (depth == 0) {
            return;
        }
        
        bool isHit = false;
        HitRecord record = getClosestHit(ray, t0, t1, &isHit);
        
        if (isHit && !record.isLight)
        {
            // refractive
            if (record.refractive_index > 0) {
                
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
                    refractRay.photon.setColor(refractRay.photon.getColor() * record.specular);
                    //                    refractRay.photon.color *= record.specular;
                    refractRay.photon.mask |= 0x2;
                    localPhotonTrace(refractRay, t0, t1, depth-1, indirect_list, caustics_list, indirect_needed, caustics_needed);
                }
                
            }
            // specular reflective
            else if (record.specular != Color3::Black())
            {
                // Pure reflective surface
                if (record.diffuse == Color3::Black()) {
                    
                    Ray reflectRay = Ray(record.position, normalize(ray.d - 2 * dot(ray.d, record.normal) * record.normal));
                    reflectRay.photon = ray.photon;
                    reflectRay.photon.setColor(reflectRay.photon.getColor() * record.specular);
                    //                    reflectRay.photon.color *= record.specular;
                    reflectRay.photon.mask |= 0x2;
                    localPhotonTrace(reflectRay, t0, t1, depth-1, indirect_list, caustics_list, indirect_needed, caustics_needed);
                    
                }
                // Hit on a surface that is both reflective and diffusive
                else
                {
                    real_t prob = random();
                    // Then there is a possibility of whether reflecting or absorbing
                    if (prob < 0.5) {
                        Ray reflectRay = Ray(record.position, normalize(ray.d - 2 * dot(ray.d, record.normal) * record.normal));
                        reflectRay.photon = ray.photon;
                        reflectRay.photon.setColor(reflectRay.photon.getColor() * record.specular);
                        //                        reflectRay.photon.color *= record.specular;
                        reflectRay.photon.mask |= 0x2;
                        localPhotonTrace(reflectRay, t0, t1, depth - 1, indirect_list, caustics_list, indirect_needed, caustics_needed);
                    }
                    else {
                        // absorb
                        if (indirect_list.size() < indirect_needed)
                        {
                            ray.photon.position = record.position;
                            setPhotonDirection(ray.photon, -ray.d);
                            indirect_list.push_back(ray.photon);
                        }
                        
                    }
                }
            }
            // diffusive
            else {
                // direct illumination, do not store
                if (ray.photon.mask == 0x0) {
                    // consider don't do direct illumination
                    // if remove this, global photons could be faster but caustics are getting far slower
                    real_t prob = random();
                    if (prob > PROB_DABSORB) {
                        Ray photonRay = Ray(record.position, uniformSampleHemisphere(record.normal));
                        photonRay.photon = ray.photon;
                        photonRay.photon.mask |= 0x1;
                        photonRay.photon.setColor(ray.photon.getColor() * record.diffuse);
 
                        localPhotonTrace(photonRay, t0, t1, depth-1, indirect_list, caustics_list, indirect_needed, caustics_needed);
                    }
                }
                // caustics
                else if (ray.photon.mask == 0x2) {
                    if (caustics_list.size() < caustics_needed) {
                        ray.photon.position = record.position;
                        setPhotonDirection(ray.photon, -ray.d);
                        caustics_list.push_back(ray.photon);
                        
                    }
                }
                // indirect illumination
                else {
                    real_t prob = random();
                    if (prob < PROB_DABSORB) {
                        // Store photon in indirect illumination map
                        if (indirect_list.size() < indirect_needed) {
                            ray.photon.position = (record.position);
                            setPhotonDirection(ray.photon, -ray.d);
                            ray.photon.setColor(ray.photon.getColor() * (real_t(1)/real_t(PROB_DABSORB)));
                            indirect_list.push_back(ray.photon);
                        }
                    }
                    else {
                        // Generate a diffusive reflect
                        Ray photonRay = Ray(record.position, uniformSampleHemisphere(record.normal));
                        photonRay.photon = ray.photon;
                        photonRay.photon.mask |= 0x1;
                        photonRay.photon.setColor(ray.photon.getColor() * record.diffuse);
//                        ray.photon.setColor(ray.photon.getColor() * (real_t(1)/real_t(1.0 - PROB_DABSORB)));
                        localPhotonTrace(photonRay, t0, t1, depth-1, indirect_list, caustics_list, indirect_needed, caustics_needed);
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
        bool isHit = false;
        HitRecord record = getClosestHit(ray, t0, t1, &isHit);
        
        if (isHit && !record.isLight) {
            return shade(ray, record, t0, t1, depth);
        }
        else if (isHit && record.isLight) {
            // TODO: lixiao debug force color to be white...
//            return scene->get_lights()[0].color * Color3(100, 100, 100);
//            return scene->get_lights()[0].color;
            return Color3(10,10,10);
        }
        else
            return scene->background_color;
        
    }
    
    // trace conterparts that enables DOF
    Color3 Raytracer::trace_DOF(HitRecord &record, Ray ray, real_t t0, real_t t1, int depth)
    {
        bool isHit = false;
        record = getClosestHit(ray, t0, t1, &isHit);
        
        if (isHit && !record.isLight) {
            return shade(ray, record, t0, t1, depth);
        }
        else if (isHit && record.isLight) {
            return Color3(10,10,10);
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
            if (record.specular != Color3::Black())
            {
                reflColor = record.specular * trace(reflectRay, t0, t1, depth - 1);
            }
            
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
        }
        
        // 4/14/2014, wrong assumption for direct illumination
//        else {
        
        // ray hits texture map
        if (record.diffuse == Color3(-1.0, -1.0, -1.0))
        {
            Color3 hdrColor = record.texture;
            hdrColor.r = hdrColor.r > 0.99 ? hdrColor.r * 3 : hdrColor.r;
            hdrColor.g = hdrColor.g > 0.99 ? hdrColor.g * 3 : hdrColor.g;
            hdrColor.b = hdrColor.b > 0.99 ? hdrColor.b * 3 : hdrColor.b;
            record.texture = hdrColor;
            radiance += hdrColor;
        }
        else if (record.diffuse != Color3::Black() && record.refractive_index == 0) {
            
            // 6/27/2014, path tracing test
            // generate a new ray, randomized along normal sphere
            
            // Trace each light source for direct illumination
            radiance += shade_direct_illumination(record, t0, t1);
            
#if ENABLE_PATH_TRACING_GI
            Vector3 dir = uniformSampleHemisphere(record.normal);
            Ray pathRay = Ray(record.position, dir);
            radiance += 0.6 * trace(pathRay, t0, t1, depth - 1);
#else
#endif
            
            
            
//            radiance += shade_direct_illumination(record, t0, t1) + 0.6 * trace(pathRay, t0, t1, depth - 1);
            
            // remove ambient light
//            radiance += scene->ambient_light * record.ambient;
            
            // add coeefficient so that darker part are brighter and brighter part get darker
            
            
//            if (radiance.r <= 0.1 && radiance.g <= 0.1 && radiance.b <= 0.1)
//                coeef *= 2;
            
            // Shading caustics
//            radiance += shade_caustics(record, 0.0001 * 1, INFINITY) * (1.0/(CAUSTICS_PHOTON_NEEDED)) * coeef; // 0.0001
            
            // Shading global illumination
//            radiance += shade_indirect_illumination(record, 0.0001 * 1, INFINITY) * (1.0/(INDIRECT_PHOTON_NEEDED)) * coeef;   // 0.0001
//            radiance += shade_photons(record, 0.0001, INFINITY) * (2.0/(CAUSTICS_PHOTON_NEEDED + INDIRECT_PHOTON_NEEDED)) * coeef;
            
            // normal
#if ENABLE_PHOTON_MAPPING
                int coeef = 25;
    #if C_PHOTON_MODE
                if (CAUSTICS_PHOTON_NEEDED + INDIRECT_PHOTON_NEEDED > 0)
                {
                    Color3 photonRadiance = shade_cphotons(record, PHOTON_QUERY_RADIUS, 0.001) * (2.0/(CAUSTICS_PHOTON_NEEDED + INDIRECT_PHOTON_NEEDED)) * coeef;
                    photonRadiance = clamp(photonRadiance, 0, 1.0);
                    radiance += photonRadiance;
                    radiance = clamp(radiance, 0, 1.0);
                }
                
    #else
                if (CAUSTICS_PHOTON_NEEDED + INDIRECT_PHOTON_NEEDED > 0)
                {
                    Color3 photonRadiance = shade_photons(record, PHOTON_QUERY_RADIUS, INFINITY) * (2.0/(CAUSTICS_PHOTON_NEEDED + INDIRECT_PHOTON_NEEDED)) * coeef;
                    photonRadiance = clamp(photonRadiance, 0, 1.0);
    //                if (photonRadiance.r > 1.0 && photonRadiance.g > 1.0 && photonRadiance.b > 1.0) {
    //                    
    //                    std::cout<<photonRadiance<<std::endl;
    //                }
                    radiance += photonRadiance;
                    
                    radiance = clamp(radiance, 0, 1.0);
                    
    //                if (radiance.r > 1.0 || radiance.g > 1.0 || radiance.b > 1.0) {
    //                    
    //                    std::cout<<radiance<<std::endl;
    //                }

                }
    #endif
#endif
            
        }
        
//            real_t progressiveCoef = real_t(TOTAL_ITERATION - num_iteration + 1)/real_t(TOTAL_ITERATION);
//            progressiveCoef *= progressiveCoef;
//            printf(" === %f\n",progressiveCoef);
            
            // Shading caustics
//            radiance += shade_caustics(record, 0.0001 * 1, INFINITY) * (1.0/(CAUSTICS_PHOTON_NEEDED)) * 100; // 0.0001
            
            // Shading global illumination
//            radiance += shade_indirect_illumination(record, 0.0001 * 1, INFINITY) * (1.0/(INDIRECT_PHOTON_NEEDED)) * 100;   // 0.0001
//        }
        
        // for raytracing ambient
        // turn amient on if needed
//        if ( record.refractive_index == 0) {
//            radiance += scene->ambient_light * record.ambient;
//        }
        
        // TODO:
        // try to add high light with blinn-phong model
        // the higher the phong, the smaller the high light
        if (record.phong > 0) {
            SphereLight light = scene->get_lights()[0];
            Vector3 lightpos = sampleLightSource(light);
            Vector3 l = normalize(lightpos - record.position);
            Ray rayToLight = Ray(record.position, l);
            bool isHit;
            getClosestHit(rayToLight, EPSILON, TMAX, &isHit);
            Color3 highLightColor = light.color * 64;  // TODO: light intensity
            if (isHit) {
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
        
        
        for (size_t i = 0; i < scene->num_lights(); i++) {
            
            // make random points on a light, this would make egdes of shadows smoother
            size_t sample_num_per_light = NUM_SAMPLE_PER_LIGHT;    // todo: add samples
            for (size_t j = 0; j < sample_num_per_light; j++)
            {
                // Get a random point on light sphere
                Vector3 sampleLightPos = sampleLightSource(scene->get_lights()[i]);
                
                // scene->get_lights()[i].position; // For hard shadows
                
                Vector3 d_shadowRay = sampleLightPos - record.position;
                Vector3 d_shadowRay_normolized = normalize(d_shadowRay);
                Ray shadowRay = Ray(record.position + EPSILON * d_shadowRay_normolized, d_shadowRay_normolized);
//                    real_t tl = d_shadowRay.x / d_shadowRay_normolized.x;
//                    real_t tl = length(d_shadowRay) / length(d_shadowRay_normolized);
                
                bool isHit = false;
//                    tl = tl < t1 ? tl : t1;
                HitRecord shadowRecord = getClosestHit(shadowRay, t0, t1, &isHit);
                
                // only when hits the light before hits anything else
                if (shadowRecord.isLight)
                {
                    // if the obejct is not in shadow and is opaque
                    if (!isHit && (record.refractive_index == 0))
                    {
                        res += (record.getPixelColorFromLight(scene->get_lights()[i])) * (real_t(1)/real_t(sample_num_per_light));
                    }
                    else
                    {
                        // in shadow
                        record.isInShadow = true;
                    }
                    
                    
                    if (record.refractive_index == 0)
                    {
                        //res = Color3::White();
                        // if the light is a sphere light source
                        if (scene->get_lights()[i].type == 1)
                        {
                            res += (record.getPixelColorFromLight(scene->get_lights()[i])) * (real_t(1)/real_t(sample_num_per_light));
                        }
                        else if (scene->get_lights()[i].type == 2)
                        {
                            // test the light normal
                            if (dot(scene->get_lights()[i].normal, d_shadowRay_normolized) < 0)
                            {
                                res += (record.getPixelColorFromLight(scene->get_lights()[i])) * (real_t(1)/real_t(sample_num_per_light));
                            }
                            else
                            {
                                // in shadow
                                record.isInShadow = true;
                            }
                        }
                        
                        
                    }
                    else
                    {
                        // in shadow
//                        record.isInShadow = true;
                    }
                }
                
                
            }
            
        }
        
        return res;
    }
    
    // Shading of caustics
    Color3 Raytracer::shade_caustics(HitRecord &record, real_t radius, size_t num_samples)
    {
        Color3 causticsColor = Color3::Black();
        std::vector<Photon> nearestCausticPhotons;
        
        if (record.refractive_index == 0)
        {
            // sample radius
//            real_t sampleSquaredRadiusCaustics = 0.04;   // 0.1 // 0.001 // 0.0225
//            real_t maxSquaredDistCaustics = 0.001;      // 0.001 // 0.001
            size_t maxPhotonsEstimate = num_samples;
            
            // gather samples
            // TODO: doing nearest neighbor search with KD Tree
            float maxSearchSquaredRadius = radius * 0.9;
            while (nearestCausticPhotons.size() == 0) {
                maxSearchSquaredRadius /= 0.9;
                locatePhotons(1, record.position, kdtree_photon_caustic_list, nearestCausticPhotons, maxSearchSquaredRadius, maxPhotonsEstimate);
            }
            
            
            // calculate radiance
            for (size_t i = 0; i < nearestCausticPhotons.size(); i++)
            {
                causticsColor +=
                record.getPhotonLambertianColor(getPhotonDirection(nearestCausticPhotons[i]), nearestCausticPhotons[i].getColor());
//                record.getPhotonLambertianColor(nearestCausticPhotons[i].direction, nearestCausticPhotons[i].getColor())
//                * getConeFilterWeight(sqrt(nearestCausticPhotons[i].squaredDistance), sqrt(sampleSquaredRadiusCaustics));
                //* getGaussianFilterWeight(nearestCausticPhotons[i].squaredDistance, maxPhotonsEstimate);
            }
            
            // color/= PI*r^2
//            causticsColor *= (real_t(1)/((real_t(1) - real_t(2)/(real_t(3) * CONE_K)) * (PI * maxSquaredDistCaustics))) * (1.0/(CAUSTICS_PHOTON_NEEDED));
            causticsColor = causticsColor * (real_t(1)/(PI * maxSearchSquaredRadius));// * 0.0001;// * (1.0/(CAUSTICS_PHOTON_NEEDED));
        }
        
        nearestCausticPhotons.clear();
        
        return causticsColor;
        
    }
    
    // Shading of Indirect lightning
    Color3 Raytracer::shade_indirect_illumination(HitRecord &record, real_t radius, size_t num_samples)
    {
        Color3 indirectColor = Color3::Black();
        std::vector<Photon> nearestIndirectPhotons;
        
        if (record.refractive_index == 0)
        {
//            real_t sampleSquaredRadiusIndirect = 0.5;       // 0.5
//            real_t maxSquaredDistIndirect = 0.001;          // 0.001
            size_t maxPhotonsEstimate = num_samples;
            float maxSearchSquaredRadius = radius*0.9;
            
            // this value is used for search pruning
            
            while (nearestIndirectPhotons.size() == 0)
            {
                maxSearchSquaredRadius /= 0.9;
                locatePhotons(1, record.position, kdtree_photon_indirect_list, nearestIndirectPhotons, maxSearchSquaredRadius, maxPhotonsEstimate);
            }
            
            // cone filter to gather radiance of indirect illumination photons
            for (size_t i = 0; i < nearestIndirectPhotons.size(); i++)
            {
                indirectColor +=
                record.getPhotonLambertianColor(getPhotonDirection(nearestIndirectPhotons[i]), nearestIndirectPhotons[i].getColor());
//                record.getPhotonLambertianColor(nearestIndirectPhotons[i].direction, nearestIndirectPhotons[i].getColor())
                //* getConeFilterWeight(sqrt(nearestIndirectPhotons[i].squaredDistance), sqrt(sampleSquaredRadiusIndirect));
                //* getGaussianFilterWeight(nearestIndirectPhotons[i].squaredDistance, maxSearchSquaredRadius);
            }
            //indirectColor *= (real_t(1)/((real_t(1) - real_t(2)/(real_t(3) * CONE_K)) * (PI * maxSquaredDistIndirect)));
            indirectColor *= (real_t(1)/(PI * maxSearchSquaredRadius));// * 0.0005;// * (1.0/(INDIRECT_PHOTON_NEEDED));
            
        }
        
        return indirectColor;
    }
    
    Color3 Raytracer::shade_photons(HitRecord &record, real_t radius, size_t num_samples)
    {
        Color3 color = Color3::Black();
        std::vector<Photon> nearestPhotons;
        
        if (record.refractive_index == 0)
        {
            size_t maxPhotonsEstimate = num_samples;
            float maxSearchSquaredRadius = radius * 0.9;
            
//            locatePhotons(1, record.position, kdtree_photon_indirect_list, nearestPhotons, maxSearchSquaredRadius, maxPhotonsEstimate);
//            locatePhotons(1, record.position, kdtree_photon_caustic_list, nearestPhotons, maxSearchSquaredRadius, maxPhotonsEstimate);
//            
//            if (!nearestPhotons.size()) {
//                return color;
            
            unsigned int start, end;
            while (nearestPhotons.size() == 0)
            {
                maxSearchSquaredRadius /= 0.9;
                start = SDL_GetTicks();
                locatePhotons(1, record.position, kdtree_photon_indirect_list, nearestPhotons, maxSearchSquaredRadius, maxPhotonsEstimate);
                end = SDL_GetTicks();
                acc_iphoton_search_time += end - start;
                
                start = SDL_GetTicks();
                locatePhotons(1, record.position, kdtree_photon_caustic_list, nearestPhotons, maxSearchSquaredRadius, maxPhotonsEstimate);
                 end = SDL_GetTicks();
                acc_cphoton_search_time += end - start;
            }
            
            // update debug varibales
            if (record.isInShadow)
            {
                radius_shadow += maxSearchSquaredRadius;
                shadow_count += 1;
            }
            else
            {
                radius_clear += maxSearchSquaredRadius;
                clear_count += 1;
            }
            
            // shade
//            printf("size = %ld\n",nearestPhotons.size());
            for (size_t i = 0; i < nearestPhotons.size(); i++)
            {
//                std::cout<<nearestPhotons[i].getColor()<<std::endl;
                color +=
                record.getPhotonLambertianColor(getPhotonDirection(nearestPhotons[i]), nearestPhotons[i].getColor());
            }
            
            color *= (real_t(1)/(PI * maxSearchSquaredRadius));
        }
        
        return color;
    }
    
    Color3 Raytracer::shade_cphotons(HitRecord &record, real_t radius, size_t num_samples)
    {
        Color3 color = Color3::Black();
        std::vector<int> nearestPhotonIndices;
        std::vector<Photon> sourcePhotons;
        
        if (record.refractive_index == 0)
        {
            size_t maxPhotonsEstimate = num_samples;
            float maxSearchSquaredRadius = radius;// * 0.9;
            
            unsigned int start, end;
            if (sourcePhotons.size() == 0)
            {
//                maxSearchSquaredRadius /= 0.9;
                start = SDL_GetTicks();
                vvhcPhotonLocate(record.position,
                                 vvh_indirect_root,
                                 ispc_cphoton_indirect_data,
                                 cphoton_indirect_data,
                                 nearestPhotonIndices,
                                 maxSearchSquaredRadius,
                                 maxPhotonsEstimate);
                
                // convert cphoton to normal photon
                for (size_t i = 0; i < nearestPhotonIndices.size(); i++)
                {
                    sourcePhotons.push_back(photon_indirect_list[nearestPhotonIndices[i]]);
                }
                nearestPhotonIndices.clear();
                end = SDL_GetTicks();
                acc_iphoton_search_time += end - start;
                
                start = SDL_GetTicks();
                vvhcPhotonLocate(record.position,
                                 vvh_caustics_root,
                                 ispc_cphoton_caustics_data,
                                 cphoton_caustics_data,
                                 nearestPhotonIndices,
                                 maxSearchSquaredRadius,
                                 maxPhotonsEstimate);
                
                // convert cphoton to normal photon
                for (size_t i = 0; i < nearestPhotonIndices.size(); i++)
                {
                    sourcePhotons.push_back(photon_caustic_list[nearestPhotonIndices[i]]);
                }
                nearestPhotonIndices.clear();
                end = SDL_GetTicks();
                acc_cphoton_search_time += end - start;
            }
            
            // update debug varibales
            if (record.isInShadow)
            {
                radius_shadow += maxSearchSquaredRadius;
                shadow_count += 1;
            }
            else
            {
                radius_clear += maxSearchSquaredRadius;
                clear_count += 1;
            }
            
//            printf("source photon = %ld\n",sourcePhotons.size());
            // shade
            for (size_t i = 0; i < sourcePhotons.size(); i++)
            {
//                std::cout<<sourcePhotons[i].getColor()<<std::endl;
                color +=
                record.getPhotonLambertianColor(getPhotonDirection(sourcePhotons[i]), sourcePhotons[i].getColor());
            }
            
            color *= (real_t(1)/(PI * maxSearchSquaredRadius));
        }
        
        return color;
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
            real_t b = _462::random_uniform() * real_t(2) - real_t(1);
            real_t c = _462::random_uniform() * real_t(2) - real_t(1);//_462::random_uniform();
//            real_t h = random() + EPSILON;
//            real_t d = _462::random_uniform();
            
            Vector3 vec1 = light.vertex1 - light.position;
            Vector3 vec2 = light.vertex2 - light.position;
//            Vector3 vec3 = normalize(cross(vec2, vec1));
            
            ran = light.position + (vec1) * b + (vec2) * c ;// + vec3 * d;
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
//            std::cout<<"light pos = "<<light.position<<std::endl;
//            std::cout<<"light nom = "<<light.normal<<std::endl;
            
//            Vector3 p = samplePointOnUnitSphere();
//            photonRay = Ray(p + light.position, p);
            
            // <-----
//            Vector3 vec1 = light.vertex1 - light.position;
//            Vector3 vec2 = light.vertex2 - light.position;
//            real_t r1 = random() * 2 - 1;
//            real_t r2 = random() * 2 - 1;
//            vec1 *= r1;
//            vec2 *= r2;
//            
//            Vector3 lpos = light.position + vec1 + vec2;
//            Vector3 lpos = sampleLightSource(light);
            Vector3 lpos = light.position;
            // revert! xiaoxiao debug
            // ------>
            
            Vector3 normal = light.normal;
            Vector3 d = uniformSampleHemisphere(normal);
//            Vector3 p = samplePointOnUnitSphere();//sampleLightSource(light);
            
//            real_t u = random();
//            real_t v = 2 * PI * random();
//            Vector3 d = -Vector3(cos(v) * sqrt(u), sin(v) * sqrt(u), sqrt(1 - u));
//            std::cout<<light.position<<std::endl;
            photonRay = Ray(lpos + d * EPSILON, d);
            
        }
        else
        {
            printf("error in light type!\n");
        }
        
        return photonRay;
    }
    
    // for higher efficiency of photon tracing pass
    Ray Raytracer::getPhotonEmissionRayForCaustics(SphereLight light)
    {
        // select a random object
        size_t index = rand() % cboxes.size();
        Box abox = cboxes[index];
        
        Vector3 endpos = abox.bounds[0];
        Vector3 startpos = Vector3::Zero();
        Vector3 dir = Vector3::Zero();
        
        Vector3 dim = abox.bounds[1] - abox.bounds[0];
        
        endpos = abox.bounds[0] + 0.5 * dim;
        Vector3 pointOnSphere = samplePointOnUnitSphere() * 0.5 * dim.x;
        endpos + pointOnSphere;
        
        if (light.type == 1)
        {
            Vector3 p = samplePointOnUnitSphere();
            startpos = p + light.position;
            dir = normalize(endpos - startpos);
        }
        // square light
        else if (light.type == 2)
        {
            Vector3 lpos = sampleLightSource(light);
            Vector3 normal = light.normal;
            Vector3 d = uniformSampleHemisphere(normal);
            
            startpos = lpos + d * EPSILON;
            dir = normalize(endpos - startpos);
            
        }
        
        return Ray(startpos, dir);
    }
    
    
    void Raytracer::generateCausticsBoxes()
    {
        // find objects that could create caustics:
        for (size_t i = 0; i < scene->num_geometries(); i++)
        {
            if (scene->get_geometries()[i]->type == eSphere)
            {
                Sphere *sphere = (Sphere *)scene->get_geometries()[i];
                if ( (sphere->material->refractive_index > 0) || (sphere->material->specular != Color3::Black()) )
                {
                    Vector3 pos = sphere->position;
                    real_t radius = sphere->radius;
                    Box aBox = Box(Vector3(pos.x - radius, pos.y-radius, pos.z-radius), Vector3(pos.x + radius, pos.y+radius, pos.z+radius));
                    cboxes.push_back(aBox);
                }
            }
            
//            if (scene->get_geometries()[i]->type == eTriangle)
//            {
//                Triangle *triangle = (Triangle *)scene->get_geometries()[i];
//                for (int j = 0; j < 3; j++)
//                {
//                    if ((triangle->vertices[j].material->refractive_index > 0) || (triangle->vertices[j].material->specular != Color3::Black()))
//                    {
//                        Box aBox = *triangle->boundingBox;
//                        cboxes.push_back(aBox);
//                    }
//                }
//            }
        }
    }
    
    
    /**
     * Retrieve the closest hit record
     * @param r             incoming ray
     * @param t0            lower limit of t
     * @param t1            upper limit of t
     * @param *isHit        indicate if there is a hit incident
     * @return HitRecord    the closest hit record
     */
    HitRecord Raytracer::getClosestHit(Ray r, real_t t0, real_t t1, bool *isHit)
    {
        HitRecord closestHitRecord;
        HitRecord tmp;
        real_t t = t1;
        
        Geometry* const* geometries = scene->get_geometries();
        *isHit = false;
        for (size_t i = 0; i < scene->num_geometries(); i++)
        {
            if (geometries[i]->hit(r, t0, t1, tmp))
            {
                if (!*isHit) {
                    *isHit = true;
                    t = tmp.t;
                    closestHitRecord = tmp;
                    
                }
                else {
                    if (tmp.t < t) {
                        t = tmp.t;
                        closestHitRecord = tmp;
                        
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
    
    Vector3 Raytracer::samplePointOnUnitSphereUniform()
    {
        real_t x = _462::random_uniform() * 2 - 1;
        real_t y = _462::random_uniform() * 2 - 1;
        real_t z = _462::random_uniform() * 2 - 1;
        
        Vector3 ran = Vector3(x, y, z);
        return normalize(ran);
    }
    
    Vector3 Raytracer::uniformSampleHemisphere(const Vector3& normal)
    {
        Vector3 newDir = samplePointOnUnitSphereUniform();
        if (dot(newDir, normal) < 0.0) {
            newDir = -newDir;
        }
//        return normalize(newDir);
        return newDir;
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
    
    void Raytracer::balance(size_t index, std::vector<Photon> &balancedKDTree, std::vector<Photon> &list)
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
        // O(N)
        Vector3 max = Vector3(-INFINITY, -INFINITY, -INFINITY);
        Vector3 min = Vector3(INFINITY, INFINITY, INFINITY);
        for (std::vector<Photon>::iterator it = list.begin(); it != list.end(); it++)
        {
            // calculate box
            max.x = it->position.x >= max.x ? it->position.x : max.x;
            max.y = it->position.y >= max.y ? it->position.y : max.y;
            max.z = it->position.z >= max.z ? it->position.z : max.z;
            
            min.x = it->position.x <= min.x ? it->position.x : min.x;
            min.y = it->position.y <= min.y ? it->position.y : min.y;
            min.z = it->position.z <= min.z ? it->position.z : min.z;
        }
        
        char splitAxis = -1;
        Vector3 diff = Vector3(max.x - min.x, max.y - min.y, max.z - min.z);
        if ((diff.x >= diff.y) && (diff.x >= diff.z))
            splitAxis = 0;
        else if ((diff.y >= diff.x) && (diff.y >= diff.z))
            splitAxis = 1;
        else if ((diff.z >= diff.x) && (diff.z >= diff.y))
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
        
        // O(NlogN)
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
    
    
    // TODO: change to c style kdtree, and use ispc for parallelism
    void Raytracer::cPhotonBalance(size_t index, std::vector<cPhoton> &balancedKDTree, std::vector<cPhoton> &list)
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
        // O(N)
        Vector3 max = Vector3(-INFINITY, -INFINITY, -INFINITY);
        Vector3 min = Vector3(INFINITY, INFINITY, INFINITY);
        // TODO: parallel
        for (std::vector<cPhoton>::iterator it = list.begin(); it != list.end(); it++)
        {
            // calculate box
            max.x = it->position[0] >= max.x ? it->position[0] : max.x;
            max.y = it->position[1] >= max.y ? it->position[1] : max.y;
            max.z = it->position[2] >= max.z ? it->position[2] : max.z;
            
            min.x = it->position[0] <= min.x ? it->position[0] : min.x;
            min.y = it->position[1] <= min.y ? it->position[1] : min.y;
            min.z = it->position[2] <= min.z ? it->position[2] : min.z;
        }
        
        char splitAxis = -1;
        Vector3 diff = Vector3(max.x - min.x, max.y - min.y, max.z - min.z);
        if ((diff.x >= diff.y) && (diff.x >= diff.z))
            splitAxis = 0;
        else if ((diff.y >= diff.x) && (diff.y >= diff.z))
            splitAxis = 1;
        else if ((diff.z >= diff.x) && (diff.z >= diff.y))
            splitAxis = 2;
        
        // Sorting the vector
        bool (*comparator)(const cPhoton &a, const cPhoton &b) = NULL;
        
        switch (splitAxis) {
            case 0:
                comparator = cPhotonComparatorX;
                break;
            case 1:
                comparator = cPhotonComparatorY;
                break;
            case 2:
                comparator = cPhotonComparatorZ;
                break;
                
            default:
                break;
        }
        
        // O(NlogN)
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
        cPhoton p = list[median];
        p.splitAxis = splitAxis;
        assert(index < balancedKDTree.size());
        balancedKDTree[index] = p;
        
//        printf("LT = %d, RT = %d\n", LT, RT);
        
        if (LT > 0) {
            std::vector<cPhoton> leftList(list.begin(), list.begin() + median);
            cPhotonBalance(2 * index, balancedKDTree, leftList);
        }
        
        if (RT > 0) {
            std::vector<cPhoton> rightList(list.begin() + median + 1, list.end());
            cPhotonBalance(2 * index + 1, balancedKDTree, rightList);
        }
        
        list.clear();
        
    }
    
    /**
     @brief locate the nearest photon list of given hit position
     @param p               balancedKDTree index, starts search from root node where p = 1
     @param position        hit position
     @param balancedKDTree  balanced kd tree, stored as a left-balanced kd tree, vector
     @param nearestPhotons  nearest photons found at hit position 
     @param sqrDist         smallest squared distance
     @param maxNum          maximum number of photons needed
     */
    void Raytracer::locatePhotons(size_t p,
                                  Vector3 position,
                                  std::vector<Photon> &balancedKDTree,
                                  std::vector<Photon> &nearestPhotons,
                                  float &sqrDist,
                                  size_t maxNum)
    {
        // examine child nodes
        Photon photon = balancedKDTree[p];
        
        if (2 * p + 1 < balancedKDTree.size())
        {
            assert(photon.splitAxis != -1);
            Vector3 diff = position - photon.position;
            real_t diffToPlane = 0.0;
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
            nearestPhotons.push_back(photon);
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
//                    sqrDist = nearestPhotons.front().squaredDistance;
//                }
//            }
        }
    }
    
    // c style locate photons
    void Raytracer::cPhotonLocate(size_t p,
                                  Vector3 position,
                                  std::vector<cPhoton> &balancedKDTree,
                                  std::vector<cPhoton> &nearestPhotons,
                                  float &sqrDist,
                                  size_t maxNum)
    {
        assert(balancedKDTree.size() > 0);
        
        // examine child nodes
        cPhoton cphoton = balancedKDTree[p];
        Vector3 cphotonpos(cphoton.position);
//        Photon photon = dataSource[cphoton.index];
        
        if (2 * p + 1 < balancedKDTree.size())
        {
            assert(cphoton.splitAxis != -1);
            Vector3 diff = position - cphotonpos;
            real_t diffToPlane = 0.0;
            switch (cphoton.splitAxis) {
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
            float sqrDiffToPlane = diffToPlane * diffToPlane;
            
            if (diffToPlane < 0)
            {
                // search left subtree
                cPhotonLocate(2 * p, position, balancedKDTree, nearestPhotons, sqrDist, maxNum);
                if (sqrDiffToPlane < sqrDist)
                {
                    // check right subtree
                    cPhotonLocate(2 * p + 1, position, balancedKDTree, nearestPhotons, sqrDist, maxNum);
                }
            }
            else
            {
                // search right subtree
                cPhotonLocate(2 * p + 1, position, balancedKDTree, nearestPhotons, sqrDist, maxNum);
                if (sqrDiffToPlane < sqrDist)
                {
                    // check left subtree
                    cPhotonLocate(2 * p, position, balancedKDTree, nearestPhotons, sqrDist, maxNum);
                }
            }
            
        }
        
        // compute true squared distance to photon
        real_t sqrDistPhoton = squared_distance(position, cphotonpos);
        if (sqrDistPhoton <= sqrDist)
        {
            nearestPhotons.push_back(cphoton);
        }
    }
    
    void Raytracer::applyGammaHDR(Color3 &color)
    {
        real_t A = (0.5);
        real_t gamma = (0.5);
        
        color.r = A * pow(color.r, gamma);
        color.g = A * pow(color.g, gamma);
        color.b = A * pow(color.b, gamma);
    }
    
    void Raytracer::setPhotonDirection(Photon &photon, Vector3 dir)
    {
        // from jensen's implementation
        int theta = int( acos(dir.z)*(256.0/M_PI) );
        if (theta > 255)
            photon.theta = 255;
        else
            photon.theta = (unsigned char)theta;
        
        int phi = int( atan2(dir.y,dir.x) * (256.0 / (2.0 * M_PI)) );
        if (phi > 255)
            photon.phi = 255;
        else if (phi < 0)
            photon.phi = (unsigned char)(phi + 256);
        else
            photon.phi = (unsigned char)phi;
    }
    
    Vector3 Raytracer::getPhotonDirection(Photon &photon)
    {
        return Vector3(sintheta[photon.theta] * cosphi[photon.phi],
                       sintheta[photon.theta] * sinphi[photon.phi],
                       costheta[photon.theta]);
    }
    
    void Raytracer::vvhKDTreePreprocess(std::vector<cPhoton>& /*source*/,
                                        std::vector<metaCPhoton>& /*sortedSource*/)
    {
//        std::vector<cPhoton> datax = source;
//        std::vector<cPhoton> datay = source;
//        std::vector<cPhoton> dataz = source;
//        
//        std::sort(datax.begin(), datax.end(), cPhotonComparatorX);
//        std::sort(datay.begin(), datay.end(), cPhotonComparatorY);
//        std::sort(dataz.begin(), dataz.end(), cPhotonComparatorZ);
        
    }
    
    /// vvh
    void Raytracer::vvhKDTreeConstruction(std::vector<cPhoton> &list, KDNode *root)
    {
        if (list.size() < 1) {
            return;
        }
        
        // VVH kdtree data
        std::vector<KDNode *> nodeList;
        std::vector<KDNode *> activeList;
        std::vector<KDNode *> smallList;
        std::vector<KDNode *> nextList;
        
        nodeList.clear();
        activeList.clear();
        smallList.clear();
        nextList.clear();
        
        // create root node
//        KDNode root;
//        root = new KDNode;
        root->isLeaf = false;
        root->head = 0;
        root->tail = list.size();
        
        activeList.push_back(root);
        
        // large node stage
        while (activeList.size() > 0)
        {
//            printf("active list size = %ld\n",activeList.size());
            vvhProcessLargeNode(activeList, smallList, nextList, list);
            nodeList.insert(nodeList.end(), activeList.begin(), activeList.end());
            
            // swap
            std::vector<KDNode *> tmp = activeList;
            activeList = nextList;
            nextList = tmp;
            
            nextList.clear();
        }
        
        // TODO
        // small node stage
        vvhPreprocessSmallNodes(smallList, nodeList, list);
        activeList = smallList;
        while (activeList.size() > 0)
        {
            vvhProcessSmallNode(activeList, nextList, list);
            nodeList.insert(nodeList.end(), activeList.begin(), activeList.end());
            
            // swap
            std::vector<KDNode *> tmp = activeList;
            activeList = nextList;
            nextList = tmp;
            
            nextList.clear();
        }
        
        // kd-tree output stage
        vvhPreorderTraversal(nodeList, list);
    }
    
    void Raytracer::vvhProcessLargeNode(std::vector<KDNode *> &activeList,
                                        std::vector<KDNode *> &smallList,
                                        std::vector<KDNode *> &nextList,
                                        std::vector<cPhoton> &list)
    {
        for (size_t i = 0; i < activeList.size(); i++)
        {
            // get the surrounding cube
            // O(N)
            Vector3 max = Vector3(-INFINITY, -INFINITY, -INFINITY);
            Vector3 min = Vector3(INFINITY, INFINITY, INFINITY);
            
            KDNode *node = activeList[i];
            int begin = node->head;
            int end = node->tail;
            
//            printf("node %d, begin = %d, end = %d\n", i, begin, end);
            for (int j = begin; j < end; j++)
            {
                // calculate box
                max.x = list[j].position[0] >= max.x ? list[j].position[0] : max.x;
                max.y = list[j].position[1] >= max.y ? list[j].position[1] : max.y;
                max.z = list[j].position[2] >= max.z ? list[j].position[2] : max.z;
                
                min.x = list[j].position[0] <= min.x ? list[j].position[0] : min.x;
                min.y = list[j].position[1] <= min.y ? list[j].position[1] : min.y;
                min.z = list[j].position[2] <= min.z ? list[j].position[2] : min.z;
            }
            
            int splitAxis = -1;
            Vector3 diff = Vector3(max.x - min.x, max.y - min.y, max.z - min.z);
            if ((diff.x >= diff.y) && (diff.x >= diff.z))
                splitAxis = 0;
            else if ((diff.y >= diff.x) && (diff.y >= diff.z))
                splitAxis = 1;
            else if ((diff.z >= diff.x) && (diff.z >= diff.y))
                splitAxis = 2;
            
            // Sorting the vector
            bool (*comparator)(const cPhoton &a, const cPhoton &b) = NULL;
            
            switch (splitAxis) {
                case 0:
                    comparator = cPhotonComparatorX;
                    break;
                case 1:
                    comparator = cPhotonComparatorY;
                    break;
                case 2:
                    comparator = cPhotonComparatorZ;
                    break;
                    
                default:
                    break;
            }
            
            // O(NlogN)
//            std::sort(list.begin() + begin, list.begin() + end, comparator);
            int median = node->head + (node->tail - node->head)/2;
            std::nth_element(list.begin() + begin, list.begin() + median, list.begin() + end, comparator);
            
            
            // Split node i at spatial median of the longest axis
            node->cphotonIndex = median;
            node->splitAxis = splitAxis;
            node->splitValue = list[node->cphotonIndex].position[splitAxis];
            list[node->cphotonIndex].splitAxis = splitAxis;
            
            KDNode *lch = new KDNode;
            KDNode *rch = new KDNode;
            
            lch->isLeaf = false;
            lch->head = begin;
            lch->tail = median;
            
            rch->isLeaf = false;
            rch->head = median + 1;
            rch->tail = end;
            
            node->left = lch;
            node->right = rch;
            
            if (lch->tail - lch->head <= SMALL_NODE_GRANULARITY)
            {
#if SIMPLE_SMALL_NODE
                lch->isLeaf = true;
#endif
                smallList.push_back(lch);
            }
            else
                nextList.push_back(lch);
            
            if (rch->tail - rch->head <= SMALL_NODE_GRANULARITY)
            {
#if SIMPLE_SMALL_NODE
                rch->isLeaf = true;
#endif
                smallList.push_back(rch);
            }
            else
                nextList.push_back(rch);
        }
    }
    
    void Raytracer::vvhPreprocessSmallNodes(std::vector<KDNode *> &smallList,
                                            std::vector<KDNode *> &nodeList,
                                            std::vector<cPhoton> & /*list*/)
    {
        // temporary treat all small nodes as leaf nodes
#if SIMPLE_SMALL_NODE
        nodeList.insert(nodeList.end(), smallList.begin(), smallList.end());
        smallList.clear();
#endif
    }
    
    void Raytracer::vvhProcessSmallNode(std::vector<KDNode *> &activelist,
                                        std::vector<KDNode *> &nextList,
                                        std::vector<cPhoton> &list)
    {
        for (size_t i = 0; i < activelist.size(); i++)
        {
            KDNode *node = activelist[i];
            float VVH0 = node->tail - node->head;
            std::vector<cPhoton> data(node->tail - node->head);
            std::copy(list.begin() + node->head, list.begin() + node->tail, data.begin());
            char splitAxis = -1;
            if (vvhComputeVVH(data, VVH0, splitAxis))
            {
                // need split
                std::copy(data.begin(), data.end(), list.begin() + node->head);
                node->splitAxis = splitAxis;
                
                int begin = node->head;
                int end = node->tail;
                int median = node->head + (node->tail - node->head)/2;
                
                // Split node i at spatial median of the longest axis
                node->cphotonIndex = median;
                node->splitAxis = splitAxis;
                node->splitValue = list[node->cphotonIndex].position[(int)splitAxis];
                list[node->cphotonIndex].splitAxis = splitAxis;
                
                KDNode *lch = new KDNode;
                KDNode *rch = new KDNode;
                
                lch->isLeaf = false;
                lch->head = begin;
                lch->tail = median;
                
                rch->isLeaf = false;
                rch->head = median + 1;
                rch->tail = end;
                
                node->left = lch;
                node->right = rch;
                
                nextList.push_back(lch);
                nextList.push_back(rch);
                
            }
            else
            {
                activelist[i]->isLeaf = true;
            }
            
        }
    }
    
    bool Raytracer::vvhComputeVVH(std::vector<cPhoton> &data, float &VVH0, char &axis)
    {
        float VVH = VVH0;
        bool ret = false;
        char splitAxis = -1;
        std::vector<cPhoton> xdata = data;
        std::vector<cPhoton> ydata = data;
        std::vector<cPhoton> zdata = data;
        std::sort(xdata.begin(), xdata.end(), cPhotonComparatorX);
        std::sort(ydata.begin(), ydata.end(), cPhotonComparatorY);
        std::sort(zdata.begin(), zdata.end(), cPhotonComparatorZ);
        
        Vector3 a(xdata[0].position[0], ydata[0].position[1], zdata[0].position[2]);
        Vector3 b(xdata[xdata.size() - 1].position[0], ydata[ydata.size() - 1].position[1], zdata[zdata.size() - 1].position[2]);
        float v = vvhComputeVolume(a, b);
        
        // x axis
        for (size_t i = 1; i < xdata.size()-1; i++)
        {
            Vector3 la(xdata[0].position);
            Vector3 lb(xdata[i-1].position);
            Vector3 ra(xdata[i+1].position);
            Vector3 rb(xdata[xdata.size()-1].position);
            float vl = vvhComputeVolume(la, lb);
            float vr = vvhComputeVolume(ra, rb);
            float vvh_i = 1 + vl/v + vr/v;
            if (vvh_i < VVH) {
                VVH = vvh_i;
                splitAxis = 0;
            }
        }
        
        // y axis
        for (size_t i = 1; i < ydata.size()-1; i++)
        {
            Vector3 la(ydata[0].position);
            Vector3 lb(ydata[i-1].position);
            Vector3 ra(ydata[i+1].position);
            Vector3 rb(ydata[ydata.size()-1].position);
            float vl = vvhComputeVolume(la, lb);
            float vr = vvhComputeVolume(ra, rb);
            float vvh_i = 1 + vl/v + vr/v;
            if (vvh_i < VVH) {
                VVH = vvh_i;
                splitAxis = 1;
            }
        }
        
        // z axis
        for (size_t i = 1; i < zdata.size()-1; i++)
        {
            Vector3 la(zdata[0].position);
            Vector3 lb(zdata[i-1].position);
            Vector3 ra(zdata[i+1].position);
            Vector3 rb(zdata[zdata.size()-1].position);
            float vl = vvhComputeVolume(la, lb);
            float vr = vvhComputeVolume(ra, rb);
            float vvh_i = 1 + vl/v + vr/v;
            if (vvh_i < VVH) {
                VVH = vvh_i;
                splitAxis = 2;
            }
        }
        
        axis = splitAxis;
        
        // need split
        if (VVH < VVH0)
        {
            switch (splitAxis) {
                case 0:
                    data = xdata;
                    ret = true;
                    break;
                case 1:
                    data = ydata;
                    ret = true;
                    break;
                case 2:
                    data = zdata;
                    ret = true;
                    break;
                    
                default:
                    break;
            }
        }
        
        return ret;
    }
    
    float Raytracer::vvhComputeVolume(Vector3 a, Vector3 b)
    {
        if (a == b)
            return SPHERE_VOLUME;
        
        Vector3 diff = b - a;
        return diff.x * diff.y * diff.z + TWO_MUL_RADIUS * (diff.x * diff.y + diff.x * diff.z + diff.y + diff.z ) + SPHERE_VOLUME;
    }
    
    void Raytracer::vvhPreorderTraversal(std::vector<KDNode *>& /*nodeList*/,
                                         std::vector<cPhoton> & /*list*/ )
    {
//        for (std::vector<KDNode *>::iterator it = nodeList.begin(); it != nodeList.end(); it++)
//        {
//            KDNode *node = *it;
//            printf("head = %d, tail = %d\n",node->head, node->tail);
//        }
//        exit(0);
    }
    
    void Raytracer::vvhcPhotonLocate(Vector3 position,
                                     KDNode *node,
                                     ispcCPhotonData &ispcCphotonData,
                                     std::vector<cPhoton> &cPhotonList,
                                     std::vector<int> &nearestPhotonsIndices,
                                     float &sqrDist,
                                     size_t maxNum)
    {
        assert(node != NULL);
        
        if (cPhotonList.size() < 1) {
            return;
        }
        
        // if the node is leaf, check children for get nearstPhotons
        if (node->isLeaf)
        {
            int begin = node->head;
            int end = node->tail;
            std::vector<cPhoton>::iterator it;

            for (int i = begin; i < end; i++)
            {
                cPhoton cphoton = cPhotonList[i];
                Vector3 cphotonpos(cphoton.position);
                
                // compute true squared distance to photon
                real_t sqrDistPhoton = squared_distance(position, cphotonpos);
                if (sqrDistPhoton <= sqrDist)
                {
                    nearestPhotonsIndices.push_back(cphoton.index);
                }
            }
        }
        // if the node is a splitting node
        else
        {
            cPhoton cphoton = cPhotonList[node->cphotonIndex];
            Vector3 cphotonpos(cphoton.position);
            
            Vector3 diff = position - cphotonpos;
            
            real_t diffToPlane = 0.0;
            switch (cphoton.splitAxis) {
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
                vvhcPhotonLocate(position,
                                 node->left,
                                 ispcCphotonData,
                                 cPhotonList,
                                 nearestPhotonsIndices,
                                 sqrDist,
                                 maxNum);
                
                if (sqrDiffToPlane < sqrDist)
                {
                    // check right subtree
                    vvhcPhotonLocate(position,
                                     node->right,
                                     ispcCphotonData,
                                     cPhotonList,
                                     nearestPhotonsIndices,
                                     sqrDist,
                                     maxNum);
                }
            }
            else
            {
                // search right subtree
                vvhcPhotonLocate(position,
                                 node->right,
                                 ispcCphotonData,
                                 cPhotonList,
                                 nearestPhotonsIndices,
                                 sqrDist,
                                 maxNum);
                
                if (sqrDiffToPlane < sqrDist)
                {
                    // check left subtree
                    vvhcPhotonLocate(position,
                                     node->left,
                                     ispcCphotonData,
                                     cPhotonList,
                                     nearestPhotonsIndices,
                                     sqrDist,
                                     maxNum);
                }
            }
            
            // compute true squared distance to photon
            real_t sqrDistPhoton = squared_distance(position, cphotonpos);
            if (sqrDistPhoton <= sqrDist)
            {
                nearestPhotonsIndices.push_back(cphoton.index);
            }
        }
    }
    
    
    
    
    
    
} /* _462 */
