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

#define MAX_DEPTH                       (5)

#define MAX_THREADS_SCATTER             (4)
#define MAX_THREADS_TRACE               (4)

#include "math/color.hpp"
#include "math/random462.hpp"
#include "math/vector.hpp"
#include "scene/scene.hpp"
#include "raytracer/Photon.hpp"
#include "raytracer/Utils.h"
#include <stack>
#include <thread>

namespace _462 {
    
    struct PhotonScatterData
    {
        std::vector<Photon> worker_photon_indirect;
        std::vector<Photon> worker_photon_caustics;
        SphereLight light;
        unsigned int indirect_needed;
        unsigned int caustics_needed;
    };
    
    
    struct TracePixelData
    {
//        Scene *scene;
        Color3 *buffer;
        int row_num;
        
    };
    
    struct ispcCPhotonData
    {
        int size;
        float *posx;
        float *posy;
        float *posz;
        
        int8_t *bitmap;
    };

    
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
        
        void perPixelRender(unsigned char* buffer);
        
        // indirect and caustics list for photons to trace
        std::vector<Photon> photon_indirect_list;
        std::vector<Photon> photon_caustic_list;
        
        // balanced kdtree
        std::vector<Photon> kdtree_photon_indirect_list;
        std::vector<Photon> kdtree_photon_caustic_list;
        
        // cPhoton source data
        std::vector<cPhoton> cphoton_indirect_data;
        std::vector<cPhoton> cphoton_caustics_data;
        
        // ispc friendly cPhoton data, for parallel computation
        ispcCPhotonData ispc_cphoton_indirect_data;
        ispcCPhotonData ispc_cphoton_caustics_data;
        
        // balanced c style photon kdtree
        std::vector<cPhoton> kdtree_cphoton_indirect;
        std::vector<cPhoton> kdtree_cphoton_caustics;
        
        KDNode *vvh_indirect_root;
        KDNode *vvh_caustics_root;
        
        // caustics bounding boxes
        // it stores bounding boxes of objects that might generate caustics
        std::vector<Box> cboxes;
        
        unsigned int num_iteration;
        
        
    private:
        
        // data used for photon direction encoding
        float costheta[256];
        float sintheta[256];
        float cosphi[256];
        float sinphi[256];
        
        // data used for measuring final gathering performance
        unsigned int master_start;  // Overall start time for raytracing
        unsigned int pass_start;    // Start time for each pass of raytracing & photon mapping
        unsigned int pass_end;      // End time for each pass
        unsigned int master_end;    // Overall end time for raytracing & photon mapping
        unsigned int acc_cphoton_search_time;
        unsigned int acc_iphoton_search_time;
        
        // data used for measuring the final gathering radius
        float radius_clear;
        float radius_shadow;
        unsigned int clear_count;
        unsigned int shadow_count;
        unsigned int acc_pass_spent;
        unsigned int acc_kdtree_cons;
        
        TracePixelData tracePixelData[MAX_THREADS_TRACE];

        Color3 *raytraceColorBuffer;
        
        Color3 trace_pixel(const Scene* scene,
                           size_t x,
                           size_t y,
                           size_t width,
                           size_t height);
        
        // Photon Scatter
        void photonScatter(const Scene* scene);
        
        // Parallel Photon Scatter
        void parallelPhotonScatter(const Scene* scene);
        
        // worker node for parallel photon scatter
        void photonScatterWorker(PhotonScatterData *data);
        
        // kdtree construction
        void kdtreeConstruction();
        
        // c style photon kdtree construction
        void cPhotonKDTreeConstruction();
        
        /////////// vvh kdtree ///////////
        // vvh kdtree preprocessing for faster construction
        void vvhKDTreePreprocess(std::vector<cPhoton>& source,
                                 std::vector<metaCPhoton>& sortedSource);
                                 
        
        // vvh kdtree construction
        void vvhKDTreeConstruction(std::vector<cPhoton> &list, KDNode *root);
        
        void vvhProcessLargeNode(std::vector<KDNode *> &activeList,
                                 std::vector<KDNode *> &smallList,
                                 std::vector<KDNode *> &nextList,
                                 std::vector<cPhoton> &list);
        
        void vvhPreprocessSmallNodes(std::vector<KDNode *> &smallList,
                                     std::vector<KDNode *> &nodeList,
                                     std::vector<cPhoton> &list);
        
        void vvhProcessSmallNode(std::vector<KDNode *> &activelist,
                                 std::vector<KDNode *> &nextList,
                                 std::vector<cPhoton> &list);
        
        void vvhPreorderTraversal(std::vector<KDNode *> &nodeList,
                                  std::vector<cPhoton> &list);
        
        void vvhcPhotonLocate(Vector3 position,
                              KDNode *node,
                              ispcCPhotonData &ispcCphotonData,
                              std::vector<cPhoton> &cPhotonList,
                              std::vector<int> &nearestPhotonsIndices,
                              float &sqrDist,
                              size_t maxNum);
        
        bool vvhComputeVVH(std::vector<cPhoton> &data, float &VVH0, char &axis);
        
        float vvhComputeVolume(Vector3 a, Vector3 b);
        
        /////////// vvh kdtree ///////////
        
        // worker node for raytracing and photon radiance gathering
        void tracePixelWorker(TracePixelData *data);
        
        // Photon Tracing for global illumination and caustics
        void photonTrace(Ray ray, real_t t0, real_t t1, int depth);
        
        // Photon Tracing function for a single thread purposes
        void localPhotonTrace(Ray ray, real_t t0, real_t t1, int depth,
                              std::vector<Photon> &indirect_list, std::vector<Photon> &caustics_list,
                              size_t indirect_needed, size_t caustics_needed);
        
        // Raytracing helper function, to decide if there is a hit on a surface to shade
        Color3 trace(Ray ray, real_t t0, real_t t1, int depth);
        
        Color3 trace_DOF(HitRecord &record, Ray ray, real_t t0, real_t t1, int depth);
        
        // Shading function, shades the hit record from a surface
        Color3 shade(Ray ray, HitRecord record, real_t t0, real_t t1, int depth);
        
        // Shading of direct illumination
        Color3 shade_direct_illumination(HitRecord &record, real_t t0, real_t t1);
        
        // Shading of caustics
        Color3 shade_caustics(HitRecord &record, real_t radius, size_t num_samples);
        
        // Shading of Indirect lightning
        Color3 shade_indirect_illumination(HitRecord &record, real_t radius, size_t num_samples);
        
        // Shade global illumination and caustics at the same time
        Color3 shade_photons(HitRecord &record, real_t radius, size_t num_samples);
        
        // Test: Shade c photons
        Color3 shade_cphotons(HitRecord &record, real_t radius, size_t num_samples);
        
        // retrieve the closest hit record
        HitRecord getClosestHit(Ray r, real_t t0, real_t t1, bool *isHit);
        
        // sample a Light source on the volumn of sphere light
        Vector3 sampleLightSource(SphereLight light);
        
        Ray getPhotonEmissionRayFromLight(SphereLight light);
        
        Ray getPhotonEmissionRayForCaustics(SphereLight light);
        
        void generateCausticsBoxes();
        
        // helper function for sampling a point on a given unit sphere
        Vector3 samplePointOnUnitSphere();
        
        Vector3 samplePointOnUnitSphereUniform();
        
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
        
        // blance kdtree, kdtree stored in compact complete binary tree format
        void balance(size_t index, std::vector<Photon> &balancedKDTree, std::vector<Photon> &list);
        
        // cPhoton balanced kdtree
        void cPhotonBalance(size_t index, std::vector<cPhoton> &balancedKDTree, std::vector<cPhoton> &list);
        
        // find nearest photons
        void locatePhotons(size_t p,
                           Vector3 position,
                           std::vector<Photon> &balancedKDTree,
                           std::vector<Photon> &nearestPhotons,
                           float &sqrDist,
                           size_t maxNum);
        
        // c style locate photons
        void cPhotonLocate(size_t p,
                           Vector3 position,
                           std::vector<cPhoton> &balancedKDTree,
                           std::vector<cPhoton> &nearestPhotons,
                           float &sqrDist,
                           size_t maxNum);
        
        void applyGammaHDR(Color3 &color);
        
//        static const unsigned char filter[];
        
        // set encoded photon direction
        void setPhotonDirection(Photon &photon, Vector3 dir);
        
        // get decoded photon direction
        Vector3 getPhotonDirection(Photon &photon);
        
    };
    
    inline real_t getGaussianFilterWeight(real_t dist_sqr, real_t radius_sqr);
    /*!
     @brief get the gaussian filter weight, jensen P67, this value is the only weight
     @param dist_x_p    distance between x and photon p
     @param dist_max    r is the maximum distance
    */
    inline real_t getConeFilterWeight(real_t dist_x_p, real_t dist_max);
    
    
    // Photon comparator for sort functions
    inline bool photonComparatorX(const Photon &a, const Photon &b)
    {
        return a.position.x < b.position.x;
    }
    
    inline bool photonComparatorY(const Photon &a, const Photon &b)
    {
        return a.position.y < b.position.y;
    }
    
    inline bool photonComparatorZ(const Photon &a, const Photon &b)
    {
        return a.position.z < b.position.z;
    }
    
    // comparator for c style photon sorting
    inline bool cPhotonComparatorX(const cPhoton &a, const cPhoton &b)
    {
        return a.position[0] < b.position[0];
    }
    
    inline bool cPhotonComparatorY(const cPhoton &a, const cPhoton &b)
    {
        return a.position[1] < b.position[1];
    }
    
    inline bool cPhotonComparatorZ(const cPhoton &a, const cPhoton &b)
    {
        return a.position[2] < b.position[2];
    }
    
} /* _462 */

#endif /* _462_RAYTRACER_HPP_ */
