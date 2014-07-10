/**
 * @file scene.hpp
 * @brief Class definitions for scenes.
 *
 */

#ifndef _462_SCENE_SCENE_HPP_
#define _462_SCENE_SCENE_HPP_

#include "math/vector.hpp"
#include "math/quaternion.hpp"
#include "math/matrix.hpp"
#include "math/camera.hpp"
#include "math/color.hpp"
#include "scene/material.hpp"
#include "scene/mesh.hpp"

//#include "scene/bvh.hpp"
#include "ray.hpp"
#include <string>
#include <vector>
#include <cfloat>

namespace _462 {
    
    enum GeometryType
    {
        eNone = 0,
        eTriangle = 1,
        eSphere = 2,
        eModel = 3,
    };
    
    /*
     @brief Changed the sphere light class to also containing parallelogram light
            type = 0 is sphere light
            type = 1 is parallelogram light
     */
    struct SphereLight
    {
        struct Attenuation
        {
            real_t constant;
            real_t linear;
            real_t quadratic;
        };
        
        SphereLight();

        // The position of the light, relative to world origin.
        Vector3 position;
        // The color of the light (both diffuse and specular)
        Color3 color;
        // attenuation
        Attenuation attenuation;
        
        // point light radius
        real_t radius;
        
        // type = 0: sphere light
        // type = 1: parallelogram light
        real_t type;
        
        // two vectors for parallelogram light
        Vector3 vertex1;
        Vector3 vertex2;
        Vector3 normal;
        real_t depth;
    };

    /**
     * @brief Data structure for intersection
     */
    struct HitRecord
    {
        Vector3 position;       // world position
        Vector3 normal;         // world normal
        
        Color3 diffuse;         // Surface diffusive color
        Color3 ambient;         // Surface ambien color
        Color3 specular;        // Surface specular color
        Color3 texture;         // Surface texture color
        
        GeometryType type;      // Type of Geometry
        
        real_t t;               // Hit time
        real_t beta, gamma;     // Barycentric parameters for triangles
        
        real_t refractive_index;    // refractive index
        
        int phong;           // p value for blinn-phong model
        
        bool isLight;
        bool isInShadow;
        
        // Constructor
        HitRecord()
        {
            position = Vector3::Zero();
            normal = Vector3::Zero();
            diffuse = Color3::White();
            ambient = Color3::White();
            specular = Color3::White();
            
            type = eNone;
            t = beta = gamma = 0;
            
            phong = 0;
            
            isLight = false;
            isInShadow = false;
        }
        
        // Get lambertian color that come from light
        Color3 getPixelColorFromLight(SphereLight light)
        {
            Vector3 l = normalize(light.position - position);
            Vector3 n = normal;
            
            real_t squared_dist = squared_distance(light.position, this->position);
            real_t dist = sqrt( squared_dist );
            
            real_t dot_nl = dot(n, l);
            
            Color3 lightColor =
            light.color * (real_t(1) / (light.attenuation.constant +
                                        dist * light.attenuation.linear +
                                        squared_dist * light.attenuation.quadratic));
            
            Color3 lambertian = this->diffuse * lightColor * (dot_nl > 0 ? dot_nl : 0);
            return lambertian;
        }
        
        Color3 getPhotonLambertianColor(Vector3 direction, Color3 photonColor)
        {
            Vector3 l = direction;//normalize(direction);
            Vector3 n = this->normal;
            
            real_t dot_nl = dot(n, l);
            
            Color3 lambertian = this->diffuse * photonColor * (dot_nl > 0 ? dot_nl : 0);
            return lambertian;
        }
    };
    
    class Geometry
    {
    public:
        Geometry();
        virtual ~Geometry();
        
        GeometryType type;
        
        /*
         World transformation are applied in the following order:
         1. Scale
         2. Orientation
         3. Position
         */
        
        // The world position of the object.
        Vector3 position;
        
        // The world orientation of the object.
        // Use Quaternion::to_matrix to get the rotation matrix.
        Quaternion orientation;
        
        // The world scale of the object.
        Vector3 scale;
        
        // Transformation matrix
        Matrix4 mat;
        // Inverse transformation matrix
        Matrix4 invMat;
        // Normal transformation matrix
        Matrix3 normMat;
        
//        float isLight;
        
        /**
         * Renders this geometry using OpenGL in the local coordinate space.
         */
        virtual void render() const = 0;
        
        /**
         * Virtual function: Create the bounding box for Geometry object
         */
        virtual void createBoundingBox() const = 0;
        
        /**
         * Virtual function: To determine if a ray hits the surface, return hit record
         */
        virtual bool hit(Ray ray, real_t t0, real_t t1, HitRecord &rec) const = 0;
        
        bool initialize();
        
        real_t isLight;
    };
    
    /**
     * @brief Bounding box class
     */
    class Box {
    public:
        
        Box() {
            bounds[0] = Vector3::Zero();
            bounds[1] = Vector3::Zero();
        }
        
        Box(const Vector3 &min, const Vector3 &max) {
            bounds[0] = min;
            bounds[1] = max;
        }
        
        // Test if a given ray hits a bounding box
        bool intersect(const Ray &r, real_t t0, real_t t1) const;
        
        // Bounds[0] is bottom left front corner,
        // Bounds[1] is top right back corner.
        Vector3 bounds[2];
    };
    
    
    
    /**
     * The container class for information used to render a scene composed of
     * Geometries.
     */
    class Scene
    {
    public:
        
        /// the camera
        Camera camera;
        /// the background color
        Color3 background_color;
        /// the amibient light of the scene
        Color3 ambient_light;
        /// the refraction index of air
        real_t refractive_index;
        
        /// Creates a new empty scene.
        Scene();
        
        /// Destroys this scene. Invokes delete on everything in geometries.
        ~Scene();
        
        bool initialize();
        
        // accessor functions
        Geometry* const* get_geometries() const;
        size_t num_geometries() const;
        const SphereLight* get_lights() const;
        size_t num_lights() const;
        Material* const* get_materials() const;
        size_t num_materials() const;
        Mesh* const* get_meshes() const;
        size_t num_meshes() const;
        
        /// Clears the scene, and invokes delete on everything in geometries.
        void reset();
        
        // functions to add things to the scene
        // all pointers are deleted by the scene upon scene deconstruction.
        void add_geometry( Geometry* g );
        void add_material( Material* m );
        void add_mesh( Mesh* m );
        void add_light( const SphereLight& l );
        

    private:
        
        typedef std::vector< SphereLight > SphereLightList;
        typedef std::vector< Material* > MaterialList;
        typedef std::vector< Mesh* > MeshList;
        typedef std::vector< Geometry* > GeometryList;
        
        // list of all lights in the scene
        SphereLightList point_lights;
        // all materials used by geometries
        MaterialList materials;
        // all meshes used by models
        MeshList meshes;
        // list of all geometries. deleted in dctor, so should be allocated on heap.
        GeometryList geometries;
        
    private:
        
        // no meaningful assignment or copy
        Scene(const Scene&);
        Scene& operator=(const Scene&);
        
    };

    /* Helper functions for common utilities */
    /**
     * @brief   Adjust the texture coordinates for an object and its surface
     * @param   tex_coord       Texture map coordinate
     * @return  Vector          Adjusted texture coordinate
     */
    inline Vector2 getAdjustTexCoord(Vector2 tex_coord)
    {
        if (tex_coord.x < 0 || tex_coord.y < 0) {
//            printf("(%f, %f)\n",tex_coord.x,tex_coord.y);
//            assert(0);
        }
        
        return Vector2(tex_coord.x > 1 ? tex_coord.x - (int)tex_coord.x : tex_coord.x, tex_coord.y > 1 ? tex_coord.y - (int)tex_coord.y : tex_coord.y);
    }
    
    /*!
     * Get refraction ray direction, text book p 305
     * @param d         Incoming ray's direction.
     * @param n         Hit instance's normal
     * @param ratio     Refractive indices' ratio
     * @param t         OUT New refraction direction
     * @return bool Is there a fraction
     */
    inline bool isRefract(Vector3 d, Vector3 n, real_t ratio, Vector3 &t)
    {
        real_t dot_dn = dot(d, n);
        real_t square = real_t(1) - (pow(ratio, 2) * (real_t(1) - pow(dot_dn, 2)));
        
        if (square < 0) {
            t = Vector3().Zero();
            return false;
        }
        
        t = (ratio * (d - (n * dot_dn))) - (n * sqrt(square));
        return true;
    }
    
    /**
     * @brief   Get the bounding volumn for triangle mesh
     * @param   A, B, C     Position of vertices of triangle
     * @return  Box         Bounding box for the triangle
     */
    inline Box getBoundingBoxForTriangle(Vector3 A, Vector3 B, Vector3 C)
    {
        real_t xmin, xmax, ymin, ymax, zmin, zmax;
        
        xmin = A.x < B.x ? A.x : B.x;
        xmin = C.x < xmin ? C.x : xmin;
        
        xmax = A.x > B.x ? A.x : B.x;
        xmax = C.x > xmax ? C.x : xmax;
        
        ymin = A.y < B.y ? A.y : B.y;
        ymin = C.y < ymin ? C.y : ymin;
        
        ymax = A.y > B.y ? A.y : B.y;
        ymax = C.y > ymax ? C.y : ymax;
        
        zmin = A.z < B.z ? A.z : B.z;
        zmin = C.z < zmin ? C.z : zmin;
        
        zmax = A.z > B.z ? A.z : B.z;
        zmax = C.z > zmax ? C.z : zmax;
        
        Vector3 min = Vector3(xmin, ymin, zmin);
        Vector3 max = Vector3(xmax, ymax, zmax);
        
        Box box = Box(min, max);
        return box;
    }
    
    /**
     * @brief   Get the barycentric paramters for a ray - triangle intersection
     * @param   r           Ray to test
     * @param   A, B, C     Vertices of triangle mesh
     * @return  Vector3     Barycentric coordinates' parameter vector
     */
    inline Vector3 getResultTriangleIntersection(Ray r, Vector3 A, Vector3 B, Vector3 C)
    {
        // Reference: http://en.wikipedia.org/wiki/Cramer's_rule
        // Text book, Page 79
        // result.x = beta, result.y = gamma, result.z = t
        Vector3 result;
        
        // ba, ca, d makes M
        Vector3 ba = A - B;
        Vector3 ca = A - C;
        Vector3 od = r.d;
        Vector3 ea = A - r.e;
        
        real_t a, b, c, d, e, f, g, h, i, j, k, l;
        a = ba.x; b = ba.y; c = ba.z;
        d = ca.x; e = ca.y; f = ca.z;
        g = od.x; h = od.y; i = od.z;
        j = ea.x; k = ea.y; l = ea.z;
        
        real_t M = a * (e * i - h * f) + b * (g * f - d * i) + c * (d * h - e * g);
        result.x = j * (e * i - h * f) + k * (g * f - d * i) + l * (d * h - e * g);
        result.y = a * (i * k - l * h) + b * (g * l - i * j) + c * (j * h - g * k);
        result.z = a * (e * l - k * f) + b * (j * f - d * l) + c * (d * k - j * e);
        
        result /= M;
        
        return result;
    }

    
} /* _462 */

#endif /* _462_SCENE_SCENE_HPP_ */
