/**
 * @file triangle.hpp
 * @brief Class definition for Triangle.
 *
 * @author Eric Butler (edbutler)
 */

#ifndef _462_SCENE_TRIANGLE_HPP_
#define _462_SCENE_TRIANGLE_HPP_

#include "scene/scene.hpp"
#include "scene/BndBox.hpp"

namespace _462 {
    
    /**
     * a triangle geometry.
     * Triangles consist of 3 vertices. Each vertex has its own position, normal,
     * texture coordinate, and material. These should all be interpolated to get
     * values in the middle of the triangle.
     * These values are all in local space, so it must still be transformed by
     * the Geometry's position, orientation, and scale.
     */
    class Triangle : public Geometry
    {
    public:
        
        struct Vertex
        {
            // note that position and normal are in local space
            Vector3 position;
            Vector3 normal;
            Vector2 tex_coord;
            const Material* material;
        };
        
        // the triangle's vertices, in CCW order
        Vertex vertices[3];
        
        // Bounding box for triangle mesh
        BndBox *bbox_local;
        BndBox *bbox_world;
        
        Triangle();
        virtual ~Triangle();
        virtual void render() const;
        
        // Creating a bounding box for triangle
        virtual void createBoundingBox() const;
        
        // Override of virtual function from Geometry
        virtual bool hit(Ray ray, real_t t0, real_t t1, HitRecord &rec) const;
    };
    
    
} /* _462 */

#endif /* _462_SCENE_TRIANGLE_HPP_ */
