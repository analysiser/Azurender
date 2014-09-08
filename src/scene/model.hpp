/**
 * @file model.hpp
 * @brief Model class
 *
 * @author Eric Butler (edbutler)
 */

#ifndef _462_SCENE_MODEL_HPP_
#define _462_SCENE_MODEL_HPP_

#include "scene/azBVHTree.hpp"
#include "scene/mesh.hpp"
#include "scene/scene.hpp"

#include <iostream>
#include <set>
#include <vector>
#include <map>
#include <list>

namespace _462 {

    class azBVHTree;
    /**
     * A mesh of triangles.
     */
    class Model : public Geometry
    {
    public:
        
        const Mesh* mesh;
        const Material* material;
        
        mutable azBVHTree *bvhTree;

        Model();

        
        virtual ~Model();

        
        virtual void render() const;

        
        //
        virtual void createBoundingBox() const;
        
        // Override of virtual function from Geometry
        virtual bool hit(Ray ray, real_t t0, real_t t1, HitRecord &rec) const;

    };
    
} /* _462 */

#endif /* _462_SCENE_MODEL_HPP_ */
