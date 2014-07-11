/**
 * @file sphere.hpp
 * @brief Class defnition for Sphere.
 *
 * @author Kristin Siu (kasiu)
 * @author Eric Butler (edbutler)
 */

#ifndef _462_SCENE_SPHERE_HPP_
#define _462_SCENE_SPHERE_HPP_

#include "scene/scene.hpp"
#include "scene/BndBox.hpp"

namespace _462 {

/**
 * A sphere, centered on its position with a certain radius.
 */
class Sphere : public Geometry
{
public:

    real_t radius;
    const Material* material;
    BndBox *bbox_local;
    BndBox *bbox_world;
//    Box *boundingBox;
//    Box *globalBBox;

    Sphere();
    virtual ~Sphere();
    virtual void render() const;
    
    virtual void createBoundingBox() const;
    
    virtual bool hit(Ray ray, real_t t0, real_t t1, HitRecord &rec) const;
    
    // pre computation, for accelerate hit test
//    Vector3 c;
    
};

} /* _462 */

#endif /* _462_SCENE_SPHERE_HPP_ */

