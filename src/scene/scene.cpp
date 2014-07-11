/**
 * @file scene.cpp
 * @brief Function definitions for scenes.
 *
 * @author Eric Butler (edbutler)
 * @author Kristin Siu (kasiu)
 */

#include "scene/scene.hpp"

namespace _462 {
    
    
    Geometry::Geometry():
    position(Vector3::Zero()),
    orientation(Quaternion::Identity()),
    scale(Vector3::Ones())
    {
        
    }
    
    Geometry::~Geometry() { }
    
    bool Geometry::initialize()
    {
        make_inverse_transformation_matrix(&invMat, position, orientation, scale);
        make_transformation_matrix(&mat, position, orientation, scale);
        make_normal_matrix(&normMat, mat);
        
        position_local = invMat.transform_point(position);
        
        return true;
    }
    
    SphereLight::SphereLight():
    position(Vector3::Zero()),
    color(Color3::White()),
	radius(real_t(0))
    {
        attenuation.constant = 1;
        attenuation.linear = 0;
        attenuation.quadratic = 0;
    }
    
    Scene::Scene()
    {
        reset();
    }
    
    Scene::~Scene()
    {
        reset();
    }
    
    bool Scene::initialize()
    {
        bool res = true;
        for (unsigned int i = 0; i < num_geometries(); i++)
            res &= geometries[i]->initialize();
        return res;
    }
    
    
    Geometry* const* Scene::get_geometries() const
    {
        return geometries.empty() ? NULL : &geometries[0];
    }
    
    size_t Scene::num_geometries() const
    {
        return geometries.size();
    }
    
    const SphereLight* Scene::get_lights() const
    {
        return point_lights.empty() ? NULL : &point_lights[0];
    }
    
    size_t Scene::num_lights() const
    {
        return point_lights.size();
    }
    
    Material* const* Scene::get_materials() const
    {
        return materials.empty() ? NULL : &materials[0];
    }
    
    size_t Scene::num_materials() const
    {
        return materials.size();
    }
    
    Mesh* const* Scene::get_meshes() const
    {
        return meshes.empty() ? NULL : &meshes[0];
    }
    
    size_t Scene::num_meshes() const
    {
        return meshes.size();
    }
    
    void Scene::reset()
    {
        for ( GeometryList::iterator i = geometries.begin(); i != geometries.end(); ++i ) {
            delete *i;
        }
        for ( MaterialList::iterator i = materials.begin(); i != materials.end(); ++i ) {
            delete *i;
        }
        for ( MeshList::iterator i = meshes.begin(); i != meshes.end(); ++i ) {
            delete *i;
        }
        
        geometries.clear();
        materials.clear();
        meshes.clear();
        point_lights.clear();
        
        camera = Camera();
        
        background_color = Color3::Black();
        ambient_light = Color3::Black();
        refractive_index = 1.0;
        
    }
    
    void Scene::add_geometry( Geometry* g )
    {
        geometries.push_back( g );
    }
    
    void Scene::add_material( Material* m )
    {
        materials.push_back( m );
    }
    
    void Scene::add_mesh( Mesh* m )
    {
        meshes.push_back( m );
    }
    
    void Scene::add_light( const SphereLight& l )
    {
        point_lights.push_back( l );
    }
    
    /**
     * @brief test if a ray hits the bounding box
     * @param   &r      reference of the ray
     * @param   t0      lower limit of t
     * @param   t1      upper limit of t
     * @return  bool    if the intersection happens
     */
    bool Box::intersect(const _462::Ray &r, real_t t0, real_t t1) const
    {
        real_t tmin, tmax, tymin, tymax, tzmin, tzmax;
        
        // x
        if (r.d.x >= 0) {
            tmin = (bounds[0].x - r.e.x) / r.d.x;
            tmax = (bounds[1].x - r.e.x) / r.d.x;
        }
        else {
            tmin = (bounds[1].x - r.e.x) / r.d.x;
            tmax = (bounds[0].x - r.e.x) / r.d.x;
        }
        
        // y
        if (r.d.y >= 0) {
            tymin = (bounds[0].y - r.e.y) / r.d.y;
            tymax = (bounds[1].y - r.e.y) / r.d.y;
        }
        else {
            tymin = (bounds[1].y - r.e.y) / r.d.y;
            tymax = (bounds[0].y - r.e.y) / r.d.y;
        }
        
        if ( (tmin > tymax) || (tymin > tmax) )
            return false;
        
        if (tymin > tmin)
            tmin = tymin;
        
        if (tymax < tmax)
            tmax = tymax;
        
        // z
        if (r.d.z >= 0) {
            tzmin = (bounds[0].z - r.e.z) / r.d.z;
            tzmax = (bounds[1].z - r.e.z) / r.d.z;
        }
        else {
            tzmin = (bounds[1].z - r.e.z) / r.d.z;
            tzmax = (bounds[0].z - r.e.z) / r.d.z;
        }
        
        if ( (tmin > tzmax) || (tzmin > tmax) )
            return false;
        if (tzmin > tmin)
            tmin = tzmin;
        if (tzmax < tmax)
            tmax = tzmax;
        
        return ( (tmin < t1) && (tmax > t0) );
        
    }
    
    
} /* _462 */

