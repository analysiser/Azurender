/**
 * @file triangle.cpp
 * @brief Function definitions for the Triangle class.
 *
 * @author Eric Butler (edbutler)
 */

#include "scene/triangle.hpp"
#include "application/opengl.hpp"

namespace _462 {

    
    Triangle::Triangle()
    {
        vertices[0].material = 0;
        vertices[1].material = 0;
        vertices[2].material = 0;
        
        bbox_local = new BndBox();
        bbox_world = new BndBox();
        
        type = eTriangle;
    }
    
    Triangle::~Triangle() {
        delete bbox_local;
        delete bbox_world;
    }
    
    void Triangle::render() const
    {
        bool materials_nonnull = true;
        for ( int i = 0; i < 3; ++i )
            materials_nonnull = materials_nonnull && vertices[i].material;
        
        // this doesn't interpolate materials. Ah well.
        if ( materials_nonnull )
            vertices[0].material->set_gl_state();
        
//        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glBegin(GL_TRIANGLES);
        
        glNormal3dv( &vertices[0].normal.x );
        glTexCoord2dv( &vertices[0].tex_coord.x );
        glVertex3dv( &vertices[0].position.x );
        
        glNormal3dv( &vertices[1].normal.x );
        glTexCoord2dv( &vertices[1].tex_coord.x );
        glVertex3dv( &vertices[1].position.x);
        
        glNormal3dv( &vertices[2].normal.x );
        glTexCoord2dv( &vertices[2].tex_coord.x );
        glVertex3dv( &vertices[2].position.x);
        
        glEnd();
        
        if ( materials_nonnull )
            vertices[0].material->reset_gl_state();
    }
    
    /**
     * @brief   Create a bounding box for triangle mesh
     */
    void Triangle::createBoundingBox() const
    {
        for (int i = 0; i < 3; i++) {
            // local bounding box
            bbox_local->include(vertices[i].position);
            bbox_world->include(mat.transform_point(vertices[i].position));
        }
    }
    
    bool Triangle::hit(Ray ray, real_t t0, real_t t1, HitRecord &rec) const
    {
        if(!bbox_world->intersect(ray, t0, t1))
            return false;
        
        // Transform ray to sphere's local space
        Ray r = Ray(invMat.transform_point(ray.e), invMat.transform_vector(ray.d));

        Vertex A = vertices[0];
        Vertex B = vertices[1];
        Vertex C = vertices[2];
        
        /// result.x = beta, result.y = gamma, result.z = t
        Vector3 result = getResultTriangleIntersection(r, A.position, B.position, C.position);
        
        if (result.z < t0 || result.z > t1) {
            return false;
        }

        if (result.y < 0 || result.y > 1) {
            return false;
        }
        
        if (result.x < 0 || result.x > 1 - result.y) {
            return false;
        }
        
        rec.type = eTriangle;
        
        rec.position = ray.e + result.z * ray.d;
        
        real_t beta = result.x;
        real_t gamma = result.y;
        real_t alpha = 1 - beta - gamma;
        
        rec.normal = normalize(alpha * (normMat * A.normal) + beta * (normMat * B.normal) + gamma * (normMat * C.normal));
        
        rec.diffuse = alpha * A.material->diffuse + beta * B.material->diffuse + gamma * C.material->diffuse;
        rec.ambient = alpha * A.material->ambient + beta * B.material->ambient + gamma * C.material->ambient;
        rec.specular = alpha * A.material->specular + beta * B.material->specular + gamma * C.material->specular;
        
        // phong intersection... weird
        rec.phong = alpha * A.material->phong + beta * B.material->phong + gamma * C.material->phong;
        
        // texture interpolation
        Vector2 tex_cood_interpolated = alpha * A.tex_coord + beta * B.tex_coord + gamma * C.tex_coord;
        tex_cood_interpolated = getAdjustTexCoord(tex_cood_interpolated);
        real_t u = tex_cood_interpolated.x, v = tex_cood_interpolated.y;
        
        int width, height;
        A.material->get_texture_size(&width, &height);
        Color3 texA = A.material->get_texture_pixel(u * width, v * height);
        B.material->get_texture_size(&width, &height);
        Color3 texB = B.material->get_texture_pixel(u * width, v * height);
        C.material->get_texture_size(&width, &height);
        Color3 texC = C.material->get_texture_pixel(u * width, v * height);
        
        
        rec.texture = alpha * texA + beta * texB + gamma * texC;
        
        rec.t = result.z;
        
        rec.refractive_index =
        alpha * A.material->refractive_index +
        beta * B.material->refractive_index +
        gamma * C.material->refractive_index;
        
//        rec.isLight = this->isLight;
        
        return true;

    }
    
    
} /* _462 */
