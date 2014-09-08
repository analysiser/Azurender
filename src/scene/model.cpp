/**
 * @file model.cpp
 * @brief Model class
 *
 * @author Eric Butler (edbutler)
 * @author Zeyang Li (zeyangl)
 */

#include "scene/model.hpp"
#include "scene/material.hpp"
#include "application/opengl.hpp"
#include "scene/triangle.hpp"
#include <iostream>
#include <cstring>
#include <string>
#include <fstream>
#include <sstream>


namespace _462 {

    Model::Model() : mesh( 0 ), material( 0 ) {
        type = eModel;
        bvhTree = nullptr;
    }
    Model::~Model() {
    }
    
    void Model::render() const
    {
        if ( !mesh )
            return;
        if ( material )
            material->set_gl_state();
        mesh->render();
        if ( material )
            material->reset_gl_state();
    }
    
    void Model::createBoundingBox() const
    {
        // if already have bounding boxes, then no need to recalculate again
        // reduced about 10 MB for each time a new ray tracing happens
//        if (this->bBox.bounds[0] == Vector3::Zero() && this->bBox.bounds[1] == Vector3::Zero()) {
//            MeshTriangle const *triangles = mesh->get_triangles();
//            
//            for (size_t i = 0; i < mesh->num_triangles(); i++) {
//                
//                Vector3 A = mesh->vertices[triangles[i].vertices[0]].position;
//                Vector3 B = mesh->vertices[triangles[i].vertices[1]].position;
//                Vector3 C = mesh->vertices[triangles[i].vertices[2]].position;
//                
//                Box boundingBox = getBoundingBoxForTriangle(A, B, C);
//                
//                BoxNode *newNode = modelBoxNodes->create(boundingBox, i);
//                modelBoxNodes->append(newNode);
//            }
//        }
        
        if (bvhTree == nullptr) {
            
            MeshTriangle const *triangles = mesh->get_triangles();
            bvhTree = new azBVHTree(mesh->num_triangles());

            for (size_t i = 0; i < mesh->num_triangles(); i++) {
                
                Vector3 A = mesh->vertices[triangles[i].vertices[0]].position;
                Vector3 B = mesh->vertices[triangles[i].vertices[1]].position;
                Vector3 C = mesh->vertices[triangles[i].vertices[2]].position;
                
                BndBox bbox(A);
                bbox.include(B);
                bbox.include(C);
                
                bvhTree->setLeaf(azBVHTree::azBVNode(bbox, i), i);
            }
            
            bvhTree->buildBVHTree();
        }
    }
    
    bool Model::hit(Ray ray, real_t t0, real_t t1, HitRecord &rec) const
    {
        Ray r = Ray(lmat.transform_point(ray.e), lmat.transform_vector(ray.d));
        
        auto rayTriangleIntersectionTest = [this](Ray rr, real_t tt0, real_t tt1, real_t &tt, INT64 triIndex){
            MeshTriangle const *triangles = mesh->get_triangles();
            
            MeshVertex A = mesh->vertices[triangles[triIndex].vertices[0]];
            MeshVertex B = mesh->vertices[triangles[triIndex].vertices[1]];
            MeshVertex C = mesh->vertices[triangles[triIndex].vertices[2]];
            
            // result.x = beta, result.y = gamma, result.z = t
            Vector3 result = getResultTriangleIntersection(rr, A.position, B.position, C.position);
            
            if (result.z < tt0 || result.z > tt1) {
                return false;
            }
            
            if (result.y < 0 || result.y > 1) {
                return false;
            }
            
            if (result.x < 0 || result.x > 1 - result.y) {
                return false;
            }
            
            tt = result.z;
            return true;
        };

        INT64 idx = -1;
        
        if (bvhTree->getFirstIntersectIndex(r, t0, t1, idx, rayTriangleIntersectionTest)) {
            
            MeshTriangle const *triangles = mesh->get_triangles();
            MeshVertex A = mesh->vertices[triangles[idx].vertices[0]];
            MeshVertex B = mesh->vertices[triangles[idx].vertices[1]];
            MeshVertex C = mesh->vertices[triangles[idx].vertices[2]];
            
            // result.x = beta, result.y = gamma, result.z = t
            Vector3 result = getResultTriangleIntersection(r, A.position, B.position, C.position);
            
            rec.type = eTriangle;
            
            rec.position = ray.e + result.z * ray.d;
            
            real_t beta = result.x;
            real_t gamma = result.y;
            real_t alpha = 1 - beta - gamma;
            
            rec.normal = normalize(alpha * (normMat * A.normal) + beta * (normMat * B.normal) + gamma * (normMat * C.normal));
            
            // For texture mapping adjustment
            A.tex_coord = getAdjustTexCoord(A.tex_coord);
            B.tex_coord = getAdjustTexCoord(B.tex_coord);
            C.tex_coord = getAdjustTexCoord(C.tex_coord);
            
            Vector2 tex_cood_interpolated =
            alpha * A.tex_coord + beta * B.tex_coord + gamma * C.tex_coord;
            
            int width = 0, height = 0;
            if (material) {
                material->get_texture_size(&width, &height);
            }
            
            rec.diffuse = material->diffuse;
            rec.ambient = material->ambient;
            rec.specular = material->specular;
            rec.phong = material->phong;
            
            rec.texture = material->get_texture_pixel(tex_cood_interpolated.x * width, tex_cood_interpolated.y * height);
            
            rec.t = result.z;
            
            rec.refractive_index = material->refractive_index;
            
            return true;
        }
        
        return false;
    
    }
    
} /* _462 */
