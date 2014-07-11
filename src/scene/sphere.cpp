/**
 * @file sphere.cpp
 * @brief Function defnitions for the Sphere class.
 *
 * @author Kristin Siu (kasiu)
 * @author Eric Butler (edbutler)
 */

#include "scene/sphere.hpp"
#include "application/opengl.hpp"

namespace _462 {

#define SPHERE_NUM_LAT 80
#define SPHERE_NUM_LON 100
    
#define SPHERE_NUM_VERTICES ( ( SPHERE_NUM_LAT + 1 ) * ( SPHERE_NUM_LON + 1 ) )
#define SPHERE_NUM_INDICES ( 6 * SPHERE_NUM_LAT * SPHERE_NUM_LON )
    // index of the x,y sphere where x is lat and y is lon
#define SINDEX(x,y) ((x) * (SPHERE_NUM_LON + 1) + (y))
#define VERTEX_SIZE 8
#define TCOORD_OFFSET 0
#define NORMAL_OFFSET 2
#define VERTEX_OFFSET 5
    
    static unsigned int Indices[SPHERE_NUM_INDICES];
    static float Vertices[VERTEX_SIZE * SPHERE_NUM_VERTICES];
    
    static void init_sphere()
    {
        static bool initialized = false;
        if ( initialized )
            return;
        
        for ( int i = 0; i <= SPHERE_NUM_LAT; i++ ) {
            for ( int j = 0; j <= SPHERE_NUM_LON; j++ ) {
                real_t lat = real_t( i ) / SPHERE_NUM_LAT;
                real_t lon = real_t( j ) / SPHERE_NUM_LON;
                float* vptr = &Vertices[VERTEX_SIZE * SINDEX(i,j)];
                
                vptr[TCOORD_OFFSET + 0] = lon;
                vptr[TCOORD_OFFSET + 1] = 1-lat;
                
                lat *= PI;
                lon *= 2 * PI;
                real_t sinlat = sin( lat );
                
                vptr[NORMAL_OFFSET + 0] = vptr[VERTEX_OFFSET + 0] = sinlat * sin( lon );
                vptr[NORMAL_OFFSET + 1] = vptr[VERTEX_OFFSET + 1] = cos( lat ),
                vptr[NORMAL_OFFSET + 2] = vptr[VERTEX_OFFSET + 2] = sinlat * cos( lon );
            }
        }
        
        for ( int i = 0; i < SPHERE_NUM_LAT; i++ ) {
            for ( int j = 0; j < SPHERE_NUM_LON; j++ ) {
                unsigned int* iptr = &Indices[6 * ( SPHERE_NUM_LON * i + j )];
                
                unsigned int i00 = SINDEX(i,  j  );
                unsigned int i10 = SINDEX(i+1,j  );
                unsigned int i11 = SINDEX(i+1,j+1);
                unsigned int i01 = SINDEX(i,  j+1);
                
                iptr[0] = i00;
                iptr[1] = i10;
                iptr[2] = i11;
                iptr[3] = i11;
                iptr[4] = i01;
                iptr[5] = i00;
            }
        }
        
        initialized = true;
    }
    
    Sphere::Sphere()
    : radius(0), material(0)
    {
        boundingBox = new Box(Vector3::Zero(), Vector3::Zero());
//        globalBBox = new Box(Vector3::Zero(), Vector3::Zero());
        type = eSphere;
//        c = invMat.transform_point(position);
    }
    
    Sphere::~Sphere()
    {
        delete boundingBox;
//        delete globalBBox;
    }
    
    void Sphere::render() const
    {
        // create geometry if we haven't already
        init_sphere();
        
        if ( material )
            material->set_gl_state();
        
        // just scale by radius and draw unit sphere
        glPushMatrix();
        glScaled( radius, radius, radius );
        glInterleavedArrays( GL_T2F_N3F_V3F, VERTEX_SIZE * sizeof Vertices[0], Vertices );
        glDrawElements( GL_TRIANGLES, SPHERE_NUM_INDICES, GL_UNSIGNED_INT, Indices );
        glPopMatrix();
        
        if ( material )
            material->reset_gl_state();
    }
    
    void Sphere::createBoundingBox() const
    {
        // Transform ray to sphere's local space
        Vector3 center = invMat.transform_point(position);
        Box box = Box(center - Vector3(radius, radius, radius), center + Vector3(radius, radius, radius));
        
        boundingBox->bounds[0] = box.bounds[0];
        boundingBox->bounds[1] = box.bounds[1];
        
        // create global bounding box
//        Vector3 gcenter = position;
//        globalBBox->bounds[0] = gcenter - Vector3(radius, radius, radius);
//        globalBBox->bounds[1] = gcenter + Vector3(radius, radius, radius);
    }

    
    bool Sphere::hit(Ray ray, real_t t0, real_t t1, HitRecord &rec) const
    {
        // Transform ray to sphere's local space
        Ray r = Ray(invMat.transform_point(ray.e), invMat.transform_vector(ray.d));
        
        if (!boundingBox->intersect(r, t0, t1)) {
            return false;
        }
        
        bool isHit = false;
        
        Vector3 e = r.e;
        Vector3 d = r.d;
        Vector3 c = invMat.transform_point(position);
        real_t R = radius;
        
        Vector3 ce = e - c;
        real_t dot_dce = dot(d, ce);
        real_t b_square = pow(dot_dce, 2);
        real_t dot_dd = dot(d, d);
        real_t ac_4 = dot_dd * (dot(ce, ce) - pow(R, 2));
        
        real_t discriminant = b_square - ac_4;
        real_t t = -1;
        
        if (discriminant < 0)
            isHit = false;
        else
            isHit = true;
        
        if (isHit) {
            
            real_t sqrt_discrim = sqrt(discriminant);
            real_t inv_dot_dd = 1.0/dot_dd;
            real_t dot_nd_ce  = dot(-d, ce);
            real_t t_1 = (dot_nd_ce + sqrt_discrim) * inv_dot_dd;
            real_t t_2 = (dot_nd_ce - sqrt_discrim) * inv_dot_dd;
            
            if(t_1 <= 0.0) {
                return false;
            }
            else if(t_2 <= 0.0)
                t = t_1;
            else
                t = t_2;

            if (t < t0 || t > t1) {
                isHit = false;
            }
            // If is hit and t is located between t0 and t1 and is ensured as the smaller value
            else {
                Vector3 localIntersectPoint = r.e + t * (r.d);
                Vector3 worldIntersectPoint = ray.e + t * (ray.d);
                Vector3 localNormal = localIntersectPoint - c;
                Vector3 worldNormal = normalize(normMat * localNormal);
                
                rec.type = eSphere;
                rec.position = worldIntersectPoint;
                rec.normal = worldNormal;
                
                rec.diffuse = material->diffuse;
                rec.ambient = material->ambient;
                rec.specular = material->specular;
                rec.phong = material->phong;
                
                int width, height;
                material->get_texture_size(&width, &height);
                
                if (width > 0 && height > 0) {
                    
                    real_t THETA = acos(rec.normal.y);
                    real_t PHI   = atan2(rec.normal.x, rec.normal.z);
                    real_t u = PHI/(2.0 * PI);
                    real_t v = (PI - THETA)/PI;
                    
                    rec.texture = material->get_texture_pixel(((int)(width*u)) % width, ((int)(height*v)) % height);
                }
                else {
                    rec.texture = Color3::White();
                }
                
                // To decide if there is not and, if there is, to calculate refractive ray
                rec.refractive_index = material->refractive_index;
                
                rec.t = t;
                
                rec.isLight = this->isLight;
                
            }
        }
        
        return isHit;
    }
    
} /* _462 */

