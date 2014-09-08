//
//  azReflection.cpp
//  Azurender
//
//  Created by Xiao Li on 8/4/14.
//
//

#ifndef __Azurender__azReflection__
#define __Azurender__azReflection__



#include "math/color.hpp"
#include "math/math.hpp"
#include "math/vector.hpp"
#include "raytracer/constants.h"
#include "raytracer/montecarlo.hpp"

namespace _462 {
    
    // PBRT P.426
    // The shading coordinate system is defined by the orthonormal basis vectors (s, t, n)
    // They lie along x, y, z axes in the coordinate system.
    // Direction vectors w(omiga) in world space are transformed into the shading coordinate system
    // before any of the BRDF or BTDF methods are called.
    
    // BSDF Inline Functions
    inline float CosTheta(const Vector3 &w) { return w.z; };
    
    inline float AbsCosTheta(const Vector3 &w) { return fabsf(w.z); }
    
    inline float SinTheta2(const Vector3 &w) {
        return std::max(0.f, 1.f - CosTheta(w) * CosTheta(w));
    }
    
    inline float SinTheta(const Vector3 &w) {
        return sqrtf(SinTheta2(w));
    }
    
    inline float CosPhi(const Vector3 &w) {
        float sintheta = SinTheta(w);
        if (sintheta == 0.f) return 1.f;
        return clamp((float)w.x / sintheta, -1.f, 1.f);
    }
    
    inline float SinPhi(const Vector3 &w) {
        float sintheta = SinTheta(w);
        if (sintheta == 0.f) return 0.f;
        return clamp((float)w.y / sintheta, -1.f, 1.f);
    }
    
    inline bool SameHemisphere(const Vector3 &w, const Vector3 &wp) {
        return w.z * wp.z > 0.f;
    }
    
    // BSDF Declarations
    enum BxDFType {
        BSDF_REFLECTION   = 1<<0,
        BSDF_TRANSMISSION = 1<<1,
        BSDF_DIFFUSE      = 1<<2,
        BSDF_GLOSSY       = 1<<3,
        BSDF_SPECULAR     = 1<<4,
        BSDF_ALL_TYPES        = BSDF_DIFFUSE |
                                BSDF_GLOSSY |
                                BSDF_SPECULAR,
        BSDF_ALL_REFLECTION   = BSDF_REFLECTION |
                                BSDF_ALL_TYPES,
        BSDF_ALL_TRANSMISSION = BSDF_TRANSMISSION |
                                BSDF_ALL_TYPES,
        BSDF_ALL              = BSDF_ALL_REFLECTION |
                                BSDF_ALL_TRANSMISSION
    };
    
    // TODO: BxDF, BRDF, BTDF, BSDF
    /*!
     @brief base class BxDF for BRDFs & BTDFs.
     @note  Not all BxDFs can be evaluated with the f() method. For example, perfectly specular
            objects like a mirror, glass, or water only scatter light from a single incident direction
            into a single outgoing direction. Such BxDFs are best described with delta distributions
            that are zero except for the single direction where light is scattered.
     */
    class BxDF {
    public:
        //BxDF Interface
        virtual ~BxDF() { }
        BxDF(BxDFType t) : type(t) { }
        
        bool MatchesFlags(BxDFType flags) const {
            return (type & flags) == type;
        }
        
        // Key method, returns the value of the distribution function for the given pair of directions
        virtual Color3 f(const Vector3 &wo, const Vector3 &wi) const = 0;
        /*!
         @brief This method is used both for handling scattering that is described by delta
                distributions as well as for randomly sampling directions from BxDFs that scatter light
                along multiple directions
         */
        virtual Color3 Sample_f(const Vector3 &wo, Vector3 *wi, float u1, float u2, float *pdf) const;
        
        /*!
         @brief hemispherical-directional reflectance is a 2D function that gives the total reflection in
                a given direction due to constant illumination over the hemisphere, or, equivalently, total
                reflection over the hemisphere due to light from a given direction. // Sec 14.5.5
         */
        virtual Color3 rho(const Vector3 &wo, int nSamples, const float *samples) const;
        
        /*!
         @brief The hemispherical-hemispherical reflectance of a surface, denoted by ρhh, is a constant
                spectral value that gives the fraction of incident light reflected by a surface when the
                incident light is the same from all directions.
         */
        virtual Color3 rho(int nSamples, const float *samples1, const float *samples2) const;
        
        virtual float  Pdf(const Vector3 &wi, const Vector3 &wo) const;
        
        // BxDF Public Data
        const BxDFType type;
    };
    
    // TODO: change ordinary surfaces to Lambertian Model, test BRDF implementation
    class Lambertian : public BxDF {
    public:
        // Lambertian public methods
        Lambertian(const Color3 &reflectance) : BxDF(BxDFType(BSDF_REFLECTION | BSDF_DIFFUSE)), R(reflectance) { }
        Color3 f(const Vector3 &wo, const Vector3 &wi) const;
        Color3 rho(const Vector3 &, int, const float *) const { return R; }
        Color3 rho(int, const float *, const float *) const { return R; }
    private:
        // Lamertian private data
        Color3 R;
    };
    
    
    class azFresnel;
    class azReflection;
    class azReflection {
        
    public:
        
        azReflection() {}
        ~azReflection() {}
        
        static Vector3 reflect(Vector3 i, Vector3 n)
        {
            return i - 2.0 * n * dot(n, i);
        }
        
        // tempoaryt
        static Vector3 refract(Vector3 d, Vector3 n, real_t ratio)
        {
            real_t dot_dn = dot(d, n);
            real_t square = real_t(1) - (pow(ratio, 2) * (real_t(1) - pow(dot_dn, 2)));
            
            if (square < 0) {
                return Vector3().Zero();
            }
            
            return (ratio * (d - (n * dot_dn))) - (n * sqrt(square));
        }
        
        // Incoming direction d
        static Vector3 refract(Vector3 d, Vector3 n, const float &etai, const float &etat)
        {
            float cosi = dot(d, n);
            cosi = clamp(cosi, -1.f, 1.f);
            
            // Compute indices of refraction for dielectric
            // From air to dieletric
            bool entering = cosi < 0.;
            float ei = etai, et = etat;
            
            // From dieletric to air
            if (!entering)
            {
                std::swap(ei, et);
                n = -n;
                cosi = -cosi;
            }
            
            float square = 1.f - (pow(ei/et, 2.f) * (1.f - pow(cosi, 2.f)));
            
            if (square < 0) {
                return Vector3().Zero();
            }
            
            return (ei/et * (d - (n * cosi))) - (n * sqrt(square));
        }
        
        static float TSMicrofacetCoeef(const Vector3 &wo, const Vector3 &wi, const float &etai, const float &etat);
        
        static float TSMicrofacetSample_f(const Vector3 &wo, const Vector3 &wi, const float &etai, const float &etat);
        
        
        static float G(const Vector3 &wo, const Vector3 &wi, const Vector3 &wh) {
            float NdotWh = AbsCosTheta(wh);
            float NdotWo = AbsCosTheta(wo);
            float NdotWi = AbsCosTheta(wi);
            float WOdotWh = fabsf(dot(wo, wh));
            return std::min(1.f, std::min((2.f * NdotWh * NdotWo / WOdotWh),
                                          (2.f * NdotWh * NdotWi / WOdotWh)));
        }
        
        static float D(const Vector3 &wh) {
            float ex = .01f, ey = .10f;
            float costhetah = AbsCosTheta(wh);
            float d = 1.f - costhetah * costhetah;
            if (d == 0.f) return 0.f;
            float e = (ex * wh.x * wh.x + ey * wh.y * wh.y) / d;
            return sqrtf((ex+2.f) * (ey+2.f)) * INV_TWOPI * powf(costhetah, e);
        }
        
    };
    
    class azFresnel {
        
    public:
        
        azFresnel() {}
        ~azFresnel() {}
        
        // PBRT version
        static float FrDiel(float cosi, float cost, const float &etai, const float &etat)
        {
            float Rparl = ((etat * cosi) - (etai * cost)) /
            ((etat * cosi) + (etai * cost));
            float Rperp = ((etai * cosi) - (etat * cost)) /
            ((etai * cosi) + (etat * cost));
            
            return (Rparl * Rparl + Rperp * Rperp) / 2.f;
        }
        
        // PBRT version
        static float FresnelDielectricEvaluate(float cosi, const float &etai, const float &etat)
        {
            // Compute Fresnel reflectance for dielectric
            cosi = clamp(cosi, -1.f, 1.f);
            
            // Compute indices of refraction for dielectric
            bool entering = cosi < 0.;
            float ei = etai, et = etat;
            
            if (!entering)  std::swap(ei, et);
            
            // Compute _sint_ using Snell's law
            float sint = ei/et * sqrtf(std::max(0.f, 1.f - cosi * cosi));
            if (sint >= 1.) {
                // Handle total internal reflection
                return 1.;
            }
            else {
                float cost = sqrtf(std::max(0.f, 1.f - sint*sint));
                return FrDiel(fabsf(cosi), cost, ei, et);
            }
        }
        
        // source: https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
        static float FresnelDielectricDielectric(float Eta, float CosTheta)
        {
            float SinTheta2 = 1 - CosTheta * CosTheta;
            
            float t0 = sqrt(1 - (SinTheta2 / (Eta * Eta)));
            float t1 = Eta * t0;
            float t2 = Eta * CosTheta;
            
            float rs = (CosTheta - t1) / (CosTheta + t1);
            float rp = (t0 - t2) / (t0 + t2);
            
            return 0.5 * (rs * rs + rp * rp);
        }
        
    };
    
    
}


#endif /* defined(__Azurender__azReflection__) */