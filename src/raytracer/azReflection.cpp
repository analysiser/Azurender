//
//  azReflection.cpp
//  Azurender
//
//  Created by Xiao Li on 8/4/14.
//
//

#include "azReflection.hpp"



namespace _462 {
    
    float azReflection::TSMicrofacetCoeef(const _462::Vector3 &wo, const _462::Vector3 &wi, const float &etai, const float &etat)
    {
        float cosThetaO = AbsCosTheta(wo);
        float cosThetaI = AbsCosTheta(wi);
        if (cosThetaI == 0.f || cosThetaO == 0.f) return (0.f);
        Vector3 wh = wi + wo;
        if (wh.x == 0. && wh.y == 0. && wh.z == 0.) return (0.f);
        
        wh = normalize(wh);
        float cosThetaH = dot(wi, wh);
        
        ///?
//        if (dot(wo, wi) < 0.f) {
//            std::cout<<"!!!"<<std::endl;
//            return 0.f;
//        }
        
        // entering, need -cosThetaH??
        float F = azFresnel::FresnelDielectricEvaluate(-cosThetaH, etai, etat);
        return 1.f * D(wh) * G(wo, wi, wh) * F / (4.f * cosThetaI * cosThetaO);
    }
    
    float azReflection::TSMicrofacetSample_f(const _462::Vector3 &wo, const _462::Vector3 &wi, const float &etai, const float &etat)
    {
        if (!SameHemisphere(wo, wi)) return (0.f);
        return TSMicrofacetCoeef(wo, wi, etai, etat);
    }
    

#pragma mark BxDF implementation
    
    // BxDF::Sample_f() takes two sample values in the
    // range [0, 1)^2 that are intended to be used by a transformation-based sampling algorithm.
    Color3 BxDF::Sample_f(const Vector3 &wo, Vector3 *wi, float u1, float u2, float *pdf) const
    {
        // Cosine-sample the hemisphere, flipping the direction if necessary
        *wi = CosineSampleHemisphere(u1, u2);
        if (wo.z < 0.) wi->z *= -1.f;
        *pdf = Pdf(wo, *wi);
        return f(wo, *wi);
    }
    
    float BxDF::Pdf(const Vector3 &wi, const  Vector3 &wo) const
    {
        return SameHemisphere(wo, wi) ? AbsCosTheta(wi) * INV_PI : 0.f;
    }
    
    
    Color3 BxDF::rho(const _462::Vector3 &w, int nSamples, const float *samples) const
    {
        // TODO: see what is going on...
        Color3 r = Color3::Black();
        for (int i = 0; i < nSamples; ++i) {
            // Estimate one term of $\rho_\roman{hd}$
            Vector3 wi;
            float pdf = 0.f;
            Color3 f = Sample_f(w, &wi, samples[2*i], samples[2*i+1], &pdf);
            if (pdf > 0.) r += f * AbsCosTheta(wi) / pdf;
        }
        return r / float(nSamples);
    }
    
    Color3 BxDF::rho(int nSamples, const float *samples1, const float *samples2) const
    {
        // TODO: see what is going on...
        Color3 r = Color3::Black();
        for (int i = 0; i < nSamples; ++i) {
            // Estimate one term of $\rho_\roman{hh}$
            Vector3 wo, wi;
            wo = UniformSampleHemisphere(samples1[2*i], samples1[2*i+1]);
            float pdf_o = INV_TWOPI, pdf_i = 0.f;
            Color3 f = Sample_f(wo, &wi, samples2[2*i], samples2[2*i+1], &pdf_i);
            if (pdf_i > 0.)
                r += f * AbsCosTheta(wi) * AbsCosTheta(wo) / (pdf_o * pdf_i);
        }
        return r / (M_PI*nSamples);
    }
    
    Color3 Lambertian::f(const Vector3 & /*wo*/, const Vector3 & /*wi*/ ) const {
        return R * INV_PI;
    }
    
    
    
}