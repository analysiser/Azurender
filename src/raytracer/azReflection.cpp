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
    
    
    Color3 BxDF::f(const Vector3 &wo, const Vector3 &wi) const
    {

    }
    
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
        
    }
    
    Color3 Lambertian::f(const Vector3 & /*wo*/, const Vector3 & /*wi*/ ) const {
        return R * INV_PI;
    }
    
}