#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    
    // Homework 1: implement this!
    // (void)reflect; // silence unuse warning, remove this when implementing hw
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    DisneyGlass glass_bsdf = DisneyGlass{bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta};
    DisneyDiffuse diffuse_bsdf = DisneyDiffuse{bsdf.base_color, bsdf.roughness, bsdf.subsurface};
    DisneyClearcoat clearcoat_bsdf = DisneyClearcoat{bsdf.clearcoat};
    DisneySheen sheen_bsdf = DisneySheen{bsdf.base_color, bsdf.sheen_tint};
    DisneyMetal metal_bsdf = DisneyMetal{bsdf.base_color, bsdf.roughness, bsdf.anisotropic};
    
    Spectrum f_glass = eval(glass_bsdf, dir_in, dir_out, vertex, texture_pool, dir);

    Spectrum f_diffuse = eval(diffuse_bsdf, dir_in, dir_out, vertex, texture_pool, dir);

    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    Real eta = bsdf.eta;
    Vector3 half_vector = normalize(dir_in + dir_out);
    Real luminance_ = luminance(base_color);
    Spectrum Ctint = luminance_>0? base_color / luminance_ : Vector3(1, 1, 1);
    Spectrum Ks = (1-specular_tint)+specular_tint * Ctint;
    Spectrum C0 = specular * pow(eta-1, 2) / pow(eta+1, 2) * (1-metallic) * Ks + metallic * base_color;
    Spectrum Fm = C0 + (1-C0)*pow(1-dot(half_vector, dir_out), 5);
    // Spectrum Fm = base_color + (1 - base_color) * pow(1 - std::fabs(dot(half_vector, dir_out)), 5);
    Real aspect = sqrt(1-0.9 * anisotropic);
    Real alpha_x = fmax(0.0001, pow(roughness, 2) / aspect);
    Real alpha_y = fmax(0.0001, pow(roughness, 2) * aspect);
    Vector3 local_half = to_local(frame, half_vector);
    Real h_x2 = local_half[0]*local_half[0];
    Real h_y2 = local_half[1]*local_half[1];
    Real h_z2 = local_half[2]*local_half[2];
    Real Dm = 1/(c_PI*alpha_x*alpha_y*pow(h_x2/alpha_x/alpha_x + h_y2/alpha_y/alpha_y + h_z2, 2));
    Real G_in = G_m(dir_in, frame, alpha_x, alpha_y);
    Real G_out = G_m(dir_out, frame, alpha_x, alpha_y);
    Real Gm = G_in * G_out;
    Spectrum f_metal = Fm * Dm * Gm / (4 * fabs(dot(dir_in, frame.n)));
    // Spectrum F_past = base_color + (1-base_color) * pow(1 - fabs(dot(half_vector, dir_out)), 5);
    // Spectrum f_metal = eval(metal_bsdf, dir_in, dir_out, vertex, texture_pool, dir) / F_past * Fm;
    
    Spectrum f_clearcoat = eval(clearcoat_bsdf, dir_in, dir_out, vertex, texture_pool, dir);

    Spectrum f_sheen = eval(sheen_bsdf, dir_in, dir_out, vertex, texture_pool, dir);

    // return f_diffuse;
    if (dot(vertex.geometric_normal, dir_in) <= 0 || !reflect) {
        return (1-metallic)*specular_transmission * f_glass;
        // return (1-metallic)*specular_transmission * (sqrt(base_color) * ((1 - F) * D * G * fabs(h_dot_out * h_dot_in)) / 
        //     (fabs(dot(frame.n, dir_in)) * sqrt_denom * sqrt_denom));
    }
    // return f_metal;
    Spectrum ret = (1-specular_transmission) * (1-metallic) * f_diffuse 
        + (1-metallic) * sheen *f_sheen
        + (1-specular_transmission * (1-metallic)) * f_metal
        + 0.25*clearcoat * f_clearcoat
        + (1-metallic) * specular_transmission * f_glass;
    if(!(ret[0]>0 && ret[1]>0 && ret[2]>0)) {
        printf("diffuse: %f, %f, %f\n", f_diffuse[0], f_diffuse[1], f_diffuse[2]);
        printf("%f, %f, %f\n", f_sheen[0], f_sheen[1], f_sheen[2]);
        printf("%f, %f, %f\n", f_metal[0], f_metal[1], f_metal[2]);
        printf("%f, %f, %f\n", f_clearcoat[0], f_clearcoat[1], f_clearcoat[2]);
        printf("%f, %f, %f\n", f_glass[0], f_glass[1], f_glass[2]);
        assert((ret[0]>0 && ret[1]>0 && ret[2]>0));
    }
    return ret;
}

Real pdf_sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    
    // Homework 1: implement this!
    // (void)reflect; // silence unuse warning, remove this when implementing hw
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    Vector3 half_vector;
    DisneyGlass glass_bsdf = DisneyGlass{bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta};
    DisneyDiffuse diffuse_bsdf = DisneyDiffuse{bsdf.base_color, bsdf.roughness, bsdf.subsurface};
    DisneyClearcoat clearcoat_bsdf = DisneyClearcoat{bsdf.clearcoat};
    DisneySheen sheen_bsdf = DisneySheen{bsdf.base_color, bsdf.sheen_tint};
    DisneyMetal metal_bsdf = DisneyMetal{bsdf.base_color, bsdf.roughness, bsdf.anisotropic};
    // return pdf_sample_bsdf(metal_bsdf, dir_in, dir_out, vertex, texture_pool, dir);
    
    
    Real s_diffuse = pdf_sample_bsdf(diffuse_bsdf, dir_in, dir_out, vertex, texture_pool, dir);
    Real s_metal = pdf_sample_bsdf(metal_bsdf, dir_in, dir_out, vertex, texture_pool, dir);
    Real s_clearcoat = pdf_sample_bsdf(clearcoat_bsdf, dir_in, dir_out, vertex, texture_pool, dir);
    Real s_glass = pdf_sample_bsdf(glass_bsdf, dir_in, dir_out, vertex, texture_pool, dir);

    // return s_metal;
    if (dot(vertex.geometric_normal, dir_in) <= 0 || !reflect) {
        return s_glass;
    }

    Real diffuse_weight = (1-specular_transmission) * (1-metallic);
    Real metal_weight = (1-specular_transmission * (1-metallic));
    Real clearcoat_weight =0.25*clearcoat;
    Real glass_weight = (1-metallic) * specular_transmission;
    Real weight_sum = diffuse_weight + metal_weight + clearcoat_weight + glass_weight;
    diffuse_weight /= weight_sum;
    metal_weight /= weight_sum;
    clearcoat_weight /= weight_sum;
    glass_weight /= weight_sum;

    // return s_diffuse;
    Real ret = diffuse_weight * s_diffuse 
        + metal_weight* s_metal
        + clearcoat_weight * s_clearcoat
        + glass_weight* s_glass;
    assert(ret!=0);
    return ret;
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyBSDF &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    
    // Homework 1: implement this!
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real metallic = eval(bsdf.metallic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_transmission = eval(bsdf.specular_transmission, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

    Real diffuse_w = (1-metallic) * (1-specular_transmission);
    Real metal_w = (1-specular_transmission * (1-metallic));
    Real glass_w = (1-metallic) * specular_transmission;
    Real clearcoat_w = 0.25 * clearcoat;
    Real sum_w = diffuse_w + metal_w + glass_w + clearcoat_w;
    diffuse_w /= sum_w;
    metal_w /= sum_w;
    glass_w /= sum_w;
    clearcoat_w /= sum_w;
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    Real alpha = roughness * roughness;

    DisneyGlass glass_bsdf = DisneyGlass{bsdf.base_color, bsdf.roughness, bsdf.anisotropic, bsdf.eta};
    DisneyDiffuse diffuse_bsdf = DisneyDiffuse{bsdf.base_color, bsdf.roughness, bsdf.subsurface};
    DisneyClearcoat clearcoat_bsdf = DisneyClearcoat{bsdf.clearcoat};
    DisneySheen sheen_bsdf = DisneySheen{bsdf.base_color, bsdf.sheen_tint};
    DisneyMetal metal_bsdf = DisneyMetal{bsdf.base_color, bsdf.roughness, bsdf.anisotropic};
    
    // test
    // test
    // Clamp roughness to avoid numerical issues.
    // // Sample a micro normal and transform it to world space -- this is our half-vector.
    // return BSDFSampleRecord{
    //         to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
    //         Real(0) 
    //     };
    
    // test
    // return sample_bsdf(diffuse_bsdf, dir_in, vertex, texture_pool, rnd_param_uv, rnd_param_w, dir);
    // return sample_bsdf(metal_bsdf, dir_in, vertex, texture_pool, rnd_param_uv, rnd_param_w, dir);
    if (dot(vertex.geometric_normal, dir_in) <= 0) {
        return sample_bsdf(glass_bsdf, dir_in, vertex, texture_pool, rnd_param_uv, rnd_param_w, dir);
    }

    if (rnd_param_w <= diffuse_w) {
        return sample_bsdf(diffuse_bsdf, dir_in, vertex, texture_pool, rnd_param_uv, rnd_param_w, dir);
    }
    else if (rnd_param_w <= diffuse_w + metal_w) {
        return sample_bsdf(metal_bsdf, dir_in, vertex, texture_pool, rnd_param_uv, rnd_param_w, dir);
    }
    else if (rnd_param_w <= diffuse_w + metal_w + glass_w) {
        Real rnd_r = (rnd_param_w - (diffuse_w + metal_w)) / glass_w;
        return sample_bsdf(glass_bsdf, dir_in, vertex, texture_pool, rnd_param_uv, rnd_r, dir);
    }
    else {
        return sample_bsdf(clearcoat_bsdf, dir_in, vertex, texture_pool, rnd_param_uv, rnd_param_w, dir);
    }
}

TextureSpectrum get_texture_op::operator()(const DisneyBSDF &bsdf) const {
    return bsdf.base_color;
}
