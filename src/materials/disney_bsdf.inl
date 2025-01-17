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
    Real subsurface_scale = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real specular_tint = eval(bsdf.specular_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real anisotropic = eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat = eval(bsdf.clearcoat, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real clearcoat_gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen = eval(bsdf.sheen, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real sheen_tint = eval(bsdf.sheen_tint, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    Real specular = eval(bsdf.specular, vertex.uv, vertex.uv_screen_size, texture_pool);
    Vector3 half_vector;
    if (dot(vertex.geometric_normal, dir_in) <= 0 || !reflect) {
        if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
            frame = -frame;
        }
        half_vector = normalize(dir_in + dir_out * eta);
        if (dot(half_vector, frame.n) < 0) {
            half_vector = -half_vector;
        }
        Real h_dot_in = dot(half_vector, dir_in);
        Real h_dot_out = dot(half_vector, dir_out);
        Real F = fresnel_dielectric(fabs(h_dot_in), fabs(h_dot_out), eta);
        Real aspect = sqrt(1-0.9 * eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool));
        Real alpha_x = fmax(0.0001, pow(roughness, 2) / aspect);
        Real alpha_y = fmax(0.0001, pow(roughness, 2) * aspect);
        Vector3 local_half = to_local(frame, half_vector);
        Real h_x2 = local_half[0]*local_half[0];
        Real h_y2 = local_half[1]*local_half[1];
        Real h_z2 = local_half[2]*local_half[2];
        Real D = 1/(c_PI*alpha_x*alpha_y*pow(h_x2/alpha_x/alpha_x + h_y2/alpha_y/alpha_y + h_z2, 2));
        Real G_in = G_m(dir_in, frame, alpha_x, alpha_y);
        Real G_out = G_m(dir_out, frame, alpha_x, alpha_y);
        Real G = G_in * G_out;
        Real sqrt_denom = h_dot_in + eta * h_dot_out;
        // Very complicated BSDF. See Walter et al.'s paper for more details.
        // "Microfacet Models for Refraction through Rough Surfaces"
        return (1-metallic)*specular_transmission * (sqrt(base_color) * ((1 - F) * D * G * fabs(h_dot_out * h_dot_in)) / 
            (fabs(dot(frame.n, dir_in)) * sqrt_denom * sqrt_denom));
    }
    if (dot(frame.n, dir_in)<0) {
        frame = -frame;
    }
    half_vector = normalize(dir_in + dir_out);
    Real h_dot_in = dot(half_vector, dir_in);
    Real h_dot_out = dot(half_vector, dir_out);
    Real Fg = fresnel_dielectric(fabs(h_dot_in), fabs(h_dot_out), eta);
    Real aspect = sqrt(1-0.9 * eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool));
    Real alpha_x = fmax(0.0001, pow(roughness, 2) / aspect);
    Real alpha_y = fmax(0.0001, pow(roughness, 2) * aspect);
    Vector3 local_half = to_local(frame, half_vector);
    Real h_x2 = local_half[0]*local_half[0];
    Real h_y2 = local_half[1]*local_half[1];
    Real h_z2 = local_half[2]*local_half[2];
    Real Dg = 1/(c_PI*alpha_x*alpha_y*pow(h_x2/alpha_x/alpha_x + h_y2/alpha_y/alpha_y + h_z2, 2));
    Real G_in = G_m(dir_in, frame, alpha_x, alpha_y);
    Real G_out = G_m(dir_out, frame, alpha_x, alpha_y);
    Real Gg = G_in * G_out;
    Spectrum f_glass = base_color * (Fg * Dg * Gg) / (4 * fabs(dot(frame.n, dir_in)));

    Real F_D90 = 0.5+ 2* eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool) * (dot(half_vector, dir_out)*dot(half_vector, dir_out));
    Real F_in = 1+(F_D90-1)*std::pow((1-dot(frame.n, dir_in)), 5);
    Real F_out = 1+(F_D90-1)*std::pow((1-dot(frame.n, dir_out)), 5);
    Spectrum base_diffuse = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool)/c_PI*F_in*F_out*std::fabs(dot(frame.n, dir_out));

    Real F_SS90 = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool) * (dot(half_vector, dir_out)*dot(half_vector, dir_out));
    Real F_SSin = 1+(F_SS90-1)*std::pow((1-std::fabs(dot(frame.n, dir_in))), 5);
    Real F_SSout = 1+(F_SS90-1)*std::pow((1-std::fabs(dot(frame.n, dir_out))), 5);
    Spectrum subsurface = 1.25*eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool)/c_PI*(F_SSin*F_SSout*(1.0/(std::fabs(dot(frame.n, dir_out))+std::fabs(dot(frame.n, dir_in)))-0.5)+0.5)*std::fabs(dot(frame.n, dir_out));
    Spectrum f_diffuse = (1-subsurface_scale)*base_diffuse+subsurface_scale*subsurface;

    Real luminance = 0.2126 * base_color[0] + 0.7152 * base_color[1] + 0.0722 * base_color[2];
    Spectrum Ctint = luminance>0? base_color / luminance : Vector3(1, 1, 1);
    Spectrum Ks = (1-specular_tint)+specular_tint * Ctint;
    Spectrum C0 = specular * pow(eta-1, 2) / pow(eta+1, 2) * (1-metallic) * base_color + metallic * base_color;
    Spectrum Fm = C0 + (1-C0)*pow(1-dot(half_vector, dir_out), 5);
    // Spectrum Fm = base_color + (1 - base_color) * pow(1 - std::fabs(dot(half_vector, dir_out)), 5);
    Real Dm = 1/(c_PI*alpha_x*alpha_y*pow(h_x2/alpha_x/alpha_x + h_y2/alpha_y/alpha_y + h_z2, 2));
    Real Gm = G_in * G_out;
    Spectrum f_metal = Fm * Dm * Gm / (4 * fabs(dot(dir_in, frame.n)));
    
    Real R0_15 = (1.5-1)*(1.5-1)/ ((1.5+1) * (1.5+1));
    Real Fc = R0_15+(1-R0_15)*pow(1-fabs(dot(half_vector, dir_out)), 5);
    Real gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha_g = (1 - gloss) * 0.1 + gloss * 0.001;
    Real Dc = (alpha_g*alpha_g-1) / (c_PI * log(alpha_g*alpha_g) * (1 + (alpha_g*alpha_g-1) * local_half[2]*local_half[2]));
    Real G_inc = G_m(dir_in, frame, 0.25, 0.25);
    Real G_outc = G_m(dir_out, frame, 0.25, 0.25);
    Real Gc = G_inc * G_outc;
    Real final = Fc * Dc * Gc / (4 * fabs(dot(dir_in, frame.n)));
    Spectrum f_clearcoat = Spectrum(final, final, final);

    Spectrum Csheen = (1-sheen_tint) * Ctint + sheen_tint;
    Spectrum f_sheen = Csheen * pow(1-fabs(dot(half_vector, dir_out)), 5) * fabs(dot(frame.n, dir_out));
    
    // return f_diffuse;
    return (1-specular_transmission) * (1-metallic) * f_diffuse 
        + (1-metallic) * sheen *f_sheen
        + (1-specular_transmission * (1-metallic)) * f_metal
        + 0.25*clearcoat * f_clearcoat
        + (1-metallic) * specular_transmission * f_glass;
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
    if (dot(vertex.geometric_normal, dir_in) <= 0 || !reflect) {
        if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
            frame = -frame;
        }
        half_vector = normalize(dir_in + dir_out * eta);
        if (dot(half_vector, frame.n) < 0) {
            half_vector = -half_vector;
        }
        // Clamp roughness to avoid numerical issues.
        Real aspect = sqrt(1-0.9 * eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool));
        roughness = std::clamp(roughness, Real(0.01), Real(1));

        // We sample the visible normals, also we use F to determine
        // whether to sample reflection or refraction
        // so PDF ~ F * D * G_in for reflection, PDF ~ (1 - F) * D * G_in for refraction.
        Real h_dot_in = dot(half_vector, dir_in);
        Real h_dot_out = dot(half_vector, dir_out);
        Real F = fresnel_dielectric(fabs(h_dot_in), fabs(h_dot_out), eta);
        Real alpha_x = fmax(0.0001, pow(roughness, 2) / aspect);
        Real alpha_y = fmax(0.0001, pow(roughness, 2) * aspect);
        Vector3 local_half = to_local(frame, half_vector);
        Real h_x2 = local_half[0]*local_half[0];
        Real h_y2 = local_half[1]*local_half[1];
        Real h_z2 = local_half[2]*local_half[2];
        Real D = 1/(c_PI*alpha_x*alpha_y*pow(h_x2/alpha_x/alpha_x + h_y2/alpha_y/alpha_y + h_z2, 2));
        Real G_in = G_m(dir_in, frame, alpha_x, alpha_y);
        Real sqrt_denom = h_dot_in + eta * h_dot_out;
        Real dh_dout = eta * eta * h_dot_out / (sqrt_denom * sqrt_denom);
        Real s_glass = (1 - F) * D * G_in * fabs(dh_dout * h_dot_in / dot(frame.n, dir_in));
        assert(s_glass!=0);
        return s_glass;
    }
    if (dot(frame.n, dir_in)<0) {
        frame = -frame;
    }
    half_vector = normalize(dir_in + dir_out);
    Real h_dot_in = dot(half_vector, dir_in);
    Real h_dot_out = dot(half_vector, dir_out);
    Real F = fresnel_dielectric(fabs(h_dot_in), fabs(h_dot_out), eta);
    Real aspect = sqrt(1-0.9 * eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool));
    Real alpha_x = fmax(0.0001, pow(roughness, 2) / aspect);
    Real alpha_y = fmax(0.0001, pow(roughness, 2) * aspect);
    Vector3 local_half = to_local(frame, half_vector);
    Real h_x2 = local_half[0]*local_half[0];
    Real h_y2 = local_half[1]*local_half[1];
    Real h_z2 = local_half[2]*local_half[2];
    Real D = 1/(c_PI*alpha_x*alpha_y*pow(h_x2/alpha_x/alpha_x + h_y2/alpha_y/alpha_y + h_z2, 2));
    Real G_in = G_m(dir_in, frame, alpha_x, alpha_y);
    Real s_glass = (F * D * G_in) / (4 * fabs(dot(frame.n, dir_in)));

    Real s_diffuse = fmax(dot(frame.n, dir_out), Real(0)) / c_PI;

    Real Dm = 1/(c_PI*alpha_x*alpha_y*pow(h_x2/alpha_x/alpha_x + h_y2/alpha_y/alpha_y + h_z2, 2));
    Real s_metal =  Dm * G_in / (4 * fabs(dot(dir_in, frame.n)));
    
    Real gloss = clearcoat_gloss;
    Real alpha_g = (1 - gloss) * 0.1 + gloss * 0.001;
    Vector3 loc_dir_in = to_local(frame, dir_in);
    Vector3 loc_dir_out = to_local(frame, dir_out);
    Real Dc = (alpha_g*alpha_g-1) / (c_PI * log(alpha_g*alpha_g) * (1 + (alpha_g*alpha_g-1) * local_half[2]*local_half[2]));
    Real s_clearcoat = Dc * fabs(dot(frame.n, half_vector)) / (4 * fabs(dot(half_vector, dir_out)));
    // return s_diffuse;
    Real ret = (1-specular_transmission) * (1-metallic) * s_diffuse 
        + (1-specular_transmission * (1-metallic)) * s_metal
        + 0.25*clearcoat * s_clearcoat
        + (1-metallic) * specular_transmission * s_glass;
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
    Real rnd = rnd_param_w;
    
    // test
    // test
    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));
    // // Sample a micro normal and transform it to world space -- this is our half-vector.
    Real alpha = roughness * roughness;
    // return BSDFSampleRecord{
    //         to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
    //         Real(0) 
    //     };
    
    // test
    if (dot(vertex.geometric_normal, dir_in) <= 0) {
        if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
            frame = -frame;
        }
        Vector3 local_dir_in = to_local(frame, dir_in);
        Vector3 local_micro_normal =
            sample_visible_normals(local_dir_in, alpha, rnd_param_uv);

        Vector3 half_vector = to_world(frame, local_micro_normal);
        // Flip half-vector if it's below surface
        if (dot(half_vector, frame.n) < 0) {
            half_vector = -half_vector;
        }
        Real h_dot_in = dot(half_vector, dir_in);
        Real h_dot_out_sq = 1 - (1 - h_dot_in * h_dot_in) / (eta * eta);
        if (h_dot_out_sq <= 0) {
            // Total internal reflection
            // This shouldn't really happen, as F will be 1 in this case.
            Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
            return BSDFSampleRecord{reflected, Real(0) /* eta */, roughness};
        }
        Real h_dot_out = sqrt(h_dot_out_sq);
        Real F = fresnel_dielectric(fabs(h_dot_in), fabs(h_dot_out), eta);

        if (rnd_param_w <= F) {
            // Reflection
            Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
            // set eta to 0 since we are not transmitting
            return BSDFSampleRecord{reflected, Real(0) /* eta */, roughness};
        } else {
            // Refraction
            // https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
            // (note that our eta is eta2 / eta1, and l = -dir_in)
            
            // flip half_vector if needed
            if (h_dot_in < 0) {
                half_vector = -half_vector;
            }
            Real h_dot_out= sqrt(h_dot_out_sq);
            Vector3 refracted = -dir_in / eta + (fabs(h_dot_in) / eta - h_dot_out) * half_vector;
            assert(dot(refracted, half_vector) * dot(dir_in, half_vector)<0);
            return BSDFSampleRecord{refracted, eta, roughness};
        }
    }

    // if (dot(frame.n, dir_in)) {
    //     frame = -frame;
    // }
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }

    Vector3 local_dir_in = to_local(frame, dir_in);
    Vector3 local_micro_normal =
        sample_visible_normals(local_dir_in, alpha, rnd_param_uv);

    Vector3 half_vector = to_world(frame, local_micro_normal);
    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    //test

        // Now we need to decide whether to reflect or refract.
        // We do this using the Fresnel term.
        // Real rnd_r = (rnd - (diffuse_w + metal_w + glass_w)) / (1-(diffuse_w + metal_w + glass_w));
        // Real h_dot_in = dot(half_vector, dir_in);
        // Real h_dot_out_sq = 1 - (1 - h_dot_in * h_dot_in) / (eta * eta);
        // if (h_dot_out_sq <= 0) {
        //     // Total internal reflection
        //     // This shouldn't really happen, as F will be 1 in this case.
        //     Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
        //     return BSDFSampleRecord{reflected, Real(0) /* eta */, roughness};
        // }
        // Real h_dot_out = sqrt(h_dot_out_sq);
        // Real F = fresnel_dielectric(fabs(h_dot_in), fabs(h_dot_out), eta);

        // if (rnd_param_w <= F) {
        //     // Reflection
        //     Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
        //     // set eta to 0 since we are not transmitting
        //     return BSDFSampleRecord{reflected, Real(0) /* eta */, roughness};
        // } else {
        //     // Refraction
        //     // https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
        //     // (note that our eta is eta2 / eta1, and l = -dir_in)
            
        //     // flip half_vector if needed
        //     if (h_dot_in < 0) {
        //         half_vector = -half_vector;
        //     }
        //     Real h_dot_out= sqrt(h_dot_out_sq);
        //     Vector3 refracted = -dir_in / eta + (fabs(h_dot_in) / eta - h_dot_out) * half_vector;
        //     assert(dot(refracted, half_vector) * dot(dir_in, half_vector)<0);
        //     return BSDFSampleRecord{refracted, eta, roughness};
        // }
    //test

    if (rnd <= diffuse_w) {
        return BSDFSampleRecord{
            to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
            Real(0) 
        };
    }
    else if (rnd <= diffuse_w + metal_w) {
        Vector3 local_dir_in = to_local(frame, dir_in);
        // Clamp roughness to avoid numerical issues.
        roughness = std::clamp(roughness, Real(0.01), Real(1));
        Real alpha = roughness * roughness;
        Vector3 local_micro_normal =
            sample_visible_normals(local_dir_in, alpha, rnd_param_uv);
        
        // Transform the micro normal to world space
        Vector3 half_vector = to_world(frame, local_micro_normal);
        // Reflect over the world space normal
        Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
        return BSDFSampleRecord{
            reflected,
            Real(0) /* eta */, roughness /* roughness */
        };
    }
    else if (rnd <= diffuse_w + metal_w + glass_w) {
        Real roughness = eval(
            bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
        // Clamp roughness to avoid numerical issues.
        roughness = std::clamp(roughness, Real(0.01), Real(1));
        // Sample a micro normal and transform it to world space -- this is our half-vector.
        Real alpha = roughness * roughness;
        Vector3 local_dir_in = to_local(frame, dir_in);
        Vector3 local_micro_normal =
            sample_visible_normals(local_dir_in, alpha, rnd_param_uv);

        Vector3 half_vector = to_world(frame, local_micro_normal);
        // Flip half-vector if it's below surface
        if (dot(half_vector, frame.n) < 0) {
            half_vector = -half_vector;
        }

        // Now we need to decide whether to reflect or refract.
        // We do this using the Fresnel term.
        Real rnd_r = (rnd - (diffuse_w + metal_w + glass_w)) / (1-(diffuse_w + metal_w + glass_w));
        Real h_dot_in = dot(half_vector, dir_in);
        Real h_dot_out_sq = 1 - (1 - h_dot_in * h_dot_in) / (eta * eta);
        if (h_dot_out_sq <= 0) {
            // Total internal reflection
            // This shouldn't really happen, as F will be 1 in this case.
            Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
            return BSDFSampleRecord{reflected, Real(0) /* eta */, roughness};
        }
        Real h_dot_out = sqrt(h_dot_out_sq);
        Real F = fresnel_dielectric(fabs(h_dot_in), fabs(h_dot_out), eta);

        if (rnd_r <= F) {
            // Reflection
            Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_vector) * half_vector);
            // set eta to 0 since we are not transmitting
            return BSDFSampleRecord{reflected, Real(0) /* eta */, roughness};
        } else {
            // Refraction
            // https://en.wikipedia.org/wiki/Snell%27s_law#Vector_form
            // (note that our eta is eta2 / eta1, and l = -dir_in)
            
            // flip half_vector if needed
            if (h_dot_in < 0) {
                half_vector = -half_vector;
            }
            Real h_dot_out= sqrt(h_dot_out_sq);
            Vector3 refracted = -dir_in / eta + (fabs(h_dot_in) / eta - h_dot_out) * half_vector;
            assert(dot(refracted, half_vector) * dot(dir_in, half_vector)<0);
            return BSDFSampleRecord{refracted, eta, roughness};
        }
    }
    else {
        Real gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
        Real alpha_g = (1 - gloss) * 0.1 + gloss * 0.001;
        Real cos_elevation = sqrt((1-pow(alpha_g*alpha_g, 1 - rnd_param_uv[0])) / (1-alpha_g*alpha_g));
        Real sin_elevation = sqrt(1-cos_elevation*cos_elevation);
        Real azimuth = 2 * c_PI * rnd_param_uv[1];
        Vector3 half = Vector3(sin_elevation*cos(azimuth), sin_elevation*sin(azimuth), cos_elevation);
        Vector3 half_world = to_world(frame, half);
        Vector3 reflected = normalize(-dir_in + 2 * dot(dir_in, half_world) * half_world);
        return BSDFSampleRecord{
            reflected,
            Real(0) /* eta */, Real(0.25) /* roughness */
        };
    }
}

TextureSpectrum get_texture_op::operator()(const DisneyBSDF &bsdf) const {
    return bsdf.base_color;
}
