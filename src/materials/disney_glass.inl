#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // (void)reflect; // silence unuse warning, remove this when implementing hw
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;

    Spectrum Ks = eval(
        bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Spectrum Kt = sqrt(Ks);
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real aspect = sqrt(1-0.9 * eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool));

    Vector3 half_vector;
    if (reflect) {
        half_vector = normalize(dir_in + dir_out);
    } else {
        // "Generalized half-vector" from Walter et al.
        // See "Microfacet Models for Refraction through Rough Surfaces"
        half_vector = normalize(dir_in + dir_out * eta);
    }

    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    // Clamp roughness to avoid numerical issues.
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    // Compute F / D / G
    // Note that we use the incoming direction
    // for evaluating the Fresnel reflection amount.
    // We can also use outgoing direction -- then we would need to
    // use 1/bsdf.eta and we will get the same result.
    // However, using the incoming direction allows
    // us to use F to decide whether to reflect or refract during sampling.
    Real h_dot_in = dot(half_vector, dir_in);
    Real h_dot_out = dot(half_vector, dir_out);
    Real F = fresnel_dielectric(fabs(h_dot_in), fabs(h_dot_out), eta);

    // hw test
    // Real F0 = ((1-eta) / (1+eta)) * ((1-eta) / (1+eta));
    // Real cos_theta_i = dot(dir_in, frame.n);
    // Real cos_theta_t_square = 1 - (1-cos_theta_i*cos_theta_i) / (eta*eta);
    // F = (cos_theta_t_square>0)?F0 + (1-F0) * pow(1-sqrt(cos_theta_t_square), 5):1;


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
    if (reflect) {
        return Ks * (F * D * G) / (4 * fabs(dot(frame.n, dir_in)));
    } else {
        // Snell-Descartes law predicts that the light will contract/expand 
        // due to the different index of refraction. So the normal BSDF needs
        // to scale with 1/eta^2. However, the "adjoint" of the BSDF does not have
        // the eta term. This is due to the non-reciprocal nature of the index of refraction:
        // f(wi -> wo) / eta_o^2 = f(wo -> wi) / eta_i^2
        // thus f(wi -> wo) = f(wo -> wi) (eta_o / eta_i)^2
        // The adjoint of a BSDF is defined as swapping the parameter, and
        // this cancels out the eta term.
        // See Chapter 5 of Eric Veach's thesis "Robust Monte Carlo Methods for Light Transport Simulation"
        // for more details.
        // Real eta_factor = dir == TransportDirection::TO_LIGHT ? (1 / (eta * eta)) : 1;
        Real sqrt_denom = h_dot_in + eta * h_dot_out;
        // Very complicated BSDF. See Walter et al.'s paper for more details.
        // "Microfacet Models for Refraction through Rough Surfaces"
        return Kt * ((1 - F) * D * G * fabs(h_dot_out * h_dot_in)) / 
            (fabs(dot(frame.n, dir_in)) * sqrt_denom * sqrt_denom);
    }
}

Real pdf_sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    bool reflect = dot(vertex.geometric_normal, dir_in) *
                   dot(vertex.geometric_normal, dir_out) > 0;
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    // (void)reflect; // silence unuse warning, remove this when implementing hw
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
    assert(eta > 0);

    Vector3 half_vector;
    if (reflect) {
        half_vector = normalize(dir_in + dir_out);
    } else {
        // "Generalized half-vector" from Walter et al.
        // See "Microfacet Models for Refraction through Rough Surfaces"
        half_vector = normalize(dir_in + dir_out * eta);
    }

    // Flip half-vector if it's below surface
    if (dot(half_vector, frame.n) < 0) {
        half_vector = -half_vector;
    }

    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    // Clamp roughness to avoid numerical issues.
    Real aspect = sqrt(1-0.9 * eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool));
    roughness = std::clamp(roughness, Real(0.01), Real(1));

    // We sample the visible normals, also we use F to determine
    // whether to sample reflection or refraction
    // so PDF ~ F * D * G_in for reflection, PDF ~ (1 - F) * D * G_in for refraction.
    Real h_dot_in = dot(half_vector, dir_in);
    Real h_dot_out = dot(half_vector, dir_out);
    Real F = fresnel_dielectric(fabs(h_dot_in), fabs(h_dot_out), eta);

    // Real F0 = ((1-eta) / (1+eta)) * ((1-eta) / (1+eta));
    // Real cos_theta_i = dot(dir_in, frame.n);
    // Real cos_theta_t_square = 1 - (1-cos_theta_i*cos_theta_i) / (eta*eta);
    // F = (cos_theta_t_square>0)?F0 + (1-F0) * pow(1-sqrt(cos_theta_t_square), 5):1;

    Real alpha_x = fmax(0.0001, pow(roughness, 2) / aspect);
    Real alpha_y = fmax(0.0001, pow(roughness, 2) * aspect);
    Vector3 local_half = to_local(frame, half_vector);
    Real h_x2 = local_half[0]*local_half[0];
    Real h_y2 = local_half[1]*local_half[1];
    Real h_z2 = local_half[2]*local_half[2];
    Real D = 1/(c_PI*alpha_x*alpha_y*pow(h_x2/alpha_x/alpha_x + h_y2/alpha_y/alpha_y + h_z2, 2));
    Real G_in = G_m(dir_in, frame, alpha_x, alpha_y);
    if (reflect) {
        return (F * D * G_in) / (4 * fabs(dot(frame.n, dir_in)));
    } else {
        Real h_dot_out = dot(half_vector, dir_out);
        Real sqrt_denom = h_dot_in + eta * h_dot_out;
        Real dh_dout = eta * eta * h_dot_out / (sqrt_denom * sqrt_denom);
        return (1 - F) * D * G_in * fabs(dh_dout * h_dot_in / dot(frame.n, dir_in));
    }
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyGlass &bsdf) const {
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) * dot(vertex.geometric_normal, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Real eta = dot(vertex.geometric_normal, dir_in) > 0 ? bsdf.eta : 1 / bsdf.eta;
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

    // Real F0 = ((1-eta) / (1+eta)) * ((1-eta) / (1+eta));
    // Real cos_theta_i = dot(dir_in, frame.n);
    // Real cos_theta_t_square = 1 - (1-cos_theta_i*cos_theta_i) / (eta*eta);
    // F = (cos_theta_t_square>0)?F0 + (1-F0) * pow(1-sqrt(cos_theta_t_square), 5):1;

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

TextureSpectrum get_texture_op::operator()(const DisneyGlass &bsdf) const {
    return bsdf.base_color;
}
