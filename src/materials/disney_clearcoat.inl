#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyClearcoat &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return make_zero_spectrum();
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
    Real R0_15 = (1.5-1)*(1.5-1)/ ((1.5+1) * (1.5+1));
    Vector3 half_vector = normalize(dir_in+dir_out);
    Real Fc = R0_15+(1-R0_15)*pow(1-fabs(dot(half_vector, dir_out)), 5);
    Real gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha_g = (1 - gloss) * 0.1 + gloss * 0.001;
    Vector3 local_half = to_local(frame, half_vector);
    Real Dc = (alpha_g*alpha_g-1) / (c_PI * log(alpha_g*alpha_g) * (1 + (alpha_g*alpha_g-1) * local_half[2]*local_half[2]));
    Real G_in = G_m(dir_in, frame, 0.25, 0.25);
    Real G_out = G_m(dir_out, frame, 0.25, 0.25);
    Real Gc = G_in * G_out;
    Real final = Fc * Dc * Gc / (4 * fabs(dot(dir_in, frame.n)));
    return Spectrum(final, final, final);
}

Real pdf_sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0 ||
            dot(vertex.geometric_normal, dir_out) < 0) {
        // No light below the surface
        return 0;
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!

    Vector3 half_vector = normalize(dir_in+dir_out);
    Real gloss = eval(bsdf.clearcoat_gloss, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha_g = (1 - gloss) * 0.1 + gloss * 0.001;
    Vector3 local_half = to_local(frame, half_vector);
    Vector3 loc_dir_in = to_local(frame, dir_in);
    Vector3 loc_dir_out = to_local(frame, dir_out);
    Real Dc = (alpha_g*alpha_g-1) / (c_PI * log(alpha_g*alpha_g) * (1 + (alpha_g*alpha_g-1) * local_half[2]*local_half[2]));
    return Dc * fabs(dot(frame.n, half_vector)) / (4 * fabs(dot(half_vector, dir_out)));
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyClearcoat &bsdf) const {
    if (dot(vertex.geometric_normal, dir_in) < 0) {
        // No light below the surface
        return {};
    }
    // Flip the shading frame if it is inconsistent with the geometry normal
    Frame frame = vertex.shading_frame;
    if (dot(frame.n, dir_in) < 0) {
        frame = -frame;
    }
    // Homework 1: implement this!
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
        Real(0) /* eta */, Real(0.25) /* roughness */};
}

TextureSpectrum get_texture_op::operator()(const DisneyClearcoat &bsdf) const {
    return make_constant_spectrum_texture(make_zero_spectrum());
}
