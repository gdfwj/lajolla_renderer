#include "../microfacet.h"

Spectrum eval_op::operator()(const DisneyMetal &bsdf) const {
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
    Vector3 half_vector = normalize(dir_in+dir_out);
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Spectrum Fm = base_color + (1 - base_color) * pow(1 - std::fabs(dot(half_vector, dir_out)), 5);
    Real aspect = sqrt(1-0.9 * eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool));
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
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
    return Fm * Dm * Gm / (4 * fabs(dot(dir_in, frame.n)));

    // return make_zero_spectrum();
}

Real pdf_sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
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
    Real aspect = sqrt(1-0.9 * eval(bsdf.anisotropic, vertex.uv, vertex.uv_screen_size, texture_pool));
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real alpha_x = fmax(0.0001, pow(roughness, 2) / aspect);
    Real alpha_y = fmax(0.0001, pow(roughness, 2) * aspect);
    Vector3 local_half = to_local(frame, half_vector);
    Vector3 loc_dir_in = to_local(frame, dir_in);
    Vector3 loc_dir_out = to_local(frame, dir_out);
    Real h_x2 = local_half[0]*local_half[0];
    Real h_y2 = local_half[1]*local_half[1];
    Real h_z2 = local_half[2]*local_half[2];
    Real Dm = 1/(c_PI*alpha_x*alpha_y*pow(h_x2/alpha_x/alpha_x + h_y2/alpha_y/alpha_y + h_z2, 2));
    Real G_in = G_m(dir_in, frame, alpha_x, alpha_y);
    return Dm * G_in / (4 * fabs(dot(dir_in, frame.n)));
}

std::optional<BSDFSampleRecord>
        sample_bsdf_op::operator()(const DisneyMetal &bsdf) const {
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

    
    // Sample from the specular lobe.

    // Convert the incoming direction to local coordinates
    Vector3 local_dir_in = to_local(frame, dir_in);
    Real roughness = eval(
        bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
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

TextureSpectrum get_texture_op::operator()(const DisneyMetal &bsdf) const {
    return bsdf.base_color;
}
