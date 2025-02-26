Spectrum eval_op::operator()(const DisneyDiffuse &bsdf) const {
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
    // return fmax(dot(frame.n, dir_out), Real(0)) *
    //        eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool) / c_PI;
    // Homework 1: implement this!
    Real roughness = eval(bsdf.roughness, vertex.uv, vertex.uv_screen_size, texture_pool);
    roughness = std::clamp(roughness, 0.01, 1.0);
    Spectrum base_color = eval(bsdf.base_color, vertex.uv, vertex.uv_screen_size, texture_pool);
    Real subsurface_scale = eval(bsdf.subsurface, vertex.uv, vertex.uv_screen_size, texture_pool);
    Vector3 half_vector = normalize(dir_in+dir_out);
    Real F_D90 = 0.5+ 2* roughness * (dot(half_vector, dir_out)*dot(half_vector, dir_out));
    Real F_in = 1+(F_D90-1)*std::pow(1-fabs(dot(frame.n, dir_in)), 5);
    Real F_out = 1+(F_D90-1)*std::pow(1-fabs(dot(frame.n, dir_out)), 5);
    Spectrum base_diffuse = base_color/c_PI*F_in*F_out*std::fabs(dot(frame.n, dir_out));

    Real F_SS90 = roughness * (dot(half_vector, dir_out)*dot(half_vector, dir_out));
    Real F_SSin = 1+(F_SS90-1)*std::pow((1-std::fabs(dot(frame.n, dir_in))), 5);
    Real F_SSout = 1+(F_SS90-1)*std::pow((1-std::fabs(dot(frame.n, dir_out))), 5);
    Spectrum subsurface = 1.25*base_color/c_PI*(F_SSin*F_SSout*(1.0/(std::fabs(dot(frame.n, dir_out))+std::fabs(dot(frame.n, dir_in)))-0.5)+0.5)*std::fabs(dot(frame.n, dir_out));
    Spectrum ret = (1-subsurface_scale)*base_diffuse+subsurface_scale*subsurface;
    if (!(ret[0]>0 && ret[1]>0 && ret[2]>0)) {
        printf("%f, %f, %f\n", base_diffuse[0], base_diffuse[1], base_diffuse[2]);
        printf("%f, %f, %f\n", subsurface[0], subsurface[1], subsurface[2]);
        assert((ret[0]>0 && ret[1]>0 && ret[2]>0));
    }
    // return base_diffuse;
    // return subsurface;
    return (1-subsurface_scale)*base_diffuse+subsurface_scale*subsurface;
}

Real pdf_sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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
    return fmax(dot(frame.n, dir_out), Real(0)) / c_PI;
}

std::optional<BSDFSampleRecord> sample_bsdf_op::operator()(const DisneyDiffuse &bsdf) const {
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
    return BSDFSampleRecord{
        to_world(frame, sample_cos_hemisphere(rnd_param_uv)),
        Real(0) /* eta */, Real(1) /* roughness */};
}

TextureSpectrum get_texture_op::operator()(const DisneyDiffuse &bsdf) const {
    return bsdf.base_color;
}
