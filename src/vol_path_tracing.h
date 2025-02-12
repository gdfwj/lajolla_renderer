#pragma once

// The simplest volumetric renderer: 
// single absorption only homogeneous volume
// only handle directly visible light sources
Spectrum vol_path_tracing_1(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    // return make_zero_spectrum();
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};

    std::optional<PathVertex> vertex = intersect(scene, ray, ray_diff);

    if (vertex) {
        Medium medium_through = scene.media[vertex->exterior_medium_id];
        Real t = distance(vertex->position, ray.org);
        Spectrum transmittance = exp(-get_sigma_a(medium_through, vertex->position) * t);
        Spectrum Le = make_zero_spectrum();

        // Spectrum C = Vector3(3.2, 0.0, 0.0);
        // Spectrum E = C / get_sigma_a(medium_through, vertex->position);

        if (is_light(scene.shapes[vertex->shape_id])) {
            Le = emission(*vertex,-ray.dir,scene);
        }
        // return E + (Le-E) * transmittance;
        return Le * transmittance;
    }
    return make_zero_spectrum();
}

// The second simplest volumetric renderer: 
// single monochromatic homogeneous volume with single scattering,
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_2(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray camera_ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    std::optional<PathVertex> isect = intersect(scene, camera_ray, ray_diff);
    int camera_medium_id = scene.camera.medium_id;
    Spectrum sigma_a = get_sigma_a(scene.media[camera_medium_id], camera_ray.org);
    Spectrum sigma_s = get_sigma_s(scene.media[camera_medium_id], camera_ray.org);
    Spectrum sigma_t = sigma_a+sigma_s;

    Real u = next_pcg32_real<Real>(rng);
    Real t = -log(1-u) / sigma_t[0];
    Real t_hit = INFINITY;
    if (isect) {
        t_hit = distance(camera_ray.org, isect->position);
    }
    if (t < t_hit) {
        Spectrum trans_pdf = exp(-sigma_t * t) * sigma_t;
        Spectrum transmittance = exp(-sigma_t * t);
        // Vector3 p = camera_ray.org + t * camera_ray.dir;
        PhaseFunction phase_fuc = get_phase_function(scene.media[camera_medium_id]);
        int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
        Light light = scene.lights[light_id];
        Vector3 now_place = camera_ray.org + t * camera_ray.dir;
        PointAndNormal point_on_light = sample_point_on_light(light, now_place, Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)), next_pcg32_real<Real>(rng), scene);
        Ray shadow_ray = Ray{now_place, normalize(point_on_light.position - now_place), get_shadow_epsilon(scene), (1-get_shadow_epsilon(scene)) * distance(point_on_light.position,now_place)};
        if (occluded(scene, shadow_ray)) {
            return make_zero_spectrum();
        }

        Vector3 dir_in = -camera_ray.dir;
        Vector3 dir_out = normalize(point_on_light.position - now_place);
        Spectrum phase = eval(phase_fuc, dir_in, dir_out);

        Spectrum Le = emission(light, -dir_out, 0, point_on_light, scene);

        Spectrum attenuate = exp(-sigma_t * distance(now_place, point_on_light.position));
        Vector3 omega = (point_on_light.position - now_place) /distance(now_place, point_on_light.position);
        Real angle = abs(dot(omega, point_on_light.normal)) / distance_squared(now_place, point_on_light.position);
        Spectrum L_s1 = Le * attenuate * phase * angle;
        Real L_s1_pdf = pdf_point_on_light(light, point_on_light, now_place, scene) * light_pmf(scene, light_id);
        return (transmittance / trans_pdf) * sigma_s * (L_s1 / L_s1_pdf);
        // L_s1_estimate, L_s1_pdf = L_s1(p, sample_point_on_light(rng))
    } else {
        Spectrum trans_pdf = exp(-sigma_t * t_hit);
        Spectrum transmittance = exp(-sigma_t * t_hit);
        Spectrum Le = make_zero_spectrum();
        if (is_light(scene.shapes[isect->shape_id])) {
            Le = emission(*isect,-camera_ray.dir,scene);
        }
        return (transmittance / trans_pdf) * Le;
    }

    return make_zero_spectrum();
}

int update_medium(const PathVertex &isect, Ray ray, int current_medium) {
    int medium = current_medium;
    if (isect.interior_medium_id != isect.exterior_medium_id) {
        // At medium transition. Update medium.
        if (dot(ray.dir, isect.geometric_normal) > 0)
            medium = isect.exterior_medium_id;
        else
            medium = isect.interior_medium_id;
    }
    return medium;
}

// The third volumetric renderer (not so simple anymore): 
// multiple monochromatic homogeneous volumes with multiple scattering
// no need to handle surface lighting, only directly visible light source
Spectrum vol_path_tracing_3(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    int current_medium = scene.camera.medium_id;
    Spectrum current_path_throughput = Vector3(1.0, 1.0, 1.0);
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;
    while(1) {
        bool scatter = false;
        std::optional<PathVertex> isect = intersect(scene, ray, ray_diff);
        Spectrum transmittance = Vector3(1.0, 1.0, 1.0);
        Spectrum trans_pdf = Vector3(1.0, 1.0, 1.0);
        if (current_medium>-1) {
            Spectrum sigma_a = get_sigma_a(scene.media[current_medium], ray.org);
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium], ray.org);
            Spectrum sigma_t = sigma_a+sigma_s;
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1-u) / sigma_t[0];
            Real t_hit = INFINITY;
            if (isect) {
                t_hit = distance(ray.org, isect->position);
            }
            if (t < t_hit) {
                scatter = true;
                trans_pdf = exp(-sigma_t * t) * sigma_t;
                transmittance = exp(-sigma_t * t);
                ray.org += ray.dir * t;
            } else {
                trans_pdf = exp(-sigma_t * t_hit);
                transmittance = exp(-sigma_t * t_hit);
                ray.org = isect->position + get_intersection_epsilon(scene) * ray.dir;
            }
        }
        current_path_throughput *= (transmittance / trans_pdf);
        if (!scatter) {
            Spectrum Le = make_zero_spectrum();
            if (is_light(scene.shapes[isect->shape_id])) {
                Le = emission(*isect,-ray.dir,scene);
            }
            radiance += current_path_throughput * Le;
        }
        if (bounces == scene.options.max_depth- 1 && scene.options.max_depth !=-1) {
            break;
        }
        if (!scatter && isect) {
            if (isect->material_id ==-1) {
                // index-matching
                current_medium = update_medium(*isect, ray, current_medium);
                bounces += 1;
                continue;
            }
        }
        if (scatter) {
            // next_dir = sample_phase_function(-ray.dir, rng)
            PhaseFunction phase_fuc = get_phase_function(scene.media[current_medium]);
            Vector3 next_dir = sample_phase_function(phase_fuc, -ray.dir, Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng))).value();
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium], ray.org);
            Spectrum dir_in = normalize(-ray.dir);
            Spectrum dir_out = normalize(next_dir);
            current_path_throughput *= (eval(phase_fuc, dir_in, dir_out) / pdf_sample_phase(phase_fuc, dir_in, dir_out)) * sigma_s;
            ray.dir = next_dir;
        } else {
            break;
        }
        Real rr_prob = 1;
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(current_path_throughput[0], 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            } else{
                current_path_throughput /= rr_prob;
            }
        }
        bounces += 1;
    }
    return radiance;
}

Spectrum next_event_estimation(Vector3 p, int current_medium, const Scene &scene, pcg32_state &rng, int bounces, Vector3 dir_in) {
    int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
    Light light = scene.lights[light_id];
    PointAndNormal p_prime = sample_point_on_light(light, p, Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)), next_pcg32_real<Real>(rng), scene);
    Vector3 dir_light = normalize(p_prime.position - p);
    // # Compute transmittance to light. Skip through index-matching shapes.
    Real T_light = 1.0;
    int shadow_medium = current_medium;
    int shadow_bounces = 0;
    Real p_trans_dir = 1.0; // # for multiple importance sampling
    Spectrum p_cache = p;
    while (1) {
        Ray shadow_ray = Ray{p, normalize(p_prime.position - p), get_shadow_epsilon(scene), (1-get_shadow_epsilon(scene)) * distance(p_prime.position,p)};
        RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
        std::optional<PathVertex> isect = intersect(scene, shadow_ray, ray_diff);
        Real next_t = distance(p, p_prime.position);
        if (isect) {
            next_t = distance(p, isect->position);
        }
        if (shadow_medium!=-1) {
            Spectrum sigma_a = get_sigma_a(scene.media[shadow_medium], shadow_ray.org);
            Spectrum sigma_s = get_sigma_s(scene.media[shadow_medium], shadow_ray.org);
            Spectrum sigma_t = sigma_a+sigma_s;
            T_light *= exp(-sigma_t[0] * next_t);
            p_trans_dir *= exp(-sigma_t[0] * next_t);
        }
        if (!isect) {
            break;
        } else {
            if (isect->material_id >=0) {
                return make_zero_spectrum();
            }
            shadow_bounces += 1;
            if(scene.options.max_depth !=-1 && bounces + shadow_bounces + 1 >= scene.options.max_depth) {
                // # Reach the max no. of vertices
                return make_zero_spectrum();
            }
            shadow_medium = update_medium(*isect, shadow_ray, shadow_medium);
            p = p + next_t * dir_light;
        }
    }
    if (T_light>0) {
        Vector3 dir_out = dir_light;
        PhaseFunction phase_fuc = get_phase_function(scene.media[current_medium]);
        Spectrum rho = eval(phase_fuc, dir_in, dir_out);
        Spectrum L = emission(light, -dir_out, 0, p_prime, scene);
        Real G = abs(dot(dir_out, p_prime.normal)) / distance_squared(p_cache, p_prime.position);
        Real pdf_nee = pdf_point_on_light(light, p_prime, p_cache, scene) * light_pmf(scene, light_id);
        Spectrum contrib = T_light * G * rho * L / pdf_nee;
        Real pdf_phase = pdf_sample_phase(phase_fuc, dir_in, dir_out) * G * p_trans_dir;
        Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);
        return w * contrib;
    }
    return make_zero_spectrum();
}

// The fourth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// still no surface lighting
Spectrum vol_path_tracing_4(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    int current_medium = scene.camera.medium_id;
    Spectrum current_path_throughput = Vector3(1.0, 1.0, 1.0);
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;
    Real dir_pdf = 0;
    Spectrum nee_p_cache = make_zero_spectrum();
    Real multi_trans_pdf = 1;
    bool never_scatter = true;
    while(1) {
        bool scatter = false;
        std::optional<PathVertex> isect = intersect(scene, ray, ray_diff);
        Spectrum transmittance = Vector3(1.0, 1.0, 1.0);
        Spectrum trans_pdf = Vector3(1.0, 1.0, 1.0);
        if (current_medium>-1) {
            // printf("in current_medium>-1\n");
            Spectrum sigma_a = get_sigma_a(scene.media[current_medium], ray.org);
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium], ray.org);
            Spectrum sigma_t = sigma_a+sigma_s;
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1-u) / sigma_t[0];
            Real t_hit = INFINITY;
            if (isect) {
                t_hit = distance(ray.org, isect->position);
            }
            if (t < t_hit) {
                scatter = true;
                trans_pdf = exp(-sigma_t * t) * sigma_t;
                transmittance = exp(-sigma_t * t);
            } else {
                trans_pdf = exp(-sigma_t * t_hit);
                transmittance = exp(-sigma_t * t_hit);
            }
            ray.org += ray.dir * t;
            // printf("out current_medium>-1\n");
        }
        multi_trans_pdf *= trans_pdf[0];
        current_path_throughput *= (transmittance / trans_pdf);
        if (!scatter && isect) {
            ray.org = isect->position + get_intersection_epsilon(scene) * ray.dir;
            if (never_scatter) {
                // printf("in never scatter\n");
                Spectrum Le = make_zero_spectrum();
                if (is_light(scene.shapes[isect->shape_id])) {
                    Le = emission(*isect,-ray.dir,scene);
                }
                radiance += current_path_throughput * Le;
                // printf("out never scatter\n");
            } else {
                // printf("in never scatter else\n");
                if (is_light(scene.shapes[isect->shape_id])) {
                    int light_id = get_area_light_id(scene.shapes[isect->shape_id]);
                    Light light = scene.lights[light_id];
                    PointAndNormal p_prime = {isect->position, isect->geometric_normal};

                    Real pdf_nee = pdf_point_on_light(light, p_prime, nee_p_cache, scene);
                    Vector3 dir_out = normalize(p_prime.position - nee_p_cache);
                    Real G = abs(dot(dir_out, p_prime.normal)) / distance_squared(nee_p_cache, p_prime.position);
                    Real dir_pdf_ = dir_pdf * multi_trans_pdf * G;
                    Real w = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);
                    Spectrum Le = emission(*isect,-ray.dir,scene);
                    radiance += current_path_throughput * Le * w;
                    // printf("out never scatter else\n");
                }
            }
        }
        if (bounces == scene.options.max_depth- 1 && scene.options.max_depth !=-1) {
            break;
        }
        if (!scatter && isect) {
            // printf("in !scatter && isect\n");
            if (isect->material_id ==-1) {
                // index-matching
                current_medium = update_medium(*isect, ray, current_medium);
                bounces += 1;
                // printf("out !scatter && isect\n");
                continue;
            }
            // printf("out !scatter && isect\n");
        }
        if (scatter) {
            multi_trans_pdf = 1;
            // printf("in scatter\n");
            PhaseFunction phase_fuc = get_phase_function(scene.media[current_medium]);
            Vector3 next_dir = sample_phase_function(phase_fuc, -ray.dir, Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng))).value();
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium], ray.org);
            Spectrum dir_in = normalize(-ray.dir);
            Spectrum dir_out = normalize(next_dir);
            dir_pdf = pdf_sample_phase(phase_fuc, dir_in, dir_out);

            never_scatter = false;
            nee_p_cache = ray.org;
            Spectrum nee = next_event_estimation(ray.org, current_medium, scene, rng, bounces, -ray.dir);
            radiance += current_path_throughput * nee * sigma_s;
            current_path_throughput *= eval(phase_fuc, dir_in, dir_out) / pdf_sample_phase(phase_fuc, dir_in, dir_out) * sigma_s;
            // current_path_throughput *= (eval(phase_fuc, dir_in, dir_out) / pdf_sample_phase(phase_fuc, dir_in, dir_out)) * sigma_s;
            ray.dir = next_dir;
            // printf("out scatter\n");
        } else {
            break;
        }
        Real rr_prob = 1;
        // printf("in rr_prob\n");
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(current_path_throughput[0], 0.95);
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            } else{
                current_path_throughput /= rr_prob;
            }
        }
        bounces += 1;
        // printf("out rr_prob\n");
    }
    return radiance;
}

Spectrum next_event_estimation_on_surface(Vector3 p, int current_medium, const Scene &scene, pcg32_state &rng, int bounces, Vector3 dir_in, std::optional<PathVertex> &isect_out) {
    const Material &mat = scene.materials[isect_out->material_id];
    int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
    Light light = scene.lights[light_id];
    PointAndNormal p_prime = sample_point_on_light(light, p, Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)), next_pcg32_real<Real>(rng), scene);
    Vector3 dir_light = normalize(p_prime.position - p);
    // # Compute transmittance to light. Skip through index-matching shapes.
    Real T_light = 1.0;
    int shadow_medium = current_medium;
    int shadow_bounces = 0;
    Real p_trans_dir = 1.0; // # for multiple importance sampling
    Spectrum p_cache = p;
    while (1) {
        Ray shadow_ray = Ray{p, normalize(p_prime.position - p), get_shadow_epsilon(scene), (1-get_shadow_epsilon(scene)) * distance(p_prime.position,p)};
        RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
        std::optional<PathVertex> isect = intersect(scene, shadow_ray, ray_diff);
        Real next_t = distance(p, p_prime.position);
        if (isect) {
            next_t = distance(p, isect->position);
        }
        if (shadow_medium!=-1) {
            Spectrum sigma_a = get_sigma_a(scene.media[shadow_medium], shadow_ray.org);
            Spectrum sigma_s = get_sigma_s(scene.media[shadow_medium], shadow_ray.org);
            Spectrum sigma_t = sigma_a+sigma_s;
            T_light *= exp(-sigma_t[0] * next_t);
            p_trans_dir *= exp(-sigma_t[0] * next_t);
        }
        if (!isect) {
            break;
        } else {
            if (isect->material_id >=0) {
                return make_zero_spectrum();
            }
            shadow_bounces += 1;
            if(scene.options.max_depth !=-1 && bounces + shadow_bounces + 1 >= scene.options.max_depth) {
                // # Reach the max no. of vertices
                return make_zero_spectrum();
            }
            shadow_medium = update_medium(*isect, shadow_ray, shadow_medium);
            p = p + next_t * dir_light + get_shadow_epsilon(scene) * dir_light;
        }
    }
    if (T_light>0) {
        Vector3 dir_out = dir_light;
        // PhaseFunction phase_fuc = get_phase_function(scene.media[shadow_medium]);
        // Spectrum rho = eval(phase_fuc, dir_in, dir_out);
        Spectrum f = eval(mat, dir_in, dir_out, *isect_out, scene.texture_pool);
        Spectrum L = emission(light, -dir_out, 0, p_prime, scene);
        Real G = abs(dot(dir_out, p_prime.normal)) / distance_squared(p_cache, p_prime.position);
        Real pdf_nee = pdf_point_on_light(light, p_prime, p_cache, scene) * light_pmf(scene, light_id);
        Spectrum contrib = T_light * G * f * L / pdf_nee;
        Real pdf_bsdf = pdf_sample_bsdf(mat, dir_in, dir_light, *isect_out, scene.texture_pool) * p_trans_dir * G;
        // Real pdf_phase = pdf_sample_phase(phase_fuc, dir_in, dir_out) * G * p_trans_dir;
        Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_bsdf * pdf_bsdf);
        return w * contrib;
    }
    return make_zero_spectrum();
}

// The fifth volumetric renderer: 
// multiple monochromatic homogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing_5(const Scene &scene,
                            int x, int y, /* pixel coordinates */
                            pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    int current_medium = scene.camera.medium_id;
    Spectrum current_path_throughput = Vector3(1.0, 1.0, 1.0);
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;
    Real dir_pdf = 0;
    Spectrum nee_p_cache = make_zero_spectrum();
    Real multi_trans_pdf = 1;
    bool never_scatter = true;
    Real eta_scale = 1.0;
    while(1) {
        bool scatter = false;
        std::optional<PathVertex> isect = intersect(scene, ray, ray_diff);
        Spectrum transmittance = Vector3(1.0, 1.0, 1.0);
        Spectrum trans_pdf = Vector3(1.0, 1.0, 1.0);
        if (current_medium>-1) {
            // printf("in current_medium>-1\n");
            Spectrum sigma_a = get_sigma_a(scene.media[current_medium], ray.org);
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium], ray.org);
            Spectrum sigma_t = sigma_a+sigma_s;
            Real u = next_pcg32_real<Real>(rng);
            Real t = -log(1-u) / sigma_t[0];
            Real t_hit = INFINITY;
            if (isect) {
                t_hit = distance(ray.org, isect->position);
            }
            if (t < t_hit) {
                scatter = true;
                trans_pdf = exp(-sigma_t * t) * sigma_t;
                transmittance = exp(-sigma_t * t);
            } else {
                trans_pdf = exp(-sigma_t * t_hit);
                transmittance = exp(-sigma_t * t_hit);
            }
            ray.org += ray.dir * t + get_intersection_epsilon(scene) * ray.dir;
            // printf("out current_medium>-1\n");
        }
        multi_trans_pdf *= trans_pdf[0];
        current_path_throughput *= (transmittance / trans_pdf);
        if (!scatter && isect) {
            ray.org = isect->position;
            if (never_scatter) {
                // printf("in never scatter\n");
                Spectrum Le = make_zero_spectrum();
                if (is_light(scene.shapes[isect->shape_id])) {
                    Le = emission(*isect,-ray.dir,scene);
                }
                radiance += current_path_throughput * Le;
                // printf("out never scatter\n");
            } else {
                // printf("in never scatter else\n");
                if (is_light(scene.shapes[isect->shape_id])) {
                    int light_id = get_area_light_id(scene.shapes[isect->shape_id]);
                    Light light = scene.lights[light_id];
                    PointAndNormal p_prime = {isect->position, isect->geometric_normal};

                    Real pdf_nee = pdf_point_on_light(light, p_prime, nee_p_cache, scene);
                    Vector3 dir_out = normalize(p_prime.position - nee_p_cache);
                    Real G = abs(dot(dir_out, p_prime.normal)) / distance_squared(nee_p_cache, p_prime.position);
                    Real dir_pdf_ = dir_pdf * multi_trans_pdf * G;
                    Real w = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);
                    Spectrum Le = emission(*isect,-ray.dir,scene);
                    radiance += current_path_throughput * Le * w;
                    // printf("out never scatter else\n");
                }
            }
        }
        if (bounces == scene.options.max_depth- 1 && scene.options.max_depth !=-1) {
            break;
        }
        if (!scatter && isect) {
            // printf("in !scatter && isect\n");
            if (isect->material_id ==-1) {
                // index-matching
                current_medium = update_medium(*isect, ray, current_medium);
                bounces += 1;
                // printf("out !scatter && isect\n");
                continue;
            } else {
                // bsdf
                nee_p_cache = ray.org;
                multi_trans_pdf = 1;
                never_scatter = false;
                ray.tnear = get_shadow_epsilon(scene);

                const Material &mat = scene.materials[isect->material_id];
                Spectrum nee_bsdf = next_event_estimation_on_surface(ray.org, current_medium, scene, rng, bounces, -ray.dir, isect);
                radiance += current_path_throughput * nee_bsdf;
                std::optional<BSDFSampleRecord> bsdf_sample_ =
                    sample_bsdf(mat,
                                -ray.dir,
                                *isect,
                                scene.texture_pool,
                                Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)),
                                next_pcg32_real<Real>(rng));
                if (!bsdf_sample_) {
                    // BSDF sampling failed. Abort the loop.
                    break;
                }
                const BSDFSampleRecord &bsdf_sample = *bsdf_sample_;
                Vector3 dir_out = bsdf_sample.dir_out;
                Real pdf_bsdf = pdf_sample_bsdf(mat, -ray.dir, dir_out, *isect, scene.texture_pool);
                current_path_throughput *= eval(mat, -ray.dir, dir_out, *isect, scene.texture_pool) / pdf_bsdf;
                ray.dir = dir_out;
                current_medium = update_medium(*isect, ray, current_medium);
                bounces += 1;

                dir_pdf = pdf_bsdf;
                if (bsdf_sample.eta == 0) {
                    // ray_diff.spread = reflect(ray_diff, isect->mean_curvature, bsdf_sample.roughness);
                } else {
                    // ray_diff.spread = refract(ray_diff, isect->mean_curvature, bsdf_sample.eta, bsdf_sample.roughness);
                    eta_scale /= (bsdf_sample.eta * bsdf_sample.eta);
                }
            }
            // printf("out !scatter && isect\n");
        }
        if (scatter) {
            multi_trans_pdf = 1;
            // printf("in scatter\n");
            PhaseFunction phase_fuc = get_phase_function(scene.media[current_medium]);
            Vector3 next_dir = sample_phase_function(phase_fuc, -ray.dir, Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng))).value();
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium], ray.org);
            Spectrum dir_in = normalize(-ray.dir);
            Spectrum dir_out = normalize(next_dir);
            dir_pdf = pdf_sample_phase(phase_fuc, dir_in, dir_out);

            never_scatter = false;
            nee_p_cache = ray.org;
            Spectrum nee = next_event_estimation(ray.org, current_medium, scene, rng, bounces, -ray.dir);
            radiance += current_path_throughput * nee * sigma_s;
            current_path_throughput *= eval(phase_fuc, dir_in, dir_out) / dir_pdf * sigma_s;
            // current_path_throughput *= (eval(phase_fuc, dir_in, dir_out) / pdf_sample_phase(phase_fuc, dir_in, dir_out)) * sigma_s;
            ray.dir = next_dir;
            // printf("out scatter\n");
        }
        Real rr_prob = 1;
        // printf("in rr_prob\n");
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(max((1 / eta_scale) * current_path_throughput), Real(0.95));
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                break;
            } else{
                current_path_throughput /= rr_prob;
            }
        }
        bounces += 1;
        // printf("out rr_prob\n");
    }
    return radiance;
}

Spectrum next_event_estimation_heterogeneous(Vector3 p, int current_medium, const Scene &scene, pcg32_state &rng, int bounces, Vector3 dir_in, const std::optional<PathVertex> &isect_out) {
    int light_id = sample_light(scene, next_pcg32_real<Real>(rng));
    Light light = scene.lights[light_id];
    PointAndNormal p_prime = sample_point_on_light(light, p, Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)), next_pcg32_real<Real>(rng), scene);
    Vector3 dir_light = normalize(p_prime.position - p);
    // # Compute transmittance to light. Skip through index-matching shapes.
    int shadow_medium = current_medium;
    int shadow_bounces = 0;
    Spectrum T_light = make_const_spectrum(1);
    Spectrum p_trans_nee = make_const_spectrum(1);
    Spectrum p_trans_dir = make_const_spectrum(1); // # for multiple importance sampling
    // auto outer_start = std::chrono::high_resolution_clock::now();
    bool debug = false;

    // must cached p for final calculation
    Vector3 p_cached = p;
    while(1) {
        if (shadow_bounces>10) {
            std::cout<<shadow_bounces<<std::endl;
        }
        // auto outer_now = std::chrono::high_resolution_clock::now();
        // std::cout<<"in loop 3"<<shadow_bounces<<std::endl;
        Ray shadow_ray = Ray{p, normalize(p_prime.position - p), get_shadow_epsilon(scene), (1-get_shadow_epsilon(scene)) * distance(p_prime.position, p)};
        
        RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
        std::optional<PathVertex> isect = intersect(scene, shadow_ray, ray_diff);
        Real next_t = distance(p, p_prime.position);
        // if (std::chrono::duration_cast<std::chrono::microseconds>(outer_now-outer_start).count()>100000) {
        //     std::cout<<"-----------------"<<std::endl;
        //     std::cout<<"loop 3 "<<shadow_bounces<<std::endl;
        //     std::cout<<p<<std::endl;
        //     std::cout<<p_prime.position<<std::endl;
        //     std::cout<<shadow_ray.dir<<std::endl;
        //     std::cout<<next_t<<std::endl;
        //     debug = true;
        // }
        if (isect) {  // have to put the integral out the isect, because we have to count the last intersection to light
            if (isect->material_id >=0) {
                return make_zero_spectrum();
            }
            next_t = distance(p, isect->position);
        } 
        if (shadow_medium!=-1) {
            Real u = next_pcg32_real<Real>(rng);
            int channel = std::clamp(int(u * 3), 0, 2);
            int iteration = 0;
            Spectrum majorant = get_majorant(scene.media[shadow_medium], shadow_ray);
            Real accum_t = 0.0;
            // auto outer_start = std::chrono::high_resolution_clock::now();
            while(1) {
                // auto outer_now = std::chrono::high_resolution_clock::now();
                // if (std::chrono::duration_cast<std::chrono::microseconds>(outer_now-outer_start).count()>100000) {
                //     std::cout<<"loop 4"<<iteration<<std::endl;
                // }
                // std::cout<<"in loop 4"<<iteration<<std::endl;
                if (majorant[channel]<=0) break;
                if (iteration > scene.options.max_null_collisions) break;
                Real t =  -log(1- next_pcg32_real<Real>(rng)) / majorant[channel];
                Real dt = next_t- accum_t;
                accum_t = min(accum_t + t, next_t);
                if (t < dt) {
                    p = p + t * dir_light;
                    Spectrum sigma_a = get_sigma_a(scene.media[shadow_medium], p);
                    Spectrum sigma_s = get_sigma_s(scene.media[shadow_medium], p);
                    Spectrum sigma_t = sigma_a+sigma_s;
                    Spectrum real_prob = sigma_t / majorant;
                    Spectrum sigma_n = majorant - sigma_t;
                    T_light *= exp(-majorant * t) * sigma_n / max(majorant);
                    p_trans_nee *= exp(-majorant * t) * majorant / max(majorant);
                    p_trans_dir *= exp(-majorant * t) * majorant * (1- real_prob) / max(majorant);
                    if (max(T_light) <= 0) {
                        break;
                    } // optimization for places where sigma_n = 0
                } else { // hit
                    T_light *= exp(-majorant * dt);
                    p_trans_nee *= exp(-majorant * dt);
                    p_trans_dir *= exp(-majorant * dt);
                    p = p + dt * dir_light;
                    break;
                }
                iteration += 1;
            }
        }
        if (!isect) {
            break;
        } else { // should put behind the break
            shadow_bounces += 1;
            if(scene.options.max_depth !=-1 && bounces + shadow_bounces + 1 >= scene.options.max_depth) {
                // # Reach the max no. of vertices
                return make_zero_spectrum();
            }
            shadow_medium = update_medium(*isect, shadow_ray, shadow_medium);
            p = isect->position;
        }
        
    }
    if (max(T_light)>0) {
        if(isect_out) {  // if we hit a surface
            const Material &mat = scene.materials[isect_out->material_id];
            Vector3 dir_out = dir_light;
            Spectrum f = eval(mat, dir_in, dir_out, *isect_out, scene.texture_pool);
            Spectrum L = emission(light, -dir_out, 0, p_prime, scene);
            Real G = max(-dot(dir_out, p_prime.normal), Real(0)) / distance_squared(p_cached, p_prime.position);
            Real pdf_nee = pdf_point_on_light(light, p_prime, p_cached, scene) * light_pmf(scene, light_id) * avg(p_trans_nee);
            Spectrum contrib = T_light * G * f * L / pdf_nee;
            Real pdf_bsdf = pdf_sample_bsdf(mat, dir_in, dir_light, *isect_out, scene.texture_pool) * avg(p_trans_dir) * G;
            // Real pdf_phase = pdf_sample_phase(phase_fuc, dir_in, dir_out) * G * p_trans_dir;
            Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_bsdf * pdf_bsdf);
            return w * contrib;
        } else {// if scattering
            Vector3 dir_out = dir_light;
            // use the medium of scattering place (outer) to get phase function
            PhaseFunction phase_fuc = get_phase_function(scene.media[current_medium]);
            Spectrum rho = eval(phase_fuc, dir_in, dir_out);
            Spectrum L = emission(light, -dir_out, 0, p_prime, scene);
            Real G = max(-dot(dir_out, p_prime.normal), Real(0)) / distance_squared(p_cached, p_prime.position);
            Real pdf_nee = pdf_point_on_light(light, p_prime, p_cached, scene) * light_pmf(scene, light_id) * avg(p_trans_nee);
            Spectrum contrib = T_light * G * rho * L / pdf_nee;
            Real pdf_phase = pdf_sample_phase(phase_fuc, dir_in, dir_out) * G * avg(p_trans_dir);
            Real w = (pdf_nee * pdf_nee) / (pdf_nee * pdf_nee + pdf_phase * pdf_phase);
            return w * contrib;
        }
    }
    return make_zero_spectrum();
}


// The final volumetric renderer: 
// multiple chromatic heterogeneous volumes with multiple scattering
// with MIS between next event estimation and phase function sampling
// with surface lighting
Spectrum vol_path_tracing(const Scene &scene,
                          int x, int y, /* pixel coordinates */
                          pcg32_state &rng) {
    // Homework 2: implememt this!
    int w = scene.camera.width, h = scene.camera.height;
    Vector2 screen_pos((x + next_pcg32_real<Real>(rng)) / w,
                       (y + next_pcg32_real<Real>(rng)) / h);
    Ray ray = sample_primary(scene.camera, screen_pos);
    RayDifferential ray_diff = RayDifferential{Real(0), Real(0)};
    int current_medium = scene.camera.medium_id;
    Spectrum current_path_throughput = Vector3(1.0, 1.0, 1.0);
    Spectrum radiance = make_zero_spectrum();
    int bounces = 0;
    Real dir_pdf = 0;
    Spectrum nee_p_cache = make_zero_spectrum();
    Spectrum multi_trans_pdf = make_const_spectrum(1);
    Spectrum multi_nee_pdf = make_const_spectrum(1);
    bool never_scatter = true;
    Real eta_scale = 1.0;
    // auto outer_start = std::chrono::high_resolution_clock::now();
    while(1) {
        // auto outer_now = std::chrono::high_resolution_clock::now();
        // if (std::chrono::duration_cast<std::chrono::microseconds>(outer_now-outer_start).count()>100000) {
        //     std::cout<<"loop 1"<<bounces<<std::endl;
        // }
        // std::cout<<"in loop 1 "<<bounces<<std::endl;
        bool scatter = false;
        ray.tnear = get_intersection_epsilon(scene);
        std::optional<PathVertex> isect = intersect(scene, ray, ray_diff);
        Spectrum transmittance = Vector3(1.0, 1.0, 1.0);
        Spectrum trans_nee_pdf = make_const_spectrum(1);
        Spectrum trans_dir_pdf = make_const_spectrum(1);
        if (current_medium>-1) {
            //edit 
            Spectrum majorant = get_majorant(scene.media[current_medium], ray);
            Real u = next_pcg32_real<Real>(rng);
            Real channel = std::clamp(int(u * 3), 0, 2);
            Real accum_t = 0;
            int iteration = 0;
            Real t_hit = INFINITY;
            if (isect) {
                t_hit = distance(ray.org, isect->position);
            }
            // auto outer_start = std::chrono::high_resolution_clock::now();
            while(1) {
                // auto outer_now = std::chrono::high_resolution_clock::now();
                // if (std::chrono::duration_cast<std::chrono::microseconds>(outer_now-outer_start).count()>100000) {
                //     std::cout<<"loop 2"<<iteration<<std::endl;
                // }
                // std::cout<<"in loop 2"<<iteration<<std::endl;
                if (majorant[channel]<=0) break;
                if (iteration > scene.options.max_null_collisions) break;
                Real t =  -log(1- next_pcg32_real<Real>(rng)) / majorant[channel];
                Real dt = t_hit- accum_t;
                accum_t = min(accum_t + t, t_hit);
                if (t < dt) {
                    ray.org = ray.org + t * ray.dir;
                    ray.tnear = get_intersection_epsilon(scene);
                    Spectrum sigma_a = get_sigma_a(scene.media[current_medium], ray.org);
                    Spectrum sigma_s = get_sigma_s(scene.media[current_medium], ray.org);
                    Spectrum sigma_t = sigma_a+sigma_s;
                    Spectrum real_prob = sigma_t / majorant;
                    Spectrum sigma_n = majorant - sigma_t;
                    if (next_pcg32_real<Real>(rng) < real_prob[channel]) { // hit real particle
                        scatter = true;
                        never_scatter = false;
                        transmittance *= exp(-majorant * t) / max(majorant);
                        trans_dir_pdf *= exp(-majorant * t) * majorant * real_prob / max(majorant);
                        break;
                    } else { // hit fake particle
                        transmittance *= exp(-majorant * t) * sigma_n / max(majorant);
                        trans_dir_pdf *= exp(-majorant * t) * majorant * (1- real_prob) / max(majorant);
                        trans_nee_pdf *= exp(-majorant * t) * majorant / max(majorant);
                    }
                } else {
                    ray.org = isect->position;
                    ray.tnear = get_intersection_epsilon(scene);
                    transmittance *= exp(-majorant * dt);
                    trans_dir_pdf *= exp(-majorant * dt);
                    trans_nee_pdf *= exp(-majorant * dt);
                    break;
                }
                iteration += 1;
            }
        }
        multi_trans_pdf *= trans_dir_pdf;
        multi_nee_pdf *= trans_nee_pdf;
        current_path_throughput *= (transmittance / avg(trans_dir_pdf));
        if (!scatter && isect) {
            ray.org = isect->position;
            ray.tnear = get_intersection_epsilon(scene);
            if (never_scatter) {
                // printf("in never scatter\n");
                Spectrum Le = make_zero_spectrum();
                if (is_light(scene.shapes[isect->shape_id])) {
                    Le = emission(*isect,-ray.dir,scene);
                }
                radiance += current_path_throughput * Le;
                // printf("out never scatter\n");
            } else {
                // printf("in never scatter else\n");
                if (is_light(scene.shapes[isect->shape_id])) {
                    int light_id = get_area_light_id(scene.shapes[isect->shape_id]);
                    Light light = scene.lights[light_id];
                    PointAndNormal p_prime = {isect->position, isect->geometric_normal};

                    Real pdf_nee = pdf_point_on_light(light, p_prime, nee_p_cache, scene) * avg(multi_nee_pdf);
                    Vector3 dir_out = normalize(p_prime.position - nee_p_cache);
                    Real G = abs(dot(dir_out, p_prime.normal)) / distance_squared(nee_p_cache, p_prime.position);
                    Real dir_pdf_ = dir_pdf * avg(multi_trans_pdf) * G;
                    Real w = (dir_pdf_ * dir_pdf_) / (dir_pdf_ * dir_pdf_ + pdf_nee * pdf_nee);
                    Spectrum Le = emission(*isect,-ray.dir,scene);
                    radiance += current_path_throughput * Le * w;
                    // printf("out never scatter else\n");
                }
            }
        }
        if (bounces == scene.options.max_depth- 1 && scene.options.max_depth !=-1) {
            break;
        }
        if (!scatter && isect) {
            // printf("in !scatter && isect\n");
            if (isect->material_id ==-1) {
                // index-matching
                current_medium = update_medium(*isect, ray, current_medium);
                bounces += 1;
                // printf("out !scatter && isect\n");
                continue;
            } else {
                // bsdf
                nee_p_cache = ray.org;
                multi_trans_pdf = make_const_spectrum(1);
                multi_nee_pdf = make_const_spectrum(1);
                ray.tnear = get_shadow_epsilon(scene);

                const Material &mat = scene.materials[isect->material_id];
                // change nee
                // Spectrum nee_bsdf = next_event_estimation_on_surface(ray.org, current_medium, scene, rng, bounces, -ray.dir, isect);
                // std::cout<<"bsdf"<<std::endl;
                Spectrum nee_bsdf = next_event_estimation_heterogeneous(ray.org, current_medium, scene, rng, bounces, -ray.dir, isect);

                radiance += current_path_throughput * nee_bsdf;
                std::optional<BSDFSampleRecord> bsdf_sample_ =
                    sample_bsdf(mat,
                                -ray.dir,
                                *isect,
                                scene.texture_pool,
                                Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng)),
                                next_pcg32_real<Real>(rng));
                if (!bsdf_sample_) {
                    // BSDF sampling failed. Abort the loop.
                    break;
                }
                const BSDFSampleRecord &bsdf_sample = *bsdf_sample_;
                Vector3 dir_out = bsdf_sample.dir_out;
                Real pdf_bsdf = pdf_sample_bsdf(mat, -ray.dir, dir_out, *isect, scene.texture_pool);
                current_path_throughput *= eval(mat, -ray.dir, dir_out, *isect, scene.texture_pool) / pdf_bsdf;
                ray.dir = dir_out;
                ray.tnear = get_intersection_epsilon(scene);
                current_medium = update_medium(*isect, ray, current_medium);
                bounces += 1;

                dir_pdf = pdf_bsdf;
                if (bsdf_sample.eta == 0) {
                    // ray_diff.spread = reflect(ray_diff, isect->mean_curvature, bsdf_sample.roughness);
                } else {
                    // ray_diff.spread = refract(ray_diff, isect->mean_curvature, bsdf_sample.eta, bsdf_sample.roughness);
                    eta_scale /= (bsdf_sample.eta * bsdf_sample.eta);
                }
            }
            // printf("out !scatter && isect\n");
        }
        if (scatter) {
            multi_trans_pdf = make_const_spectrum(1);
            multi_nee_pdf = make_const_spectrum(1);
            // printf("in scatter\n");
            PhaseFunction phase_fuc = get_phase_function(scene.media[current_medium]);
            Vector3 next_dir = sample_phase_function(phase_fuc, -ray.dir, Vector2(next_pcg32_real<Real>(rng), next_pcg32_real<Real>(rng))).value();
            Spectrum sigma_s = get_sigma_s(scene.media[current_medium], ray.org);
            Spectrum dir_in = normalize(-ray.dir);
            Spectrum dir_out = normalize(next_dir);
            dir_pdf = pdf_sample_phase(phase_fuc, dir_in, dir_out);

            nee_p_cache = ray.org;
            // change nee
            // Spectrum nee = next_event_estimation(ray.org, current_medium, scene, rng, bounces, -ray.dir);
            // std::cout<<"scatter"<<std::endl;
            // Spectrum nee = next_event_estimation(ray.org, current_medium, scene, rng, bounces, -ray.dir);

            Spectrum nee = next_event_estimation_heterogeneous(ray.org, current_medium, scene, rng, bounces, -ray.dir, std::nullopt);

            radiance += current_path_throughput * nee * sigma_s;
            current_path_throughput *= eval(phase_fuc, dir_in, dir_out) / dir_pdf * sigma_s;
            // current_path_throughput *= (eval(phase_fuc, dir_in, dir_out) / pdf_sample_phase(phase_fuc, dir_in, dir_out)) * sigma_s;
            ray.dir = next_dir;
            // printf("out scatter\n");
        }
        Real rr_prob = 1;
        // printf("in rr_prob\n");
        if (bounces >= scene.options.rr_depth) {
            rr_prob = min(max((1 / eta_scale) * current_path_throughput), Real(0.95));
            // std::cout<<rr_prob<<std::endl;
            if (next_pcg32_real<Real>(rng) > rr_prob) {
                // std::cout<<"break good";
                break;
            } else{
                current_path_throughput /= rr_prob;
            }
        }
        bounces += 1;
        never_scatter = false;
        // printf("out rr_prob\n");
    }
    return radiance;
}
