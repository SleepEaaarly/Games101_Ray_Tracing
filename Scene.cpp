//
// Created by Göksu Güvendiren on 2019-05-14.
//

#include "Scene.hpp"


void Scene::buildBVH() {
    printf(" - Generating BVH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::NAIVE);
}

void Scene::buildSAH() {
    printf(" - Generating SAH...\n\n");
    this->bvh = new BVHAccel(objects, 1, BVHAccel::SplitMethod::SAH);
}

Intersection Scene::intersect(const Ray &ray) const
{
    return this->bvh->Intersect(ray);
}

void Scene::sampleLight(Intersection &pos, float &pdf) const
{
    float emit_area_sum = 0;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
        }
    }
    float p = get_random_float() * emit_area_sum;  // 随机均匀采样的系数
    emit_area_sum = 0; // 置0是为了下一步uniform采样
    for (uint32_t k = 0; k < objects.size(); ++k) {
        if (objects[k]->hasEmit()){
            emit_area_sum += objects[k]->getArea();
            if (p <= emit_area_sum){
                objects[k]->Sample(pos, pdf);
                break; // 只采样一个点
            }
        }
    }
}

bool Scene::trace(
        const Ray &ray,
        const std::vector<Object*> &objects,
        float &tNear, uint32_t &index, Object **hitObject)
{
    *hitObject = nullptr;
    for (uint32_t k = 0; k < objects.size(); ++k) {
        float tNearK = kInfinity;
        uint32_t indexK;
        Vector2f uvK;
        if (objects[k]->intersect(ray, tNearK, indexK) && tNearK < tNear) {
            *hitObject = objects[k];
            tNear = tNearK;
            index = indexK;
        }
    }


    return (*hitObject != nullptr);
}

// Implementation of Path Tracing
Vector3f Scene::castRay(const Ray &ray, int depth) const
{
    // TO DO Implement Path Tracing Algorithm here
    Intersection inter = intersect(ray); // 入射光线求交
    if (!inter.happened) return this->backgroundColor; // 如果没有交点
    Vector3f l_dir(0, 0, 0), l_indir(0, 0, 0), l_emit(0, 0, 0);
    if (inter.m->hasEmission())         // 如果物体本身发光，要加上发光项
        l_emit = inter.m->getEmission();

    // Contribution from the light source.
    Intersection light_point;
    float pdf_light;
    sampleLight(light_point, pdf_light);
    Vector3f light_coords = light_point.coords, ws = (light_point.coords - inter.coords).normalized(), 
             light_normal = light_point.normal, emit = light_point.emit;
    Vector3f wo = ray.direction;

    Intersection block_check = intersect(Ray(inter.coords, ws));
    if ((light_coords - inter.coords).norm() - block_check.distance < EPSILON) {
        l_dir = emit * inter.m->eval(wo, ws, inter.normal) * dotProduct(ws, inter.normal) * dotProduct(-ws, light_normal) 
        / dotProduct(light_coords, inter.coords) / pdf_light; 
    }

    // Contribution from other reflectors.
    if (get_random_float() > RussianRoulette) {
        return l_dir + l_emit;
    }
    Vector3f wi = inter.m->sample(wo, inter.normal);
    Ray r(inter.coords, wi);
    Intersection non_emit_check = intersect(r);
    if (non_emit_check.happened && !non_emit_check.m->hasEmission()) {
        float pdf = inter.m->pdf(wo, wi, inter.normal);     // 0的处理
        if (pdf > EPSILON) {
            l_indir = castRay(r, depth+1) * inter.m->eval(wo, wi, inter.normal) * dotProduct(wi, inter.normal) / pdf / RussianRoulette;
        }
    }

    return l_dir + l_indir + l_emit;
}