#include <algorithm>
#include <cassert>
#include "BVH.hpp"

BVHAccel::BVHAccel(std::vector<Object*> p, int maxPrimsInNode,
                   SplitMethod splitMethod)
    : maxPrimsInNode(std::min(255, maxPrimsInNode)), splitMethod(splitMethod),
      primitives(std::move(p))
{
    time_t start, stop;
    time(&start);
    if (primitives.empty())
        return;

    if (splitMethod == SplitMethod::NAIVE) {
        root = recursiveBuild(primitives);
    } else {
        root = recursiveBuildSAH(primitives);
    }

    time(&stop);
    double diff = difftime(stop, start);
    int hrs = (int)diff / 3600;
    int mins = ((int)diff / 60) - (hrs * 60);
    int secs = (int)diff - (hrs * 3600) - (mins * 60);

    if (splitMethod == SplitMethod::NAIVE) {
        printf(
            "\rBVH Generation complete: \nTime Taken: %i hrs, %i mins, %i secs\n\n",
            hrs, mins, secs);
    } else {
        printf(
            "\rSAH Generation complete: \nTime Taken: %i hrs, %i mins, %i secs\n\n",
            hrs, mins, secs);
    }
}

BVHBuildNode* BVHAccel::recursiveBuild(std::vector<Object*> objects)
{
    BVHBuildNode* node = new BVHBuildNode();

    // Compute bounds of all primitives in BVH node
    Bounds3 bounds;
    for (int i = 0; i < objects.size(); ++i)
        bounds = Union(bounds, objects[i]->getBounds());
    if (objects.size() == 1) {
        // Create leaf _BVHBuildNode_
        node->bounds = objects[0]->getBounds();
        node->object = objects[0];
        node->left = nullptr;
        node->right = nullptr;
        node->area = objects[0]->getArea();
        return node;
    }
    else if (objects.size() == 2) {
        node->left = recursiveBuild(std::vector{objects[0]});
        node->right = recursiveBuild(std::vector{objects[1]});

        node->bounds = Union(node->left->bounds, node->right->bounds);
        node->area = node->left->area + node->right->area;
        return node;
    }
    else {
        Bounds3 centroidBounds;
        for (int i = 0; i < objects.size(); ++i)
            centroidBounds =
                Union(centroidBounds, objects[i]->getBounds().Centroid());
        int dim = centroidBounds.maxExtent();
        switch (dim) {
        case 0:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().x <
                       f2->getBounds().Centroid().x;
            });
            break;
        case 1:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().y <
                       f2->getBounds().Centroid().y;
            });
            break;
        case 2:
            std::sort(objects.begin(), objects.end(), [](auto f1, auto f2) {
                return f1->getBounds().Centroid().z <
                       f2->getBounds().Centroid().z;
            });
            break;
        }

        auto beginning = objects.begin();
        auto middling = objects.begin() + (objects.size() / 2);
        auto ending = objects.end();

        auto leftshapes = std::vector<Object*>(beginning, middling);
        auto rightshapes = std::vector<Object*>(middling, ending);

        assert(objects.size() == (leftshapes.size() + rightshapes.size()));

        node->left = recursiveBuild(leftshapes);
        node->right = recursiveBuild(rightshapes);

        node->bounds = Union(node->left->bounds, node->right->bounds);
        node->area = node->left->area + node->right->area;
    }

    return node;
}

BVHBuildNode* BVHAccel::recursiveBuildSAH(std::vector<Object*> objects)
{   
    // 物体数小于12时用BVH划分
    if (objects.size() < 12) {
        return recursiveBuild(objects);
    }
    BVHBuildNode* node = new BVHBuildNode();

    Bounds3 bounds;
    for (int i = 0; i < objects.size(); ++i)
        bounds = Union(bounds, objects[i]->getBounds());
    
    int dim, div_ind;
    double min_cost = std::numeric_limits<double>::max();
    float div[] = {1./6, 2./6, 3./6, 4./6, 5./6};

    for (int i = 0; i < 3; i++) {
        std::sort(objects.begin(), objects.end(), cmp[i]);
        for (int j = 0; j < 5; j++) {
            int interval = (int)(div[j] * objects.size());
            auto left_objs = std::vector<Object*>(objects.begin(), objects.begin()+interval);
            auto right_objs = std::vector<Object*>(objects.begin()+interval, objects.end());
            // 光线碰到盒的概率是外盒与内盒表面积之比
            Bounds3 left_bds, right_bds;
            for (auto obj : left_objs) {
                left_bds = Union(left_bds, obj->getBounds());
            }
            for (auto obj : right_objs) {
                right_bds = Union(right_bds, obj->getBounds());
            }
            double cost = 120.f + (left_bds.SurfaceArea()*left_objs.size() + right_bds.SurfaceArea()*right_objs.size()) / bounds.SurfaceArea();
            if (cost < min_cost) {
                min_cost = cost;
                dim = i;
                div_ind = j;
            }
        }
    }
    std::sort(objects.begin(), objects.end(), cmp[dim]);
    int mid = div[div_ind] * objects.size();
    
    auto leftshapes = std::vector<Object*>(objects.begin(), objects.begin() + mid);
    auto rightshapes = std::vector<Object*>(objects.begin() + mid, objects.end());

    node->left = recursiveBuildSAH(leftshapes);
    node->right = recursiveBuildSAH(rightshapes);

    node->bounds = Union(node->left->bounds, node->right->bounds);

    return node;
}

Intersection BVHAccel::Intersect(const Ray& ray) const
{
    Intersection isect;
    if (!root)
        return isect;
    isect = BVHAccel::getIntersection(root, ray);
    return isect;
}

Intersection BVHAccel::getIntersection(BVHBuildNode* node, const Ray& ray) const
{
    // TODO Traverse the BVH to find intersection
    std::array<int, 3> dirIsNeg{ray.direction.x > 0, ray.direction.y > 0, ray.direction.z > 0};
    if (!node->bounds.IntersectP(ray, ray.direction_inv, dirIsNeg)) {
        return Intersection();
    }
    // 叶子节点内仅有一个物体
    /* 
        注意: MeshTriangle作为一个object储存在node->object中，而其中包含物体的Meshes
            在调用MeshTriangle的getIntersection时，调用其bounding box的求交函数Intersect
            (即紧邻上面的函数)，此时BVHAccel为MeshTriangle内的对象，其中的root也不同与scene
            的root，再调用本函数时，即与MeshTriangle中的各个mesh求交。
    */
    if (node->left == nullptr && node->right == nullptr) {
        return node->object->getIntersection(ray);
    }
    
    // 非叶子节点一定左右节点都包含
    Intersection inter_left, inter_right;
    inter_left = getIntersection(node->left, ray);
    inter_right = getIntersection(node->right, ray);

    if (!inter_left.happened && !inter_right.happened) {
        return Intersection();
    }
    if (inter_left.distance < inter_right.distance) {
        return inter_left;
    }
    return inter_right;
}


void BVHAccel::getSample(BVHBuildNode* node, float p, Intersection &pos, float &pdf){
    if(node->left == nullptr || node->right == nullptr){
        node->object->Sample(pos, pdf);
        pdf *= node->area;
        return;
    }
    if(p < node->left->area) getSample(node->left, p, pos, pdf);
    else getSample(node->right, p - node->left->area, pos, pdf);
}

// 叶子节点面积除以根节点面积
void BVHAccel::Sample(Intersection &pos, float &pdf){
    float p = std::sqrt(get_random_float()) * root->area;
    getSample(root, p, pos, pdf);
    pdf /= root->area;
}