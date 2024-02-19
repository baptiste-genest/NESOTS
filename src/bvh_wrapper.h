#ifndef BVH_WRAPPER_H
#define BVH_WRAPPER_H


#include <bvh/v2/bvh.h>
#include <bvh/v2/vec.h>
#include <bvh/v2/ray.h>
#include <bvh/v2/node.h>
#include <bvh/v2/default_builder.h>
#include <bvh/v2/thread_pool.h>
#include <bvh/v2/executor.h>
#include <bvh/v2/stack.h>
#include <bvh/v2/tri.h>

#include <iostream>

namespace BVH_WRAPPER {

using Scalar  = double;
using Vec3    = bvh::v2::Vec<Scalar, 3>;
using BBox    = bvh::v2::BBox<Scalar, 3>;
using Tri     = bvh::v2::Tri<Scalar, 3>;
using Node    = bvh::v2::Node<Scalar, 3>;
using Bvh     = bvh::v2::Bvh<Node>;
using Ray     = bvh::v2::Ray<Scalar, 3>;

using PrecomputedTri = bvh::v2::PrecomputedTri<Scalar>;

class BVH{
public:
    BVH(){}
    BVH(const std::vector<Tri>& tris);
    struct intersection {
        size_t id;
        Scalar u,v,distance;
        bool valid = true;
    };
    intersection get_intersection(const Vec3& origin,const Vec3& dir);
    intersection get_intersection_check_unique(const Vec3& origin,const Vec3& dir,size_t& other_intersect);
    static constexpr size_t invalid_id = std::numeric_limits<size_t>::max();
private:
    // Permuting the primitive data allows to remove indirections during traversal, which makes it faster.
    static constexpr bool should_permute = false;
    std::vector<PrecomputedTri> precomputed_tris;
    Bvh bvh;

};


}

#endif // BVH_WRAPPER_H
