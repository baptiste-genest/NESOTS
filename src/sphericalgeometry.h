#ifndef SPHERICALGEOMETRY_H
#define SPHERICALGEOMETRY_H

#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#include "utils.h"
#include "slicedgeometry.h"

#include <iostream>
#include <fstream>

class SphericalGeometry : public SlicedGeometry<vec,vec,vec>
{
public:
    SphericalGeometry() {}

    scalar R = 1;

    vec getSlice() const override;
    vec getAverage(const vecs &X) const ;
    vec advect(const vec& x, const vec &a, const vec &b, const vec &s) const override;
    vec projectOnManifold(const vec &p) const ;
    vec projectOnSlice(const vec &p, const vec &s) const override;
    vec interpolate(const vec &a, const vec &b, scalar t) const ;
    scalar curvilinearAbscissaAlongSlice(const point &x, const slice &s) const override;
    scalar distanceAlongSlice(scalar t1, scalar t2) const override;
    vecs samples(int N) const override;
    vec Log(const vec &x, const vec& y) const override;
    vec Exp(const vec &x, const vec& y,scalar t = 1) const override;
    scalar distance(const vec &a, const vec &b) const override;


    inline static scalar circle_distance(scalar t1,scalar t2){
        scalar d = std::abs(t1-t2);
        return std::min(d,1-d);
    }

    inline static mat bracket(const vec &n)
    {
        mat brack(3,3);
        brack << 0.0 , -n(2), n(1),
                n(2), 0.0 , -n(0),
                -n(1) , n(0),0.0 ;
        return brack;
    }

    static inline mat transition(const vec &a, const vec &b){
        scalar c = a.dot(b);
        if (std::abs( c + 1.0) < 0.00001)
            return -mat::Identity();
        auto vv = a.cross(b);
        mat skew = bracket(vv);
        return mat::Identity() + skew +
                1.0 / (1.0 + c) * skew * skew;
    }
    vector projectOnTangentSpace(const point &x, const vector &v) const override{
        return make_ortho_proj(x)*v;
    }

    vector parallelTransport(const point &x, const point &y, const vector &v) const override {
        return transition(x,y)*v;
    }

    int computeOptimalCut(const labeled<scalar> &a, const labeled<scalar> &b, int p) const override{
        return computeOptimalCut2(projB(a),projB(b));
    }
};


#endif // SPHERICALGEOMETRY_H
