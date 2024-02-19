#ifndef CYLINDRICALGEOMETRY_H
#define CYLINDRICALGEOMETRY_H

#include "slicedgeometry.h"
#include <iostream>
#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include <cmath>
#include <iostream>

class CylindricalGeometry : public SlicedGeometry<vec,vec,vec>
{
public:
    vec getSlice() const override;
    vec getAverage(const vecs &X) const ;
    vec advect(const vec& x, const vec &a, const vec &b,const vec& s) const override;
    vec projectOnManifold(const vec &p) const ;
    vec projectOnSlice(const vec &p, const vec &s) const override;
    vec interpolate(const vec &a, const vec &b, scalar t) const ;
    labeled<scalar> orderAlongSlice(const vecs &X, const vec &s) const override;
    scalar distanceAlongSlice(scalar x,scalar y) const override;
    vecs samples(int N) const override;

    inline static std::pair<vec,vec> getEllipseGenerators(const vec& s) {
        auto Z = vec::UnitZ();
        vec steepest = (Z - s(2)*s);
        vec lowest = Eigen::AngleAxisd(M_PI_2,s)*steepest;
        auto a = steepest/xy(steepest).norm();
        auto b = lowest/xy(lowest).norm();
        return {a,b};
    }


    inline static vec projectOnCircle(const vec& x) {
        auto f = std::sqrt(x(0)*x(0)+x(1)*x(1));
        return x/f;
    }

    inline static mat rotZ(scalar th){
        mat R;R.setIdentity();
        R(0,0) = cos(th);
        R(0,1) = -sin(th);
        R(1,0) = sin(th);
        R(1,1) = cos(th);
        return R;
    }

    inline static vec xy(const vec& x){
        return vec(x(0),x(1),0);
    }

    inline static scalar arg(const vec& x){
        return std::atan2(x(1),x(0));
    }

    inline static scalar circle_distance(scalar t1,scalar t2){
        scalar d = std::abs(t1-t2);
        return std::min(d,1-d);
    }



    inline static scalar circle_cross(const vec& x,const vec& y){
        return x(0)*y(1) - x(1)*y(0);
    }

    inline static scalar circle_dot(const vec& x,const vec& y){
        return x(0)*y(0) + x(1)*y(1);
    }
    inline static scalar angle_between(const vec& x,const vec& y){
        //return std::abs((2-(xy(x-y)).squaredNorm())/2);
        if (circle_cross(x,y) > 0)
            return std::acos(circle_dot(x,y));
        else
            return -std::acos(circle_dot(x,y));
    }

    inline static vec slide_along_z(const vec& x,const vec& d,scalar c = 0){
        // <x + t*z,d> = c
        // => t = (c-<x,d>)/<z,d>
        return x + (c- x.dot(d))/d(2)*vec::UnitZ();
    }

    inline vec param(scalar th,scalar z) const {
        return vec(cos(th),sin(th),z);
    }


    // SlicedGeometry interface
public:
    scalar curvilinearAbscissaAlongSlice(const point &x, const slice &s) const override {return 0;}
};

#endif // CYLINDRICALGEOMETRY_H
