#include "cylindricalgeometry.h"

vec CylindricalGeometry::getSlice() const
{
    //return uniform_sphere();
    auto th = uniform_01()*2*M_PI;
    return param(th,uniform()*10000).normalized();
    vec d = uniform_sphere();
    if (std::abs(d(2)) < 0.1)
        return getSlice();
    return d.normalized();
}

vec CylindricalGeometry::getAverage(const vecs &X) const
{
    scalar sz = 0;
    vec s = vec::Zero();
    for (const auto& x : X){
        s += xy(x);
        sz += x(2);
    }
    s.normalize();
    return vec(s(0),s(1),sz/X.size());
}

vec CylindricalGeometry::advect(const vec& x, const vec &a, const vec &b, const vec &s) const
{
    auto th = angle_between(a,b);
    return rotZ(th)*x + vec::UnitZ()*(b(2)-a(2));
}

vec CylindricalGeometry::projectOnManifold(const vec &x) const
{
    auto f = std::sqrt(x(0)*x(0)+x(1)*x(1));
    return vec(x(0)/f,x(1)/f,x(2));
}

vec CylindricalGeometry::projectOnSlice(const vec &p, const vec &s) const
{
    return projectOnCircle(p - s*s.dot(p));
}

vec CylindricalGeometry::interpolate(const vec &a, const vec &b, scalar t) const
{
    auto th = angle_between(a,b);
    vec r = rotZ(th*t)*a;
    r(2) = (1-t)*a(2) + t*b(2);
    return r;
}

labeled<scalar> CylindricalGeometry::orderAlongSlice(const vecs &X, const vec &s) const
{
    auto [av,bv] = getEllipseGenerators(s);
    auto a = av.norm();
    auto b = bv.norm();
    auto m = 1 - b*b/a/a;
    auto Em = std::comp_ellint_2(m);
    labeled<scalar> Y(X.size());
    auto L = a*Em*4;
    auto k = a/b;
    for (int i = 0;i<X.size();i++){
        auto th = angle_between(av,X[i]);
        if (th < 0)
            th = 2*M_PI+th;
        auto T = std::atan2(k*std::sin(th),std::cos(th));
        Y[i] = {i, a*(Em + std::ellint_2(m,T - M_PI_2))/L+0.5};
    }
    return Y;
}

scalar CylindricalGeometry::distanceAlongSlice(scalar x, scalar y) const {
    return circle_distance(x,y);
}

vecs CylindricalGeometry::samples(int N) const
{
    vecs X(N);
    for (auto& x : X)
        x = param(uniform_01()*2*M_PI,uniform_01()*4-2);
    return X;
}
