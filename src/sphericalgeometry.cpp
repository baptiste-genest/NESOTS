#include "sphericalgeometry.h"

SphericalGeometry::vector SphericalGeometry::getSlice() const
{
    return uniform_sphere();
}

SphericalGeometry::vector SphericalGeometry::getAverage(const vecs &X) const
{
    SphericalGeometry::vector s = vec::Zero();
    for (const auto& x : X)
        s += x;
    return projectOnManifold(s);
}

SphericalGeometry::vector SphericalGeometry::advect(const vec& x, const SphericalGeometry::vector &a, const SphericalGeometry::vector &b, const SphericalGeometry::vector &) const
{
    return transition(a,b)*x;
}

SphericalGeometry::vector SphericalGeometry::projectOnManifold(const SphericalGeometry::vector &p) const
{
    return p.normalized();
}

SphericalGeometry::vector SphericalGeometry::projectOnSlice(const SphericalGeometry::vector &p, const SphericalGeometry::vector &s) const
{
    return (p-s*s.dot(p)).normalized();
}

SphericalGeometry::vector SphericalGeometry::interpolate(const SphericalGeometry::vector &a, const SphericalGeometry::vector &b, scalar t) const
{
    auto th = distance(a,b);
    auto s = sin(th);
    if (std::abs(s) < 1e-10)
        return a;
    vector r = (sin((1-t)*th)*a + sin(t*th)*b)/s;
    if (isNan(r))
        assert(false);
    return r.normalized();
}

scalar SphericalGeometry::curvilinearAbscissaAlongSlice(const point &x, const slice &s) const
{
    auto p = projectOnSlice(x,s);
    return (M_PI + std::atan2(-p(1),-p(0)))/(2*M_PI);
}

scalar SphericalGeometry::distanceAlongSlice(scalar t1, scalar t2) const
{
    return circle_distance(t1,t2);
}

vecs SphericalGeometry::samples(int N) const
{
    vecs S(N);
    for (auto& x : S)
        x = uniform_sphere();
    return S;
}

SphericalGeometry::vector SphericalGeometry::Log(const SphericalGeometry::vector &a, const SphericalGeometry::vector &b) const
{
    return ortho_proj_b_against_a(a,(b-a)).normalized()*distance(a,b); // Log_a(b)
}

SphericalGeometry::vector SphericalGeometry::Exp(const SphericalGeometry::vector &p, const SphericalGeometry::vector &v, scalar t) const
{
    if (std::abs(p.dot(v)) > 1e-6){
        std::cerr << "v is not in the tangent plane of p" << std::endl;
        std::cerr <<"[p] " <<  p.transpose() << std::endl;
        std::cerr <<"[p norm] " <<  p.norm() << std::endl;
        std::cerr <<"[v] " << v.transpose() << std::endl;
    }
    t *= v.norm();
    return p*cos(t) + sin(t)*v.normalized();
}

scalar SphericalGeometry::distance(const SphericalGeometry::vector &a, const SphericalGeometry::vector &b) const{
    auto d = a.dot(b);
    d = std::clamp<scalar>(d,-1.,1.);
    return std::acos(d);
}
