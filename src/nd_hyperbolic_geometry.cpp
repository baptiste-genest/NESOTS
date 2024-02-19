#include "nd_hyperbolic_geometry.h"

NDHyperbolicGeometry::NDHyperbolicGeometry(int d) : dim(d)
{
    origin = Vec::Zero(d);origin(d-1) = 1;

}

NDHyperbolicGeometry::points NDHyperbolicGeometry::samples(int N) const
{
    points S(N,Vec(dim));
    for (int i= 0;i<N;i++){
        auto t = sqrt(uniform_01());
        S[i] = cosh(t)*origin + sinh(t)*getSlice();
    }
    return S;
}

NDHyperbolicGeometry::slice NDHyperbolicGeometry::getSlice() const
{
    auto v = NSphericalGeometry(dim-1).sampleSphere();
    slice s(dim);
    s << v , 0;
    return s;
}

NDHyperbolicGeometry::vector NDHyperbolicGeometry::projectOnTangentSpace(const point &x, const vector &v) const
{
    return v;
}

NDHyperbolicGeometry::vector NDHyperbolicGeometry::parallelTransport(const point &x, const point &y, const vector &v) const
{
    return v;
}

NDHyperbolicGeometry::point NDHyperbolicGeometry::advect(const point &x, const point &a, const point &b, const slice &) const
{
    if ((b-a).norm() < 1e-6)
        return x;
    const auto& e2 = origin;
    const auto& e1 = orthoproj_b_against_a(e2,b-a).normalized();
    point xO = orthoproj_b_against_a(e1,x);
    xO = orthoproj_b_against_a(e2,xO);
    vector xr = x - xO;
    auto th = distance(a,b);
    auto c = std::cosh(th);
    auto s = std::sinh(th);
    scalar ca = xr.dot(e1);
    scalar cb = xr.dot(e2);
    point rslt = xO + e1*(c*ca + s*cb) + e2*(s*ca + c*cb);
    return rslt;
}


NDHyperbolicGeometry::point NDHyperbolicGeometry::projectOnSlice(const point &x, const slice &v) const
{
    point p = projectOnSpan(v,x);
    return p/sqrt(std::abs(-dotL(p,p)));
}

scalar NDHyperbolicGeometry::curvilinearAbscissaAlongSlice(const point &x, const slice &s) const
{
    return sgn(x.dot(s))*distance(origin,x);
}

scalar NDHyperbolicGeometry::distanceAlongSlice(scalar x, scalar y) const
{
    return std::abs(x-y);
}

scalar NDHyperbolicGeometry::distance(const point &x, const point &y) const
{
    auto d = std::max(1.,-dotL(x,y));
    return std::acosh(d);
}

/*
NDHyperbolicGeometry::vector NDHyperbolicGeometry::Log(const point &p1, const point &p2) const
{
    vector x = sheetToPoincare(p1);
    vector y = sheetToPoincare(p2);

    vector s = MobiusSum(-x,y);
    auto n = std::clamp(s.norm(),1e-8,1.);
    auto l = conformalFactor(x);
    return 2/l * std::atanh(n)*s/n;
}

NDHyperbolicGeometry::point NDHyperbolicGeometry::Exp(const point &x, const vector &v, scalar t) const
{
    point z = sheetToPoincare(x);
    auto l = conformalFactor(z);
    vector d = std::tanh(l*v.norm()*0.5*t)*v.normalized();
    return poincareToSheet(MobiusSum(z,d));
}
*/

NDHyperbolicGeometry::vector NDHyperbolicGeometry::Log(const point &p1, const point &p2) const
{
    scalar d= dotL(p1,p2);
    scalar L = std::acosh(std::max(1.,-d))/std::sqrt(std::max(1e-5,d*d-1));
    return L*(p2 + d*p1);
}

NDHyperbolicGeometry::point NDHyperbolicGeometry::Exp(const point &x, const vector &v, scalar t) const
{
    auto norml = std::sqrt(std::max(1e-5,dotL(v,v)));
    return std::cosh(norml*t)*x + std::sinh(norml*t)*v/norml;
}
