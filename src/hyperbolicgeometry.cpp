#include "hyperbolicgeometry.h"

vecs HyperbolicLorentzGeometry::samples(int N) const
{
    vecs X(N);
    for (int i = 0;i<N;i++){
        X[i] = param(uniform_tau(),std::sqrt(uniform_01()));
    }
    return X;
}

vec HyperbolicLorentzGeometry::getSlice() const
{
    auto th = uniform_tau();
    return vec(cos(th),sin(th),0);

}

// Einstein midpoint
// https://openaccess.thecvf.com/content_CVPR_2020/papers/Khrulkov_Hyperbolic_Image_Embeddings_CVPR_2020_paper.pdf
vec HyperbolicLorentzGeometry::getAverage(const vecs &X) const
{
    auto D = apply<vec,complex>(sheetToPoincare,X);
    auto K = apply<complex,complex>(poincareToKlein,D);
    auto L = apply<complex,scalar>(kleinLorentzFactor,K);
    complex avg = 0;
    scalar S = 0;
    for (int i = 0;i<X.size();i++){
        avg += K[i]*L[i];
        S += L[i];
    }
    return poincareToSheet(kleinToPoincare(avg/S));
}

vec HyperbolicLorentzGeometry::advect(const vec &x, const vec &a, const vec &b, const vec &) const
{
    auto za = sheetToPlane(a);
    auto zb = sheetToPlane(b);
    auto d = HyperbolicLorentzGeometry::distance(a,b);
    auto z = sheetToPlane(x);
    int sgn;
    auto E = getGeodesicExtremities(za,zb,sgn);
    if (sgn == 0)
        return x;
    auto M = TranslateAlongGeodesic(
                E,
                d*sgn
                );
    return planeToSheet(M*z);
}

vec HyperbolicLorentzGeometry::projectOnManifold(const vec &p) const
{
    return normalizedL(p);
}

vec HyperbolicLorentzGeometry::projectOnSlice(const vec &x, const vec &v) const
{
    auto p = projectOnSpan(v,x);
    return p/sqrt(std::abs(-dotL(p,p)));
}

vec HyperbolicLorentzGeometry::interpolate(const vec &a, const vec &b, scalar t) const
{
    auto za = sheetToPlane(a);
    auto zb = sheetToPlane(b);
    auto d = HyperbolicLorentzGeometry::distance(a,b);
    int sgn;
    auto E = getGeodesicExtremities(za,zb,sgn);
    auto M = TranslateAlongGeodesic(
                E,
                d*sgn*t
                );
    return planeToSheet(M*za);
}

scalar HyperbolicLorentzGeometry::distanceAlongSlice(scalar t1, scalar t2) const
{
    return std::abs(t1-t2);
}

complex HyperbolicLorentzGeometry::Log(const vec &p1, const vec &p2) const
{
    auto x = sheetToPoincare(p1);
    auto y = sheetToPoincare(p2);

    auto s = MobiusSum(-x,y);
    auto n = std::min(std::abs(s),1.);
    auto l = conformalFactor(x);
    return 4/l * std::atanh(n)*s/n;
}

vec HyperbolicLorentzGeometry::Exp(const vec &x, const complex &v, scalar t) const
{
    auto z = sheetToPoincare(x);
    auto l = conformalFactor(z);
    auto d = std::tanh(l*std::abs(v)*0.5)*v/std::abs(v);
    return poincareToSheet(MobiusSum(z,d));
}

scalar HyperbolicLorentzGeometry::curvilinearAbscissaAlongSlice(const point &x, const slice &s) const
{
    return sgn(x.dot(s))*distance(origin(),x);
}

twins<> HyperbolicLorentzGeometry::getGeodesicExtremities(const complex &z1, const complex &z2,int& sgn){
    const auto& x1 = z1.real();
    const auto& y1 = z1.imag();
    const auto& x2 = z2.real();
    const auto& y2 = z2.imag();
    if (std::abs(x1-x2) < 1e-8){
        sgn = 0;
        return {0,0};
    }
    auto x = 0.5*(x1*x1 - x2*x2 + y1*y1 - y2*y2)/(x1-x2);
    auto r = std::abs(z1-complex(x,0));
    sgn = (x1 < x2) ? 1 : -1;
    /// FIND A BETTER WAY TO REMOVE NANS
    if (r > 1e6)
        sgn = 0;
    return {x-r,x+r};
}

HyperbolicLorentzGeometry::MobiusMap<> HyperbolicLorentzGeometry::TranslateAlongGeodesic(const twins<> &E, scalar lambda)
{
    auto [x1,x2] = E;
    auto c = std::cosh(lambda*0.5);
    auto s = std::sinh(lambda*0.5);

    mat2 M;
    M(0,0) = (x2-x1)*c + (x2+x1)*s;
    M(0,1) = -2*x1*x2*s;
    M(1,0) = 2*s;
    M(1,1) = - (x2+x1)*s + (x2-x1)*c;
    return {M};
}
