#ifndef HYPERBOLICGEOMETRY_H
#define HYPERBOLICGEOMETRY_H

#include "slicedgeometry.h"

class HyperbolicLorentzGeometry : public SlicedGeometry<vec,vec,complex>
{
    // SlicedGeometry interface
public:
    HyperbolicLorentzGeometry() : SlicedGeometry() {}
    vecs samples(int N) const override;
    vec getSlice() const override;
    vec getAverage(const vecs &X) const;
    vec advect(const vec &x, const vec &a, const vec &b, const vec &s) const override;
    vec projectOnManifold(const vec &p) const ;
    vec projectOnSlice(const vec &p, const vec &s) const override;
    vec interpolate(const vec &a, const vec &b, scalar t) const;
    scalar distanceAlongSlice(scalar t1, scalar t2) const override;
    virtual int computeOptimalCut(const labeled<scalar> &a, const labeled<scalar> &b, int p) const override {return 0;}
    complex Log(const vec &x, const vec &y) const override;
    vec Exp(const vec &x, const complex &y, scalar t = 1.) const override;
    scalar curvilinearAbscissaAlongSlice(const point &x, const slice &s) const override;

    static point orthoproj_b_against_a(const point &a, const point &b){
        return b - a*a.dot(b)/a.squaredNorm();
    }

    template<class T = scalar>
    struct MobiusMap {
        Eigen::Matrix<T,2,2> M;
        complex operator*(const complex& z) const {
            return (M(0,0)*z + M(0,1))/(M(1,0)*z + M(1,1));
        }
    };


    static inline vec geodesic(const vec& x,const vec& v,scalar t=1) {
        return std::cosh(t)*x + std::sinh(t)*v;
    }

    static inline vec param(scalar th,scalar t) {
        return std::cosh(t)*origin() + std::sinh(t)*(vec(cos(th),sin(th),0));
    }

    static inline vec origin() {return vec::UnitZ();}

    inline scalar normL(const vec& x) const {
        return distance(origin(),x);
    }

    inline vec normalizedL(const vec& x) const {
        return x/std::sqrt(normL(x));
    }

    static inline vec projectOnSpan(const vec& v,const vec& x) {
        return make_ortho_proj(origin().cross(v))*x;
    }

    static inline scalar dotL(const vec& x,const vec& y) {
        return x(0)*y(0) + x(1)*y(1) - x(2)*y(2);
    }
    inline virtual scalar distance(const vec& x,const vec& y) const override {
        auto d = std::max(1.,-dotL(x,y));
        return std::acosh(d);
    }

    static twins<> getGeodesicExtremities(const complex& z1,const complex& z2,int& sgn);

    static MobiusMap<> TranslateAlongGeodesic(const twins<>& E,scalar lambda);

    static inline complex sheetToPoincare(const vec& x) {
        return complex(x(0)/(1+x(2)),x(1)/(1+x(2)));
    }
    static inline complex poincareToPlane(const complex& z) {
        return (-1.+I*z)/(1.+I*z)*(-I);
    }
    static inline complex sheetToKlein(const vec& x) {
        return poincareToKlein(sheetToPoincare(x));
    }

    static inline complex poincareToKlein(const complex& z) {
        return 2.*z/(1+std::norm(z));
    }
    static inline complex kleinToPoincare(const complex& z) {
        return z/(1.+std::sqrt(std::max(0.,1-std::norm(z))));
    }

    static inline complex sheetToPlane(const vec& x){
        return poincareToPlane(sheetToPoincare(x));
    }

    static inline scalar kleinLorentzFactor(const complex& x){
        return 1/std::sqrt(1-std::norm(x));
    }

    static inline complexs sheetToPlane(const vecs& X){
        complexs H(X.size());
        for (int i = 0;i<X.size();i++)
            H[i] = sheetToPlane(X[i]);
        return H;
    }

    static inline scalar conformalFactor(const complex& z) {
        return 1/(1-std::norm(z));
    }

    static inline complex planeToPoincare(const complex& z){
        return (1.+I*z)/(1.-I*z)*(-I);
    }
    static inline vec poincareToSheet(const complex& z){
        auto sq = std::norm(z);
        return vec(2*z.real(),2*z.imag(),1+sq)/(1-sq);
    }
    static inline vec planeToSheet(const complex& z){
        return poincareToSheet(planeToPoincare(z));
    }

    static inline complex MobiusSum(const complex& x,const complex& y) {
        auto dot = x.real()*y.real() + x.imag()*y.imag();
        auto nx = std::norm(x);
        auto ny = std::norm(y);
        return ((1+ 2*dot + ny)*x + (1-nx)*y)/( 1 + 2*dot + nx*ny);
    }


    complex convertToRiemannianGradient(const vec &x, const complex &v) const{
        return v;
        auto l = conformalFactor(sheetToPoincare(x));
        return 1/(1-l*l)*v;
    }

    vector projectOnTangentSpace(const point &x, const vector &v) const override {return vector();}
    vector parallelTransport(const point &x, const point &y, const vector &v) const override {return vector();}
};

#endif // HYPERBOLICGEOMETRY_H
