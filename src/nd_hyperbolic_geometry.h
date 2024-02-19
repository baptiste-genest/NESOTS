#ifndef NDHYPERBOLICGEOMETRY_H
#define NDHYPERBOLICGEOMETRY_H

#include "slicedgeometry.h"
#include "n_spherical_geometry.h"
#include "hyperbolicgeometry.h"

class NDHyperbolicGeometry : public SlicedGeometry<Vec,Vec,Vec>
{
public:
    NDHyperbolicGeometry(int d);

    inline scalar dotL(const point& x,const point& y) const{
        scalar D = 0;
        for (int i = 0;i<dim-1;i++)
            D += x(i)*y(i);
        return D - x(dim-1)*y(dim-1);
    }


    inline point sheetToPoincare(const point& x) const{
        return point(x.head(dim-1))/(1+x(dim-1));
    }

    inline point sheetToKlein(const point& x) const{
        return point(x.head(dim-1))/x(dim-1);
    }

    inline point poincareToSheet(const vector& v) const {
        scalar sq = v.squaredNorm();
        vector S(dim);
        S << 2*v,1+sq;
        S /= 1-sq;
        return S;
    }

    static point orthoproj_b_against_a(const point &a, const point &b){
        return b - a*a.dot(b)/a.squaredNorm();
    }

    ///
    /// \brief projectOnSpan project x on span(xO,v)
    /// \return projected
    ///
    inline point projectOnSpan(const point& v,const point& x) const {
        auto e1 = origin;
        auto e2 = orthoproj_b_against_a(e1,v).normalized();
        return e1*e1.dot(x) + e2*e2.dot(x);
    }

    vector MobiusSum(const vector& x,const vector& y) const {
        auto dot = x.dot(y);
        auto nx = x.squaredNorm();
        auto ny = y.squaredNorm();
        return ((1+ 2*dot + ny)*x + (1-nx)*y)/( 1 + 2*dot + nx*ny);
    }

    static inline scalar conformalFactor(const vector& z) {
        return 1/(1-z.squaredNorm());
    }

    const point& Origin() const {
        return origin;
    }

private:
    int dim;
    point origin;

    // SlicedGeometry interface
public:
    points samples(int N) const override;
    slice getSlice() const override;
    vector projectOnTangentSpace(const point &x, const vector &v) const override;
    vector parallelTransport(const point &x, const point &y, const vector &v) const override;
    point advect(const point &x, const point &a, const point &b, const slice &s) const override;
    point projectOnSlice(const point &p, const slice &s) const override;
    scalar curvilinearAbscissaAlongSlice(const point &x, const slice &s) const override;
    scalar distanceAlongSlice(scalar x, scalar y) const override;
    scalar distance(const point &x, const point &y) const override;
    vector Log(const point &x, const point &y) const override;
    point Exp(const point &x, const vector &y, scalar t=1.) const override;
};

#endif // NDHYPERBOLICGEOMETRY_H
