#ifndef NSPHERICALGEOMETRY_H
#define NSPHERICALGEOMETRY_H

#include "slicedgeometry.h"

class NSphericalGeometry : public SlicedGeometry<Vec,twins<Vec>,Vec>
{
public:
    NSphericalGeometry(int d) : d(d){};

    point sampleSphere() const;
    static point orthoproj_b_against_a(const point& a,const point& b);
    static point proj_b_along_a(const point& a,const point& b);

    int d;

    static bool isOrthogonal(const point& x,const point& y) {
        return std::abs(x.dot(y)) < 1e-8;
    }

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
    point Exp(const point &x, const vector&y, scalar t = 1) const override;
    std::pair<points,scalar> computeOTSliceDirection(const points &mu, const points &nu, scalar epsilon, int nb_batches,bool useMedian = true,const std::function<scalar(const point& p)>& step_functionnal =
            [](const point& p) {return 1e8;});
    int computeOptimalCut(const labeled<scalar> &a, const labeled<scalar> &b, int p) const override{
        return computeOptimalCut2(projB(a),projB(b));
    }

    scalar computeOTEnergytoUniform(const points &mu, int nb_batches) const;
};

#endif // NSPHERICALGEOMETRY_H
