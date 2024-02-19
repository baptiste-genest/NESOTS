#ifndef EUCLIDEANGEOMETRY_H
#define EUCLIDEANGEOMETRY_H

#include "slicedgeometry.h"

class EuclideanGeometry : public SlicedGeometry<Vec,Vec,Vec>
{
public:
    EuclideanGeometry(int dim);

    // SlicedGeometry interface
public:
    points samples(int N) const override;
    slice getSlice() const override;
    point getAverage(const points &X) const;
    point advect(const point &x, const point &a, const point &b, const slice &s) const override;
    point projectOnManifold(const point &p) const;
    point projectOnSlice(const point &p, const slice &s) const override;
    point interpolate(const point &a, const point &b, scalar t) const;
    scalar distanceAlongSlice(scalar x, scalar y) const override;
    scalar distance(const point &x, const point &y) const override;
    vector Log(const point &x, const point &y) const override;
    point Exp(const point &x, const vector &y, scalar t= 1) const override;
    scalar curvilinearAbscissaAlongSlice(const point &x, const slice &s) const override;
    vector projectOnTangentSpace(const point &x, const vector &v) const override {return v;}
    vector parallelTransport(const point &x, const point &y, const vector &v) const override {return v;}

private:
    int dim;
};

#endif // EUCLIDEANGEOMETRY_H
