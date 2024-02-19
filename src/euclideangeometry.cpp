#include "euclideangeometry.h"

EuclideanGeometry::EuclideanGeometry(int dim) : SlicedGeometry(),dim(dim)
{

}

EuclideanGeometry::points EuclideanGeometry::samples(int N) const
{
    points P(N);
    for (int i = 0;i<N;i++)
        P[i] = vector::Random(dim);
    return P;
}

EuclideanGeometry::slice EuclideanGeometry::getSlice() const
{
    static std::mt19937 gen;
    static std::normal_distribution<scalar> dist{0.0,1.0};
    Vec dir(dim);
    for (int i = 0;i<dim;i++)
        dir(i) = dist(gen);
    return dir.normalized();
}

EuclideanGeometry::point EuclideanGeometry::getAverage(const points &X) const
{
    return vector();
}

EuclideanGeometry::point EuclideanGeometry::advect(const point &x, const point &a, const point &b, const slice &s) const
{
    return x + (b-a);
}

EuclideanGeometry::point EuclideanGeometry::projectOnManifold(const point &p) const
{
    return p;
}

EuclideanGeometry::point EuclideanGeometry::projectOnSlice(const point &p, const slice &s) const
{
    return s*s.dot(p);
}

EuclideanGeometry::point EuclideanGeometry::interpolate(const point &a, const point &b, scalar t) const
{
    return (1-t)*b + t*a;
}

scalar EuclideanGeometry::distanceAlongSlice(scalar x, scalar y) const
{
    return std::abs(x-y);
}

scalar EuclideanGeometry::distance(const point &x, const point &y) const
{
    return (x-y).norm();
}

EuclideanGeometry::vector EuclideanGeometry::Log(const point &x, const point &y) const
{
    return y - x;
}

EuclideanGeometry::point EuclideanGeometry::Exp(const point &x, const vector &y, scalar t) const
{
    return x + y*t;

}

scalar EuclideanGeometry::curvilinearAbscissaAlongSlice(const point &x, const slice &s) const
{
    return x.dot(s);
}
