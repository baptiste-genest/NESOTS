#ifndef CONFORMALLAYOUT_H
#define CONFORMALLAYOUT_H

#include "hyperbolicgeometry.h"
#include "yamabeflow.h"
#include <queue>

struct HyperbolicLayout
{
    using triangle = std::array<int,3>;
    using meshTopo = std::vector<triangle>;
    using FaceMapping = std::vector<Face>;
    HyperbolicLorentzGeometry G;

    DiscreteMetric L;

    vecs pos;
    FaceMapping faceMap;
    meshTopo triangles;

    void clear() {
        faceMap.clear();
        pos.clear();
        triangles.clear();
    }

    HyperbolicLayout(){}
    HyperbolicLayout(const DiscreteMetric& L) : L(L) {}

    bool are_neighboors(size_t a,size_t b) {
        auto f1 = faceMap[a];
        auto f2 = faceMap[b];
        for (auto n : f1.adjacentFaces())
            if (n == f2)
                return true;
        return false;
    }

    struct EmbeddedHalfedge {
        EmbeddedHalfedge(Halfedge h,const vec& tail,const vec& tip,int r) : h(h),tail(tail),tip(tip),round(r) {}
        Halfedge h;
        vec tail,tip;
        int round;
    };

    // a = h.tail
    // b = h.tip
    vec getThirdPoint(const EmbeddedHalfedge& h,bool& skipped) const;

    //if max_z < 0, then no max z, note that on the hyperboloid, z starts at 1
    std::vector<int> build(const Mesh& M,int v0 = 0,scalar max_z = -1);

    vec hyperbolicBarycentricCoordinates(const vec& p, const triangle &T) const;

    vec posFromBarycentric(const vec& b,const triangle& T) const {
        vec p = vec::Zero();
        int i = 0;
        for (auto v : T)
            p += pos[v]*b(i++);
        return p;
    }

    vec BarycentricCoordinates(const vec& p,const triangle& T) const {
        const auto& v0 = pos[T[0]];
        const auto& v1 = pos[T[1]];
        const auto& v2 = pos[T[2]];
        vec B;
        // compute the plane's normal
        vec v0v1 = v1 - v0;
        vec v0v2 = v2 - v0;
        // no need to normalize
        vec N = v0v1.cross(v0v2); // N
        scalar denom = N.squaredNorm();

        // edge 1
        vec edge1 = v2 - v1;
        vec vp1 = p - v1;
        vec C = edge1.cross(vp1);
        B(0) = N.dot(C)/denom;
        // edge 2
        vec edge2 = v0 - v2;
        vec vp2 = p - v2;
        C = edge2.cross(vp2);
        B(1) = N.dot(C)/denom; // P is on the right side;
        B(2) = 1-B(0)-B(1);
        return B;
    }

    static scalar angle(scalar a,scalar b,scalar c) {
        using namespace std;
        auto l = (cosh(b)*cosh(c) - cosh(a))/(sinh(b)*sinh(c));
        if (std::abs(l) > 1)
            std::cout << l << std::endl;
        return acos(std::clamp(l,-1.,1.));
    }
    static scalar area(scalar a,scalar b,scalar c) {
        auto th0 = angle(a,b,c);
        auto th1 = angle(c,a,b);
        auto th2 = angle(b,c,a);
        return M_PI - th0 - th1 - th2;
    }


};

#endif // CONFORMALLAYOUT_H
