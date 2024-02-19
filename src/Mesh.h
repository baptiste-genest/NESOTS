#ifndef MESH_H
#define MESH_H

//#include "sphericalgeometry.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/surface_mesh.h"
#include "utils.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

inline vec toVec(const Vector3& x){
    return vec(x.x,x.y,x.z);
}

template<class T>
using mesh_map = std::map<int,T>;

using scalar_mesh_map = mesh_map<scalar>;
using vec_mesh_map = mesh_map<vec>;

template<class T>
T& set(mesh_map<T>& M,int i){
    int k = 0;
    for (auto& kv : M){
        if (i == k)
            return kv.second;
        k++;
    }
}

template<class T>
T get(const mesh_map<T>& M,int i){
    int k = 0;
    for (const auto& kv : M){
        if (i == k)
            return kv.second;
        k++;
    }
    std::cerr << "map out index" << std::endl;
    return T();
}

template<class T>
int get_key(const mesh_map<T>& M,int i){
    int k = 0;
    for (const auto& kv : M){
        if (i == k)
            return kv.first;
        k++;
    }
    return -1;
}



struct Mesh {
    std::unique_ptr<ManifoldSurfaceMesh> topology;
    std::unique_ptr<VertexPositionGeometry> geometry;

    scalars CurvatureBasedDensity(scalar t = 1.) const {
        scalars A(topology->nFaces(),0);
        geometry->requireVertexMeanCurvatures();
        for (auto f : topology->faces()){
            A[f.getIndex()] = geometry->faceArea(f);
            for (auto v : f.adjacentVertices())
                A[f.getIndex()]  += t*std::abs(geometry->vertexMeanCurvature(v));
        }
        return A;
    }


    scalars faceAreas() const {
        scalars A(topology->nFaces(),0);
        for (auto f : topology->faces())
            A[f.getIndex()] = geometry->faceArea(f);
        return A;
    }

    mat tangentBasis(int face) const {
        geometry->requireVertexPositions();
        mat M;
        auto F = topology->face(face);
        auto V = F.adjacentVertices();
        const auto& pos = geometry->vertexPositions;
        std::vector<vec> P;
        for (auto v : V)
            P.push_back(toVec(pos[v]));

        M.col(0) = P[1]-P[0];
        M.col(1) = P[2]-P[0];
        M.col(2) = M.col(0).cross(M.col(1)).normalized();
        return M;
    }

    inline vec vertex(const Vertex& v) const {
        return toVec(geometry->vertexPositions[v]);
    }


    vec_mesh_map vertexInFaceMap(const Face& face) const {
        vec_mesh_map V;
        for (const auto& v : face.adjacentVertices())
            V[v.getIndex()] = vertex(v);
        return V;
    }

    std::array<vec,3> vertexInFace(const Face& face) const {
        std::array<vec,3> V;
        int i = 0;
        for (const auto& v : face.adjacentVertices())
            V[i++] = vertex(v);
        return V;
    }

    scalar distanceToTriangle(const vec& p,const Face& face) const{
        auto V = vertexInFace(face);
        const vec& a = V[0];
        const vec& b = V[1];
        const vec& c = V[2];
        static const auto cross = [] (const vec& a,const vec& b) {
            return a.cross(b);};
        static const auto dot = [] (const vec& a,const vec& b) {
            return a.dot(b);};
        static const auto dot2 = [] (const vec& a) {
            return a.dot(a);};
        static const auto sign = [] (const scalar& a) {
            return (a < 0) ? -1 : 1;};
        vec ba = b - a; vec pa = p - a;
        vec cb = c - b; vec pb = p - b;
        vec ac = a - c; vec pc = p - c;
        vec nor = cross( ba, ac );

          return sqrt(
            (sign(dot(cross(ba,nor),pa)) +
             sign(dot(cross(cb,nor),pb)) +
             sign(dot(cross(ac,nor),pc))<2.0)
             ?
             std::min( std::min(
             dot2(ba*clamp<scalar>(dot(ba,pa)/dot2(ba),0.0,1.0)-pa),
             dot2(cb*clamp<scalar>(dot(cb,pb)/dot2(cb),0.0,1.0)-pb) ),
             dot2(ac*clamp<scalar>(dot(ac,pc)/dot2(ac),0.0,1.0)-pc) )
             :
             dot(nor,pa)*dot(nor,pa)/dot2(nor) );
    }

    Face closestFace(const vec& p) const {
        Face F;
        scalar best = 1e7;
        for (const auto& f : topology->faces()){
            auto d = distanceToTriangle(p,f);
            if (d < best) {
                best = d;
                F = f;
            }
        }
        return F;
    }

    Face containingFace(const vec& p) {
        scalar scale = 100;
        for (const auto& f : topology->faces()){
            auto V = vertexInFace(f);
            if (rayTriangleIntersection(p*scale,V[0]*scale,V[1]*scale,V[2]*scale))
                return f;
        }
        std::cout << "no containing triangle" << std::endl;
        return Face();
    }

    scalar_mesh_map distanceWeights(const vec& p,Face f) const{
        auto V = vertexInFaceMap(f);
        scalar_mesh_map W;
        scalar s = 0;
        for (const auto& m : V){
            W[m.first] = (p-m.second).squaredNorm();
            s += W[m.first];
        }
        for (auto& v : W)
            v.second /= s;
        return W;
    }

    vec Barycentric(const vec& p, Face f) const
    {
        auto V = vertexInFace(f);
        vec v0 = V[1] - V[0], v1 = V[2] - V[0], v2 = p - V[0];
        //v2 = v2 - proj_b_along_a(v0.cross(v1),v2);
        scalar d00 = v0.dot(v0);
        scalar d01 = v0.dot(v1);
        scalar d11 = v1.dot(v1);
        scalar d20 = v2.dot(v0);
        scalar d21 = v2.dot(v1);
        scalar denom = d00 * d11 - d01 * d01;
        vec L;
        L(0) = (d11 * d20 - d01 * d21) / denom;
        L(1) = (d00 * d21 - d01 * d20) / denom;
        L(2) = 1.0f - L(0) - L(1);
        return L;
    }

    vec Barycentric2(const vec& p,Face f) const {
        auto V = vertexInFace(f);
        const auto& v0 = V[0];
        const auto& v1 = V[1];
        const auto& v2 = V[2];
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


    scalar_mesh_map BarycentricMap(const vec& p, Face f) const
    {
        auto VM = vertexInFaceMap(f);
        vecs V;
        for (const auto& v : VM )
            V.push_back(v.second);
        vec v0 = V[1] - V[0], v1 = V[2] - V[0], v2 = p - V[0];
        v2 = v2 - proj_b_along_a(v0.cross(v1),v2);
        scalar d00 = v0.dot(v0);
        scalar d01 = v0.dot(v1);
        scalar d11 = v1.dot(v1);
        scalar d20 = v2.dot(v0);
        scalar d21 = v2.dot(v1);
        scalar denom = d00 * d11 - d01 * d01;
        scalar_mesh_map L;
        L[get_key(VM,0)] = (d11 * d20 - d01 * d21) / denom;
        L[get_key(VM,1)] = (d00 * d21 - d01 * d20) / denom;
        L[get_key(VM,2)] = 1.0f - get(L,0) - get(L,1);
        return L;
    }

    scalar diameter(Face f) const {
        scalar M = 0;
        for (auto e : f.adjacentEdges())
            M = std::max(M,geometry->edgeLength(e));
        return M;
    }


    /*
    vec sphericalBarycentric(const vec& p, Face f) const
    {
        auto V = vertexInFace(f);
        auto denom = sphericalTriangleArea(V[0],V[1],V[2]);
        vec L;
        L(0) =  sphericalTriangleArea(p,V[1],V[2])/denom;
        L(1) =  sphericalTriangleArea(V[0],p,V[2])/denom;
        L(2) =  1 - L(0) - L(1);
        return L;
    }
    */

    bool rayTriangleIntersection(const vec& dir,const vec& v0,const vec& v1,const vec& v2) const {  // compute the plane's normal
        vec v0v1 = v1 - v0;
        vec v0v2 = v2 - v0;
        // no need to normalize
        vec N = v0v1.cross(v0v2); // N
        scalar area2 = N.norm();
        vec orig = vec::Zero();
        scalar t;

        // Step 1: finding P

        // check if the ray and plane are parallel.
        scalar NdotRayDirection = N.dot(dir);
        if (std::abs(NdotRayDirection) < 1e-8) // almost 0
            return false; // they are parallel, so they don't intersect!

        // compute d parameter using equation 2
        scalar d = -N.dot(v0);

        // compute t (equation 3)
        t = -(N.dot(orig) + d) / NdotRayDirection;

        // check if the triangle is behind the ray
        if (t < 0) return false; // the triangle is behind

        // compute the intersection point using equation 1
        vec P = orig + t * dir;

        // Step 2: inside-outside test
        vec C; // vector perpendicular to triangle's plane

        // edge 0
        vec edge0 = v1 - v0;
        vec vp0 = P - v0;
        C = edge0.cross(vp0);
        if (N.dot(C) < 0) return false; // P is on the right side

        // edge 1
        vec edge1 = v2 - v1;
        vec vp1 = P - v1;
        C = edge1.cross(vp1);
        if (N.dot(C) < 0)  return false; // P is on the right side

        // edge 2
        vec edge2 = v0 - v2;
        vec vp2 = P - v2;
        C = edge2.cross(vp2);
        if (N.dot(C) < 0) return false; // P is on the right side;

        return true; // this ray hits the triangle
    }

    vec barycenter() const {
        vec avg = vec::Zero();
        for (const auto& v : topology->vertices())
            avg += vertex(v);
        return avg/topology->nVertices();
    }

    scalar area() const {
        scalar A = 0;
        for (const auto& f : topology->faces())
            A += geometry->faceArea(f);
        return A;
    }

    scalar totalEdgeLength() const {
        scalar EL = 0;
        for (const auto& e : topology->edges())
            EL += geometry->edgeLength(e);
        return EL;
    }



    vec posFromWeights(const scalar_mesh_map& B,Face f) const {
        vec p = vec::Zero();
        for (auto v : f.adjacentVertices())
            p += toVec(geometry->vertexPositions[v])*B.at(v.getIndex());
        return p;
    }

    vec posFromWeights(const vec& B,Face f) const {
        vec p = vec::Zero();
        auto V = vertexInFace(f);
        for (int i = 0;i<3;i++)
            p += V[i]*B(i);
        return p;
    }

    inline Vector3 toVector3(const vec& x) const {
        Vector3 rslt;
        rslt.x = x(0);
        rslt.y = x(1);
        rslt.z = x(2);
        return rslt;
    }

    vec faceBarycenter(Face f) const {
        return posFromWeights(vec::Ones()/3,f);
    }

    vec bary = vec::Zero();
    scalar A = 1;

    vec MapOnNotNormalized(vec x) const {
        x = x*A + bary;
        return x;
    }

    vecs MapOnNotNormalized(vecs X) const {
        for (auto& x : X)
            x = MapOnNotNormalized(x);
        return X;
    }

    void normalize(scalar f = 1) {
        bary = barycenter();
        A = sqrt(area());
        for (const auto& v : topology->vertices())
            geometry->vertexPositions[v] = toVector3(vertex(v)-bary)/A*f;
    }


};


#endif // MESH_H
