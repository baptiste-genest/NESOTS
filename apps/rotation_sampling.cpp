#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#ifdef __APPLE__
#include "Eigen/Dense"
#else
#include "eigen3/Eigen/Dense"
#endif
#include "../src/utils.h"
#include "../src/n_spherical_geometry.h"
#include "../src/sampling.h"
#include "../src/ssw.h"
#include "../src/fluid_sim_sdot.h"
#include "../src/point_cloud_io.h"
#include "../src/bvh_wrapper.h"
#include "../src/halton.hpp"

#include <fenv.h>

#include <iostream>
#include <fstream>

using Geometry = NSphericalGeometry;
using point = NSphericalGeometry::point;
using points = NSphericalGeometry::points;

int N = 16;
int d = 4;
Geometry G(d);
Geometry::points MU,NU;
int batch_size = 64;
float epsilon = 0.4;

polyscope::PointCloud* PC;
polyscope::PointCloud* MU_PC;

std::vector<float> SW;

float z_pole = 1.1;

bool toggle = false;
int id = 0;

vecs stereoproj(const Geometry::points& X){
    vecs P(X.size());
    for (int i = 0;i<X.size();i++){
        const auto& x = X[i];
        P[i] = vec(x(0),x(1),x(2))/(z_pole-x(3));
    }
    return P;
}

using Quat = Eigen::Matrix<scalar,4,4>;

Quat toQuat(const point& x) {
    const auto& a = x(0);
    const auto& b = x(1);
    const auto& c = x(2);
    const auto& d = x(3);
    Quat M;
    M << a,-b,-c,-d,
         b,a,-d,c,
         c,d,a,-b,
         d,-c,b,a;
    return M;
}

Quat vecToQuat(const vec& x) {
    point X = point::Zero(4);
    for (int i = 0;i<3;i++)
        X(1+i) = x(i);
    return toQuat(X);
}

Quat toInverseQuat(point X) {
    for (int i = 0;i<3;i++)
        X(i+1) *= -1;
    return toQuat(X/X.squaredNorm());
}
vec pureQuatToVec(const Quat& X) {
    Eigen::Vector<scalar,4> c = X.col(0);
    if (std::abs(c(0)) > 1e-6)
        std::cerr << "[non pure quaternion] " << c(0) << std::endl;
    return vec(c(1),c(2),c(3));
}


Geometry::points super_fibo(int n) {
    Geometry::points SF(n,Vec::Zero(4));
    static auto phi = std::sqrt(2);
    static auto psi = 1.533751168755204288118041;
    for (int i = 0;i<n;i++){
        auto s = i + 0.5;
        auto t = s/n;
        auto d = 2*M_PI*s;
        auto r = sqrt(t);
        auto R = sqrt(std::abs(1.-t));
        auto alpha = d/phi;
        auto beta = d/psi;
        SF[i](0) = r*sin(alpha);
        SF[i](1) = r*cos(alpha);
        SF[i](2) = R*sin(beta);
        SF[i](3) = R*cos(beta);
    }
    return SF;
}


polyscope::SurfaceMesh* INPUT;
std::vector<polyscope::SurfaceMesh*> copies;
std::string meshname;
vecs meshPos;
Mesh input_mesh;
void loadMesh() {
    std::tie(input_mesh.topology, input_mesh.geometry) = readManifoldSurfaceMesh(meshname);
    input_mesh.normalize();
    //INPUT = polyscope::registerSurfaceMesh("input_mesh",
    //                       input_mesh.geometry->vertexPositions, input_mesh.topology->getFaceVertexList());
    auto& T = *input_mesh.topology;
    auto& G = *input_mesh.geometry;
    meshPos.resize(T.nVertices());
    for (auto v : T.vertices())
        meshPos[v.getIndex()] = toVec(G.inputVertexPositions[v]);
}


void rotateMesh(const point& X,int k,const vec& post_offset,const vec& pre_offset) {
    auto q = toQuat(X);
    auto qi = toInverseQuat(X);
    auto P = meshPos;
    for (int i = 0;i<meshPos.size();i++){
        auto v = vecToQuat(P[i]+pre_offset);
        P[i] = pureQuatToVec(q*v*qi) + post_offset;
    }
    polyscope::SurfaceMesh* pc = polyscope::registerSurfaceMesh("rotated mesh " + std::to_string(k),
                           P, input_mesh.topology->getFaceVertexList());
    copies.push_back(pc);
}

Vec x0;
Mat B;
int coord = 0;
Mat getOrthoBasis(const Vec& x) {
    Mat M = Mat::Identity(4,4);
    M.col(coord) = x;
    Mat Q = M.householderQr().householderQ();
    return Q;
}

Vec allBut(const Vec& X,int i){
    Vec Y(X.size()-1);
    for (int j = 0;j<Y.size();j++)
        Y(j) = X[j + (j >= i)];
    return Y;
}

void plotRots(const points& P) {
    int n = sqrt(N);
    int i = 0;
    for (const auto& x : P){
        int sign = ((x0-x).squaredNorm() < (x0+x).squaredNorm())?1:-1;
        Vec l = G.Log(x0,x);
        Vec off = (B.transpose()*l*sign).tail(3);
        //rotateMesh(MU[j*n+i],j*n+i,vec(i*pad,j*pad,0),vec(0,0,0));
        rotateMesh(x,i++,vec(off),vec(0,0,0));
        //rotateMesh(x,i++,vec(-off),vec(0,0,0));
    }
}

void myCallBack() {
    if (ImGui::Button("toggle")) {
        toggle = toggle ^ true;
        if (!toggle)
            std::cout << id << std::endl;
    }
    if (toggle || ImGui::Button("Iteration")){
        for (int i = 0;i<10;i++)
            MU = G.computeOTSliceDirection(MU,NU,epsilon,batch_size,true).first;
        x0 = MU[0];
        B = getOrthoBasis(x0);
        plotRots(MU);
    }
    if (ImGui::Button("reset")){
        MU = Geometry::points(N,Vec::Ones(d).normalized());
        plotRots(MU);
    }
    ImGui::SliderFloat("epsilon",&epsilon,0,1);
}

void init () {
    loadMesh();
    MU = G.samples(N);
    int M = 500'000;
    NU = G.samples(M);
    x0 = MU[0];
    B = getOrthoBasis(x0);
    plotRots(MU);
}


pcg32 GLOBALRNG;
int main(int argc,char** argv) {
    if (argc >= 2){
        meshname = argv[1];
    }
    polyscope::init();
    init();

    polyscope::state::userCallback = myCallBack;
    polyscope::show();
    return 0;
}

/*
 *
void plotRots() {
    int n = sqrt(N);
    for (int j = 0;j<n;j++)
        for (int i = 0;i<n;i++){
            int sign = ((x0-MU[j*n+i]).squaredNorm() < (x0+MU[j*n+i]).squaredNorm())?1:-1;
            Vec l = G.Log(x0,MU[j*n+i]);
            Vec off = (B.transpose()*l*sign).tail(3);
            //rotateMesh(MU[j*n+i],j*n+i,vec(i*pad,j*pad,0),vec(0,0,0));
            rotateMesh(MU[j*n+i],j*n+i,vec(off),vec(0,0,0));
        }
}

*/
