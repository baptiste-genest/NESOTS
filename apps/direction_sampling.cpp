#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#ifdef __APPLE__
#include "Eigen/Dense"
#else
#include "eigen3/Eigen/Dense"
#endif
#include "../src/nesots.h"

#include <fenv.h>

#include <iostream>
#include <fstream>

#include "../deps/CLI11.hpp"

using Geometry = NSphericalGeometry;

int N = 16;
int d = 3;
Geometry G(d);
using points = Geometry::points;
points MU,NU;
int batch_size = 64;
float epsilon = 0.4;

bool toggle = false;

int id = 0;

vec get_ortho(const vec& x) {
    return vec(
        std::copysign(x(2),x(0)),
        std::copysign(x(2),x(1)),
        -std::copysign(std::abs(x(0)) + std::abs(x(1)),x(2))
        ).normalized();
}

void plotPlanes(const points& X,glm::vec3 c,std::string label,vec offset = vec::Zero(),bool centered = false){
    std::vector<std::vector<int>> F(1);
    int M = 50;
    F[0].resize(M);
    std::iota(F[0].begin(),F[0].end(),0);

    if (!centered)
        polyscope::registerPointCloud(label,vecs{offset})->setPointRadius(0.5,false);

    for (int i = 0;i< N;i++){
        vec e1 = get_ortho(X[i])*0.2;
        vec e2 = e1.cross((vec)X[i]);

        vecs P;
        for (float t = 0;t<2*M_PI;t+= 2*M_PI/M)
            P.push_back(e1*cos(t) + e2*sin(t) + offset);
        if (!centered) {
            auto Pp = P;
            for (auto& x : Pp)
                x += X[i];
            polyscope::registerSurfaceMesh(label+std::to_string(i) + "p",Pp,F)->setSurfaceColor(c);
            auto Pn = P;
            for (auto& x : Pn)
                x -= X[i];
            polyscope::registerSurfaceMesh(label+std::to_string(i) + "n",Pn,F)->setSurfaceColor(c);
        }
        else
            polyscope::registerSurfaceMesh(label+std::to_string(i),P,F)->setSurfaceColor(c);

    }
}

void plotLines(const points& X,glm::vec3 c,std::string label,vec offset = vec::Zero()) {
    polyscope::registerPointCloud("center "+label,vecs{offset})->setPointRadius(1,false);
    //auto color = C->getPointColor();
    //std::cout << color.x << " " << color.y << " " << color.z << std::endl;
    for (int i = 0;i< N;i++){
        vecs line = {-X[i]*2 + offset,X[i]*2 + offset};
        polyscope::registerCurveNetworkLine(label+std::to_string(i),line)->setColor(c);
    }
}

points toPoints(const vecs& X){
    points rslt(X.size());
    for (int i = 0;i<N;i++)
        rslt[i] = (vec)X[i];
    return rslt;
}

points F = toPoints(FibonacciSphere(N));

void plot2DLines(const points& X,glm::vec3 c,std::string label,vec2 offset) {
    auto rot = [](Vec x) {
        std::swap(x(0),x(1));
        x(0)*=-1;
        return x;
    };

    polyscope::registerPointCloud(label,X);
    /*
    auto Y = X;
    for (auto& y : Y) y *= -1;
    polyscope::registerPointCloud("neg",Y);
    */
    polyscope::registerPointCloud("center",vecs{vec::Zero()})->setPointRadius(1,false);
    for (int i = 0;i< X.size();i++){
        scalar d= -X[i](2);
        vec2 n = X[i].head(2);
        vec2 l = rot(n).normalized();
        vec2 x = l + n*d/n.squaredNorm();
        auto check = std::abs(x.dot(n) - d);
        if (check > 1e-3)
            std::cout << "not on line " << check << std::endl;
        std::vector<vec2> line = {l + n*d/n.squaredNorm() + offset,-l + n*d/n.squaredNorm() + offset};
        polyscope::registerCurveNetworkLine2D(label+std::to_string(i),line)->setColor(c)->setRadius(0.01,false);
    }
}

points Transform(points X,scalar s,vec offset){
    for (auto& x : X){
        x *= s;
        x += offset;
    }
    return X;
}

polyscope::PointCloud* plotClouds(const points& X,glm::vec3 color,std::string label,vec offset) {
    auto C = polyscope::registerPointCloud("centerC "+label,vecs{offset})->setPointRadius(1,false);
    polyscope::registerPointCloud(label,Transform(X,1,offset))->setPointColor(color)->setPointRadius(0.05,false);
    color.x = 1-color.x;
    color.y = 1-color.y;
    color.z = 1-color.z;
    polyscope::registerPointCloud("opposite " + label,Transform(X,-1,offset))->setPointColor(color)->setPointRadius(0.05,false);
    return C;
}

bool lines3D = false;
bool lines2D = false;
bool clouds = false;
bool planes = false;

void myCallBack() {
    if (ImGui::Button("toggle")) {
        toggle = toggle ^ true;
        if (!toggle)
            std::cout << id << std::endl;
    }
    if (toggle || ImGui::Button("Iteration")){
        auto rslt = G.computeOTSliceDirection(MU,NU,epsilon,batch_size,true);
        MU = rslt.first;
        id ++;
        ImGui::SliderFloat("epsilon",&epsilon,0,1);
        if (lines3D) {
            plotLines(MU,glm::vec3(1,0,0),"ours",vec(3,-3,0));
            plotLines(toPoints(FibonacciSphere(N)),glm::vec3(0,0,1),"fibo",vec(0,-3,0));
        }
        else if (lines2D) {
            plot2DLines(MU,glm::vec3(1,0,0),"ours",vec2(2,0));
            plot2DLines(toPoints(FibonacciSphere(N)),glm::vec3(0,0,1),"fibo",vec2(0,0));
        }
        else if (planes) {
            plotPlanes(MU,glm::vec3(1,0,0),"ours",vec(2,0,0));
            plotPlanes(toPoints(FibonacciSphere(N)),glm::vec3(0,0,1),"fibo",vec(0,0,0));
        }
        if (true) {
            auto c1 = plotClouds(MU,glm::vec3(1,0,0),"ours",vec(3,0,0));
            auto c2 = plotClouds(toPoints(FibonacciSphere(N)),glm::vec3(0,0,1),"fibo",vec(0,0,0));
            c1->setPointColor(glm::vec3(0.109,0.388,0.89));
            c2->setPointColor(glm::vec3(0.89,0.758,0.109));
        }
    }
}

points loadHoughTransform() {
    std::ifstream file("/tmp/hough_transform.lines");
    points X;
    double x,y,z;
    while (file >> x >> y >> z){
        X.push_back(vec(x,y,z).normalized());
        X.push_back(vec(-X.back()));
    }
    return X;
}


void init () {
    /*
    NU = PointCloudIO::read_point_cloud("/tmp/nu.pts",3);
    auto nu_size = NU.size();
    for (int i  = 0;i<nu_size;i++){
        auto x = NU[i];
        x *= -1;
        NU.push_back(x);
    }
*/
    NU = G.samples(500'000);
    //NU = loadHoughTransform();
    //plot2DLines(NU,glm::vec3(0,0,1),"hough",vec2(3,0));
    /*
    for (auto& x : NU){
        x(2) *= 8;
        x.normalize();
    }
*/

    //plotLines()
   // plot2DLines(G.sub_sample(NU,200),glm::vec3(0,0,1),"hough",vec2(3,0));
    MU = G.sub_sample(NU,N);

}

pcg32 GLOBALRNG;

int main(int argc,char** argv) {
    CLI::App app("Projective Plane Sampling");

    app.add_option("--sample_size", N, "Number of samples to generate");
    app.add_flag("--lines3D", lines3D, "3D lines");
    app.add_flag("--lines2D", lines2D, "affine 2D lines");
    app.add_flag("--planes", planes, "3D planes");
    app.add_flag("--clouds", clouds, "3D clouds");

    CLI11_PARSE(app, argc, argv);
    polyscope::init();
    init();

    polyscope::state::userCallback = myCallBack;
    polyscope::show();
    return 0;
}
