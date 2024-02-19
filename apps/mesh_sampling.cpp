#include "geometrycentral/surface/remeshing.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#ifdef __APPLE__
#include "Eigen/Dense"
#else
#include "eigen3/Eigen/Dense"
#endif
#include "../src/utils.h"
#include "../src/sphericalgeometry.h"
#include "../src/sampling.h"
#include "../src/ssw.h"
#include "../src/cylindricalgeometry.h"
#include "../src/hyperbolicgeometry.h"
#include "../src/euclideangeometry.h"
#include "../src/fluid_sim_sdot.h"
#include "../src/point_cloud_io.h"
#include "../src/bvh_wrapper.h"
#include "../src/halton.hpp"

 #include <fenv.h>

#include <iostream>
#include <fstream>

#include "../deps/CLI11.hpp"

BVH_WRAPPER::Vec3 toVec3(const vec& x){
    return BVH_WRAPPER::Vec3(x(0),x(1),x(2));
}


using Geometry = SphericalGeometry;

int N = 2048;
Geometry G;
Geometry::points MU,NU;
int batch_size = 64;
float epsilon = 0.4;

BVH_WRAPPER::BVH BVH;

Mesh input_mesh;
Mesh spherical_mesh;

std::string meshname;
std::string spherical_meshname;
std::string meshname_other = "";

bool other_geometry = false;
Mesh other_geometry_mesh;

bool viz = false;

void create_BVH() {
    std::vector<BVH_WRAPPER::Tri> tris;
    const auto& G = spherical_mesh.geometry->vertexPositions;
    int id = 0;
    for (auto T : spherical_mesh.topology->getFaceVertexList()){
        BVH_WRAPPER::Tri tri;
        for (int i = 0;i<3;i++)
        {
            auto x = G[T[i]];
            tri(i) = BVH_WRAPPER::Vec3(x.x,x.y,x.z);
        }
        tris.push_back(tri);
    }
    BVH = BVH_WRAPPER::BVH(tris);
}

void loadMesh() {
    std::tie(input_mesh.topology, input_mesh.geometry) = readManifoldSurfaceMesh(meshname);
    if (other_geometry)
        std::tie(other_geometry_mesh.topology, other_geometry_mesh.geometry) = readManifoldSurfaceMesh(meshname_other);
    //input_mesh.normalize();
    //other_geometry_mesh.normalize();
    if (viz){
        polyscope::registerSurfaceMesh("input_mesh",
                                               input_mesh.geometry->vertexPositions, input_mesh.topology->getFaceVertexList());
        if (other_geometry)
        polyscope::registerSurfaceMesh("other mesh",
                                               other_geometry_mesh.geometry->vertexPositions, other_geometry_mesh.topology->getFaceVertexList());
    }
    std::tie(spherical_mesh.topology, spherical_mesh.geometry) = readManifoldSurfaceMesh(spherical_meshname);

    return;
}

vecs mapToMesh(const vecs& X,const Mesh& M){
    vecs V;
    auto F = spherical_mesh.topology->getFaceIndices();
    std::vector<int> faceid;
    int i = 0;
    for (const auto& v : X){
        auto I = BVH.get_intersection(toVec3(vec::Zero()),toVec3(v));
        if (!I.valid)
            continue;
        auto f = spherical_mesh.topology->face(I.id);
        if (f.getIndex() == geometrycentral::INVALID_IND)
            continue;
        auto b= spherical_mesh.Barycentric2(v,f);
        V.push_back(M.posFromWeights(b,M.topology->face(f.getIndex())));
    }
    return V;
}

bool toggle = false;
int id = 0;
scalar step = 0.99;
scalar stept = 0.7;

void iter() {
    MU = G.computeOTSlice(MU,NU,epsilon,batch_size,true,[](const vec& v){
              auto I = BVH.get_intersection(toVec3(vec::Zero()),toVec3(v));
              auto f = spherical_mesh.topology->face(I.id);
              scalar sa = std::max(1e-7,spherical_mesh.diameter(f));
              return sa;
          }).first;
    stept *= step;
    epsilon = stept;
    auto map1 = mapToMesh(MU,input_mesh);
    if (viz){
        polyscope::registerPointCloud("map input",map1);
        if (other_geometry)
            polyscope::registerPointCloud("map other geom",mapToMesh(MU,other_geometry_mesh));
    }
    id ++;
}

void myCallBack() {
    if (ImGui::Button("toggle")) {
        toggle = toggle ^ true;
        if (!toggle)
            std::cout << id << std::endl;
    }
    if (ImGui::Button("show spherical mesh")) {
        polyscope::registerSurfaceMesh("spherical mesh",spherical_mesh.geometry->vertexPositions, spherical_mesh.topology->getFaceVertexList());
    }
    if (toggle || ImGui::Button("Iteration")){
        iter();
    }
    ImGui::SliderFloat("epsilon",&epsilon,0,1);
}

void init () {
    loadMesh();
    create_BVH();
    auto FA = input_mesh.faceAreas();
    NU = sampleMesh(spherical_mesh,100000,FA);
    MU = G.sub_sample_distinct(NU,N);
}


pcg32 GLOBALRNG;

int main(int argc,char** argv) {
    CLI::App app("Spherical Mesh Sampling") ;

    app.add_option("--input_mesh", meshname, "Input Mesh")->required();
    app.add_option("--sphere_mesh", spherical_meshname, "Spherical version of input_mesh, must have same topology")->required();

    app.add_option("--sample_size", N, "Number of samples to generate");

    app.add_option("--batch_size", batch_size, "Number of batches per slice (parallel)");

    int nb_iter = 1000;
    app.add_option("--iter", nb_iter, "Number of iterations of NESOT");

    std::string outfile = "/tmp/out.pts";
    app.add_option("--output_pts", outfile, "Output points file");

    app.add_flag("--viz", viz, "interactive interface with polyscope viz");

    
    app.add_option("--other_geometry_mesh", meshname_other, "maps points also to this mesh, must have same topology");

    std::string outfile_OG = "/tmp/out_other_geometry.pts";
    app.add_option("--output_other_geometry_pts", outfile_OG, "Output points file on other geometry");

    CLI11_PARSE(app, argc, argv);

    if (!meshname_other.empty())
        other_geometry = true;

    if (viz){
        polyscope::init();
        init();
        polyscope::view::upDir = polyscope::view::UpDir::ZUp;

        polyscope::state::userCallback = myCallBack;
        polyscope::show();
    }
    else {
        init();
        for (int i = 0;i<nb_iter;i++){
            if (i % 20 == 0)
                std::cout << "[ Running NESOTS ] " << i << "/" << nb_iter << std::endl;
            MU = G.computeOTSlice(MU,NU,stept,batch_size,true,[](const vec& v){
                      auto I = BVH.get_intersection(toVec3(vec::Zero()),toVec3(v));
                      auto f = spherical_mesh.topology->face(I.id);
                      scalar sa = std::max(1e-7,spherical_mesh.diameter(f));
                      return sa;
                  }).first;
            stept *= step;
        }
        PointCloudIO::write_point_cloud(outfile,mapToMesh(MU,input_mesh));
        std::cout << "[done] exported in " << outfile << std::endl;
        if (other_geometry){
            PointCloudIO::write_point_cloud(outfile_OG,mapToMesh(MU,other_geometry_mesh));
            std::cout << "[done] exported in " << outfile_OG << std::endl;
        }
    }
    return 0;
}
