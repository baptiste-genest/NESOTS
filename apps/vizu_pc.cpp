#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"
#ifdef __APPLE__
#include "Eigen/Dense"
#else
#include "eigen3/Eigen/Dense"
#endif

#include <fenv.h>

#include <iostream>
#include <fstream>

#include "geometrycentral/surface/surface_mesh.h"
#include "geometrycentral/surface/vertex_position_geometry.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

#include "../deps/CLI11.hpp"

polyscope::SurfaceMesh* INPUT;

using scalar = float;
using vec = Eigen::Vector<scalar,3>;
using points = std::vector<vec>;

void readMesh(std::string file) {
    std::unique_ptr<ManifoldSurfaceMesh> topology;
    std::unique_ptr<VertexPositionGeometry> geometry;
    std::tie(topology, geometry) = readManifoldSurfaceMesh(file);
    INPUT = polyscope::registerSurfaceMesh("input mesh",
                                           geometry->vertexPositions, topology->getFaceVertexList());
}

void readPts(const std::string& filename) {
    std::ifstream file(filename);
    points P;
    scalar x,y,z;
    while (file) {
        file >> x>>y>>z;
        P.push_back(vec(x,y,z));
    }
    polyscope::registerPointCloud("input point cloud",P);
}

int main(int argc,char** argv) {
    CLI::App app("Mesh and points vizualisation") ;

    std::string meshname;
    app.add_option("--input_mesh", meshname, "Input Mesh");

    std::string pcname;
    app.add_option("--input_points", pcname, "Input Points");

    polyscope::options::autoscaleStructures = false;
    polyscope::options::autocenterStructures = false;

    CLI11_PARSE(app, argc, argv);

    polyscope::init();
    polyscope::view::upDir = polyscope::view::UpDir::YUp;

    if (!meshname.empty())
        readMesh(meshname);
    if (!pcname.empty())
        readPts(pcname);

    polyscope::show();
    return 0;
}
