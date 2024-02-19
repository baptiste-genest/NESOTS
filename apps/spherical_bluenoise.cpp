#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include "polyscope/point_cloud.h"

#include "../src/nesots.h"
#include <fenv.h>

#include <iostream>
#include <fstream>

#include "../deps/CLI11.hpp"


using Geometry = NSphericalGeometry;

int N = 2048;
int dim = 3;
std::unique_ptr<Geometry> G;
Geometry::points MU,NU;
int batch_size = 32;


polyscope::PointCloud* PC;
polyscope::PointCloud* MU_PC;

std::vector<float> SW;

bool toggle = false;
int id = 0;

std::string outfile = "/tmp/out.pts";

bool median = true;
bool nonunif= true;
float epsilon = 0.4;
scalar step = 0.999;
scalar stept = 0.9;

void iter() {
  auto rslt = G->computeOTSlice(MU,NU,epsilon,batch_size,median);
  MU = rslt.first;
  SW.push_back(std::log(rslt.second));
  id ++;
  stept *= step;
  epsilon = stept;

  // PRINT ENERGY
  //static std::ofstream file("/tmp/energy.data");
  //file <<G->computeOTEnergy(MU,NU,1000) <<std::endl;
}

void myCallBack() {
    if (ImGui::Button("toggle")) {
        toggle = toggle ^ true;
        if (!toggle)
            std::cout << id << std::endl;
    }
    if (toggle || ImGui::Button("Iteration")){
        iter();
        polyscope::registerPointCloud("MU",MU);
    }
    ImGui::Checkbox("use median descent",&median);
    if (ImGui::Button("reset")){
        MU = G->samples(N);
        stept = 0.9;
    }
    if (ImGui::Button("export")){
        PointCloudIO::write_point_cloud("/tmp/mu.pts",MU);
    }
    if (ImGui::Button("show nu"))
        polyscope::registerPointCloud("NU",NU);
    ImGui::PlotLines("SW", SW.data(), SW.size());
    ImGui::SliderFloat("epsilon",&epsilon,0,1);
}

void init_polyscope_data () {
    polyscope::registerPointCloud("MU",MU);
    Geometry::points center = {Geometry::vector::Zero(3)};
    polyscope::registerPointCloud("CENTER",center)->setPointRadius(1,false);
}


pcg32 GLOBALRNG;

int main(int argc,char** argv) {
    
    GLOBALRNG.seed(time(0));

    CLI::App app("Spherical Bluenoise Sampling") ;

    app.add_option("--dim", dim, "Dimensions of points (default 3)");

    std::string nu_file;
    app.add_option("--target_measure", nu_file, "Target measure to sample, must be on the sphere, if void then uniform");

    int target_size = 500'000;
    app.add_option("--target_size", target_size, "Number of samples to generate if uniform target density");

    app.add_option("--sample_size", N, "Number of samples to generate (default 2048)");

    app.add_option("--batch_size", batch_size, "Number of batches per slice (parallel)");

    app.add_option("--use_median", median, "Use geometric median for descent (better for non uniform target)");


    int nb_iter = 1000;
    app.add_option("--iter", nb_iter, "Number of iterations of NESOTS");

    app.add_option("--output_pts", outfile, "Output points file");

    bool viz = false;
    app.add_flag("--viz", viz, "interactive interface with polyscope viz (only for dim = 3)");

    CLI11_PARSE(app, argc, argv);

    G = std::make_unique<Geometry>(dim);

    if (nu_file.empty()){
        NU = G->samples(target_size);
        //you can modify NU here to sample non uniform densities
    }
    else
        NU = PointCloudIO::read_point_cloud(nu_file,dim);

    MU = Geometry::sub_sample(NU,N);

    if (viz) {
        polyscope::init();
        init_polyscope_data();

        polyscope::state::userCallback = myCallBack;
        polyscope::show();
    }
    else {
        for (int i = 0; i < nb_iter; i++){
            if (i % 20 == 0)
                std::cout << "[computing NESOT] " << i << "/" << nb_iter << std::endl;
            iter();
        }
        std::cout << "[done] exported in " << outfile << std::endl;
        PointCloudIO::write_point_cloud(outfile,MU);
    }
    return 0;
}
