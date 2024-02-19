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
#include <queue>
#include <cmath>
#include <time.h>
#include "../deps/CLI11.hpp"

using Geometry = NDHyperbolicGeometry;

polyscope::SurfaceMesh* INPUT;

YamabeFlow Y;
polyscope::SurfaceMesh* Layout;
BVH_WRAPPER::BVH BVH;

using points = Geometry::points;

bool nonunif=false;
bool viz = false;

BVH_WRAPPER::Vec3 toKlein(const vec& z) {
    return BVH_WRAPPER::Vec3(z(0),z(1),z(2));
    static Geometry G(3);
    auto x = G.sheetToKlein(z);
    return BVH_WRAPPER::Vec3(5*x(0),5*x(1),5);
}

scalar klein_z = 5;
BVH_WRAPPER::BVH create_BVH(const HyperbolicLayout& L) {
    std::vector<BVH_WRAPPER::Tri> tris;
    int id = 0;
    for (auto T : L.triangles){
        BVH_WRAPPER::Tri tri;
        for (int i = 0;i<3;i++){
            tri(i) = toKlein(L.pos[T[i]]);
        }
        tris.push_back(tri);
    }
    return BVH_WRAPPER::BVH(tris);
}



std::string meshname;
std::string meshname_other = "";

bool other_geometry = false;
Mesh other_geometry_mesh;



void loadMesh(bool ps = true) {
    std::tie(Y.mesh.topology, Y.mesh.geometry) = readManifoldSurfaceMesh(meshname);
    Y.mesh.normalize();
    if (ps)
        INPUT = polyscope::registerSurfaceMesh("input_mesh",
                                               Y.mesh.geometry->vertexPositions, Y.mesh.topology->getFaceVertexList());
    if (other_geometry)
        std::tie(other_geometry_mesh.topology, other_geometry_mesh.geometry) = readManifoldSurfaceMesh(meshname_other);
}

int vid = 0;
Geometry G(3);

BVH_WRAPPER::Vec3 toVec3(const vec& x){
    return BVH_WRAPPER::Vec3(x(0),x(1),x(2));
}

MeshSamples mu;
StaticMeshSamples nu;
HyperbolicLayout layout;

ints visit_counts;
ints visit_queue;

int iter_per_patch = 10;
int batch_size = 32;
float max_z = 1.5;

float epsilon = 0.2;
scalar step = 0.999;
scalar stept = 0.4;

vec meshToLayout(const Mesh& M,const vec& x,const Face& f,int lfid) {
    auto b = M.Barycentric2(x,f);
    return layout.posFromBarycentric(b,layout.triangles[lfid]);
}

using int2s = std::vector<std::pair<int,int>>;


points subSampleThenMap(const StaticMeshSamples& S,int n,const int2s& F,std::discrete_distribution<int>& d) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    points samples(n);
    for (int i = 0;i<n;i++){
        auto f = F[d(gen)];
        const auto& sf = S.faceToSamples.at(f.second);
        samples[i] = meshToLayout(Y.mesh,S.samples[sf[rand()%sf.size()]],Y.mesh.topology->face(f.second),f.first);
    }
    return samples;
}

points computeOTSlice(const Geometry::points &mu, const StaticMeshSamples &nu,const int2s& faceRestriction, scalar epsilon, int nb_batches)
{
    std::vector<Geometry::vectors> gradients(nb_batches);
    std::vector<int> batches(nb_batches);
    std::iota(batches.begin(),batches.end(),0);

    scalars weights(faceRestriction.size());
    for (int i = 0;i<faceRestriction.size();i++)
        weights[i] = nu.faceToSamples.at(faceRestriction[i].second).size();
    std::discrete_distribution<int> d(weights.begin(),weights.end());
    points sub_nu[nb_batches];
    for(auto k=0; k < nb_batches; ++k)
        sub_nu[k] = subSampleThenMap(nu,mu.size(),faceRestriction,d);

#pragma omp parallel for
    for(auto k=0; k < nb_batches; ++k)
    {
        vec gamma = G.getSlice();
        scalar c = 0;
        gradients[k] = G.computeBalencedSWGradients(mu,sub_nu[k],gamma,c);
    }

    auto advected = mu;
#pragma omp parallel for
    for (int i = 0;i<mu.size();i++){
        Geometry::vectors stochastic_gradients(nb_batches);
        for (int k = 0;k<nb_batches;k++)
            stochastic_gradients[k] = gradients[k][i];
        auto tangent_descent = geoalgo::geometric_median(stochastic_gradients);
        //auto tangent_descent = geoalgo::mean(stochastic_gradients);
        advected[i] = G.Exp(mu[i],tangent_descent*epsilon);
    }
    return advected;
}

bool exportmu = false;

size_t points_per_patch = 0;
bool mu_on_layout = false;
void LocalOptim(int v,int iter,MeshSamples& mu,const StaticMeshSamples& nu,const Mesh& M) {
    layout.clear();
    layout.build(M,v,max_z);
    //if the given vertex germ gives an empty layer we won't consider it as a germ in the future anymore
    if (layout.triangles.empty()){
        visit_counts[v] = 1000000;
        std::sort(visit_queue.begin(),visit_queue.end(),[](int a,int b){
            return visit_counts[a] < visit_counts[b];
        });
        return;
    }

    points sub_mu,sub_nu,original_pos;
    ints sub_mu_id,original_face;

    //auto old_mu = mu;
    int2s layout_faces;

    std::set<int> visitedVertex;
    for (int lfid = 0;lfid < layout.faceMap.size();lfid++){
        auto f = layout.faceMap[lfid];
        for (auto v : f.adjacentVertices())
            visitedVertex.insert(v.getIndex());
        if (mu.faceToSamples.contains(f.getIndex()))
            for (const auto& s : mu.faceToSamples[f.getIndex()]){
                sub_mu.push_back(meshToLayout(M,mu.samples[s],f,lfid));
                original_pos.push_back(mu.samples[s]);
                sub_mu_id.push_back(s);
                original_face.push_back(f.getIndex());
            }
        mu.faceToSamples[f.getIndex()].clear();
        if (nu.faceToSamples.contains(f.getIndex()))
            layout_faces.push_back({lfid,f.getIndex()});
    }
    points_per_patch += sub_mu.size();

    if (exportmu){
        vecs tmp(sub_mu.size()),tmp2(sub_mu.size());
        for (int i = 0;i<sub_mu.size();i++){
            tmp[i] = sub_mu[i];
            tmp2[i] = mu.samples[sub_mu_id[i]];
        }
        PointCloudIO::write_point_cloud("/tmp/mu_layout.pts",tmp);
        PointCloudIO::write_point_cloud("/tmp/mu.pts",Y.mesh.MapOnNotNormalized(tmp2));
        std::cout << "export" << std::endl;
        exportmu = false;
        return;
    }

    if (mu_on_layout)
        polyscope::registerPointCloud("submu",sub_mu);

    for (auto v : visitedVertex)
        visit_counts[v]++;


    for (int i = 0;i<iter;i++)
        sub_mu = computeOTSlice(sub_mu,nu,layout_faces,epsilon,batch_size);

    auto bvh = create_BVH(layout);

    for (int i = 0;i<sub_mu.size();i++){
        //auto p = toVec3(sub_mu[i]);
        auto I = bvh.get_intersection(toVec3(vec::Zero()),toKlein(sub_mu[i]));
        if (!I.valid){
            mu.faceToSamples[original_face[i]].insert(sub_mu_id[i]);
            continue;
        }
        vec b = layout.BarycentricCoordinates(sub_mu[i],layout.triangles[I.id]);
        bool ok = true;
        for (int i = 0;i<2;i++)
            if (b(i) < 0 || b(i) > 1)  {
                ok = false;
                break;
            }
        if (!ok){
            mu.faceToSamples[original_face[i]].insert(sub_mu_id[i]);
            continue;
        }
        auto fid = layout.faceMap[I.id];
        mu.samples[sub_mu_id[i]] = M.posFromWeights(b,fid);
        mu.faceToSamples[fid.getIndex()].insert(sub_mu_id[i]);
    }
    std::sort(visit_queue.begin(),visit_queue.end(),[](int a,int b){
        return visit_counts[a] < visit_counts[b];
    });
}

int N = 4096;
float layout_max_error = 1e-5;

Vec VertexDensity() {
    Y.mesh.geometry->requireVertexMeanCurvatures();
    Y.mesh.geometry->requireVertexDualAreas();
    Vec H(Y.mesh.topology->nVertices());
    for (auto v : Y.mesh.topology->vertices())
        H[v.getIndex()] = std::abs(Y.mesh.geometry->vertexMeanCurvature(v));
    {
    std::ofstream file("/tmp/density_onlyH.data");
    file << H << std::endl;file.close();
    }
    for (auto v : Y.mesh.topology->vertices())
        H[v.getIndex()] += Y.mesh.geometry->vertexDualArea(v);
    {
    std::ofstream file("/tmp/density.data");
    file << H << std::endl;file.close();
    }
    return H;
}

void init () {
    loadMesh();

    Y.init();
    visit_counts.resize(Y.mesh.topology->nVertices(),0);
    visit_queue.resize(Y.mesh.topology->nVertices(),0);
    std::iota(visit_queue.begin(),visit_queue.end(),0);
    while (!Y.done)
        Y.flow(true,layout_max_error);
    layout = HyperbolicLayout(Y.length);
    auto face_density = Y.mesh.faceAreas();
    INPUT->addFaceScalarQuantity("density",face_density);
    auto D = VertexDensity();
    INPUT->addVertexScalarQuantity("density",D);
    nu = labeledStaticSampleMesh(Y.mesh,500'000,face_density);
    mu = labeledSampleMesh(Y.mesh,N,face_density);
}

bool toggle = false;
bool layoutmode = false;

points mapMuBetweenMesh(const Mesh& a,const Mesh& b) {
    points mapped;
    for (auto & F : mu.faceToSamples){
        auto f = a.topology->face(F.first);
        for (auto v : F.second){
            auto bary = a.Barycentric2(mu.samples[v],f);
            mapped.push_back(b.posFromWeights(bary,b.topology->face(f.getIndex())));
        }
    }
    return mapped;
}

int record_OT_energy_slices = 0;
TimeTypeSec cumulated_OT_time = 0;
void global_OT_energy() {
    static bool init = false;
    static HyperbolicLayout L(Y.length);
    static std::ofstream file("/tmp/energy.data");
    static points global_nu;
    static auto meshToLayout = [] (const vec& x,const Face& f,int lfid) {
        auto b = Y.mesh.Barycentric2(x,f);
        return L.posFromBarycentric(b,L.triangles[lfid]);
    };
    if (!init){
        L.build(Y.mesh);
        for (int lfid = 0;lfid < L.faceMap.size();lfid++){
            auto f = L.faceMap[lfid];
            if (nu.faceToSamples.contains(f.getIndex()))
                for (const auto& s : nu.faceToSamples.at(f.getIndex()))
                    global_nu.push_back(meshToLayout(nu.samples[s],f,lfid));
        }
        init = true;
    }

    points sub_mu;
    for (int lfid = 0;lfid < L.faceMap.size();lfid++){
        auto f = L.faceMap[lfid];
        if (mu.faceToSamples.contains(f.getIndex()))
            for (const auto& s : mu.faceToSamples[f.getIndex()]){
                sub_mu.push_back(meshToLayout(mu.samples[s],f,lfid));
            }
    }
    file << cumulated_OT_time << " " << G.computeOTEnergy(sub_mu,global_nu,record_OT_energy_slices) << std::endl;
}


int nb_iter = 0;
void iter(){
    auto t1 = Time::now();
    LocalOptim(visit_queue[0],iter_per_patch,mu,nu,Y.mesh);
    cumulated_OT_time += TimeFrom(t1);
    if (record_OT_energy_slices)
        global_OT_energy();
    stept *= step;
    epsilon = stept;
    nb_iter++;

    if (viz) {
        polyscope::registerPointCloud("map mu",mu.samples);
        Layout = polyscope::registerSurfaceMesh("layout",layout.pos, layout.triangles);
        scalars onlayout(Y.mesh.topology->nFaces(),0);
        for (auto& f : layout.faceMap)
            onlayout[f.getIndex()] = 1;
        std::ofstream file("/tmp/onlayout.dat");
        for(auto i=0; i < onlayout.size(); ++i)
          file<<onlayout[i]<<std::endl;
        file.close();
        INPUT->addFaceScalarQuantity("on layout",onlayout);
        
        if (other_geometry) {
            polyscope::registerSurfaceMesh("other geometry",other_geometry_mesh.geometry->vertexPositions,other_geometry_mesh.topology->getFaceVertexList());
            polyscope::registerPointCloud("map mu other",mapMuBetweenMesh(Y.mesh,other_geometry_mesh));
        }
    }
}

void export_layout() {
    std::ofstream file("/tmp/layout.obj");
    for (auto& x : layout.pos)
        file << "v " << x.transpose() << std::endl;
    for (auto& t : layout.triangles){
        file << "f ";
        for (int i = 0;i<3;i++)
            file << t[i]+1 << " ";
        file << std::endl;
    }
    file.close();
}




void myCallBack() {
    if (ImGui::Button("iteration")){
        iter();
    }
    ImGui::SliderInt("iter per patch",&iter_per_patch,1,100);
    ImGui::SliderFloat("epsilon",&epsilon,0.1,1);
    if (ImGui::Button("toggle")){
        toggle ^= true;
        if (!toggle)
            std::cout << nb_iter << std::endl;
    }
    ImGui::SliderFloat("max_z",&max_z,1.01,1.5);
    if (toggle){
        iter();
    }
    if (ImGui::Button("show nu")){
        polyscope::registerPointCloud("nu",nu.samples);
    }
    if (ImGui::Button("export")){
        //exportmu = true;
        //iter();
        //export_layout();
        PointCloudIO::write_point_cloud("/tmp/mu.pts",Y.mesh.MapOnNotNormalized(mu.samples));
    }
}

pcg32 GLOBALRNG;

int main(int argc,char** argv) {
    CLI::App app("Local Hyperbolic Mesh Sampling") ;

    app.add_option("--input_mesh", meshname, "Input Mesh")->required();

    app.add_option("--sample_size", N, "Number of samples to generate");

    app.add_option("--max_z_layout", max_z, "Max z coordinate for layout construction");

    app.add_option("--yamabe_max_error", layout_max_error, "Convergence threshold of Yamabe flow");

    app.add_option("--batch_size", batch_size, "(L) Number of batches per slice (parallel)");

    app.add_option("--iter_per_patch", iter_per_patch, "(K) nb of NESOT iter per layout patch");

    app.add_option("--other_geometry_mesh", meshname_other, "maps points also to this mesh, must have same topology");

    app.add_option("--record_OT_energy_slices", record_OT_energy_slices, " if >0, will print OT energy evolution at each iteration, much slower, in /tmp/energy.data");

    app.add_flag("--nonunif", nonunif, "hardcoded non unif sampling");
    app.add_flag("--mu_on_layout", mu_on_layout, "display sub mu on layout");

    size_t seed=0;
    app.add_option("--seed", seed, "RNG seed");

    int nb_iter = 500;
    app.add_option("--iter", nb_iter, "Number of iterations of NESOT");

    std::string outfile = "/tmp/out.pts";
    app.add_option("--output_pts", outfile, "Output points file");
    std::string outOtherFile = "/tmp/out_other.pts";
    app.add_option("--output_other_geometry_pts", outOtherFile, "Output points file on the other_geometry_mesh");

    app.add_flag("--viz", viz, "interactive interface with polyscope viz");

    CLI11_PARSE(app, argc, argv);

    if (seed != 0)
      GLOBALRNG.seed(seed);
    else
      GLOBALRNG.seed(time(0));


  
    if (viz) {

        if (!meshname_other.empty()){
            other_geometry = true;
        }

        polyscope::init();
        init();
        polyscope::view::upDir = polyscope::view::UpDir::YUp;
        polyscope::state::userCallback = myCallBack;
        polyscope::show();
    }
    else {
        if (!meshname_other.empty()){
            other_geometry = true;
        }
        loadMesh(false);
        Y.init();
        visit_counts.resize(Y.mesh.topology->nVertices(),0);
        visit_queue.resize(Y.mesh.topology->nVertices(),0);
        std::iota(visit_queue.begin(),visit_queue.end(),0);
        auto t_yamabe = Time::now();
        while (!Y.done)
            Y.flow(true,layout_max_error);
        std::cout << "[ Yamabe flow timing ] " << TimeFrom(t_yamabe) << "ms" << std::endl;
        layout = HyperbolicLayout(Y.length);
        auto FA = Y.mesh.faceAreas();
        
        
        nu = labeledStaticSampleMesh(Y.mesh,500'000,FA);
        mu = labeledSampleMesh(Y.mesh,N,FA);

        auto t_iters = Time::now();
        for (int i = 0;i<nb_iter;i++){
            if (i % 20 == 0)
                std::cout << "[ Running NESOTS ] " << i << "/" << nb_iter << std::endl;
            iter();
        }
        std::cout << "[ NESOTS timing ] " << TimeFrom(t_iters) << "ms" << std::endl;

        std::ofstream ppp("/tmp/ppp.data",std::fstream::app);
        ppp << scalar(points_per_patch)/nb_iter << std::endl;

        PointCloudIO::write_point_cloud(outfile,Y.mesh.MapOnNotNormalized(mu.samples));
        std::cout << "[done] exported in " << outfile << std::endl;

        if (other_geometry)
        {
            auto mubis = mapMuBetweenMesh(Y.mesh,other_geometry_mesh);
            PointCloudIO::write_point_cloud(outOtherFile,mubis);
            std::cout << "[done] exported in " << outOtherFile << std::endl;
        }
    }


    return 0;
}

/*
 *
        //Non unif hack
        if (nonunif)
        {
            std::ofstream out("/tmp/vert-nonunif.dat");
            vec orig = {-5,0.5,-1.1};
            auto i=0;
            double rescale;
            for(auto f: Y.mesh.topology->faces())
            {
                auto p = Y.mesh.vertexInFace(f);
                rescale = 1.0/( std::pow((p[0] - orig).norm(), 3.0));
                FA[i] *= rescale;
                ++i;
            }
            i=0;
            for(auto v: Y.mesh.topology->vertices())
            {
                auto p = Y.mesh.geometry->vertexPositions[v];
                vec pp = {p.x,p.y,p.z};
                rescale = 1.0/( std::pow((pp - orig).norm(), 3.0));
                out << rescale <<std::endl;
                ++i;
            }
            out.close();
        }
 */
