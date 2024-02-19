#ifndef SAMPLING_H
#define SAMPLING_H

#include "utils.h"
#include <iostream>
#include <random>
#include <set>
#include "Mesh.h"
#include "pcg32.h"

//Random generator init to draw random line directions

extern pcg32 GLOBALRNG;

inline vec uniform_sphere__() {
    static std::mt19937 gen;
    static std::normal_distribution<scalar> dist{0.0,1.0};
    vec dir;
    dir(0) = dist(gen);
    dir(1) = dist(gen);
    dir(2) = dist(gen);
    return dir.normalized();
}

// Function to convert a random number in [0, 1) to a normal distribution
inline double randomToNormal(double mean, double stddev)
{
  double u1 = GLOBALRNG.nextDouble();
  double u2 = GLOBALRNG.nextDouble();
  double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
  return mean + stddev * z0;
}

inline vec uniform_sphere() {
  vec dir;
  dir(0) = randomToNormal(0.0,1.0);
  dir(1) = randomToNormal(0.0,1.0);
  dir(2) = randomToNormal(0.0,1.0);
  return dir.normalized();
}

inline Vec uniform_sphere(int d) {
  Vec dir(d);
  for (int i = 0;i<d;i++)
    dir(i) = randomToNormal(0.0,1.0);
  return dir.normalized();
}


inline scalar uniform_01__() {
    static std::mt19937 gen;
    static std::uniform_real_distribution<scalar> dist{0.0,1.0};
    return std::clamp(dist(gen),0.,1.);
}

inline scalar uniform_01() {
  return GLOBALRNG.nextDouble();
}

inline scalars uniform_01(int N) {
    scalars U(N);
    for (auto& x : U)
        x = uniform_01();
    return U;
}


inline scalar uniform_tau() {
    return uniform_01()*2*M_PI;
}

inline scalar uniform() {
    return uniform_01()*2-1;
}


inline mat uniform_proj() {
    return make_ortho_proj(uniform_sphere());
}


inline std::vector<vec> FibonacciSphere(int n)
{
    static scalar goldenRatio = (1 + std::sqrt(5.))/2.;
    std::vector<vec> FS(n);
    for (int i = 0;i<n;i++){
        scalar theta = 2 * M_PI * i / goldenRatio;
        scalar phi = std::acos(1 - 2*(i+0.5)/n);
        FS[i] = vec(cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi)).normalized();
    }
    return FS;
}

inline std::vector<vec2> grid(int w,int h,scalar W = 1,scalar H = 1) {
    std::vector<vec2> G(w*h);
    int k = 0;
    for (int i = 0;i<w;i++)
        for (int j = 0;j<h;j++){
            auto x = vec2(scalar(i)/(w-1),scalar(j)/(h-1));
            G[k] = (x - vec2::Ones()*0.5)*2;
            G[k](0) *= W;
            G[k](1) *= H;
            k++;
        }
    return G;
}
/*
inline vecs samples(int N) {
    vecs S(N);
    for (auto& x : S)
        x = uniform_sphere();
    return S;
}
*/

struct MeshSamples {
    std::vector<vec> samples;
    std::map<int,std::set<int>> faceToSamples;
};
struct StaticMeshSamples {
    std::vector<vec> samples;
    using idList = std::vector<int>;
    std::map<int,idList> faceToSamples;
};


vecs sampleMesh(const Mesh& M,int sampleNum,const scalars& face_weights);
MeshSamples labeledSampleMesh(const Mesh&,int sampleNum,const scalars& face_weights);
StaticMeshSamples labeledStaticSampleMesh(const Mesh& M,int sampleNum,const scalars& face_weights);

struct HyperbolicLayout;
vecs sampleLayout(const HyperbolicLayout& M,int sampleNum,const scalars& face_weights);
MeshSamples labeledSampleLayout(const HyperbolicLayout& M,int sampleNum,const scalars& face_weights);

inline scalar gaussian(scalar x,scalar sigma){
    return std::exp(-x*x/sigma/2);
}

inline std::set<int> randomDistinctNumbers(int n,int max) {
    std::set<int> res;
    int count = max + 1;
    for (int i = count - n; i < count; i++) {
        auto item = rand()%(i+1);
        if (res.contains(item))
            res.insert(i);
        else
            res.insert(item);
    }
    if (res.size() != n)
        std::cerr << "nope" << std::endl;
    return res;
}


#endif // SAMPLING_H
