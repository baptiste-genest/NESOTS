#ifndef SAMPLER_H
#define SAMPLER_H
#include <vector>
#ifdef __APPLE__
#include "Eigen/Dense"
#else
#include "eigen3/Eigen/Dense"
#endif

template<int dim>
class Sampler_
{
public:
    using vec = Eigen::Vector<double,dim>;
    using Samples = std::vector<vec>;
    Sampler_() {}
    virtual Samples getSamples(size_t n);
    virtual vec operator()() = 0;
    inline Samples operator()(size_t n) {return getSamples(n);}
};

#endif // SAMPLER_H
