#ifndef HALTONSAMPLER_H
#define HALTONSAMPLER_H
#include <string>
#include "sampler.h"
double *halton ( int i, int m );
double *halton_base ( int i, int m, int b[] );
int halton_inverse ( double r[], int m );
double *halton_sequence ( int i1, int i2, int m );
int i4vec_sum ( int n, int a[] );
int prime ( int n );
double r8_mod ( double x, double y );
void r8mat_print ( int m, int n, double a[], std::string title );
void r8mat_print_some ( int m, int n, double a[], int ilo, int jlo, int ihi,int jhi, std::string title );
void timestamp ( );

template<int dim>
class HaltonSampler : public Sampler_<dim>
{
private:
    int count = 0;
public:
    HaltonSampler() {}

    // Sampler interface
public:
    using vec = Eigen::Vector<double,dim>;
    using Samples = std::vector<vec>;
    Samples getSamples(size_t n) override {
        double* v = halton_sequence(count,count+n,dim);
        count += n;
        Samples S(n,vec::Zero(dim));
        for (int i = 0;i<n;i++)
            for (int j = 0;j<dim;j++)
                S[i](j) = (v[dim*i+j]-0.5)*2;
        delete[] v;
        return S;
    }

    vec operator ()() override {
        double* v = halton(++count,dim);
        vec S(dim);
        for (int i = 0;i<dim;i++)
            S[i] = v[i];
        delete[] v;
        return S;
    }
};


        #endif
