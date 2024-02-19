#ifndef GEOMETRICALGORITHMS_H
#define GEOMETRICALGORITHMS_H

#include "utils.h"

namespace geoalgo {

template<class T>
T mean(const std::vector<T>& X){
    using vector = T;
    auto n = X.size();
    vector y = X[0];
    for (int i = 1;i<n;i++)
        y += X[i];
    y /= scalar(n);
    return y;
}


template<class T>
inline T zero() {return T::Zero();}

template<>
inline complex zero<complex>() {return complex(0.,0.);}

inline scalar distance(const complex& a,const complex& b) {return std::abs(a-b);}

template<int N>
inline scalar distance(const Eigen::Vector<scalar,N>& a,const Eigen::Vector<scalar,N>& b) {return (a-b).norm();}

//https://en.wikipedia.org/wiki/Geometric_median
template<class T>
T geometric_median(const std::vector<T>& X){
    using vector = T;
    auto n = X.size();
    auto y = mean(X);

    scalar eps = 1e-7;
    scalar jump = 2*eps;
    while (jump > eps){
        auto d = 1e-7 + distance(y,X[0]);
        scalar w = 1./d;
        vector next = X[0]/d;
        for (int i = 1;i<n;i++){
            d = 1e-7 + distance(y,X[i]);
            w += 1./d;
            next += X[i]/d;
        }
        next /= w;
        jump = distance(y,next);
        y = next;
    }
    return y;
}

}

#endif // GEOMETRICALGORITHMS_H
