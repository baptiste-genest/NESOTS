#ifndef UTILS_H
#define UTILS_H
#ifdef __APPLE__
#include "Eigen/Dense"
#include "Eigen/Sparse"
#else
#include "eigen3/Eigen/Dense"
#include "eigen3/Eigen/Sparse"
#endif

#include <vector>
#include <algorithm>
#include <execution>



using scalar = double;
using vec = Eigen::Vector<scalar,3>;
using vec4 = Eigen::Vector<scalar,4>;
using vec2 = Eigen::Vector<scalar,2>;
using mat = Eigen::Matrix<scalar,3,3>;
using mat2 = Eigen::Matrix<scalar,2,2>;
using Mat = Eigen::Matrix<double,-1,-1>;

template<class T = scalar>
using twins = std::pair<T,T>;

using complex = std::complex<scalar>;
using complexs = std::vector<complex>;
static const complex I = complex(0,1);

using vecs = std::vector<vec>;
using scalars = std::vector<scalar>;
using ints = std::vector<int>;

using Vec = Eigen::Vector<scalar,-1>;
using smat = Eigen::SparseMatrix<scalar>;
using smatd = Eigen::SparseMatrix<double>;

template<class T>
using labeled = std::vector<std::pair<int,T>>;

template<class T>
struct const_zip {

    const std::vector<T>& ref;
    const_zip(const std::vector<T>& r) : ref(r) {}

    struct const_zip_iterator {
        const std::vector<T>& ref;
        int id;

        const_zip_iterator(const std::vector<T>& ref,int id) : ref(ref),id(id){}

        bool operator!=(const const_zip_iterator& o) {return id!=o.id;}
        std::pair<int,const T&> operator*(){return {id,ref[id]};};
        void operator++(){id++;}
    };

    const_zip_iterator begin() {
        return const_zip_iterator(ref,0);
    }
    const_zip_iterator end() {
        return const_zip_iterator(ref,ref.size());
    }

};

template<class T>
struct zip {

    std::vector<T>& ref;
    zip(std::vector<T>& r) : ref(r) {}

    struct zip_iterator {
        std::vector<T>& ref;
        int id;

        zip_iterator(std::vector<T>& ref,int id) : ref(ref),id(id){}

        bool operator!=(const zip_iterator& o) {return id!=o.id;}
        std::pair<int,T&> operator*(){return {id,ref[id]};};
        void operator++(){id++;}
    };

    zip_iterator begin() {
        return zip_iterator(ref,0);
    }
    zip_iterator end() {
        return zip_iterator(ref,ref.size());
    }

};

template<class A,class B>
std::vector<A> projA(const std::vector<std::pair<A,B>>& X){
    std::vector<A> P(X.size());
    for (auto&& [id,p] : const_zip(X))
        P[id] = p.first;
    return P;
}
template<class A,class B>
std::vector<B> projB(const std::vector<std::pair<A,B>>& X){
    std::vector<B> P(X.size());
    for (auto&& [id,p] : const_zip(X))
        P[id] = p.second;
    return P;
}

inline int mod(int a,int b) {
    while (a < 0)
        a +=b;
    while (a >= b)
        a -= b;
    return a;
}
#include <functional>

template<class T1,class T2>
std::vector<T2> apply(const std::function<T2(const T1&)>& f,const std::vector<T1>& X){
    std::vector<T2> rslt(X.size());
    for (int i = 0;i<X.size();i++)
        rslt[i] = f(X[i]);
    return rslt;
}

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

inline vec proj_b_along_a(const vec& a,const vec& b){
    return a*a.dot(b)/a.squaredNorm();
}


inline mat make_ortho_proj(const vec& x) {
    return mat::Identity() - x*x.transpose()/x.squaredNorm();
}

inline vec ortho_proj_b_against_a(const vec& a,const vec& b) {
    return b - a*a.dot(b)/a.squaredNorm();
}

inline bool isNan(const vec& x){
    if (std::isnan(x(0)) ||std::isnan(x(1)) || std::isnan(x(2))){
        assert(false);
        return true;
    }
    return false;
}

inline complex cap_norm(const complex& x,scalar l) {
    auto n = std::abs(x);
    if (n < l)
        return x;
    return x*l/n;
}

template<class T>
inline T cap_norm(const T& x,scalar l) {
    auto n = x.norm();
    if (n < l)
        return x;
    return x.normalized()*l;
}

inline scalar offset(scalar x) {
    return sgn(x)*1e-10 + x;
}


#endif // UTILS_H
