#ifndef FLUID_SIM_SDOT_H
#define FLUID_SIM_SDOT_H

#include "sphericalgeometry.h"

template<class GeometryType>
class FluidSimSDOT
{
public:

    using point = typename GeometryType::point;
    using points = typename GeometryType::points;
    using slice = typename GeometryType::slice;
    using vector = typename GeometryType::vector;
    using vectors = typename GeometryType::vectors;

    scalar c0;

    FluidSimSDOT(){}
    FluidSimSDOT(const points& init,const points& target_measure,int n,const scalar &epsilon, const scalar &tau,const scalars& rho,const vector& gravity) :
        eps2(epsilon*epsilon),
        tau(tau),
        n(n),
        rho(rho),
        M(init),
        nu(target_measure),
        gravity(gravity)
    {
        V.resize(n,geoalgo::zero<vector>());
        std::tie(M,c0) = geometry.computeOTSlice(M,nu,0.2,100);
    }

    points timestep() {
        auto A = M;
        scalar c = c0*2;
        for (int i = 0;i<3;i++){
            std::tie(A,c) = geometry.computeOTSlice(A,nu,1,16);
        }
        for (int i = 0;i<n;i++){
            V[i] += tau * (
                        geometry.Log(M[i],A[i])/eps2 +
                        geometry.projectOnTangentSpace(M[i],gravity)*rho[i]);
            auto old = M[i];
            M[i] = geometry.Exp(M[i],V[i]*tau);
            V[i] = geometry.projectOnTangentSpace(M[i],V[i]);
            //V[i] = geometry.parallelTransport(old,M[i],V[i]);
        }
        return A;
    }

    const points& getPos() const {
        return M;
    }

    points M,nu;

private:
    GeometryType geometry;
    scalar eps2,tau;
    int n;
    scalars rho;
    vectors V;
    vector gravity;


};

#endif // FLUID_SIM_SDOT_H
