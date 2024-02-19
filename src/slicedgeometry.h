#ifndef SLICEDGEOMETRY_H
#define SLICEDGEOMETRY_H

#include "utils.h"
#include "sampling.h"
#include "geometric_algorithms.h"
#include "ssw.h"

#include <iostream>
#include <set>

template<class pointType,class sliceType,class vectorType>
class SlicedGeometry
{
public:
    using point = pointType;
    using slice = sliceType;
    using vector = vectorType;
    using points = std::vector<point>;
    using vectors = std::vector<vector>;
    int N;
    std::vector<scalar> cost;

    SlicedGeometry() {
    }

    inline static points sub_sample(const points& X,int n){
        points sub_nu(n);
        for (int i = 0;i<n;i++)
            sub_nu[i] = X[rand()%X.size()];
        return sub_nu;
    }

    //    https://stackoverflow.com/questions/3722430/most-efficient-way-of-randomly-choosing-a-set-of-distinct-integers
    inline points sub_sample_distinct(const points& X,int n) const {
        auto D = randomDistinctNumbers(n,X.size());
        points sub(n);
        int i = 0;
        for (auto x : D)
            sub[i++] = X[x];
        return sub;
    }

    virtual points samples(int N) const = 0;

    virtual slice getSlice() const = 0;

    virtual vector projectOnTangentSpace(const point& x,const vector& v) const = 0;

    virtual vector parallelTransport(const point& x,const point& y,const vector& v) const = 0;

    virtual point advect(const point& x,const point& a,const point& b,const slice& s) const = 0;

    virtual point projectOnSlice(const point& p,const slice& s) const = 0;

    virtual points projectOnSlice(const points& p,const slice& s) const {
        points P = p;
        for (auto& x : P)
            x = projectOnSlice(x,s);
        return P;
    }

    virtual scalar curvilinearAbscissaAlongSlice(const point& x,const slice& s) const = 0;

    virtual scalar distanceAlongSlice(scalar x,scalar y) const  = 0;

    virtual scalar distance(const point& x,const point& y) const = 0;

    virtual vector Log(const point& x,const point& y) const = 0;
    virtual point Exp(const point& x,const vector& y,scalar t = 1) const = 0;

    virtual labeled<scalar> orderAlongSlice(const points& X,const slice& s) const {
        labeled<scalar> P(X.size());
        for (int i = 0;i<X.size();i++){
            P[i] = {i,curvilinearAbscissaAlongSlice(X[i],s)};
        }
        return P;
    }


    virtual scalar transport1D(const labeled<scalar> &a, const labeled<scalar> &b, int offset, int p) const {
        scalar cost = 0;
        auto n = a.size();
        auto m = b.size();
        /*
        auto k = b.size()/a.size();
        for (int i = 0;i<n;i++)
            for (int j = 0;j<k;j++)
                cost += std::pow(std::abs(distanceAlongSlice(a[i].second,b[(k*i + j +offset)%m].second)),p);
        */
        for (int i = 0;i<n;i++)
                cost += std::pow(std::abs(distanceAlongSlice(a[i].second,b[(i+offset)%m].second)),p);
        return cost/n;
    }

    virtual int computeOptimalCut(const labeled<scalar> &a, const labeled<scalar> &b, int p) const {
        return 0;
        // naïve version n^2 cost
        auto n = b.size();
        scalar best_cost = transport1D(a,b,0,p);
        int best = 0;
        for (int i = 1;i<n;i++){
            scalar d = transport1D(a,b,i,p);
            if (d < best_cost){
                best_cost = d;
                best = i;
            }
        }
        return best;
    }

    points projectAlongSlice(points X,const slice& s) const{
        for (auto& x : X)
            x = projectOnSlice(x,s);
        return X;
    }

    vectors computeBalencedSWGradients(const points &mu, const points &nu,const slice& gamma,scalar& cost) const
    {
        auto n = mu.size();

        auto Pi_MU = projectAlongSlice(mu,gamma);
        points Pi_NU = projectAlongSlice(nu,gamma);

        static const auto cmp = [](const std::pair<int,scalar>& a,const std::pair<int,scalar>& b){return a.second < b.second;};

        auto coords_mu = orderAlongSlice(Pi_MU,gamma);
        std::sort(coords_mu.begin(),coords_mu.end(),cmp);

        auto coords_nu = orderAlongSlice(Pi_NU,gamma);
        std::sort(coords_nu.begin(),coords_nu.end(),cmp);

        int cut = computeOptimalCut(coords_mu,coords_nu,1);

        cost = 0;
        vectors gradients(mu.size());
        for (auto&& [id,x] : const_zip(coords_mu)){
            cost += distanceAlongSlice(coords_mu[id].second,coords_nu[(id+cut)%n].second);
            gradients[x.first] = Log(mu[x.first],advect(mu[x.first],Pi_MU[x.first],Pi_NU[coords_nu[(id+cut)%n].first],gamma));
        }

        return gradients;
    }

    virtual scalar computeOTEnergy(const points &mu, const points &nu, int nb_batches)
    {
        auto n = mu.size();
        scalar cost = 0;
#pragma omp parallel for reduction(+:cost)
        for(auto i =0; i < nb_batches; ++i)
        {
            slice gamma = getSlice();
            auto Pi_MU = projectAlongSlice(mu,gamma);
            points Pi_NU = projectAlongSlice(sub_sample(nu,mu.size()),gamma);
            //points Pi_NU = projectAlongSlice(nu,gamma);

            static const auto cmp = [](const std::pair<int,scalar>& a,const std::pair<int,scalar>& b){return a.second < b.second;};

            auto coords_mu = orderAlongSlice(Pi_MU,gamma);
            std::sort(coords_mu.begin(),coords_mu.end(),cmp);

            auto coords_nu = orderAlongSlice(Pi_NU,gamma);
            std::sort(coords_nu.begin(),coords_nu.end(),cmp);

            int cut = computeOptimalCut(coords_mu,coords_nu,1);
            for (auto&& [id,x] : const_zip(coords_mu)){
                cost += distanceAlongSlice(coords_mu[id].second,coords_nu[(id+cut)%n].second);
            }
        }
        return cost/nb_batches/mu.size();
    }

    virtual std::pair<points,scalar> computeOTSlice(const points &mu, const points &nu, scalar epsilon, int nb_batches,bool useMedian = true,const std::function<scalar(const point& p)>& step_functionnal =
                                                                                                   [](const point& p) {return 1e8;})
    {
        std::vector<vectors> gradients(nb_batches);
        scalar SW=0;
#pragma omp parallel for reduction(+:SW)
        for (uint32_t i =0; i < nb_batches; ++i)
        {
            slice gamma = getSlice();
            scalar c = 0;
            gradients[i] = computeBalencedSWGradients(mu,sub_sample(nu,mu.size()),gamma,c);
            SW += c;
        }
        auto cost = scalar(SW)/nb_batches/mu.size();

        auto advected = mu;
        for (int i = 0;i<mu.size();i++){
            vectors stochastic_gradients(nb_batches);
            for (int k = 0;k<nb_batches;k++)
                stochastic_gradients[k] = gradients[k][i];
            vector G;
            if (useMedian)
                G = geoalgo::geometric_median(stochastic_gradients);
            else
                G = geoalgo::mean(stochastic_gradients);
            advected[i] = Exp(mu[i],cap_norm(G,step_functionnal(mu[i]))*epsilon);
        }
        return {advected,cost};
    }
};
#endif // SLICEDGEOMETRY_H
