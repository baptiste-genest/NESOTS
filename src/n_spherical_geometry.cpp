#include "n_spherical_geometry.h"

NSphericalGeometry::point NSphericalGeometry::sampleSphere() const {
    return uniform_sphere(d);
}

NSphericalGeometry::point NSphericalGeometry::proj_b_along_a(const point &a, const point &b)
{
    return a*a.dot(b)/a.squaredNorm();
}

NSphericalGeometry::point NSphericalGeometry::orthoproj_b_against_a(const point &a, const point &b){
    return b - proj_b_along_a(a,b);
}


NSphericalGeometry::points NSphericalGeometry::samples(int N) const
{
    points X(N);
    for (int j = 0;j<N;j++)
        X[j] = sampleSphere();
    return X;
}

NSphericalGeometry::slice NSphericalGeometry::getSlice() const
{
    twins<Vec> B;
    B.first = sampleSphere();
    B.second = sampleSphere();
    while (std::abs(B.first.dot(B.second)) > 0.999)
        B.second = sampleSphere();
    B.second = orthoproj_b_against_a(B.first,B.second).normalized();
    return B;
}

NSphericalGeometry::vector NSphericalGeometry::projectOnTangentSpace(const point &x, const vector &v) const
{
    return orthoproj_b_against_a(x,v);
}

NSphericalGeometry::vector NSphericalGeometry::parallelTransport(const point &x, const point &y, const vector &v) const
{
    return v;
}

NSphericalGeometry::point NSphericalGeometry::advect(const point &x, const point &a, const point &b, const slice &) const
{
    if ((b-a).norm() < 1e-6)
        return x;
    const auto& e1 = a;
    const auto& e2 = orthoproj_b_against_a(e1,b).normalized();
    auto xO = orthoproj_b_against_a(e1,x);
    xO = orthoproj_b_against_a(e2,xO);
    if (!isOrthogonal(xO,e1) || !isOrthogonal(xO,e2)){
        return x;
        std::cerr <<"[wrong basis]" <<std::endl;
    }
    vector xr = x - xO;
    auto th = distance(a,b);
    auto c = std::cos(th);
    auto s = std::sin(th);
    auto ca = xr.dot(e1);
    auto cb = xr.dot(e2);
    return xO + e1*(c*ca - s*cb) + e2*(s*ca + c*cb);
}

NSphericalGeometry::point NSphericalGeometry::projectOnSlice(const point &p, const slice &s) const
{
    vector x = s.first*s.first.dot(p) + s.second*s.second.dot(p);
    return x.normalized();
}

scalar NSphericalGeometry::curvilinearAbscissaAlongSlice(const point &x, const slice &s) const
{
    return (M_PI + std::atan2(x.dot(s.second),x.dot(s.first)))/(2*M_PI);
}

scalar NSphericalGeometry::distanceAlongSlice(scalar x, scalar y) const
{
    auto d = std::abs(x-y);
    return std::min(d,1-d);
}

scalar NSphericalGeometry::distance(const point &x, const point &y) const
{
    return std::acos(std::clamp(x.dot(y),-1.,1.));
}

NSphericalGeometry::vector NSphericalGeometry::Log(const point &x, const point &y) const
{
    return projectOnTangentSpace(x,y-x).normalized()*distance(x,y);
}

NSphericalGeometry::point NSphericalGeometry::Exp(const point &p, const vector &v, scalar t) const
{
    if (std::abs(p.dot(v)) > 1e-6){
        std::cerr << "v is not in the tangent plane of p" << std::endl;
        std::cerr <<"[p] " <<  p.transpose() << std::endl;
        std::cerr <<"[p norm] " <<  p.norm() << std::endl;
        std::cerr <<"[v] " << v.transpose() << std::endl;
    }
    t *= v.norm();
    return p*cos(t) + sin(t)*v.normalized();
}

std::pair<NSphericalGeometry::points, scalar> NSphericalGeometry::computeOTSliceDirection(const points &mu0, const points &nu, scalar epsilon, int nb_batches, bool useMedian, const std::function<scalar (const point &)> &step_functionnal)
{
    int n = mu0.size();
    points mu(2*n);
    for (int i = 0;i<n;i++){
        mu[i] = mu0[i];
        mu[i+n] = -mu0[i];
    }
    std::vector<vectors> gradients(nb_batches);
    std::vector<int> batches(nb_batches);

    scalar SW=0;
#pragma omp parallel for reduction(+:SW)
    for (int32_t i = 0; i < nb_batches; ++i)
    {
        slice gamma = getSlice();
        scalar c = 0;
        gradients[i] = computeBalencedSWGradients(mu,sub_sample(nu,mu.size()),gamma,c);
        SW += c;
    }
    auto cost = scalar(SW)/nb_batches/mu.size();

    auto advected = mu0;
    for (int i = 0;i<n;i++){
        vectors stochastic_gradients(2*nb_batches);
        for (int k = 0;k<nb_batches;k++){
            stochastic_gradients[k] = gradients[k][i];
            stochastic_gradients[k+nb_batches] = -gradients[k][i+n];
        }
        vector G;
        if (useMedian)
            G = geoalgo::geometric_median(stochastic_gradients);
        else
            G = geoalgo::mean(stochastic_gradients);
        advected[i] = Exp(mu[i],cap_norm(G,step_functionnal(mu[i]))*epsilon);
    }
    return {advected,cost};
}

scalar NSphericalGeometry::computeOTEnergytoUniform(const points &mu, int nb_batches) const
{
    static const auto cmp = [](const std::pair<int,scalar>& a,const std::pair<int,scalar>& b){return a.second < b.second;};
    auto n = mu.size();

    scalar cost = 0;
#pragma omp parallel for reduction(+:cost)
    for (int i = 0;i<nb_batches;i++){
        auto gamma = getSlice();
        auto Pi_MU = projectAlongSlice(mu,gamma);
        points Pi_NU(mu.size());

        scalar t = M_PI*2./n;
        for (int j = 0;j<n;j++){
            Pi_NU[j] = cos(j*t)*gamma.first + sin(j*t)*gamma.second;
        }

        auto coords_mu = orderAlongSlice(Pi_MU,gamma);
        std::sort(coords_mu.begin(),coords_mu.end(),cmp);

        auto coords_nu = orderAlongSlice(Pi_NU,gamma);
        std::sort(coords_nu.begin(),coords_nu.end(),cmp);

        auto&& cut = computeOptimalCut(coords_mu,coords_nu,1);
        for (auto&& [id,x] : const_zip(coords_mu)){
            cost += distanceAlongSlice(coords_mu[id].second,coords_nu[(id+cut)%n].second);
        }
    }
    return cost/nb_batches/n;
}
