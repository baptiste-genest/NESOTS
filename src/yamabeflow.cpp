#include "yamabeflow.h"

void YamabeFlow::init() {
    mesh.geometry->requireEdgeLengths();
    for (auto e : mesh.topology->edges())
        original_length[e] = mesh.geometry->edgeLength(e);
    //u = Vec::Zero(mesh.topology->nVertices());
    u = Vec::Zero(mesh.topology->nVertices());
    updateEdgeLength();
    for (const auto& f : mesh.topology->faces()){
        if (violateTriangleInequality(f))
            std::cout << f << " violates TI" << std::endl;
    }
    computeGoalAngles();
    auto K = compute_angle_deflect();
    previous_norm = K.lpNorm<Eigen::Infinity>();
}

void YamabeFlow::computeGoalAngles() {
    goal = Vec::Ones(mesh.topology->nVertices())*2*M_PI;
    int i = 0;
    for (const auto& bl : mesh.topology->boundaryLoops()){
        for (const auto& v : bl.adjacentVertices()){
            i++;
            goal(v.getIndex()) = M_PI;
        }
    }
    std::cout << "nb of boundary vertices " << i << std::endl;
}

void YamabeFlow::updateEdgeLength() {
    for (auto e : mesh.topology->edges()){
        auto i = e.firstVertex().getIndex();
        auto j = e.secondVertex().getIndex();
        length[e] = 2*std::asinh(std::exp(0.5*(u(i)+u(j)))*
                     std::sinh(original_length[e]*0.5));
    }
}

bool YamabeFlow::violateTriangleInequality(const Halfedge &h) const {
    return length.at(h.next().next().edge()) + length.at(h.next().edge()) < length.at(h.edge()) - 1e-6;
}

bool YamabeFlow::violateTriangleInequality(const Face &f) const {
    for (const auto& e : f.adjacentHalfedges())
        if (violateTriangleInequality(e)){
            //std::cout << "TIV" << std::endl;
            return true;
        }
    return false;
}

scalar YamabeFlow::getAngle(const Vertex &v, const Halfedge &h) const {
    auto e = h.edge();
    if (violateTriangleInequality(h.face())){
        if (violateTriangleInequality(h))
            return M_PI;
        return 0;
    }
    assert(v != e.firstVertex());
    assert(v != e.secondVertex());
    auto j = mesh.topology->connectingEdge(v,e.firstVertex());
    auto k = mesh.topology->connectingEdge(v,e.secondVertex());
    const auto& li = length.at(e);
    const auto& lj = length.at(j);
    const auto& lk = length.at(k);
    auto x = (std::cosh(lj)*std::cosh(lk)-std::cosh(li))/offset(std::sinh(lj)*std::sinh(lk));
    return std::acos(std::clamp(x,-1.,1.));
}

std::array<Vertex, 3> YamabeFlow::getVertices(Face f) const {
    std::array<Vertex,3> V;
    auto i = 0;
    for (auto v : f.adjacentVertices())
        V[i++] = v;
    return V;
}

Vec YamabeFlow::compute_angle_deflect() {
    auto& T = *mesh.topology;
    Vec K = goal;
    for (const auto& v : T.vertices())
        for (const auto& h : v.outgoingHalfedges()){
            if (!h.isInterior())
                continue;
            K(v.getIndex()) -= getAngle(v,h.next());
        }
    return K;
}

scalar YamabeFlow::faceAngleDeflect(Halfedge h) {
    scalar A = M_PI;
    for (int i = 0;i<3;i++){
        A += getAngle(h.next().tipVertex(),h)*(i==0 ? 1 : -1);
        h = h.next();
    }
    return A;
}

scalar YamabeFlow::cotan(scalar x) {
    return std::clamp(cos(offset(x))/sin(offset(x)),-10e8,10e8);
}

scalar YamabeFlow::computeWij(const Edge &e) {
    scalar wij = 0;
    for (const auto& h : e.adjacentHalfedges()){
        if (!h.isInterior())
            continue;
        if (!violateTriangleInequality(h.face()))
            wij += cotan(0.5*faceAngleDeflect(h));
    }
    return 0.5*wij;
}

smat YamabeFlow::buildHessian() {
    auto& T = *mesh.topology;
    auto N  = T.nVertices();
    smat H(N,N);
    H.reserve(Eigen::VectorXi::Constant(N,1+2*10));
    for (size_t i = 0;i<N;i++)
        H.coeffRef(i,i) = 0;
    for (const auto& e: T.edges()) {
        auto i = e.firstVertex().getIndex();
        auto j = e.secondVertex().getIndex();
        auto w = computeWij(e);
        auto t = std::tanh(length[e]*0.5);
        auto diag = 0.5*w*(1 + t*t);
        auto cross = 0.5*w*(-1 + t*t);
        H.coeffRef(i,j) = cross;
        H.coeffRef(j,i) = cross;
        H.coeffRef(i,i) += diag;
        H.coeffRef(j,j) += diag;
    }
    H.makeCompressed();
    //std::cout << H << std::endl;
    return H;
}

void YamabeFlow::flow(bool newton, scalar tol, scalar stuck_tol) {
    if (done)
        return;

    auto K = compute_angle_deflect();
    //previous_norm = K.norm();
    scalar step = 10;
    Vec d;
    if (newton){
        smat H = buildHessian();
        Eigen::SimplicialLDLT<smat> S(H);
        if (S.info() != Eigen::Success){
            d = K;
            //std::cerr << "non invertible system" << std::endl;
        }
        else
            d = S.solve(K);
    }
    else
        d = K;

    auto old_u = u;


    u -= step*d;
    updateEdgeLength();

    auto Kn = compute_angle_deflect();
    scalar gn = Kn.lpNorm<Eigen::Infinity>();

    while((gn > previous_norm) || (std::isnan(gn))) {
        step *= 0.5;

        u = old_u - step*d;
        updateEdgeLength();

        Kn = compute_angle_deflect();
        gn = Kn.lpNorm<Eigen::Infinity>();
    }
    std::cout << "[Yamabe-Flow Backtrack] " << previous_norm << " -> " << gn << std::endl;
    if (std::abs(previous_norm-gn) < stuck_tol){
        std::cout << "stuck" << std::endl;
        done = true;
    }
    previous_norm = gn;
    if (previous_norm < tol){
        std::cout << "done" << std::endl;
        done = true;
    }
}
