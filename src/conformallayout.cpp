#include "conformallayout.h"

#include "sampling.h"
vecs sampleLayout(const HyperbolicLayout& M,int sampleNum,const scalars& face_weights){
    vecs sampleList;

    int triNum = M.triangles.size();
    double currentArea;

    srand(time(NULL));

    //generate a vector of cumulative areas to remove size bias from face selection
    scalars cumulative(triNum, 0.0);
    cumulative[0]=face_weights[0];

    //calculate each triangle's area and add to cumulative
    for (int i = 1;i<face_weights.size();i++) {
        cumulative[i] = cumulative[i - 1] + face_weights[i];
    }

    for (int j = 0; j < sampleNum; j++) {
        double sample;
        int faceidx;

        //random face selection
        sample = (cumulative.back() - 0) * ((double)rand() / (double)RAND_MAX) + 0;
        for (faceidx = 0; cumulative[faceidx] < sample && faceidx < cumulative.size(); faceidx++);

        //random point generation within previously selected face area
        vecs V;
        for (auto v : M.triangles[faceidx]){
            V.push_back(M.pos[v]);
        }

        scalar alpha, beta;
        alpha = (1 - 0) * ((double)rand() / (double)RAND_MAX) + 0;
        beta = (1 - 0) * ((double)rand() / (double)RAND_MAX) + 0;

        scalar a, b, c;
        a = 1 - sqrt(beta);
        b = (sqrt(beta)) * (1 - alpha);
        c = sqrt(beta) * alpha;

        //resulting sample
        vec P = a*V[0] + b*V[1] + c*V[2];
        sampleList.push_back(P);
    }

    return sampleList;
}

MeshSamples labeledSampleLayout(const HyperbolicLayout& M,int sampleNum,const scalars& face_weights){
    MeshSamples MS;

    MS.samples.resize(sampleNum);

    int triNum = M.triangles.size();
    double currentArea;

    srand(time(NULL));

    //generate a vector of cumulative areas to remove size bias from face selection
    scalars cumulative(triNum, 0.0);
    cumulative[0]=face_weights[0];

    //calculate each triangle's area and add to cumulative
    for (int i = 1;i<face_weights.size();i++) {
        cumulative[i] = cumulative[i - 1] + face_weights[i];
    }

    for (int j = 0; j < sampleNum; j++) {
        double sample;
        int faceidx;

        //random face selection
        sample = (cumulative.back() - 0) * ((double)rand() / (double)RAND_MAX) + 0;
        for (faceidx = 0; cumulative[faceidx] < sample && faceidx < cumulative.size(); faceidx++);

        //random point generation within previously selected face area
        vecs V;
        for (auto v : M.triangles[faceidx]){
            V.push_back(M.pos[v]);
        }

        scalar alpha, beta;
        alpha = (1 - 0) * ((double)rand() / (double)RAND_MAX) + 0;
        beta = (1 - 0) * ((double)rand() / (double)RAND_MAX) + 0;

        scalar a, b, c;
        a = 1 - sqrt(beta);
        b = (sqrt(beta)) * (1 - alpha);
        c = sqrt(beta) * alpha;

        //resulting sample
        vec P = a*V[0] + b*V[1] + c*V[2];
        MS.samples[j] = P;
        if (!MS.faceToSamples.contains(faceidx))
            MS.faceToSamples[faceidx] = {};
        MS.faceToSamples[faceidx].insert(j);
    }

    return MS;
}


vec HyperbolicLayout::getThirdPoint(const EmbeddedHalfedge &h, bool &skipped) const {
    const auto& a = h.tail;
    const auto& b = h.tip;
    vec na = vec(-a(0),-a(1),a(2));
    vec nb = vec(-b(0),-b(1),b(2));
    vec n = na.cross(nb);
    vec nap = n.cross(na);
    vec nbp = n.cross(nb);
    auto da = std::cosh(L.at(h.h.next().next().edge()));
    auto db = std::cosh(L.at(h.h.next().edge()));
    vec p0 = db*nap/nb.dot(nap) + da*nbp/na.dot(nbp);
    auto c0 = 1+G.dotL(p0,p0);
    auto c1 = 2*G.dotL(n,p0);
    auto c2 = G.dotL(n,n);
    skipped = false;
    if (c1*c1 < 4*c2*c0){
        skipped = true;
    }
    auto t = (-c1 + sqrt(std::max(0.,c1*c1 - 4*c2*c0)))/offset(2*c2);
    return p0 + t*n;
}

std::vector<int> HyperbolicLayout::build(const Mesh &M, int vstart, scalar max_z) {
    auto v0 = M.topology->vertex(vstart);
    if (v0.isBoundary())
        return ints();
    vec x0 = G.param(0,0);
    Halfedge h0 = *v0.outgoingHalfedges().begin();
    auto l0 = L[h0.edge()];
    auto c = std::cosh(l0);
    vec x1 = vec(sqrt(c*c-1),0,c);
    std::queue<EmbeddedHalfedge> queue;
    queue.push(EmbeddedHalfedge(h0,x0,x1,0));
    std::set<Face> done;

    int visited = 0;

    std::vector<int> visit_order(M.topology->nFaces(),-1);


    while (!queue.empty()){
        auto gh = queue.front();queue.pop();
        auto r = gh.round+1;
        auto f = gh.h.face();
        if (done.contains(f))// || (max_radius > 0 && r > max_radius))
            continue;
        visited++;
        bool skipped;
        vec x = getThirdPoint(gh,skipped);
        if (skipped || (max_z > 0 && x(2) > max_z))
            continue;
        done.insert(f);
        visit_order[f.getIndex()] = done.size();
        std::map<Vertex,vec> ordermap;
        ordermap[gh.h.tailVertex()] = gh.tail;
        ordermap[gh.h.tipVertex()] = gh.tip;
        ordermap[gh.h.next().tipVertex()] = x;
        for (auto v : f.adjacentVertices())
            pos.push_back(ordermap[v]);
        int id = pos.size();
        triangles.push_back({id-3,id-2,id-1});
        faceMap.push_back(f);
        queue.push(EmbeddedHalfedge(gh.h.next().next().twin(),gh.tail,x,r));
        queue.push(EmbeddedHalfedge(gh.h.twin(),gh.tip,gh.tail,r));
        queue.push(EmbeddedHalfedge(gh.h.next().twin(),x,gh.tip,r));
    }
    return visit_order;
}

vec HyperbolicLayout::hyperbolicBarycentricCoordinates(const vec &p, const triangle& T) const
{
    auto l0 = G.distance(p,pos[T[0]]);
    auto l1 = G.distance(p,pos[T[1]]);
    auto l2 = G.distance(p,pos[T[2]]);
    auto l01 = G.distance(pos[T[0]],pos[T[1]]);
    auto l02 = G.distance(pos[T[0]],pos[T[2]]);
    auto l12 = G.distance(pos[T[1]],pos[T[2]]);
    vec B;
    auto A = area(l01,l02,l12);
    B(0) = area(l1,l2,l12)/A;
    B(1) = area(l0,l2,l02)/A;
    B(2) = 1. - B(0) - B(1);
    return B;
}
