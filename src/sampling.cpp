#include "sampling.h"


vecs sampleMesh(const Mesh &M, int sampleNum, const scalars &face_weights) {
    vecs sampleList(sampleNum);

    const auto& pos = M.geometry->vertexPositions;

    int triNum = M.topology->nFaces();
    double currentArea;

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
        sample = (cumulative.back() - 0) * GLOBALRNG.nextDouble() + 0;
        for (faceidx = 0; cumulative[faceidx] < sample && faceidx < cumulative.size(); faceidx++);

        //random point generation within previously selected face area
        vecs V;
        for (auto v : M.topology->face(faceidx).adjacentVertices()){
            V.push_back(toVec(pos[v]));
        }

        scalar alpha, beta;
        alpha = (1 - 0) * GLOBALRNG.nextDouble() + 0;
        beta = (1 - 0) * GLOBALRNG.nextDouble() + 0;

        scalar a, b, c;
        a = 1 - sqrt(beta);
        b = (sqrt(beta)) * (1 - alpha);
        c = sqrt(beta) * alpha;

        //resulting sample
        vec P = a*V[0] + b*V[1] + c*V[2];
        sampleList[j] = P;
    }

    return sampleList;
}

MeshSamples labeledSampleMesh(const Mesh& M,int sampleNum,const scalars& face_weights){
    MeshSamples MS;

    MS.samples.resize(sampleNum);

    const auto& pos = M.geometry->vertexPositions;

    int triNum = M.topology->nFaces();
    double currentArea;


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
        sample = (cumulative.back() - 0) * GLOBALRNG.nextDouble() + 0;
        for (faceidx = 0; cumulative[faceidx] < sample && faceidx < cumulative.size(); faceidx++);

        //random point generation within previously selected face area
        vecs V;
        for (auto v : M.topology->face(faceidx).adjacentVertices()){
            V.push_back(toVec(pos[v]));
        }

        scalar alpha, beta;
        alpha = (1 - 0) * GLOBALRNG.nextDouble() + 0;
        beta = (1 - 0) * GLOBALRNG.nextDouble() + 0;

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

StaticMeshSamples labeledStaticSampleMesh(const Mesh &M, int sampleNum, const scalars &face_weights)
{
    StaticMeshSamples MS;

    MS.samples.resize(sampleNum);

    const auto& pos = M.geometry->vertexPositions;

    int triNum = M.topology->nFaces();
    double currentArea;


    //generate a vector of cumulative areas to remove size bias from face selection
    scalars cumulative(triNum, 0.0);
    cumulative[0]=face_weights[0];

    //calculate each triangle's area and add to cumulative
    for (int i = 1;i<face_weights.size();i++) {
        cumulative[i] = cumulative[i - 1] + face_weights[i];
    }

    ints sample_per_face(M.topology->nFaces(),0);

    for (int j = 0; j < sampleNum; j++) {
        double sample;
        int faceidx;

        //random face selection
        sample = (cumulative.back() - 0) * GLOBALRNG.nextDouble() + 0;
        for (faceidx = 0; cumulative[faceidx] < sample && faceidx < cumulative.size(); faceidx++);
        sample_per_face[faceidx]++;
    }
    int k = 0;//sample id
    for (int faceidx = 0;faceidx < M.topology->nFaces();faceidx++){
        //random point generation within previously selected face area
        vec V[3];
        int i = 0;
        for (auto v : M.topology->face(faceidx).adjacentVertices()){
            V[i++] = toVec(pos[v]);
        }
        MS.faceToSamples[faceidx] = ints(sample_per_face[faceidx],0);
        for (int i = 0;i<sample_per_face[faceidx];i++){
            scalar alpha, beta;
            alpha = (1 - 0) * GLOBALRNG.nextDouble() + 0;
            beta = (1 - 0) * GLOBALRNG.nextDouble() + 0;

            scalar a, b, c;
            a = 1 - sqrt(beta);
            b = (sqrt(beta)) * (1 - alpha);
            c = sqrt(beta) * alpha;

            //resulting sample
            vec P = a*V[0] + b*V[1] + c*V[2];
            MS.samples[k] = P;
            MS.faceToSamples[faceidx][i] = k;
            k++;
        }
    }

    return MS;
}
