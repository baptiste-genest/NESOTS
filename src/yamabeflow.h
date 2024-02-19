#ifndef YAMABEFLOW_H
#define YAMABEFLOW_H

#include "utils.h"
#include "Mesh.h"

using DiscreteMetric = std::map<Edge,scalar>;

struct YamabeFlow {

    scalar previous_norm = 1000;
    void init();

    Vec goal;

    void computeGoalAngles();

    void updateEdgeLength();

    bool violateTriangleInequality(const Halfedge& h) const;

    bool violateTriangleInequality(const Face& f) const;

    scalar getAngle(const Vertex& v,const Halfedge& h) const;

    std::array<Vertex,3> getVertices(Face f) const;

    Vec compute_angle_deflect();

    scalar faceAngleDeflect(Halfedge h);

    static scalar cotan(scalar x);

    scalar computeWij(const Edge& e);

    smat buildHessian();

    bool done = false;

    void flow(bool newton = true,scalar tol = 1e-6,scalar stuck_tol = 1e-10);

    Mesh mesh;
    DiscreteMetric length,original_length;
    Vec u;
};

#endif // YAMABEFLOW_H
