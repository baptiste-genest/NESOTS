#ifndef POINTCLOUDIO_H
#define POINTCLOUDIO_H

#ifdef __APPLE__
#include "Eigen/Dense"
#else
#include "eigen3/Eigen/Dense"
#endif
#include <iostream>
#include <fstream>
#include <vector>

namespace PointCloudIO
{

using vec = Eigen::Vector3d;
using vecs = std::vector<vec>;
using Vec = Eigen::VectorXd;
using Vecs = std::vector<Vec>;

void write_point_cloud(std::string filename,const vecs& X);
void write_point_cloud(std::string filename,const Vecs& X);

vecs read_3D_point_cloud(std::string filename);
Vecs read_point_cloud(std::string filename,int d);

}

#endif // POINTCLOUDIO_H
