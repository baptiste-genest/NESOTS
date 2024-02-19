#include "point_cloud_io.h"

void PointCloudIO::write_point_cloud(std::string filename, const vecs &X)
{
    std::ofstream file(filename);
    for (const auto& x : X)
        file << x.transpose() << std::endl;
    file.close();
}

void PointCloudIO::write_point_cloud(std::string filename, const Vecs &X)
{
    std::ofstream file(filename);
    for (const auto& x : X)
        file << x.transpose() << std::endl;
    file.close();
}


PointCloudIO::vecs PointCloudIO::read_3D_point_cloud(std::string filename)
{
    vecs X;
    std::ifstream file(filename);
    double x,y,z;
    while (file >> x >> y >> z)
        X.push_back(vec(x,y,z));
    return X;
}

PointCloudIO::Vecs PointCloudIO::read_point_cloud(std::string filename,int d)
{
    Vecs X;
    std::ifstream file(filename);
    Vec input(d);
    while (file){
        for (int i = 0;i<d;i++)
            file >> input(i);
        X.push_back(input);
    }
    return X;
}
