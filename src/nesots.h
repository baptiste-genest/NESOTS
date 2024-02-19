#ifndef NESOTS_H
#define NESOTS_H

#ifdef __APPLE__
#include "Eigen/Dense"
#else
#include "eigen3/Eigen/Dense"
#endif
#include "sphericalgeometry.h"
#include "n_spherical_geometry.h"
#include "hyperbolicgeometry.h"
#include "nd_hyperbolic_geometry.h"
#include "euclideangeometry.h"
#include "point_cloud_io.h"
#include "bvh_wrapper.h"
#include "yamabeflow.h"
#include "conformallayout.h"
#include <chrono>

using Time = std::chrono::high_resolution_clock;
using TimeStamp = std::chrono::time_point<std::chrono::high_resolution_clock>;
using TimeTypeSec = float;
using DurationSec = std::chrono::duration<TimeTypeSec,std::milli>;

inline TimeTypeSec TimeFrom(const TimeStamp& A){
    return DurationSec(Time::now()-A).count();
}

#endif // NESOTS_H
