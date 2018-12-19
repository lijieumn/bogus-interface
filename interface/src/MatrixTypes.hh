#ifndef ARGUS_MATRIX_TYPES_HH
#define ARGUS_MATRIX_TYPES_HH

#include <Eigen/Core>

namespace argus {

typedef Eigen::VectorXd DynVec ;                              // Resizable vector
typedef Eigen::Matrix<double, 3, Eigen::Dynamic >  DynMat3  ; // 3xn matrix
typedef Eigen::Matrix<   int, 3, Eigen::Dynamic >  DynMat3i ; // 3xn matrix of indices
typedef Eigen::Matrix3d Mat3 ;                                // 3x3 matrix



} //argus

#endif
