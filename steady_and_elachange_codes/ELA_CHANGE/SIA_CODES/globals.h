#ifndef GLOBALS_H
#define GLOBALS_H

#define EIGEN_NO_DEBUG              // turn off assertions in eigen library to save time
#include <math.h>                   // the always needed math lib, stdio and stdlib
#include <stdio.h>
#include <stdlib.h>
#include <iostream>                 // for standard input output scheme in c++
#include <Eigen/Sparse>             // include the eigen sparse header files
#include <vector>                   // vector template class of c++ to ease the vector operations ...needed with eigen library
#include <fstream>                  // handle input output of files
#include <omp.h>                    // needed if you are using openMP routines
#include <chrono>
#include <Eigen/Core>               // eigen core headers with datastructures 

#include <Eigen/IterativeLinearSolvers>                    // eigen headers for iterative methods like CG and BiCG   
#include <unsupported/Eigen/IterativeSolvers>              // methods like housholders algo https://eigen.tuxfamily.org/dox/unsupported/group__IterativeSolvers__Module.html
#include <Python.h>
#include <string.h>
#define inv_2_dx_dx  (1/(double)(2*DX*DX))	          // constants as required in the discretized equation 
#define inv_2_dx_sq  (1/(double)(2*2*DX*DX))              // constants as required in the discretized equation

//#define DT_inv_2_dx_dx  (DT*inv_2_dx_dx)                  // constants as required in the discretized equation
#define C_inv_2_dx_sq  (C*inv_2_dx_sq)                    // constants as required in the discretized equation

using Eigen::VectorXi;                  // define eigen vector for integer data             
using Eigen::VectorXd;                  // define eigen vector for doubles data
using Eigen::SparseMatrix;              // define eigen sparse matrix
using Eigen::ConjugateGradient;         // define eigen iterative solver to be used, i.e., CG

using Eigen::MatrixXi;                  // define eigen matrix for integer data
using Eigen::MatrixXd;                  // define eigen matrix for double data
using std::cout;                        // to print output
using std::endl;                        // to print on new line

const Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "];"); //Stores a set of parameters controlling the way matrices are printed

typedef SparseMatrix<double> SpMat;       //Define a sparse matrix which will be used later as A matrix

#define PI 3.14159265

#define xmax 46
#define ymax 49
#define n0 388

#define grid_size xmax*ymax

#define h_init 0
#define DX  100.0
#define DY  100.0
#define tmax 1000

#define path1 "./first/"

#define C 0.0000111
#define ela_rate 0.005
extern double ela_change; 
#define a0 1.0

#define mbgrad 0.006

#define layers 6

#endif
