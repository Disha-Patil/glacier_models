#ifndef CALC_MATRICES_H
#define CALC_MATRICES_H

#include "globals.h"

void fill_d(MatrixXd& d, MatrixXd& h, MatrixXd& bed, MatrixXi& m,VectorXd& neigh1,VectorXd& neigh2,VectorXd& neigh3,VectorXd& neigh4);
void fill_mb(MatrixXd& mb, MatrixXd& bed, MatrixXd& ela, MatrixXd& debris, MatrixXd& h, std::vector<std::array<int,2>> domain);
void fill_b(VectorXd& b,VectorXd& flux,  MatrixXd& d, MatrixXd& mb, MatrixXi& m,MatrixXd& bed, MatrixXd& h,double DT,double DT_inv_2_dx_dx, std::vector<std::array<int,2>> domain);
void fill_A(SpMat& A, MatrixXd& d,MatrixXi& m,MatrixXi& m_idx, double DT_inv_2_dx_dx, std::vector<std::array<int,2>> domain);
void calc_areaf(MatrixXd& h,MatrixXi& ice_mask, double& areaf, double& area);
#endif
