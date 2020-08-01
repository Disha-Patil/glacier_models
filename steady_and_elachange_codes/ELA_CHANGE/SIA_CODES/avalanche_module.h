#ifndef AVALANCHE_MODULE_H
#define AVALANCHE_MODULE_H

#include "globals.h"


void gen_ij_and_fracs(std::vector<std::vector<int>> &sorted_ij_, std::vector<std::vector<double>> &frac_);
void gen_catchment_and_c_frac(std::vector<std::vector<int>> &catchment_, std::vector<std::vector<double>> &c_frac_);
void avalanche(Eigen::MatrixXd& bed_, std::vector<std::vector<int>> &sorted_ij_, std::vector<std::vector<double>> &frac_, Eigen::MatrixXd& h, 
 Eigen::MatrixXd& av_height_,double avalanche_frac);
void spread(Eigen::MatrixXd& h, Eigen::MatrixXd& av_height_, std::vector<std::vector<int>> &catchment_, 
std::vector<std::vector<double>> &c_frac_, int layer_mask_[xmax][ymax], std::vector<std::vector<int>> &layer_array_, std::vector<std::array<int,2>> domain );
void gen_layer_mask_and_array(int layer_mask_[xmax][ymax], std::vector<std::vector<int>> &layer_array_);
void smoothen(Eigen::MatrixXd& h,std::vector<std::vector<int>> &catchment_, Eigen::MatrixXd& bed_, std::vector<std::array<int,2>> domain, Eigen::MatrixXi& m, Eigen::MatrixXd& av_height);
#endif
