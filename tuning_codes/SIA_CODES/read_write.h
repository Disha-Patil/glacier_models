#ifndef READ_WRITE_H
#define READ_WRITE_H

#include "globals.h"

void read_bed_and_mask(Eigen::MatrixXd& bed, Eigen::MatrixXi& mask, Eigen::MatrixXd& bed0);
void read_mask(Eigen::MatrixXi& ice_mask);
void read_ela(Eigen::MatrixXd& ela);
void read_ela0(Eigen::MatrixXd& ela);
void read_debris(Eigen::MatrixXd& debris);
void output_current(std::vector<std::array<int,2>> domain, Eigen::MatrixXd& h, Eigen::MatrixXd& mb, double t);
void output_time(double t);
void overwrite_surface(std::vector<std::array<int,2>> domain, Eigen::MatrixXd& h, Eigen::MatrixXd& bed0 );
void read_cluster(Eigen::MatrixXi& cluster);
void read_basins(Eigen::MatrixXi& basins);
void output_patchy(Eigen::MatrixXd& ela, std::vector<std::array<int,2>> domain);
void output_correction(std::vector<int> corrections, int t);
void output_h_init(std::vector<std::array<int,2>> domain, Eigen::MatrixXd& h);
void output_slope(double slope, double areaf, double tau);
void delete_ice();
#endif
