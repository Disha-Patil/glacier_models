#include "globals.h"


void gen_ij_and_fracs(std::vector<std::vector<int>> &sorted_ij_, std::vector<std::vector<double>> &frac_){
    sorted_ij_.clear();
    frac_.clear();
    int i; int j;
    std::vector<int> temp(2,0);
    std::ifstream ij_file("ij_sorted.txt");
    while(ij_file >> i >> j){
        temp[0] = i; temp[1] = j;
        sorted_ij_.push_back(temp);
    }
    ij_file.close();


    double f1, f2, f3, f4;
    std::vector<double> tmp(4,0);
    std::ifstream frac_file("fractions.txt");
    while(frac_file >> f1 >> f2 >> f3 >> f4){
        tmp[0] = f1; tmp[1] = f2; tmp[2] = f3; tmp[3] = f4;
        frac_.push_back(tmp);
    }
    frac_file.close();
}


void gen_catchment_and_c_frac(std::vector<std::vector<int>> &catchment_, std::vector<std::vector<double>> &c_frac_){
    catchment_.clear();
    c_frac_.clear();
    int i; int j;
    std::vector<int> temp(2,0);
    std::ifstream c_ij_file("catchment.txt");
    while(c_ij_file >> i >> j){
        temp[0] = i; temp[1] = j;
        catchment_.push_back(temp);
    }
    c_ij_file.close();

    double f1, f2, f3, f4;
    std::vector<double> tmp(4,0);
    std::ifstream c_frac_file("c_fractions.txt");
    while(c_frac_file >> f1 >> f2 >> f3 >> f4){
        tmp[0] = f1; tmp[1] = f2; tmp[2] = f3; tmp[3] = f4;
        c_frac_.push_back(tmp);
    }
    c_frac_file.close();

}


void gen_layer_mask_and_array(int layer_mask_[xmax][ymax], std::vector<std::vector<int>> &layer_array_){
    layer_array_.clear();
    memset(layer_mask_, 0, xmax*ymax);
    int i, j, k;
    std::vector<int> temp(3,0);
    std::ifstream lm_file("layer_array.txt");
    while(lm_file >> i >> j >> k){
        temp[0] = i, temp[1] = j; temp[2] = k;
        layer_mask_[i][j] = k;
        layer_array_.push_back(temp);
    }
    lm_file.close();
}

// check logic
void avalanche(Eigen::MatrixXd& bed_, std::vector<std::vector<int>> &sorted_ij_, std::vector<std::vector<double>> &frac_, Eigen::MatrixXd& h, 
 Eigen::MatrixXd& av_height_,double avalanche_frac){
    int l = sorted_ij_.size();
    int i , j, i1, j1;
    std::vector<double> current_fracs;
    double current_height;
    int iv[4] = {-1,0,0,1};
    int jv[4] = {0,-1,1,0};

    for(int s = 0; s < l ; s++){
        i = sorted_ij_[s][0];
        j = sorted_ij_[s][1];
        current_fracs = frac_[s];

        current_height = avalanche_frac*h(i,j);
        if(current_height > 0.05){current_height = 0.05;}
        
        for(int k = 0; k < 4; k++){
            i1 = i + iv[k];
            j1 = j + jv[k];
            if(current_fracs[k]){
                av_height_(i1,j1) += current_fracs[k]*current_height + current_fracs[k]*av_height_(i,j);
                av_height_(i,j) -= current_fracs[k]*av_height_(i,j);
                h(i,j) -= current_fracs[k]*current_height;
                
            }
        }
    }
}

//check logic
void spread(Eigen::MatrixXd& h, Eigen::MatrixXd& av_height_, std::vector<std::vector<int>> &catchment_, 
std::vector<std::vector<double>> &c_frac_, int layer_mask_[xmax][ymax], std::vector<std::vector<int>> &layer_array_, std::vector<std::array<int,2>> domain ){
    int i, j, layer_no;
    std::vector<double> current_fracs;
    double current_height;
    double portion = 1.0;
    int iv[4] = {-1,0,0,1};
    int jv[4] = {0,-1,1,0};
    int l = layer_array_.size();
    
    int i1, j1;
    int next_layer_no;
    for(int s = 0; s < l ; s++){
        i = layer_array_[s][0];
        j = layer_array_[s][1];
        current_fracs = c_frac_[s];
        layer_no = layer_mask_[i][j];
        
        portion = (layer_no - 1)/(double)layer_no ;

        current_height = portion*av_height_(i,j);
        

        for(int k = 0; k < 4; k++){
            i1 = i+iv[k];
            j1 = j+jv[k];
            next_layer_no = layer_mask_[i1][j1];
            
            av_height_(i1, j1) += current_fracs[k]*current_height;
            av_height_(i,j) -= current_fracs[k]*current_height;
        }
    }

    /*// UPDATE H
    for(int ij=0;ij<n;ij++){
        i = domain[ij][0];
		j = domain[ij][1];
        h(i,j)+=av_height_(i,j);
    
    }*/

    

}

void smoothen(Eigen::MatrixXd& h,std::vector<std::vector<int>> &catchment_, Eigen::MatrixXd& bed_, std::vector<std::array<int,2>> domain, Eigen::MatrixXi& mask, Eigen::MatrixXd& av_height){

    int iv[8] = {-1,0,0,1,-1,1,-1,1};
    int jv[8] = {0,-1,1,0,-1,-1,1,1};
    double maxgrad = 0.577;
    double h0;
    int i , j, i1, j1;
    int l = catchment_.size();
    
    int temp;
    double min_h;
    double current_height;
    for(int s = 0; s < l ; s++){
        i = catchment_[s][0];
        j = catchment_[s][1];

        h0 = av_height(i,j);

        temp = 0;
        min_h = h0;
        for(int m = 0; m < 4; m++){
            i1 = i + iv[m];
            j1 = j + jv[m];
            if(av_height(i1,j1) < min_h){min_h = av_height(i1, j1);}
        }

        if(h0 < (DX*maxgrad + min_h)){continue;}

        current_height = av_height(i,j)-min_h;

        int total = 2;
        for(int m = 0; m < 4; m++){if(mask(i+iv[m], j+jv[m])){ total += 3;}}
        for(int m = 4; m < 8; m++){if(mask(i+iv[m], j+jv[m])){ total += 2;}}

        av_height(i,j) -= ((total-2)/(double)total)*current_height;

        for(int m = 0; m < 4; m++){
            i1 = i + iv[m];
            j1 = j + jv[m];
            if(mask(i1,j1)){av_height(i1,j1) += current_height*(3/(double)total);}
        }

        for(int m = 4; m < 8; m++){
            i1 = i + iv[m];
            j1 = j + jv[m];
            if(mask(i1,j1)){av_height(i1,j1) += current_height*(2/(double)total);}
        }
    }

    /*for(int ij=0;ij<n;ij++){
        i = domain[ij][0];
		j = domain[ij][1];
        h(i,j) += av_height(i,j);
        av_height(i,j) = 0;
    }*/

}
