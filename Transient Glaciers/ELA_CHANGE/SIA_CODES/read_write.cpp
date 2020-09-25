#include "globals.h"


void read_bed_and_mask(Eigen::MatrixXd& bed, Eigen::MatrixXi& mask, Eigen::MatrixXd& bed0){
    int i,j;
    double bed_h;

    std::ifstream bed_file("smoothbed.txt");
    while(bed_file >> i >> j >> bed_h ){
        bed(i,j) = bed_h; 
        if(bed_h){mask(i,j) = 1;}
    }
    bed_file.close();

    std::ifstream bed0_file("bed.txt");
    while(bed0_file >> i >> j >> bed_h){ bed0(i,j) = bed_h; }
    bed0_file.close();
}

void read_ela(Eigen::MatrixXd& ela){
    double value;
    int i,j;
    std::ifstream ela_file("ela.txt");
    while(ela_file >> i >> j >> value ){ela(i,j) = value;}
    ela_file.close();
}

void read_ela0(Eigen::MatrixXd& ela){
    double value;
    int i,j;
    std::ifstream ela_file("ela0.txt");
    while(ela_file >> i >> j >> value ){ela(i,j) = value;}
    ela_file.close();
}

void read_debris(Eigen::MatrixXd& debris){
    double value;
    int i,j,k;
    std::ifstream db_file("debris.txt");
    while(db_file >> i >> j >> k){debris(i,j) = k;}
    db_file.close();
}

void read_basins(Eigen::MatrixXi& basins){
    int value;
    int i,j;
    std::ifstream ela_file("basins");
    while(ela_file >> i >> j >> value ){basins(i,j) = value;}
    ela_file.close();
}

void read_cluster(Eigen::MatrixXi& cluster){
    double value;
    int i,j;
    std::ifstream ela_file("glabtop_ij");
    while(ela_file >> i >> j >> value ){if(value > 80){cluster(i,j) = value;}}
    ela_file.close();
}

void read_mask(Eigen::MatrixXi& ice_mask){
    double value;
    int i,j;
    std::ifstream ela_file("ice_mask");
    while(ela_file >> i >> j >> value ){ice_mask(i,j) = value;}
    ela_file.close();
}

void output_current(std::vector<std::array<int,2>> domain, Eigen::MatrixXd& h, Eigen::MatrixXd& mb, double t){
    int i,j;
    char filename1[64];
    sprintf(filename1, "h_steady.txt" ); 
	FILE *test1 = fopen(filename1,"w");
	for(int ij=0;ij<n0;ij++){
		i = domain[ij][0];
		j = domain[ij][1];
		fprintf(test1,"%d %d %lf %lf\n",i,j,h(i,j),mb(i,j));
	}
    
    fclose(test1);
}

void output_va(std::vector<std::array<int,2>> domain,Eigen::MatrixXd& h,double t){
    int i,j;
    double total_h = 0;
    double area = 0;
    char filename1[64];
    sprintf(filename1, "va"); 
	FILE *test1 = fopen(filename1,"a");
	for(int ij=0;ij<n0;ij++){
		i = domain[ij][0];
		j = domain[ij][1];
        if(h(i,j)>0){
            total_h+= h(i,j)/100000;
            area+= 0.01;
        }
    }
    fprintf(test1,"%.2lf %.6lf %.2lf\n",t,total_h,area);
    fclose(test1);
}

void output_h_init(std::vector<std::array<int,2>> domain, Eigen::MatrixXd& h){
    int i,j;
    char filename1[64];
    sprintf(filename1, "h_init_file.txt" ); 
	FILE *test1 = fopen(filename1,"w");
	for(int ij=0;ij<n0;ij++){
		i = domain[ij][0];
		j = domain[ij][1];
		fprintf(test1,"%d %d %lf\n",i,j,h(i,j));
	}
    fclose(test1);
}

/*void output_time(double t, auto start){
    char filename2[64];
    auto stop = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> timeelapsed = stop - start;
	double tt = timeelapsed.count();
	sprintf(filename2, "time_time.txt"); 
	FILE *test2 = fopen(filename2,"a");                                                         // file is appended
	fprintf(test2,"%lf %lf\n",t, tt );
	fclose(test2);
}*/

void overwrite_surface(std::vector<std::array<int,2>> domain, Eigen::MatrixXd& h, Eigen::MatrixXd& bed0 ){
    char filename[64];
    int i,j;
    sprintf(filename, "surface.txt" );
	FILE *test = fopen(filename,"w");
	for(int ij=0;ij<n0;ij++){
		i = domain[ij][0];
		j = domain[ij][1];
		fprintf(test,"%d %d %lf\n", i,j, h(i,j)+bed0(i,j));
	}
	fclose(test);		
}

void output_patchy(Eigen::MatrixXd& ela, std::vector<std::array<int,2>> domain){
    char filename[64];
    int i,j;
    sprintf(filename, "patchy_ela" );
	FILE *test = fopen(filename,"w");
	for(int ij=0;ij<n0;ij++){
		i = domain[ij][0];
		j = domain[ij][1];
		fprintf(test,"%d %d %lf\n", i,j,ela(i,j));
	}
	fclose(test);		
}

void output_correction(std::vector<int> corrections, int t){
    char filename[64];
    int i,j;
    sprintf(filename, "correction%d", t );
	FILE *test = fopen(filename,"w");
	for(int i = 1; i < corrections.size(); i++){
        fprintf(test,"%d\n",corrections[i]);
    }
	fclose(test);		
}

void output_slope(double slope, double areaf, double tau){
    char filename[64];
    sprintf(filename, "steady_slope");
	FILE *test = fopen(filename,"w");
    fprintf(test,"%.12f %f %f\n",slope,areaf, tau);
	fclose(test);		
}

void delete_ice(){
    char filename[64];
    sprintf(filename, "ice.txt");
	FILE *test = fopen(filename,"w");
    fclose(test);		
}

void delete_va(){
    char filename[64];
    sprintf(filename, "va");
	FILE *test = fopen(filename,"w");
    fclose(test);		
}