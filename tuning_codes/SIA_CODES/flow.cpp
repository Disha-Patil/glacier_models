#include "read_write.h"
#include "calc_matrices.h"
#include "avalanche_module.h"

int main(){

	
	int i, j;
  	FILE *test, *test1 ,*test2;                           // be very careful when writing the files using these pointer, you might get segmentation fault for using wrong pointers
  	char filename[64],filename1[64],filename2[64];        // keep the name of files within 64 characters
  	VectorXd x(n0), b(n0),flux(n0), x_guess(n0);                          // initialise vectors : x= solution vector, b=RHS, flux= divergence of flux minus h(i,j); size n0  
  	SpMat A(n0,n0);                                         // the sparse matrix of sigen is defined as SpMat in begining of this code. Set it to A of size nbyn
  	MatrixXd  mb(xmax,ymax);                              // mass balance matrix size xmax by ymax
  	MatrixXd h = MatrixXd::Zero(xmax,ymax) ;             // ice thickness matrix size xmax by ymax; initialised to 0 
  	MatrixXd d = MatrixXd::Zero(xmax,ymax) ;             // diffusion matrix size xmax by ymax; initialised to 0

	//fill mask and bed
	MatrixXi m = MatrixXi::Zero(xmax,ymax) ;
	MatrixXd bed = MatrixXd::Zero(xmax,ymax);
	MatrixXd bed0 = MatrixXd::Zero(xmax,ymax) ;
	MatrixXi ice_mask = MatrixXi::Zero(xmax,ymax) ;
	read_bed_and_mask(bed,m,bed0);
	read_mask(ice_mask);


	// add extra layer around bed
	MatrixXd bed1 = MatrixXd::Zero(xmax,ymax ) ;
	for (int i = 1; i < xmax-1; i++){
    	for (int j = 1; j < ymax-1; j++){
       		if(!m(i,j)){ // access points outside domain and average over neighbours
	   			bed1(i, j) =  ( (4*(bed(i,j) * m(i,j)) )  +  
			       (2* ( (bed(i+1,j) * m(i+1,j)) + (bed(i-1,j) * m(i-1,j)) + (bed(i,j+1) * m(i,j+1))+(bed(i,j-1) * m(i,j-1)) )     ) +
			       (1* ( (bed(i+1,j+1) * m(i+1,j+1)) + (bed(i-1,j-1) * m(i-1,j-1)) + (bed(i-1,j+1) * m(i-1,j+1))+ (bed(i+1,j-1) * m(i+1,j-1)) )  )
			       ) / ( (4*(m(i,j)) )  +  (2* ( ( m(i+1,j)) + (m(i-1,j)) + (m(i,j+1))+ (m(i,j-1)) ) ) + 
                                                       (1* ( (m(i+1,j+1)) + (m(i-1,j-1)) + (m(i-1,j+1))+ (m(i+1,j-1)) ) )  )  ;
	 		}
       		else {bed1(i,j)=bed(i,j);} //inside domain do nothing

     	}
    }

	for (int row = 0; row < xmax; row++){
    	for (int col = 0; col < ymax; col++){
       		bed(row,col)=bed1(row,col);
		}
	}


	// collect domain points
	std::vector<std::array<int,2>> domain;
	for(int i = 0; i < xmax; i++){for(int j = 0; j < ymax; j++){if(m(i,j)){domain.push_back({i,j});}}}

	MatrixXd ela = MatrixXd::Zero(xmax,ymax) ;
	read_ela(ela);

	MatrixXd debris = MatrixXd::Zero(xmax,ymax) ;
	read_debris(debris);

	// fill numbered mask
	int k=0;
	MatrixXi m_idx = MatrixXi::Zero(xmax,ymax ) ;
	for(int ij=0;ij<n0;ij++){
		i = domain[ij][0];
		j = domain[ij][1];
		m_idx(i,j)=k; 
		k++;
	}

	// initialize the vectors for neighbour values calculation in diffusion matrix filling function fill_d

	//sizes are 4*n0 because every point has four neighbours
	VectorXd neigh1(grid_size);// i+1 neighbour
	VectorXd neigh2(grid_size);// j+1 nieghbour
	VectorXd neigh3(grid_size);// i-1 nieghbour
	VectorXd neigh4(grid_size);// j-1 nieghbour

	// Iitialize variables required in the time loop

 	int t1 = 0;                              // the counter used for printing every few iterations
	int t2 = 0;
	int t3 = 0;
	int ti = 0;
	double ice_arr[2] = {};
	double mb_arr[2] = {};
	double areaf_arr[2] = {};
	double areaf = 0;
	double area = 0;
	double tau = 0;
 	double ice_gain=0;                     // ice gain cummulative   
 	double ice_melt=0;                     // ice melt cummulative
 	double ice_b_melt=0 ;                  // ice lost at boundary cummulative
 	double totalice=0;                     // total ice in domain at any point
 

	// Initialize the iterative solver used, in thos case conjugate gradient
	// you can set tolerance or maximum iteration for CG or both
	Eigen::ConjugateGradient<SparseMatrix<double>, Eigen::Lower|Eigen::Upper, Eigen::IdentityPreconditioner > cg;
	cg.setTolerance(1.0e-9);

	int hi1, hi2;
	double h_init1;
	std::ifstream init_file("h_init_file.txt");
	while(init_file >> hi1 >> hi2 >> h_init1){ 
	h(hi1,hi2) = h_init1; }
	init_file.close();
	int count = 0;
	/*
	//output surface = bed
	overwrite_surface(domain,h,bed0);

	//run sortbyheight.py
	Py_Initialize();
	FILE *fd = fopen("sortbyheight.py", "r");
	PyRun_SimpleFileEx(fd, "sortbyheight.py", 1);
	Py_Finalize();

	std::vector<std::vector<int>> sorted_ij;
	std::vector<std::vector<double>> frac;
	gen_ij_and_fracs(sorted_ij, frac);

	std::vector<std::vector<int>> catchment;
	std::vector<std::vector<double>> c_frac;
	gen_catchment_and_c_frac(catchment, c_frac);

	int layer_mask[xmax][ymax];
	std::vector<std::vector<int>> layer_array;
	gen_layer_mask_and_array(layer_mask, layer_array);
	//

	Eigen::MatrixXd av_height = MatrixXd::Zero(xmax, ymax);
	*/
	double DT = 0.01;
	double DT_inv_2_dx_dx = DT*inv_2_dx_dx;
	double t = 0.00;
	auto start = std::chrono::high_resolution_clock::now();

	double prev_gain = 0;
	double prev_melt = 0;
	double prev_b_melt = 0;
	double prev_totalice = 0;
	//remove("ice_minus%d.txt", ela_change);
	double total1 =0;
	double total2 = 0;
	double steady_slope;

	delete_ice();
	//double avalanche_loss = 0;
	//av_height = MatrixXd::Zero(xmax, ymax)
	int ice_count = 0;
	double net_mb = 0.0;
	int inst_area = 0;
	
	while (t <= tmax){
		count += 1;
		ice_count++;

		
		calc_areaf(h,ice_mask,areaf,area);
		

		if(areaf >= 0.8 && tau == 0){ 
			tau = (int(t/25)+1)*25;
		}

		// CALCULATING SLOPE EVERY 10 YEARS
		

		// BREAK CONDITIONS

		if(t == 5000){
			tau = 1500;
			output_slope(steady_slope,areaf,tau);
			cout << "break after 5 tau" << endl;
			break;
		}

		

		if(t > 200 && areaf_arr[0]==areaf_arr[1] && areaf < 0.8 && steady_slope < 0.000001){
			tau = t;
			output_slope(steady_slope,areaf,tau);
			cout << "never reached 80% " << endl;
			break;
		}

		/*if(tau > 0 && t > 3*tau && areaf > 0.95){
			if(steady_slope > 0.001){output_slope(steady_slope,areaf,tau);cout << "break after 3 tau" << endl;break;}
		}*/

		/*if(tau > 0 && t > 3*tau){
			steady_slope = mb_arr[1]-mb_arr[0];
			if(steady_slope > 0.00001){output_slope(steady_slope,areaf,tau);cout << "5555555555" << endl;break;}
		}*/

		 


		// Output current state of basin
		/*if(t1 == 5000){ 
			t1=0; 
			//printf("%lf\n",t);
			output_current(domain, h,mb,t);
		}*/

		int p = 0;
		/*
		// Avalanche
		avalanche(bed, sorted_ij, frac, h, av_height,1);
		spread(h, av_height, catchment, c_frac, layer_mask, layer_array, domain);
		for(int i = 0; i < 4; i++){smoothen(h, catchment, bed, domain, m, av_height);}
		for(int ij=0;ij<n0;ij++){
			i = domain[ij][0];
			j = domain[ij][1];
			h(i,j) += av_height(i,j);
			av_height(i,j) = 0;
		}

		


		// UPDATING X
		for(int ij=0;ij<n0;ij++){
			i = domain[ij][0];
			j = domain[ij][1];
			x(p) = h(i,j);p+=1;
		}
		*/

		
		                             								// print current iteration 
		fill_mb(mb,bed,ela, debris, h, domain);                         						// update mass balance
		
		fill_d(d,h,bed,m,neigh1,neigh2,neigh3,neigh4); 						// update diffusion
		
		
		fill_b( b,flux, d,mb,m, bed, h, DT, DT_inv_2_dx_dx, domain);               	// update b vector
		
		fill_A( A,d,m,m_idx, DT_inv_2_dx_dx, domain);                          		// update A matrix
		
		cg.compute(A);                                 								// compute A matrix with CG

		x = cg.solveWithGuess(b,x);                               					// solve for Ax=b
		//std::cout << "#iterations:     " << cg.iterations() << std::endl;

		p = 0;                                    								// counter for accessing the solution vector (x) values 


		// fill h matrix
		for(int ij=0;ij<n0;ij++){
			i = domain[ij][0];
			j = domain[ij][1];

			/* calculate total ice gain and ice melt. Gain sum over positive mass balance points and melt is sum over negative points.
			When, solution shows ice thicness as negative, error at boundaries is the reason. 
			Calculate this extra melt at boundary by summing over negative  massbalance there and add the negative ice. 
			The result is extra melt at these points.*/

			if ((x(p)) > 0){
				if (mb(i,j)>0) {ice_gain +=  DT*mb(i,j);}
				else {ice_melt+= DT*mb(i,j);}
			}
			else {ice_b_melt += (-DT*mb(i,j) + x(p));}
	
			// update h matrix with new solution
			h(i,j) = x(p);
			// increase counter for x vector
			p++;       
				
			// negative ice thickness is set to zero
			if ((h(i,j) < 0)){ h(i,j) = 0;}   
		}

		
		t2++;


		


		// look at ice currently in glacier domain by summing over all h values
		
		totalice=0;
		net_mb = 0.0;
		inst_area = 0;                         
		for(int ij=0;ij<n0;ij++){
			i = domain[ij][0];
			j = domain[ij][1];
			totalice += h(i,j);
			if(h(i,j)>0){
				net_mb+=mb(i,j);
				inst_area += 1;
				}		
		}
		//avalanche_loss -= totalice;
		double conservation_sum = 0;
 
		conservation_sum += (ice_gain+ice_melt-ice_b_melt-totalice);
		conservation_sum -= (prev_gain+prev_melt-prev_b_melt-prev_totalice);

		//if(fabs(conservation_sum) > 10){std::cout << "ICE CONSERVATION ERROR" << std::endl;break;}
		

		// print for ice conservation check every iteration 
		
		

		if(count == 20000){
			count = 0;
			areaf_arr[0] = areaf_arr[1];
			areaf_arr[1] = areaf;

			steady_slope = (totalice-prev_totalice)/(totalice*200);

			if(areaf_arr[0]==areaf_arr[1] && steady_slope > 0.001){
				tau = t;
				output_slope(1.0,areaf,tau);
				cout << "TOO MUCH ICE" << endl;
				break;

			}

			if(areaf_arr[0]==areaf_arr[1] && steady_slope < 0.000001){
			tau = t;
			output_slope(steady_slope,areaf,tau);
			cout << "break after steady condition" << endl;
			break;
			}
		}

		if(ice_count==20000){
			ice_count = 0;
			sprintf( filename2, "ice.txt"); 
			test2=fopen(filename2,"a");                                                         // file is appended
			fprintf(test2,"%lf %lf %lf %lf %lf %lf %lf\n",t, ice_gain,ice_melt, ice_b_melt, totalice, conservation_sum,steady_slope);
			fclose(test2);
			prev_totalice = totalice;
		}

		prev_gain = ice_gain;
		prev_melt = ice_melt;
		prev_b_melt = ice_b_melt;
		


		

		t1++;
		ti++;
		t3++;
		/*
		// UPDATE AVALANCHE	
		if(t2 == 500){
			t2 = 0;
			//output surface
			overwrite_surface(domain,h,bed0);

			//run sortbyheight.py
			Py_Initialize();
			FILE *fd = fopen("sortbyheight.py", "r");
			PyRun_SimpleFileEx(fd, "sortbyheight.py", 1);
			Py_Finalize();

			//read in updated files
			gen_ij_and_fracs(sorted_ij, frac);
			gen_catchment_and_c_frac(catchment, c_frac);
			gen_layer_mask_and_array(layer_mask, layer_array);
			
		}
		*/
		
		t += DT;

	}

 	auto stop = std::chrono::high_resolution_clock::now();
 	std::chrono::duration<double> timeelapsed = stop - start;
 	//cout << "Time taken  -  " << timeelapsed.count() << endl;

	output_current(domain, h, mb, t);
}
