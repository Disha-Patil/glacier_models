#include "globals.h"

void fill_d(MatrixXd& d, MatrixXd& h, MatrixXd& bed, MatrixXi& m,VectorXd& neigh1,VectorXd& neigh2,VectorXd& neigh3,VectorXd& neigh4){

VectorXd neigh_m1(grid_size); // i+1 point mask
VectorXd neigh_m2(grid_size); // j+1 point mask
VectorXd neigh_m3(grid_size); // i-1 point mask
VectorXd neigh_m4(grid_size); // j-1 point mask


int ng1=0;int ng2=0;int ng3=0;int ng4=0;  // index counter

// run the loop so that all points from 0 to xmax-1 are accessible                    
for(int i=1;i<(xmax-1);i++){
    for(int j=1;j<(ymax-1);j++){

        if(m(i+1,j)){neigh1(ng1)=bed(i+1,j)+h(i+1,j);neigh_m1(ng1)=1;} // if neighbour present central difference, hence 1 is stored in mask value 
	else {neigh1(ng1)=(bed(i,j)+h(i,j));neigh_m1(ng1)=2;}          // else forward/backward difference, hence to adjust for denominator 2 is stored in mask value
	ng1 = ng1+1;

	if(m(i,j+1)){neigh2(ng2)=bed(i,j+1)+h(i,j+1);neigh_m2(ng2)=1;} // if neighbour present central difference, hence 1 is stored in mask value
	else if(!m(i,j+1)){neigh2(ng2)=(bed(i,j)+h(i,j));neigh_m2(ng2)=2;} // else forward/backward difference, hence to adjust for denominator 2 is stored in mask value
	ng2 = ng2+1;

	if(m(i-1,j)){neigh3(ng3)=bed(i-1,j)+h(i-1,j);neigh_m3(ng3)=1;} // if neighbour present central difference, hence 1 is stored in mask value
	else if(!m(i-1,j)){neigh3(ng3)=(bed(i,j)+h(i,j));neigh_m3(ng3)=2;}// else forward/backward difference, hence to adjust for denominator 2 is stored in mask value
	ng3 = ng3+1;

	if(m(i,j-1)){neigh4(ng4)=bed(i,j-1)+h(i,j-1);neigh_m4(ng4)=1;} // if neighbour present central difference, hence 1 is stored in mask value
	else if(!m(i,j-1)){neigh4(ng4)=(bed(i,j)+h(i,j));neigh_m4(ng4)=2;}// else forward/backward difference, hence to adjust for denominator 2 is stored in mask value
	ng4 = ng4+1;

    }
    }	



int ngg=0;
for(int i=1;i<(xmax-1);i++)
  {
    for(int j=1;j<(ymax-1);j++)
      {
	d(i,j) =   (C_inv_2_dx_sq) *  (h(i,j)*h(i,j)*h(i,j)*h(i,j)*h(i,j))   *  
				(
				(neigh_m1(ngg)*neigh_m3(ngg)*(neigh1(ngg) - neigh3(ngg))) * 
				(neigh_m1(ngg)*neigh_m3(ngg)*(neigh1(ngg) - neigh3(ngg)))   +    
				(
				(neigh_m2(ngg)*neigh_m4(ngg)*(neigh2(ngg) - neigh4(ngg))) *
				(neigh_m2(ngg)*neigh_m4(ngg)*(neigh2(ngg) - neigh4(ngg)))
				)
				 );
	ngg=ngg+1;
	}
  }
  

	
 return;
}


void fill_mb(MatrixXd& mb, MatrixXd& bed, MatrixXd& ela, MatrixXd& debris, MatrixXd& h, std::vector<std::array<int,2>> domain) {
  //double R2,pabs; //switch on for lemeur bed
  double q;   // for neatness of calculation
	int i, j;
  for(int ij=0;ij<n0;ij++){
    i = domain[ij][0];
		j = domain[ij][1];
	  //linear bed with accumulation at centre and rest is ablation
	  /*if(i>8 && i<12 && j>3 && j<7 ) {  mb(i,j)=1;  }
	  else { mb(i,j)= -(1);  }*/
	  
	  //mass balance as described in lemeur's paper
	  /*R2 = ((1750 - ((i+0.5)*(DX)))  *  (1750 - ((i+0.5)*(DX)))) + (((j)*(DY) - 2000) * ((j)*(DY) - 2000));
	  pabs = abs ( ((Ra1)*(Ra1)) - R2 );
	  q = ((Ra1)*(Ra1)) - R2;
	  mb(i,j) = (a0 * (pabs/q) * (sqrt(pabs)/(Ra1)));*/

	  q = mbgrad*(bed(i,j)+h(i,j) - (ela(i,j)+(double)ela_change));  // calculate mass balance with ice thickness feedback
	  mb(i,j)=q;
		if(q>=1.0){mb(i,j)=1.0;}                     // restrict the accumulation to 1 as maximum 
		//if(debris(i,j) && q <= -2){mb(i,j) = -2.0;}
  }
 return;
}

void calc_areaf(MatrixXd& h,MatrixXi& ice_mask, double& areaf, double& area){
	double nu = 0;
	double den = 0;
	area = 0;
	for(int i = 0; i < xmax; i++){
		for(int j = 0; j<ymax;j++){
			if(h(i,j)>0 && ice_mask(i,j)>0){nu += 1;}
			if(ice_mask(i,j)>0){den += 1;}
			if(h(i,j)>0){area += 1;}
		}
	}
	areaf = nu/den;
	return;
}


void fill_b(VectorXd& b,VectorXd& flux,  MatrixXd& d, MatrixXd& mb, MatrixXi& m,MatrixXd& bed, MatrixXd& h,double DT,double DT_inv_2_dx_dx, std::vector<std::array<int,2>> domain){

  // the constants from discretised equation
  double alpha, beta, gamma, eps;


  //std::cout << " Filling b " << std::endl;  // use to print b column vector
  int k=0;
  int i, j;
  for(int ij=0;ij<n0;ij++){
  i = domain[ij][0];
	j = domain[ij][1];
	/*Define all the constant terms which appear in the matrix A and bed vector, only for points for which maskk is non zero.
	Loop runs only on non zero mask values, this ensures we are finding solution of only the points in glacier domain and not outside it. This also matches the size of b vector which is n0 (NNZ)
	*/
	alpha = ((d(i+1,j) + d(i,j)));
	beta = ((d(i,j) + d(i-1,j))) ;
	gamma = ((d(i,j+1) + d(i,j)));
	eps  = ((d(i,j) + d(i,j-1))) ;
	//c1 = 1 + ((DT_inv_2_dx_dx*(alpha*m[i+1][j] + bedeta*m[i-1][j] + gamma*m[i][j+1] + eps*m[i][j-1]))*.5); //needed in A matrix calculations
	   
	// calculating the following value is optional. It is not the flux, but divergence of flux minus h(i,j)
  flux(k) = DT_inv_2_dx_dx*(
				     (alpha*(bed(i+1,j) - bed(i,j))*(m(i+1,j)))
				     - (beta*(bed(i,j) - bed(i-1,j))*(m(i-1,j)))
				     + (gamma*(bed(i,j+1) - bed(i,j))*(m(i,j+1)))
				     - (eps*(bed(i,j) - bed(i,j-1))*(m(i,j-1)))
				     + (   ( (alpha*(h(i+1,j) - h(i,j))*(m(i+1,j)))
                                     - (beta*(h(i,j)  - h(i-1,j))*(m(i-1,j)))
                                    +(gamma*(h(i,j+1) - h(i,j))*(m(i,j+1)))
                                     - (eps*(h(i,j) - h(i,j-1))*(m(i,j-1)))
					     *0.5))
                                     );

	   // calculating the b vector
	b(k)=(h(i,j) +  ((DT)*mb(i,j)) + DT_inv_2_dx_dx*(
                                     (alpha*(bed(i+1,j) - bed(i,j))*(m(i+1,j)))
                                   - (beta*(bed(i,j) - bed(i-1,j))*(m(i-1,j)))
                                   + (gamma*(bed(i,j+1) - bed(i,j))*(m(i,j+1)))
                                   - (eps*(bed(i,j) - bed(i,j-1))*(m(i,j-1)))
                                 + (((alpha*(h(i+1,j) - h(i,j))*(m(i+1,j)))
                               - (beta*(h(i,j)  - h(i-1,j))*(m(i-1,j)))
                                + (gamma*(h(i,j+1) - h(i,j))*(m(i,j+1)))
                                - (eps*(h(i,j) - h(i,j-1))*(m(i,j-1))))*0.5)));
	
	   //printf("%d %d %lf \n",i,j,b(k));				
	   k++;
    //printf("%i \n", i);
    }
			    
    // std::cout << " Filled b " << std::endl << b << std::endl;   // use to print b column vector
 
 return;
}


void fill_A(SpMat& A, MatrixXd& d,MatrixXi& m,MatrixXi& m_idx, double DT_inv_2_dx_dx, std::vector<std::array<int,2>> domain) {
    double alpha, beta, gamma, eps, c1;

    // initialization to fill a sparse eigen matrix
 typedef Eigen::Triplet<double> T;  // stores i j indices and value to be stored at i,j
 std::vector<T> aa;                 // the triplet builds a vector aa

    //std::cout << " Filling A " << std::endl;   // use for printing 

 int l=0;                           //counter to fill the matrix
    // int k=0;
    int i, j;
  for(int ij=0;ij<n0;ij++){
    i = domain[ij][0];
		j = domain[ij][1];
	/*m_idx stores sequential index for all the non zero elements in the domain. 
          So neighbours at i,j point are directly accessed by the sequential index number as it appears in m_idx matrix. 
	  For more look at how m_idx is filled in main().
	  */
	//printf("done\n");
	//Define all the constant terms which appear in the matrix A and b vector
	
	  alpha = ((d(i+1,j) + d(i,j)));
	  beta = ((d(i,j) + d(i-1,j))) ;
	  gamma = ((d(i,j+1) + d(i,j)));
	  eps  = ((d(i,j) + d(i,j-1))) ;
	  c1 = 1 + ((DT_inv_2_dx_dx*(alpha*m(i+1,j) + beta*m(i-1,j) + gamma*m(i,j+1) + eps*m(i,j-1)))*.5);

	  
	   aa.push_back( T(l,l , c1) );
	  //aa.push_back( T(k,k, c1));//do not use this
	  if(m(i+1,j))
	    {
	      aa.push_back( T(l,m_idx(i+1,j) , ((-DT_inv_2_dx_dx*alpha*m(i+1,j))*.5)) );
	      //aa.push_back( T(k,((i)*(xmax-2)) + j-1, ((-DT_inv_2_dx_dx*alpha*m(i+1,j))*.5))); //risky to use such indexing method, but if you want to go down this path, be my guest!
	    }
	  if(m(i-1,j))
	    {
	      aa.push_back( T(l,m_idx(i-1,j) , ((-DT_inv_2_dx_dx*beta*m(i-1,j))*.5) ) );
	      //aa.push_back( T(k,((i-2)*(xmax-2)) + j-1, ((-DT_inv_2_dx_dx*beta*m(i-1,j))*.5))); //risky to use such indexing method, but if you want to go down this path, be my guest!
	    }

	  if(m(i,j+1))
	    {
	      aa.push_back(T( l,m_idx(i,j+1), ((-DT_inv_2_dx_dx*gamma*m(i,j+1))*.5) ) );
	      //aa.push_back( T(k,((i-1)*(xmax-2)) + j , ((-DT_inv_2_dx_dx*gamma*m(i,j+1))*.5))); //risky to use such indexing method, but if you want to go down this path, be my guest!
	    }
	  if(m(i,j-1))
	    {
	      aa.push_back(T(l,m_idx(i,j-1), ((-DT_inv_2_dx_dx*eps*m(i,j-1))*.5) ) );
	      //aa.push_back( T(k,((i-1)*(xmax-2)) + j - 2 , ((-DT_inv_2_dx_dx*eps*m(i,j-1))*.5))); //risky to use such indexing method, but if you want to go down this path, be my guest!
	    }
	  l++;
	  //k++;
	
  
    //printf("%i \n", i);
  }     

 A.setFromTriplets(aa.begin(), aa.end());       // form the sparse matrix from triplets vector aa
    //printf("done\n");
 //std::cout << " Filled A " << std::endl ;//<< A << std::endl; // use for printing
 
 return;
}
