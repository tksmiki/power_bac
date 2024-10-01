//automatic parallelization
// -parallel -par-report3 -par-threshold100 

//for spatially constant sinking rate, copied from PDE_pp_ratio03new.c at 2015/01/08

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "../new_util/nrutil_new.h"
#include "../new_util/forode.h"
#include "../new_util/nrutil_new.c"
#include "../new_util/rk4_new.c"

#include "./global_params_pp2021_pow_approx.c"
//#include "./pp_ratio2021di_pow_approx.c"

//For OPENMP, -fopenmp (gcc OSX) or -openmp (icc Linux)

#define NUM_R 400000
#define NUM_C 2

#define DELTAZERO 1.0e-2  //0.01 for AGG but 0.01 for MONODZERO
#define EPSS 1.0e-6
#define ZERO 1.0e-20

#define NUM_LAYER 297  //number of layer
#define MIXED_LAYER  50.0 //[m]
//#define min_app_N 200

#define INTERVAL_DAY 10000 //interval for printf

void differential(double time, double in[], double out[]);
double dRkdt(double time, double vr[], int k);			//for resource at k-th layer
double dPkdt(double time, double vr[], int k);			//for predator (free-living consumer) at k-th layer
double consumption(double resource, double predator);
double growthP(double resource, double predator, int k);
double suppl(int k); //depth specific supply rate
double diff_coeff(int k); //depth specific diffusion coefficient

//for approximately calculate the average value
double av_cont(double **discont_data, double interval, double lamd);			
//for approximately calculate the differential coefficient
double av_cont_diff(double **av_data, double interval_av, double interval_diff, double lamd); 


int jP[NUM_LAYER + 1], jR[NUM_LAYER +1];

double **matrix_lambda; //need to be global to avoid using them as parameters in differential()
double interval_discontinuous, interval_differential; //need to be global to avoid using them as parameters in differential()

double depth_zero_boundary_R, depth_zero_boundary_P;  //for zero_flux at depth z = 0
double depth_bottom_boundary_R, depth_bottom_boundary_P; //for zero_flux via diffusion and no sinking
double depth_zero_boundary_PA, depth_bottom_boundary_PA; //for particle-attached bacteria

int main(int arg, char *argv[])
{	
  double t, tzero, end_time, deltat;
  double *v;
  double *vn;
  double *dfdt;
  
  long n_variable;
  long i, j, k, m;
  long write_index;
  
  long ran_seed;
  unsigned long ran_seedmt;
  double *average;
  double *cv;
  
  double *totalPA;
  double lammda;
  
  long appNk;
  
  double totalB, totalR;
  
  double depth_integ_PP; //depth integrated primary production   [mgC/m2/d]
  double sinking_rate_forPlot, POC_flux_forPlot;
  
  FILE *fp1;
  matrix_lambda = d_matrix(NUM_R, NUM_C); //allocation of memory
  
  fp1 = fopen("av_lam_data_alf=2.0e-3_v2.csv", "r");
  //fp1 = fopen("av_lam_data_alf=2.0e-2_v3.csv", "r");
    if (fp1 == NULL) message_error("File cannot be opened.");
  
  for(i = 1; i <= NUM_R; i++) {
    for(j = 1; j <= NUM_C; j++) {
      fscanf(fp1, "%lf,", &matrix_lambda[i][j]);
    }
  }
 
  interval_discontinuous = 0.1;
  interval_differential = 0.01;
  
  tzero = 0.0;
  end_time = 365.0*20.0;
  //end_time = 365.0;

  
  //assign POC supply rate 
  s0 = primary_prod/(bac_carbon*cPperCB*MIXED_LAYER);
  
  deltat = DELTAZERO;
  
  diff0 /= (delta_x*delta_x);
  sinking_R /= delta_x;
  //sinking_P /= delta_x;
  
  n_variable = 2*NUM_LAYER;
  v = d_vector(n_variable);
  vn = d_vector(n_variable);
  dfdt = d_vector(n_variable);
  average = d_vector(n_variable);
  cv = d_vector(n_variable); 
  totalPA = d_vector(NUM_LAYER);
  
  write_index = 1;

  alfa = 0.5*(d1 - a1);		//for density_dependent behavior
  
  for(k = 1; k <= NUM_LAYER; k++) {
    jR[k] = 2*k - 1;
    jP[k] = 2*k;
  }
  
  //assignment of initial conditions	
  
  for(k = 1; k <= NUM_LAYER; k++) {
    v[jR[k]] = 1.0;
    v[jP[k]] = 1.0;
  }
#ifdef MONOD_ZERO
  for(k = 1; k <= NUM_LAYER; k++) {
    v[jP[k]] = (1.0 + (a0/d0)*exp(-0.5*(a1 + d1))*v[jR[k]])*v[jP[k]];
  }	
#endif	
  //setting global parameters
  
  
#ifdef TIMESERIES	
  
  //for(k = 1; k <= NUM_LAYER; k++) printf("%.10lf\t%lf\t%.10lf\t%.10lf\n", t, k*delta_x, v[jR[k]], v[jP[k]]);
  
  for(t = tzero; t <= end_time +ZERO;) {
    //setting boundary conditions (for reflecting boundary////
    //depth_zero_boundary_R = v[jR[2]];  //for reflecting boundary
    //depth_zero_boundary_P = v[jP[2]];  //for reflecting boundary
    depth_bottom_boundary_R = v[jR[NUM_LAYER - 1]]; //for reflecting boundary
    depth_bottom_boundary_P = v[jP[NUM_LAYER - 1]]; //for reflecting boundary
    
#ifdef AGG		
    for(k = 1; k <= NUM_LAYER; k++) {
      lammda = (a0/d0)*exp(-0.5*(a1 + d1))*v[jP[k]];
      totalPA[k] = av_cont(matrix_lambda, interval_discontinuous, lammda)*v[jR[k]];
    }
#endif		
    
    differential(t, v, dfdt);
    rk4(v, dfdt, n_variable, t, deltat, v, differential);
    t += deltat;
    //v[jP[1]] = 10.0;  //for fixed boundary condition
    //v[jR[1]] = 1.0;  //for fixed boundary condition
    
    //for(k = 1; k <= n_variable; k++) v[k] = v[k];  //renewal of vector
    
  // printf("%.5lf\n", t);
   
    if(write_index < INTERVAL_DAY/deltat) write_index++;
    else {
      for(k = 1; k <= NUM_LAYER; k++) {
#ifdef AGG	
        lammda = (a0/d0)*exp(-0.5*(a1 + d1))*v[jP[k]];
        totalPA[k] = av_cont(matrix_lambda, interval_discontinuous, lammda)*v[jR[k]];
        //printf("%.10lf\n", av_cont(matrix_lambda, interval_discontinuous, lammda));
#endif	
      }// end of for k
      
      printf("%.10lf\t%lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n", t, 1*delta_x, v[jR[1]], v[jP[1]], totalPA[1], v[jP[1]]+totalPA[1], diff_coeff(1));
      printf("%.10lf\t%lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n", t, 5*delta_x, v[jR[1]], v[jP[1]], totalPA[1], v[jP[1]]+totalPA[1], diff_coeff(1));
      
      for(k = 2; k <= NUM_LAYER; k++) {
        printf("%.10lf\t%lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n", t, (k - 1 + MIXED_LAYER/delta_x)*delta_x, v[jR[k]], v[jP[k]], totalPA[k], v[jP[k]]+totalPA[k], diff_coeff(k));
      }//end of for k

      printf("\n");	
      write_index = 1;
    }
        			
  }//end of for t
//printf("End of t loop\n");
#endif //end of #ifdef for TIMESERIES
  for(k = 1; k <= NUM_LAYER; k++) {
#ifdef AGG	
    lammda = (a0/d0)*exp(-0.5*(a1 + d1))*v[jP[k]];
    totalPA[k] = av_cont(matrix_lambda, interval_discontinuous, lammda)*v[jR[k]];
#endif
  }		
  depth_integ_PP = bac_carbon*cPperCB*s0*MIXED_LAYER/1000.0; //change the unit from micro to mili
  
  
  sinking_rate_forPlot = sinking_R*delta_x;
  printf("Depth_integrated_PP = %lf\n", depth_integ_PP);
  printf("t\tdepth\tsinking_rate\tPOC_flux\tPOC\tFL_B\tPA_B\ttotalB\tav_PAB\tdiff\n");
  
  POC_flux_forPlot = sinking_rate_forPlot*bac_carbon*cPperCB*v[jR[1]]/1000.0;
  printf("%.5lf\t%lf\t%.5lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n", t, 1*delta_x, sinking_rate_forPlot, POC_flux_forPlot, v[jR[1]], v[jP[1]], totalPA[1], v[jP[1]]+totalPA[1], totalPA[1]/v[jR[1]], diff_coeff(1));

  //for plot Pi distribution//
  //for(j = 0; j <= app_N; j++) printf("%.10lf\t%ld\t%.10lf\n", 1*delta_x, j, Pi(v[jP[1]], j, app_N));
  //printf("\n");
  //for plot Pi distribution//
	printf("%.5lf\t%lf\t%.5lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n", t, 5*delta_x, sinking_rate_forPlot, POC_flux_forPlot, v[jR[1]], v[jP[1]], totalPA[1], v[jP[1]]+totalPA[1], totalPA[1]/v[jR[1]], diff_coeff(1));
 
  
  //for plot Pi distribution//
  //for(j = 0; j <= app_N; j++) printf("%.10lf\t%ld\t%.10lf\n", 5*delta_x, j, Pi(v[jP[1]], j, app_N));
  //printf("\n");
  //for plot Pi distribution//

  for(k = 2; k <= NUM_LAYER; k++) {
  	POC_flux_forPlot = sinking_rate_forPlot*bac_carbon*cPperCB*v[jR[k]]/1000.0;	
    printf("%.5lf\t%lf\t%.5lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n", t, (k - 1 + MIXED_LAYER/delta_x)*delta_x, sinking_rate_forPlot, POC_flux_forPlot, v[jR[k]], v[jP[k]], totalPA[k], v[jP[k]]+totalPA[k], totalPA[k]/v[jR[k]], diff_coeff(k));

    //for plot Pi distribution//
    //for(j = 0; j <= app_N; j++) printf("%.10lf\t%ld\t%.10lf\n", (k -1 + MIXED_LAYER/delta_x)*delta_x, j, Pi(v[jP[k]], j, app_N));
    //printf("\n");
    //for plot Pi distribution//
  }//end of for k

  free_d_matrix(matrix_lambda); //must release memory space!!
  fclose(fp1);
  
  return 0;
}


void differential(double time, double in[], double out[])
{
  int j;
 #ifdef _OPENMP
 #pragma omp parallel for private (j)
 #endif
  for(j = 1; j <= NUM_LAYER; j++) {
    out[jR[j]] = dRkdt(time, in, j);
    out[jP[j]] = dPkdt(time, in, j);
  }
}

double dRkdt(double time, double vr[], int k)
{
  double temp1 = 0.0;
  //#ifdef AGG
  
  //double lammda_k;
  //double diffR;
  
  //if(k*delta_x <= MIXED_LAYER) diffR = diff0*1.0;
  //else diffR = diff0;
  //diffR = diff0;
  
  //lammda_k = (a0/d0)*exp(-0.5*(a1 + d1))*vr[jP[k]];
  
  
  //PA_predator_k = vr[jR[k]]*lammda_k*fAA(vr[jP[k]], appr_N_k); 
  
  if(k == 1) {
    temp1 = suppl(1) - consumption(vr[jR[1]], vr[jP[1]]);
    temp1 -= (diff_coeff(1) + sinking_R)*(delta_x/MIXED_LAYER)*vr[jR[1]];  //adjusted by volume
    temp1 += diff_coeff(1)*(delta_x/MIXED_LAYER)*vr[jR[2]];   //adjusted by volume
  }
  else if (k == NUM_LAYER) {
    temp1 = -consumption(vr[jR[NUM_LAYER]], vr[jP[NUM_LAYER]]);
    temp1 += (diff_coeff(NUM_LAYER - 1) + sinking_R)*vr[jR[NUM_LAYER -1]];
    temp1 -= sinking_R*vr[jR[NUM_LAYER]]; //no sinking at bottom or??
    temp1 += diff_coeff(NUM_LAYER)*depth_bottom_boundary_R;
    temp1 -= (diff_coeff(NUM_LAYER) + diff_coeff(NUM_LAYER - 1))*vr[jR[NUM_LAYER]]; //notice!!
  }
  else{
    temp1 = suppl(k);
    temp1 -= consumption(vr[jR[k]], vr[jP[k]]);
    temp1 -= (diff_coeff(k - 1) + diff_coeff(k) + sinking_R)*vr[jR[k]];
    temp1 += diff_coeff(k)*vr[jR[k + 1]];
    temp1 += (diff_coeff(k - 1) + sinking_R)*vr[jR[k - 1]];
  }
  
  //PA_predator_k = vr[jR[k]]*lammda_k*fAA(vr[jP[k]], appr_N_k); 
  //temp1 += recycling_C*mB*vr[jP[k]] /cPperCB;   //recycled POC from predator mortality
  
  //#endif
  
  return temp1;
}

double dPkdt(double time, double vr[], int k)
{
  double temp1;
  //double diffP;
  double g_mma;
  
  //if(k*delta_x <= MIXED_LAYER) diffP = diff0*1.0;
  //else diffP = diff0;
  //diffP = diff0;
  
#ifdef MONOD_ZERO
 
#endif
  
#ifdef AGG
  double lammda_k, lammda_km1, lammda_kp1;
  double myu_k, myu_km1, myu_kp1;

  double diff0_temp, sinking_P_temp;
  
  lammda_k = (a0/d0)*exp(-0.5*(a1 + d1))*vr[jP[k]];
  if(k > 1) lammda_km1 = (a0/d0)*exp(-0.5*(a1 + d1))*vr[jP[k - 1]];
  if(k < NUM_LAYER) lammda_kp1 = (a0/d0)*exp(-0.5*(a1 + d1))*vr[jP[k + 1]];
  
  myu_k = av_cont(matrix_lambda, interval_discontinuous, lammda_k);
  if(k > 1) myu_km1 = av_cont(matrix_lambda, interval_discontinuous, lammda_km1);
  if(k < NUM_LAYER) myu_kp1 = av_cont(matrix_lambda, interval_discontinuous, lammda_kp1);
  

  if(k == 1) {
    temp1 = 0.0;
    temp1 += 0.0;
    temp1 += diff_coeff(k)*(vr[jP[k + 1]] - vr[jP[k]] + (myu_kp1 - myu_k)*vr[jR[k + 1]]);
    
  }
  else if (k == NUM_LAYER) {
    temp1 = sinking_R*(myu_km1 - myu_k)*vr[jR[k - 1]];
    temp1 += diff_coeff(k - 1)*(vr[jP[k - 1]] - vr[jP[k]] + (myu_km1 - myu_k)*vr[jR[k - 1]]);
    temp1 += diff_coeff(k)*(depth_bottom_boundary_P - vr[jP[k]] + (myu_kp1 - myu_k)*depth_bottom_boundary_R);
  }
  else{
    temp1 = sinking_R*(myu_km1 - myu_k)*vr[jR[k - 1]];
    temp1 += diff_coeff(k - 1)*(vr[jP[k - 1]] - vr[jP[k]] + (myu_km1 - myu_k)*vr[jR[k - 1]]);
    temp1 += diff_coeff(k)*(vr[jP[k + 1]] - vr[jP[k]] + (myu_kp1 - myu_k)*vr[jR[k + 1]]);
  }
  
  //normalizing spatial components only
  temp1 /= (1.0 + (a0/d0)*exp(-0.5*(a1 + d1))*vr[jR[k]]*av_cont_diff(matrix_lambda, interval_discontinuous, interval_differential, lammda_k)); 
  
  temp1 += growthP(vr[jR[k]], vr[jP[k]], k);
  
#endif		
  
  return temp1;
}

double consumption(double resource, double predator)
{
  double temp, lambda_ap;
  
  //normal MONOD equation for B_T and A_T
#ifdef  MONOD_ZERO
  temp = h1*(a0/d0)*exp(-0.5*(a1 + d1))*resource*predator/(1.0 + (a0/d0)*exp(-0.5*(a1 + d1))*resource) + mR*resource;
#endif	
  
#ifdef MONOD	
  temp = Vmax*resource*predator/(Kh + resource) + mR*resource;
#endif
  
#ifdef AGG
  lambda_ap = (a0/d0)*exp(-0.5*(a1 + d1))*predator;
  temp = h1*av_cont(matrix_lambda, interval_discontinuous, lambda_ap)*resource + mR*resource;
#endif
  
  return temp;
}

double growthP(double resource, double predator, int k)
{
  double temp, lambda_ap;
  double gB, myu;
  
  lambda_ap = (a0/d0)*exp(-0.5*(a1 + d1))*predator;
  myu = av_cont(matrix_lambda, interval_discontinuous, lambda_ap);
  gB = (cPperCB*bge + myu)*h1;
  
  //normal MONOD equation for B_T and A_T
#ifdef MONOD_ZERO
  temp = (gB - mR)*(a0/d0)*exp(-0.5*(a1 + d1))*resource/(1.0 + (a0/d0)*exp(-0.5*(a1 + d1))*resource) - mB*pow(predator, dellta2);
  temp *= predator;
#endif	
  
#ifdef MONOD	
  temp = bge*cPperCB*Vmax*resource*predator/(Kh + resource) - mB*predator;
#endif
  
#ifdef AGG
  temp = 1.0/(1.0 + (a0/d0)*exp(-0.5*(a1 + d1))*resource*av_cont_diff(matrix_lambda, interval_discontinuous, interval_differential, lambda_ap));
  temp *= (gB - mB)*myu*resource - (mB*predator + myu*suppl(k));
  
  

  //when density-dependent mortality of attached bacteria is considered
  if(dellta > 0.0) {
    //temp = (a0/d0)*exp(-0.5*(a1 + d1))*f_AA*((gB + h1*l_amda(predator)*f_AA)*resource - suppl(k))*predator;
    //temp -= mB*(pow(predator, 1.0 + dellta2) + fMB(predator, appN)*resource);//mortality of FL and PA bacteria
  }

  else{
    //temp = (a0/d0)*exp(-0.5*(a1 + d1))*f_AA*((gB + h1*l_amda(predator)*f_AA - mB)*resource - suppl(k)) - mB*pow(predator, dellta2);
    //temp *=predator;
  }
  //do not normalized by inverse term
#endif	
  
  return temp;	
}

double suppl(int k)
{
  double temp;
  if (k == 1) temp = s0;  //note that actual NPP/m2 is s0*MIXED_LAYER!!
  else temp = 0.0;
  //temp = s0;	//homogeneous
  return temp;
}

double diff_coeff(int k)
{
  double temp;
  
  if(k == 1) temp = diff0*10.0;
  else if((k + 4)*delta_x < MIXED_LAYER + 50.0) temp = diff0*(10.0 - 1.8*(k + 4  - MIXED_LAYER/delta_x));
  else temp = diff0*1.0;  //from 100m 1.0cm2/s
  
  //temp = diff0*10.0; //homogeneous
  
  return temp;
}

double av_cont(double **discont_data, double interval, double lamd) {
  int index_lam;
  double lam_down, lam_up, lam_diff;
  double av_lam_down, av_lam_up, av_lam_approx;
  
  if(lamd > NUM_R*interval) return(-1.0);
  
  index_lam = (int) lamd/interval;
  lam_down = index_lam*interval;
  lam_diff = lamd - lam_down;
  lam_up = lam_down + interval;
  
  if(index_lam == 0) {
    av_lam_down = 0.0;
    av_lam_up = discont_data[index_lam + 1][2];
  } else {
    av_lam_down = discont_data[index_lam][2];
    av_lam_up = discont_data[index_lam + 1][2];
  }
  
  //linear interpolation
  av_lam_approx = av_lam_down+ (av_lam_up - av_lam_down)*lam_diff/interval;
  
  return(av_lam_approx);
  
}

double av_cont_diff(double **av_data, double interval_av, double interval_diff, double lamd){
  double av_lam_diff_approx;
  av_lam_diff_approx = (av_cont(av_data, interval_av, lamd + interval_diff) - av_cont(av_data, interval_av, lamd))/interval_diff;
  
  return(av_lam_diff_approx);
  
} 

	
	
