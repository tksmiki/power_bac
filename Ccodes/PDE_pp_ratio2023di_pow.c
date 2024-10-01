//automatic parallelization
// -parallel -par-report3 -par-threshold100 

//for spatiall constant sinking rate, copied from PDE_pp_ratio03new.c at 2015/01/08

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#include "../new_util/nrutil_new.h"
#include "../new_util/forode.h"
#include "../new_util/nrutil_new.c"
#include "../new_util/rk4_new.c"

#include "./global_params_pp2020.c"
#include "./pp_ratio2020di_pow.c"

//For OPENMP, -fopenmp (gcc OSX) or -openmp (icc Linux)

#define DELTAZERO 1.0e-1  //0.1 for AGG but 0.01 for MONODZERO
#define EPSS 1.0e-6
#define ZERO 1.0e-20

#define NUM_LAYER 297  //number of layer
#define MIXED_LAYER  50.0 //[m]
//#define min_app_N 200

#define INTERVAL_DAY 30 //interval for printf

void differential(double time, double in[], double out[]);
double dRkdt(double time, double vr[], int k);			//for resource at k-th layer
double dPkdt(double time, double vr[], int k);			//for predator (free-living consumer) at k-th layer
double consumption(double resource, double predator, long appN);
double growthP(double resource, double predator, long appN, int k);
//double suppl(int k); //depth specific supply rate
double diff_coeff(int k); //depth specific diffusion coefficient

int jP[NUM_LAYER + 1], jR[NUM_LAYER +1];

long app_N;
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
  
  
  
  tzero = 0.0;
  end_time = 365.0*20.0;
  //end_time = 365*2.0;
  //end_time = 2.0;

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
  
  //a1 = 0.0;  //for density_independent behavior
  //d1 = 0.0;  //for density_independent behavior
  //d1 = 1.0;
  //a1 = 0.9;
  alfa = 0.5*(d1 - a1);		//for density_dependent behavior
  
  //setting for global parameters
  min_app_N = 20;  //this is critical parameter for numerical accuracy, at least for MONOD_ZERO
  //assignment of vector index
  
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
    
    //dynamical allocation of app_N value at each time step, which is common between depth (not used)
    lammda = 0.0;
    for(k = 1; k <= NUM_LAYER; k++) {
      if(lammda < (a0/d0)*exp(-0.5*(a1 + d1))*v[jP[k]]) lammda = (a0/d0)*exp(-0.5*(a1 + d1))*v[jP[k]];
    }
    if(lammda*2.0 > min_app_N) app_N = lammda*2.0;
    else app_N = min_app_N;
    
    
    //setting boundary conditions (for reflecting boundary////
    //depth_zero_boundary_R = v[jR[2]];  //for reflecting boundary
    //depth_zero_boundary_P = v[jP[2]];  //for reflecting boundary
    depth_bottom_boundary_R = v[jR[NUM_LAYER - 1]]; //for reflecting boundary
    depth_bottom_boundary_P = v[jP[NUM_LAYER - 1]]; //for reflecting boundary
    
#ifdef AGG		
    lammda = l_amda(v[jP[2]]);
    if(lammda*2.0 > min_app_N) appNk = lammda*2.0;
    else appNk = min_app_N;			
    totalPA[2] = (a0/d0)*exp(-0.5*(a1 + d1))*v[jP[2]]*v[jR[2]]*fAA(v[jP[2]], appNk);
    lammda = l_amda(v[jP[NUM_LAYER - 1]]);
    if(lammda*2.0 > min_app_N) appNk = lammda*2.0;
    else appNk = min_app_N;			
    totalPA[NUM_LAYER - 1] = (a0/d0)*exp(-0.5*(a1 + d1))*v[jP[NUM_LAYER - 1]]*v[jR[NUM_LAYER - 1]]*fAA(v[jP[NUM_LAYER - 1]], appNk);
    
    //depth_zero_boundary_PA = totalPA[2];
    depth_bottom_boundary_PA = totalPA[NUM_LAYER - 1];	
#endif		
    //end of setting boundary conditions/////	
    
    
    differential(t, v, dfdt);
    rk4(v, dfdt, n_variable, t, deltat, v, differential);
    t += deltat;
    //v[jP[1]] = 10.0;  //for fixed boundary condition
    //v[jR[1]] = 1.0;  //for fixed boundary condition
    
    //for(k = 1; k <= n_variable; k++) v[k] = v[k];  //renewal of vector
    
 
    if(write_index < INTERVAL_DAY/deltat) write_index++;
    else {
      totalB = 0.0;
      totalR = 0.0;
      for(k = 1; k <= NUM_LAYER; k++) {
#ifdef AGG	
	lammda = l_amda(v[jP[k]]);
	if(lammda*2.0 > min_app_N) appNk = lammda*2.0;
	else appNk = min_app_N;			
	totalPA[k] = (a0/d0)*exp(-0.5*(a1 + d1))*v[jP[k]]*v[jR[k]]*fAA(v[jP[k]], appNk);
	totalB += v[jP[k]] + totalPA[k];
	totalR += v[jR[k]];
#endif	
      }// end of for k
      
      printf("%.10lf\t%lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n", t, 1*delta_x, v[jR[1]], v[jP[1]], totalPA[1], v[jP[1]]+totalPA[1], diff_coeff(1));
      printf("%.10lf\t%lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n", t, 5*delta_x, v[jR[1]], v[jP[1]], totalPA[1], v[jP[1]]+totalPA[1], diff_coeff(1));
      
      for(k = 2; k <= NUM_LAYER; k++) {
	printf("%.10lf\t%lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n", t, (k - 1 + MIXED_LAYER/delta_x)*delta_x, v[jR[k]], v[jP[k]], totalPA[k], v[jP[k]]+totalPA[k], diff_coeff(k));
      }//end of for k

      //plot of statial integrated values
      //printf("%.10lf\t%.10lf\t%.10lf\n", t, totalB, totalR); 
      //printf("%lf\n", t);
      printf("\n");	
      write_index = 1;
    }
        			
  }//end of for t
  
#endif
  for(k = 1; k <= NUM_LAYER; k++) {
#ifdef AGG	
	lammda = l_amda(v[jP[k]]);
	if(lammda*2.0 > min_app_N) appNk = lammda*2.0;
	else appNk = min_app_N;			
	totalPA[k] = (a0/d0)*exp(-0.5*(a1 + d1))*v[jP[k]]*v[jR[k]]*fAA(v[jP[k]], appNk);	
#endif
  }		
  depth_integ_PP = bac_carbon*cPperCB*s0*MIXED_LAYER/1000.0; //change the unit from micro to milLi
  
  
  sinking_rate_forPlot = sinking_R*delta_x;
  printf("Depth_integrated_PP = %lf\n", depth_integ_PP);
  printf("t\tdepth\tsinking_rate\tPOC_flux\tPOC\tFL_B\tPA_B\ttotalB\tdiff\n");
  
  POC_flux_forPlot = sinking_rate_forPlot*bac_carbon*cPperCB*v[jR[1]]/1000.0;
  printf("%.10lf\t%lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n", t, 1*delta_x, sinking_rate_forPlot, POC_flux_forPlot, v[jR[1]], v[jP[1]], totalPA[1], v[jP[1]]+totalPA[1], diff_coeff(1));

  //for plot Pi distribution//
  //for(j = 0; j <= app_N; j++) printf("%.10lf\t%ld\t%.10lf\n", 1*delta_x, j, Pi(v[jP[1]], j, app_N));
  //printf("\n");
  //for plot Pi distribution//
	printf("%.10lf\t%lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n", t, 5*delta_x, sinking_rate_forPlot, POC_flux_forPlot, v[jR[1]], v[jP[1]], totalPA[1], v[jP[1]]+totalPA[1], diff_coeff(1));
 
  
  //for plot Pi distribution//
  //for(j = 0; j <= app_N; j++) printf("%.10lf\t%ld\t%.10lf\n", 5*delta_x, j, Pi(v[jP[1]], j, app_N));
  //printf("\n");
  //for plot Pi distribution//

  for(k = 2; k <= NUM_LAYER; k++) {
  	POC_flux_forPlot = sinking_rate_forPlot*bac_carbon*cPperCB*v[jR[k]]/1000.0;	
    printf("%.10lf\t%lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n", t, (k - 1 + MIXED_LAYER/delta_x)*delta_x, sinking_rate_forPlot, POC_flux_forPlot, v[jR[k]], v[jP[k]], totalPA[k], v[jP[k]]+totalPA[k], diff_coeff(k));

    //for plot Pi distritbution//
    //for(j = 0; j <= app_N; j++) printf("%.10lf\t%ld\t%.10lf\n", (k -1 + MIXED_LAYER/delta_x)*delta_x, j, Pi(v[jP[k]], j, app_N));
    //printf("\n");
    //for plot Pi distribution//
  }//end of for k

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
  
  double lammda_k;
  long appr_N_k;
  
  double diffR;
  
  double PA_predator_k;
  
  //if(k*delta_x <= MIXED_LAYER) diffR = diff0*1.0;
  //else diffR = diff0;
  diffR = diff0;
  
  lammda_k = l_amda(vr[jP[k]]);
  
  if(lammda_k*2.0 > min_app_N) appr_N_k = lammda_k*2.0;
  else appr_N_k = min_app_N;
  
  //PA_predator_k = vr[jR[k]]*lammda_k*fAA(vr[jP[k]], appr_N_k); 
  
  if(k == 1) {
    temp1 = suppl(1) - consumption(vr[jR[1]], vr[jP[1]], appr_N_k);
    temp1 -= (diff_coeff(1) + sinking_R)*(delta_x/MIXED_LAYER)*vr[jR[1]];  //adjusted by volume
    temp1 += diff_coeff(1)*(delta_x/MIXED_LAYER)*vr[jR[2]];   //adjusted by volume
  }
  else if (k == NUM_LAYER) {
    temp1 = -consumption(vr[jR[NUM_LAYER]], vr[jP[NUM_LAYER]], appr_N_k);
    temp1 += (diff_coeff(NUM_LAYER - 1) + sinking_R)*vr[jR[NUM_LAYER -1]];
    temp1 -= sinking_R*vr[jR[NUM_LAYER]]; //no sinking at bottom or??
    temp1 += diff_coeff(NUM_LAYER)*depth_bottom_boundary_R;
    temp1 -= (diff_coeff(NUM_LAYER) + diff_coeff(NUM_LAYER - 1))*vr[jR[NUM_LAYER]]; //notice!!
  }
  else{
    temp1 = suppl(k);
    temp1 -= consumption(vr[jR[k]], vr[jP[k]], appr_N_k);
    temp1 -= (diff_coeff(k - 1) + diff_coeff(k) + sinking_R)*vr[jR[k]];
    temp1 += diff_coeff(k)*vr[jR[k + 1]];
    temp1 += (diff_coeff(k - 1) + sinking_R)*vr[jR[k - 1]];
  }
  
  //PA_predator_k = vr[jR[k]]*lammda_k*fAA(vr[jP[k]], appr_N_k); 
  temp1 += recycling_C*mB*vr[jP[k]] /cPperCB;   //recycled POC from predator mortality
  
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
  g_mma = (a0/d0)*exp(-0.5*(a1 + d1));	
  if(k == 1) {
    temp1 =  growthP(vr[jR[1]], vr[jP[1]], 0, k);
    temp1 -= diff_coeff(1)*(delta_x/MIXED_LAYER)*vr[jP[1]];  //adjusted by volume
    temp1 += diff_coeff(1)*(delta_x/MIXED_LAYER)*vr[jP[2]];  //adjusted by volume
    temp1 -= sinking_R*(delta_x/MIXED_LAYER)*vr[jP[1]]*g_mma*vr[jR[1]]/(1.0 + g_mma*vr[jR[1]]);	//adjusted by volume & only attached fraction sinks
  }
  else if (k == NUM_LAYER) {
    temp1 = growthP(vr[jR[NUM_LAYER]], vr[jP[NUM_LAYER]], 0, k);
    temp1 += diff_coeff(NUM_LAYER - 1)*vr[jP[NUM_LAYER - 1]];
    temp1 += diff_coeff(NUM_LAYER)*depth_bottom_boundary_P;
    temp1 -= (diff_coeff(NUM_LAYER) + diff_coeff(NUM_LAYER - 1))*vr[jP[NUM_LAYER]];
    temp1 -= sinking_R*vr[jP[NUM_LAYER]]*g_mma*vr[jR[NUM_LAYER]]/(1.0 + g_mma*vr[jR[NUM_LAYER]]); //adjusted
    temp1 += sinking_R*vr[jP[NUM_LAYER - 1]]*g_mma*vr[jR[NUM_LAYER - 1]]/(1.0 + g_mma*vr[jR[NUM_LAYER - 1]]);
    
  }
  else{
    temp1  = growthP(vr[jR[k]], vr[jP[k]], 0, k);
    temp1 -= (diff_coeff(k) + diff_coeff(k - 1))*vr[jP[k]];
    temp1 += diff_coeff(k)*vr[jP[k + 1]];
    temp1 += diff_coeff(k - 1)*vr[jP[k - 1]];
    temp1 -= sinking_R*vr[jP[k]]*g_mma*vr[jR[k]]/(1.0 + g_mma*vr[jR[k]]);
    temp1 += sinking_R*vr[jP[k - 1]]*g_mma*vr[jR[k - 1]]/(1.0 + g_mma*vr[jR[k - 1]]);
  }
#endif
  
#ifdef AGG
  double lammda_k, lammda_2, lammda_km1, lammda_kp1;
  double PA_predator_k, PA_predator_2, PA_predator_km1, PA_predator_kp1;
  long appr_N_k, appr_N_2, appr_N_km1, appr_N_kp1;
	
  double fAA_k;
  
  double diff0_temp, sinking_P_temp;
  
  lammda_k = (a0/d0)*exp(-0.5*(a1 + d1))*vr[jP[k]];
  //lammda_2 = (a0/d0)*exp(-0.5*(a1 + d1))*vr[jP[2]];
  if(k > 1) lammda_km1 = (a0/d0)*exp(-0.5*(a1 + d1))*vr[jP[k - 1]];
  if(k < NUM_LAYER) lammda_kp1 = (a0/d0)*exp(-0.5*(a1 + d1))*vr[jP[k + 1]];
  
  if(lammda_k*2.0 > min_app_N) appr_N_k = lammda_k*2.0;
  else appr_N_k = min_app_N;
  
  //if(lammda_2*2.0 > min_app_N) appr_N_2 = lammda_2*2.0;
  //else appr_N_2 = min_app_N;
  
  if(k > 1) {
  if(lammda_km1*2.0 > min_app_N) appr_N_km1 = lammda_km1*2.0;
  else appr_N_km1 = min_app_N;
  }
  
  if(k < NUM_LAYER) {
  if(lammda_kp1*2.0 > min_app_N) appr_N_kp1 = lammda_kp1*2.0;
  else appr_N_kp1 = min_app_N;
  }
  fAA_k = fAA(vr[jP[k]], appr_N_k);
  PA_predator_k = vr[jR[k]]*lammda_k*fAA_k;
  //PA_predator_2 = vr[jR[2]]*lammda_2*fAA(vr[jP[2]], appr_N_2);
  if(k > 1) PA_predator_km1 = vr[jR[k - 1]]*lammda_km1*fAA(vr[jP[k - 1]], appr_N_km1);
  if(k < NUM_LAYER) PA_predator_kp1 = vr[jR[k + 1]]*lammda_kp1*fAA(vr[jP[k + 1]], appr_N_kp1);
  
  if(k == 1) {
    temp1 =  growthP(vr[jR[1]], vr[jP[1]], appr_N_k, k);
    temp1 -= diff_coeff(1)*(delta_x/MIXED_LAYER)*vr[jP[1]]; //adjusted by volume
    temp1 += diff_coeff(1)*(delta_x/MIXED_LAYER)*vr[jP[2]]; //adjusted by volume	
    temp1 += diff_coeff(k)*(delta_x/MIXED_LAYER)*(PA_predator_kp1 - vr[jR[k + 1]]*fAA_k*lammda_k); //adjusted by volume	
    
  }
  else if (k == NUM_LAYER) {
    temp1 = growthP(vr[jR[NUM_LAYER]], vr[jP[NUM_LAYER]], appr_N_k, k);
    temp1 += diff_coeff(NUM_LAYER - 1)*vr[jP[NUM_LAYER - 1]];
    temp1 += diff_coeff(NUM_LAYER)*depth_bottom_boundary_P;
    temp1 -= (diff_coeff(NUM_LAYER) + diff_coeff(NUM_LAYER - 1))*vr[jP[NUM_LAYER]];
    temp1 += sinking_R*(PA_predator_km1 - vr[jR[k - 1]]*fAA_k*lammda_k);
    temp1 += diff_coeff(k - 1)*(PA_predator_km1 - vr[jR[k - 1]]*fAA_k*lammda_k);	
    temp1 += diff_coeff(k)*(depth_bottom_boundary_PA - depth_bottom_boundary_R*fAA_k*lammda_k);	
  }
  else{
    temp1  = growthP(vr[jR[k]], vr[jP[k]], appr_N_k, k);
    temp1 -= (diff_coeff(k - 1) + diff_coeff(k))*vr[jP[k]];
    temp1 += diff_coeff(k)*vr[jP[k + 1]];
    temp1 += diff_coeff(k - 1)*vr[jP[k - 1]];
    temp1 += sinking_R*(PA_predator_km1 - vr[jR[k - 1]]*fAA_k*lammda_k);
    temp1 += diff_coeff(k - 1)*(PA_predator_km1 - vr[jR[k - 1]]*fAA_k*lammda_k);	
    temp1 += diff_coeff(k)*(PA_predator_kp1 - vr[jR[k + 1]]*fAA_k*lammda_k);	
  }
  
  temp1 *= fB(vr[jR[k]], vr[jP[k]], appr_N_k);//normalized by fB
#endif		
  
  return temp1;
}

double consumption(double resource, double predator, long appN)
{
  double temp;
  
  //normal MONOD equation for B_T and A_T
#ifdef  MONOD_ZERO
  temp = h1*(a0/d0)*exp(-0.5*(a1 + d1))*resource*predator/(1.0 + (a0/d0)*exp(-0.5*(a1 + d1))*resource) + mR*resource;
#endif	
  
#ifdef MONOD	
  temp = Vmax*resource*predator/(Kh + resource) + mR*resource;
#endif
  
#ifdef AGG
  temp = h1*l_amda(predator)*fAA(predator, appN)*resource + mR*resource;
#endif
  
  return temp;
}

double growthP(double resource, double predator, long appN, int k)
{
  double temp;
  double gB;
  double f_AA;
  
  gB = h1*cPperCB*bge;
  
  //normal MONOD equation for B_T and A_T
#ifdef MONOD_ZERO
  temp = (gB - mR)*(a0/d0)*exp(-0.5*(a1 + d1))*resource/(1.0 + (a0/d0)*exp(-0.5*(a1 + d1))*resource) - mB*pow(predator, dellta2);
  temp *= predator;
#endif	
  
#ifdef MONOD	
  temp = bge*cPperCB*Vmax*resource*predator/(Kh + resource) - mB*predator;
#endif
  
#ifdef AGG
  f_AA = fAA(predator, appN);

  //when density-dependent mortality of attached bacteria is considered
  if(dellta > 0.0) {
    temp = (a0/d0)*exp(-0.5*(a1 + d1))*f_AA*((gB + h1*l_amda(predator)*f_AA)*resource - suppl(k))*predator;
    temp -= mB*(pow(predator, 1.0 + dellta2) + fMB(predator, appN)*resource);//mortality of FL and PA bacteria
  }

  else{
    temp = (a0/d0)*exp(-0.5*(a1 + d1))*f_AA*((gB + h1*l_amda(predator)*f_AA - mB)*resource - suppl(k)) - mB*pow(predator, dellta2);
    temp *=predator;
  }
  //do not normalized by inverse term
#endif	
  
  return temp;	
}

double suppl(int k)
{
  double temp;
  if (k == 1) temp = s0;  //note that actural NPP/m2 is s0*MIXED_LAYER!!
  else temp = 0.0;
  //temp = s0;	//homogenenous
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


	
	
