//from 2015/01/08, copied from pp_ratio03.c

#ifdef CHECK_FUNC

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
//#include "./global_params_pp2020.c"
#endif

#define ZERO 1.0e-20

double Pzero(double free, long n);
double Pi(double free, long i, long n);
double theta_n(double free, long n);
double phi_n(double free, long n);
double psi_n(double free, long n);
double one_p_rBArBF(double host, double free, long n);

double l_amda(double free);

double fMB(double free, long n); //density dependent mortality of attached bacteria

double suppl(int k); //supply rate 
//double dfAdlammda(double free, int n);
//double dPAdF(double resource, double free, int n);


#ifdef CHECK_FUNC
int main(int arg, char *argv[])
{	
  int k;
  int m = 5;
  int max_n = 50;
  double hostT = 10.0;
  double lammda;
  double free_b = 0.20;
  
  d1 = 0.05*1.0;
  a1 = 0.01*1.0;
  //alfa = 0.5*(d1 - a1);
  
  lammda = (a0/d0)*free_b;
  
  if(lammda*2.0 > 50) max_n = lammda*2.0;
  else max_n = 50;
  
  //for(k = 0; k <= 1000; k++) printf("%d\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n", k, Pzero(4.0, k), fA(4.0, k), dfAdlammda(4.0, k), dPAdF(1.0, 4.0, k));
  //for(k = 0; k <= 10000; k++) printf("%.10lf\t%.10lf\n", 1.0*k, dPAdF(1.0, 1.0*k, 100));
  for(k = 0; k <= max_n; k++) printf("%d\t%.10lf\n", k, Pi(free_b, k, max_n));
  printf("Pzero is\t%.50lf\n", Pzero(free_b, max_n));
  printf("theta_n is\t%.10lf\n", theta_n(free_b, max_n));
  printf("phi_n is\t%.10lf\n", phi_n(free_b, max_n));
  printf("psi_n is\t%.10lf\n", psi_n(free_b, max_n));
  printf("one_p_rBArBF is\t%.10lf\n", one_p_rBArBF(hostT, free_b, max_n));
  return 0;
  
}
#endif

double Pzero(double free, long n)
{	//numerical approximation for Po(alpha, lamda) with order n
  double temp;
  
  
  temp = 1.0/(1.0 + psi_n(free, n));
  return temp;
}

double Pi(double free, long i, long n)
{	//numerical approximaiton for Pi
  double temp1, temp2;
  double lammda;
  long j, k;
  
  lammda = (a0/d0)*free;
  if (i > n) temp1 = 0.0;
  
  else {
    temp1 = Pzero(free, n);
//#ifdef _OPENMP
//#pragma omp parallel for private (j)
//#endif
    for(j = 1; j <= i; j++) {
      temp1 *= (lammda/j)*(1.0 + (a1/a0)*(j - 1))/(1.0 + (d1/d0)*j);
    }//end of for j
  }//end of else
  
  return temp1;
  
}

double theta_n(double free, long n)
{	//numerical approximation for Po(alpha, lamda) with order n
  double temp1, temp2, temp3=0.0;
  double lammda;
  long i, j;
  
  lammda = (a0/d0)*free;
  
  for(i = 1; i <= n; i++) {
    temp1 = 1.0;
//#ifdef _OPENMP
//#pragma omp parallel for private (j)
//#endif
    for(j = 1; j <= i; j++) {
      temp2 = (lammda/j)*(1.0 + (a1/a0)*(j - 1))/(1.0 + (d1/d0)*j);
      temp1 *= temp2;
    }//end of for j
    temp1 *= i;
    temp3 += temp1;
  }//end of for i

  return temp3;
}

double phi_n(double free, long n)
{	//numerical approximation for Po(alpha, lamda) with order n
  double temp1, temp2, temp3=0.0;
  double lammda;
  long i, j;
  
  lammda = (a0/d0)*free;

  for(i = 1; i <= n; i++) {
    temp1 = 1.0;

    for(j = 1; j <= i; j++) {
      temp2 = (lammda/j)*(1.0 + (a1/a0)*(j - 1))/(1.0 + (d1/d0)*j);
      temp1 *= temp2;
    }//end of for j
    temp1 *= i*i;
    temp1 /= lammda;
    temp3 += temp1;
  }//end of for i
  
  return temp3;
}

double psi_n(double free, long n)
{	//numerical approximation for Po(alpha, lamda) with order n
  double temp1, temp2, temp3=0.0;
  double lammda;
  long i, j;

  lammda = (a0/d0)*free;
  
  for(i = 1; i <= n; i++) {
    temp1 = 1.0;

    for(j = 1; j <= i; j++) {
      temp2 = (lammda/j)*(1.0 + (a1/a0)*(j - 1))/(1.0 + (d1/d0)*j);
      temp1 *= temp2;
    }//end of for j
    
    temp3 += temp1;
  }//end of for i
  
  return temp3;
}

double one_p_rBArBF(double host, double free, long n)
{
  double temp;
  double lammda;
  
  lammda = (a0/d0)*free;
  
  temp = 1.0 - host*(a0/d0)*(1.0/lammda)*(theta_n(free, n)/(1.0 + psi_n(free,n)))*(theta_n(free, n)/(1.0 + psi_n(free,n)));
  temp += host*(a0/d0)*phi_n(free,n)/(1.0 + psi_n(free,n));
                                                                            
  return temp;  
}



double l_amda(double free)
{
  return (a0/d0)*free;
}


