//from 2015/01/08, copied from pp_ratio03.c

#ifdef CHECK_FUNC

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
//#include "./global_params_pp2015.c"
#endif

#define ZERO 1.0e-20

double fa(long n); //density dependent attachment rate
double fd(long n); //density dependent detachment rate
double Pzero(double free, long n);
double Psub(double free, long n);
double Pi(double free, long i, long n);
double av_n(double free, long n); //average of n
double av_n2(double free, long n);//average of n2
double rBArBF(double free, long n); //partial derivative of B_A with respect to B_F
double rBArAT(double free, long n);//

double fA(double free, long n);  //in fact not necessary
double fB(double resource, double free, long n);
double fF(double resource, double free, long n);  //in fact not necessary
double fAA(double free, long n);
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
  int max_n = 500;
  double lammda;
  double free_b;
  
  d1 = 0.05;
  a1 = 0.01;
  alfa = 0.5*(d1 - a1);
  free_b = 2.0;
  
  lammda = (a0/d0)*exp(-0.5*(a1 + d1))*free_b;
  
  if(lammda*2.0 > 50) max_n = lammda*2.0;
  else max_n = 50;
  
  //for(k = 0; k <= 1000; k++) printf("%d\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n", k, Pzero(4.0, k), fA(4.0, k), dfAdlammda(4.0, k), dPAdF(1.0, 4.0, k));
  //for(k = 0; k <= 10000; k++) printf("%.10lf\t%.10lf\n", 1.0*k, dPAdF(1.0, 1.0*k, 100));
  for(k = 0; k <= max_n; k++) printf("%d\t%.10lf\n", k, Pi(free_b, k, max_n));
  
  return 0;
  
}
#endif

double Pzero(double free, long n)
{	//numerical approximation for Po(alpha, lamda) with order n
  double temp1, temp2, temp3;
  double lammda;
  long i, j;
  
  lammda = (a0/d0)*exp(-0.5*(a1 + d1))*free;
  
  if (lammda < ZERO) return 1.0;		//special case with free = 0
  
  temp2 = 1.0;	//when i = 0, i! = 1
  for(i = 1; i <= n; i++) {
    temp1 = 1.0;
    for(j = 1; j <= i; j++) {
      temp1 *= (lammda*exp(-alfa*i)/j);
    }//end of for j
    temp2 += temp1;
  }//end of for i
  
  temp3 = 1.0/temp2;
  return temp3;
}

double Pi(double free, long i, long n)
{	//numerical approximaiton for Pi
  double temp1, temp2;
  double lammda;
  long j, k;
  
  lammda = (a0/d0)*exp(-0.5*(a1 + d1))*free;
  if (i > n) temp1 = 0.0;
  
  else {
    temp1 = Pzero(free, n);
    for(j = 1; j <= i; j++) {
      temp1 *= (lammda*exp(-alfa*i)/j);
    }//end of for j
  }//end of else
  
  return temp1;
  
}

double l_amda(double free)
{
  return (a0/d0)*exp(-0.5*(a1 + d1))*free;
}

double fAA(double free, long n)
{                    //numerical approximation ofr fAA(BF) with order n
  double temp1, temp2;
  double lammda, p0;
  long i, j;

  lammda = (a0/d0)*exp(-0.5*(a1 + d1))*free;
  p0 = Pzero(free, n);
  temp1 = 1.0*exp(-alfa)*p0; //for i = 0, i! = 1

  //#ifdef _OPENMP
  //#pragma omp parallel for private (i, j, temp2)
  //#endif
  for(i = 1; i <= n; i++) {
    temp2 = 1.0;
    for(j = 1; j <= i; j++) temp2 *= (lammda*exp(-alfa*i-2.0*alfa)/j);
    temp2 *= (p0*exp(-alfa));
    temp1 += temp2;
  }//end of for i
  //temp1 *= Pzero(free, n);

  return temp1;
}

double fB(double resource, double free, long n)
{
  double temp1, temp2, temp3, temp4;
  double lammda;
  long i,j;
  double p_zero, f_AA;

  lammda = (a0/d0)*exp(-0.5*(a1 + d1))*free;
  p_zero = Pzero(free, n);
  f_AA = fAA(free, n);

  temp1 = 1.0*exp(-alfa)*p_zero; //for i = 0 i! = 0

  //#ifdef _OPENMP
  //#pragma omp parallel for private (i, j, temp2)
  //#endif
  for(i = 1; i <= n; i++) {
    temp2 = 1.0;
    for(j = 1; j <= i; j++) temp2 *= (lammda*exp(-alfa*i-2.0*alfa)/j);
    temp2 *= (exp(-alfa)*(i + 1.0));
    temp1 += (temp2*p_zero);
  } //end of for i

  temp3 = temp1 - lammda*f_AA*f_AA;

  temp4 = 1.0 + resource*(a0/d0)*exp(-0.5*(a1 + d1))*temp3;

  return 1.0/temp4;
}
  
double fMB(double free, long n)
{
  //density dependent mortality of attached bacteria
  double temp1, temp2;
  double lammda, p_zero;
  long i, j;

  lammda = l_amda(free);
  p_zero = Pzero(free, n);

  temp1 = p_zero;
  for(i = 1; i <= n; i++) {
    temp2 = 1.0;
    for(j = 1; j <= i; j++) temp2*= (lammda*exp(-alfa*i + dellta - 2.0*alfa)/j);
    temp1 += (temp2*p_zero);
  }//end of for i
  temp1 *= (lammda*exp(dellta - alfa));

  return temp1;
}
