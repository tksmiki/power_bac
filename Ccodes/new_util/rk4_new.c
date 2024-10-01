#include "nrutil_new.h"

void rk4(double y[], double k1[], int n, double t, double h, double yout[], void (*derivs) (double, double [], double []))
{
	int i;
	double th, hh, h6, *k2, *k3, *k4, *yt;
	
	k2 = d_vector(n);
	k3 = d_vector(n);
	k4 =d_vector(n);
	yt = d_vector(n);
	hh = h*0.5;
	h6 = h/6.0;
	th = t + hh;
	for (i=1;i<=n;i++) yt[i] = y[i] + hh*k1[i];
	(*derivs)(th, yt, k2);
	for (i=1;i<=n;i++) yt[i] = y[i] + hh*k2[i];
	(*derivs)(th, yt, k3);
	for (i=1;i<=n;i++) yt[i] = y[i] + h*k3[i];
	(*derivs)(t+h, yt, k4);
	for (i=1;i<=n;i++) yout[i] = y[i] + h6*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
	

	free_d_vector(yt);
	free_d_vector(k2);
	free_d_vector(k3);
	free_d_vector(k4);
}

