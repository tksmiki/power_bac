//list of biological parameters
//from 2015/01/08, copied from global_params_pp2015.c

double bac_carbon = 20.0; // micro gC per 1.0e+9 cells/m3 unit?
double primary_prod = 100.0*1000.0*10; //micro gC /m2 /day

//list of global parameters

double a0 = 0.1;
double d0 = 0.02;
double a1 = 0.1;
double d1 = 0.1 + 4.0e-3;
double alfa;
double beeta = 0.0; //density_dependent sinking rate
double dellta = 0.0; //density dependent mortality rate coefficient for attached bacteria
double dellta2 = 0.0; //density dependent mortality rate exponent for free-living bacteria

double bge = 0.5;			//bacterial growth efficiency
double cPperCB = 1.0e+6;			//carbon content per particle
double h1 = 5.0e-5;			//hydrolysis rate
double mR = 0.0;		//constant loss rate of resource (mA in the manuscript)
double mB = 0.015; 		//mortality of bacteria (attached)
//double mF = 0.01;

double min_app_N;    //minimum number of summation for approximation
double s0;		//supply rate of resource

double Vmax = 0.2;			//maximum uptake rate
double Kh = 1.0;		//half saturation constant

double recycling_C = 0.0; //recycling ratio of dead C of predator to POC

double diff0 = 8.64e+0; 		//diffusion rate 1.0 cm2 s-1 at deep layer,
double sinking_R = 10.0; 		//sinking rate of resource and predator
//double sinking_P = 10.0;			//sinking rate of predator (10m/d)

double delta_x = 10.0;			//spatial discretaization (10 m)
                    
                    
                    
                    