//this file is modified at 2011/05/03
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <complex.h>
#include <sys/time.h>

#ifdef LIN
#include <sys/resource.h>
#endif

//#include "nrutil_new.h"
#define NR_END 1
#define FREE_ARG char*
#define PAI_A 3.14159265358979 

void message_error(char error_text[])
/*Standard error handler*/
{
  printf("There are some errors...\n");
  printf("%s\n", error_text);
  printf("...now existing to system...\n");
  exit(1);
}


double *d_vector(long size)	//for generating double vector[1]...vector[size]
{
  double *x;		//pointer
  x = (double *) malloc((size_t) ((size + 1)*sizeof(double)));	//allocation of memory space to stock (size + 1) double variables & malloc returns the address of the top of the allocated memory space to x
  //the reason to prepare (size + 1) is just to use from x[1] to x[size]
  if(x == NULL) message_error("allocation failure in d_vector()");//if memory allocation was failed, malloc returns NULL and x becomes NULL
  return x;		//return the address of the top of the allocated memory space to x
}

int *i_vector(long size)	//for generating double vector[1]...vector[size]
{
  int *x;		//pointer
  x = (int *) malloc((size_t) ((size + 1)*sizeof(int)));	//allocation of memory space to stock (size + 1) double variables & malloc returns the address of the top of the allocated memory space to x
  //the reason to prepare (size + 1) is just to use from x[1] to x[size]
  if(x == NULL) message_error("allocation failure in d_vector()");//if memory allocation was failed, malloc returns NULL and x becomes NULL
  return x;		//return the address of the top of the allocated memory space to x
}

double **d_matrix(long size_row, long size_column) 	//for generating double matrix[1][1]...matrix[size_row][size_column]
{
	double **x;			//pointer to 'pointer valuable', valuable to stock the address of the pointer valuable
	long i;
	long size_row_P = size_row + 1;  //technical (not necessary) statement just to start from [1][1]
	long size_column_P = size_column + 1;  //technical (not necessary) statement just to start from [1][1]
	
	x = (double **) malloc((size_t)  (size_row_P*sizeof(double *)));			//allocation of memory spate to stock (size_row) pointer valuables to double, & malloc returns the address of the top of the allocated memory space to x
	if(x == NULL) message_error("allocation failure in d_vector()");//if memory allocation was failed, malloc returns NULL and x becomes NULL
	
	x[0] = (double *) malloc((size_t) (size_row_P*size_column_P*sizeof(double)));
	//allocation of memory scape to stock (size_row*zsize_column) doubles variables & malloc returns the address of the top of the allocated memory space to x[0]
	//Note that x[0] (== *x) is the value (with type pointer to double) of the pointer valuable , pointed by x
	
	if(x[0] == NULL) message_error("allocation failure in d_vector()");	//if memory allocation was failed, malloc returns NULL and x becomes NULL		
	
	for(i = 1; i < size_row_P; i++) x[i] = x[0] + i*size_column_P;	//operating on pointer
	
	return x;
	
}


int **i_matrix(long size_row, long size_column) 	//for generating double matrix[1][1]...matrix[size_row][size_column]
{
	int **x;			//pointer to 'pointer valuable', valuable to stock the address of the pointer valuable
	long i;
	long size_row_P = size_row + 1;  //technical (not necessary) statement just to start from [1][1]
	long size_column_P = size_column + 1;  //technical (not necessary) statement just to start from [1][1]
	
	x = (int **) malloc((size_t)  (size_row_P*sizeof(int *)));			//allocation of memory spate to stock (size_row) pointer valuables to double, & malloc returns the address of the top of the allocated memory space to x
	if(x == NULL) message_error("allocation failure in d_vector()");//if memory allocation was failed, malloc returns NULL and x becomes NULL
	
	x[0] = (int *) malloc((size_t) (size_row_P*size_column_P*sizeof(int)));
	//allocation of memory scape to stock (size_row*zsize_column) doubles variables & malloc returns the address of the top of the allocated memory space to x[0]
	//Note that x[0] (== *x) is the value (with type pointer to double) of the pointer valuable , pointed by x
	
	if(x[0] == NULL) message_error("allocation failure in d_vector()");	//if memory allocation was failed, malloc returns NULL and x becomes NULL		
	
	for(i = 1; i < size_row_P; i++) x[i] = x[0] + i*size_column_P;	//operating on pointer
	
	return x;
	
}

void free_d_matrix(double **x)
{
	free(x[0]);
	free(x);
}

void free_i_matrix(int **x)
{

	free(x[0]);
	free(x);
}

void free_d_vector(double *x)
{
	free(x);
}

void free_i_vector(int *x)
{
	free(x);
}

#ifdef LIN
double getETime(void)
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + (double) tv.tv_usec*1.0e-6;
}

double getCPUTime(void)
{
	struct rusage RU;
	getrusage(RUSAGE_SELF, &RU);
	return RU.ru_utime.tv_sec + (double) RU.ru_utime.tv_usec*1.0e-6;
}
#endif



