//this file is modified at 2011/05/03
void message_error(char error_text[]);
double *d_vector(long size);
int *i_vector(long size);
double **d_matrix(long size_row, long size_column);
int **i_matrix(long size_row, long size_column);
void free_d_vector(double *x);
void free_i_vector(int *x);
void free_d_matrix(double **x);
void free_i_matrix(int **x);

double getETime(void);
double getCPUTime(void);

float ran2(long *idum);
float gasdev(long *idum);
double gasdev2();
void ar1(double *a, double p, long number,  double series_min, unsigned long seed);
void ar2(double *a, double p, long number,  double series_min, unsigned long seed);
void pink(double *a, double p, long number,  double series_min, unsigned long seed);




