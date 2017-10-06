#ifndef SOLUCIONADORES_H
#define SOLUCIONADORES_H
double *diagonal(double *d, double *b, int n);
double *tinferior(double **L, double *b, int n);
double *tsuperior(double **U, double *b, int n);
double *factLU(double **A, double *b, int n);

#endif 
