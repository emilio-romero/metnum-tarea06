#ifndef EIGENV_H
#define EIGENV_H
#include <stdio.h>
#include <stdlib.h> 
#include <math.h> 
#include "lectura.h" 
#include "solucionadores.h"
double *CopiaVec(double *v, int n);
double *InitVec(int n);
double **InitMat(int n);
double *matxvec(double **A, double *x, int n);
double *NormalizeVec(double *v, int n);
double ProdPunto(double *a, double *b, int n);
double kEigenValue(double **A, double *v, int n);
double Norma2V(double *v, int n);
double ErrorVec(double **A, double *v, double lambda, int n);
double EigenValue(double **A, int n, int iter, double tol);
double InversePower(double **A, double dlta, int n, int iter, double tol);
double normaInf(double **A, int n);
void paresEigen(double **A, int n, int iter, double tol);
void encontrarMax(double **A, int n, int *mi, int *mj);
double sgn(double x);
double **Givens(int n, int mi, int mj, double c, double s);
double **mulAG(double **A, int mi, int mj, int n,double c, double s);
double **mulGA(double **A, int mi, int mj, int n, double c, double s);
double *Jacobi(double **A, int n, int iter, double tol);
#endif


