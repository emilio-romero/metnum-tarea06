#include "eigenv.h" 
int main(int argc, char *argv[]){
char archivo[30]; 
int iter; double tol; 
int nr, nc; 
double **A,**v; 

readParams(argc,argv,archivo,&iter,&tol);
A=readMatrix(archivo,&nr,&nc);
printf("%d %d %lf \n",nr,iter,tol);
v=Jacobi(A,nr,iter,sqrt(DBL_EPSILON));
//printf("La norma de A es: %lf \n",normaInf(A,nr));
//InversePower(A,0.0,nr,iter,tol);
//paresEigen(A,nr,iter,tol,50);
freeMatrix(A);
freeMatrix(v);
printf("\nSu programa ha terminado excelso caballero/a\n");
return 0;}
