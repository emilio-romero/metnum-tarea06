#include "eigenv.h" 
int main(int argc, char *argv[]){
char archivo[30]; 
int iter; double tol; 
int nr, nc,nrb,ncb; 
double **A,**v; 
double **B;
readParams(argc,argv,archivo,&iter,&tol);
A=readMatrix(archivo,&nr,&nc);
B=readMatrix("mat005.bin",&nrb,&ncb);
v=createMatrix(nr,nc);
printf("==========================\n");
printf("Ejecucion de la potencia inversa\n");
printf("==========================\n");
paresEigen(B,nrb,iter,tol,50);
printf("==========================\n");
printf("Ejecucion de Jacobi\n");
printf("==========================\n");
Jacobi(v,A,nr,iter,sqrt(DBL_EPSILON));


freeMatrix(A);
freeMatrix(v);
printf("\nSu programa ha terminado excelso caballero/a\n");
return 0;}
