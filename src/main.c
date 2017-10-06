#include "eigenv.h" 
int main(int argc, char *argv[]){
char archivo[30]; 
int iter; double tol; 
int nr, nc; 
double **A; 

readParams(argc,argv,archivo,&iter,&tol);
A=readMatrix(archivo,&nr,&nc);
printf("%d %d %lf \n",nr,iter,tol);
Jacobi(A,nr,iter,tol);
printf("La norma de A es: %lf \n",normaInf(A,nr));
InversePower(A,0.0,nr,iter,tol);
printf("\n anterior calculo de uno\n");
//paresEigen(A,nr,iter,tol);
printf("\nSu programa ha terminado excelso caballero/a\n");
return 0;}
