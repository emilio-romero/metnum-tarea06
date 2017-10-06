#include <stdio.h>
#include <stdlib.h> 
#include "lectura.h"
#include "solucionadores.h"
#include <math.h> 

double *diagonal(double *d,double *b,int n){
double *vec=(double*)malloc(n*sizeof(double)); 
double error=0.0;
int i;
for(i=0;i<n;i++){
  if(d[i]==0.0){
  printf("\nNo hay solución a este sistema\n");
  return(NULL);
  }
} 
for(i=0;i<n;i++){
vec[i]=b[i]/d[i]; 
} 
printf("\nEl tamanio de la matriz es: %d",n);
for(i=0;i<n;++i)
  error+=(d[i]*vec[i]-b[i])*(d[i]*vec[i]-b[i]);
printf("\nEl error es: %g\n",sqrt(error));

return(vec);}
//===========================================
double *tinferior(double **L, double *b,int n){
  double *vec=(double*)malloc(n*sizeof(double)); 
  int i,j,k;
  for(i=0;i<n;i++) vec[i]=0; 
  double tol=10E-10;
  double suma; 
for(i=0;i<n;i++){
 if(fabs(L[i][i])<tol){ 
  printf("\nEl sistema no tiene soluciones\n");
  return(NULL);}}
 
vec[0]=b[0]/L[0][0]; 
for(k=1;k<n;++k){
  suma=0;
  for(j=0;j<k;j++){
    suma+=L[k][j]*vec[j];
  }
  vec[k]=(b[k]-suma)/L[k][k]; 
}
//Agregar el tamaño y el error
/*double error=0; 
double *vaux=(double*)malloc(n*sizeof(double)); 
for(i=0;i<n;i++){
 for(j=0;j<n;j++){
 vaux[i]+=L[i][j]*vec[j];}
error+=(vaux[i]-b[i])*(vaux[i]-b[i]);}  
*///printf("\nEl error de L es: %g \n",sqrt(error));
return(vec);}
//===========================================
double *tsuperior(double **U, double *b,int n){
  double *vec=(double*)malloc(n*sizeof(double)); 
  int i,j,k;
  for(i=0;i<n;i++) vec[i]=0; 
  double tol=10E-10, suma; 
for(i=0;i<n;i++){
 if(fabs(U[i][i])<tol){ 
  printf("\nEl sistema no tiene soluciones\n");
  return(NULL);}} 

vec[n-1]=b[n-1]/U[n-1][n-1]; 
for(i=n-2;i>=0;i--){
  suma=0;
  for(j=i;j<n;j++){
    suma+=U[i][j]*vec[j];
  }
  vec[i]=(b[i]-suma)/U[i][i]; 
}
//Agregar el tamaño y el error 
//printf("\nEl tamaño de la matriz es de: %dx%d\n",n,n);
/*double error=0; 
double *vaux=(double*)malloc(n*sizeof(double)); 
for(i=0;i<n;i++){
 for(j=0;j<n;j++){
 vaux[i]+=U[i][j]*vec[j];   
}
error+=(vaux[i]-b[i])*(vaux[i]-b[i]);
}*/
//printf("\nEl error de U es: %g \n",sqrt(error));
return(vec);}
//============================================
double *factLU(double **A,double *b, int n){
int j,k; 
double *vec=(double*)malloc(n*sizeof(double));
double **L, **U;
double tol=1E-10, sumal, sumau; 
L=createMatrix(n,n);
U=createMatrix(n,n);
for(int i=0;i<n;i++){
 for(j=0;j<n;j++){L[i][j]=0; U[i][j]=0;}
  U[i][i]=1.0;
  L[i][0]=A[i][0];
  U[0][i]=A[0][i]/L[0][0];
}
//printf("\nentras a LU\n");
for(int i=1;i<n;i++){
  for(j=1;j<n;j++){
    sumal=0; sumau=0;
    if(i>=j){
    for(k=0;k<j;k++) sumal+=L[i][k]*U[k][j]; 
    L[i][j]=A[i][j]-sumal; 
   } 
    else if(j>i){
     for(k=0;k<i;k++) sumau+=L[i][k]*U[k][j];
     U[i][j]=(A[i][j]-sumau)/L[i][i];     
   }
  }
}
double *vaux=(double*)malloc(n*sizeof(double));
vaux=tinferior(L,b,n);
vec=tsuperior(U,vaux,n);
//====== calculo del error ||Ax-b||
double error=0; 
double *vaux2=(double*)malloc(n*sizeof(double));
for(int i=0;i<n;i++){
 vaux2[i]=0; 
 for(j=0;j<n;j++){
  vaux2[i]+=A[i][j]*vec[j];  
 }
error+=(vaux2[i]-b[i])*(vaux2[i]-b[i]); 
} 
//printf("\nUna medida del error de LU es: %g \n",sqrt(error));
//=====
free(vaux); free(vaux2); 
freeMatrix(L); freeMatrix(U);
return(vec);}











