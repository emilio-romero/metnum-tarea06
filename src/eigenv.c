#include "eigenv.h" 
void CopiaVec(double *copia,double *v, int n){
  for(int i=0;i<n;i++) copia[i]=v[i];
}

void InitVec(double *aux,int n){
  for(int i=0;i<n;i++) aux[i]=1.0;
}

void InitMat(double **mat,int n){
  for(int i=0;i<n;i++){
   for(int j=0;j<n;j++){
    mat[i][j]=0.0;
   }
   mat[i][i]=1.0;
  }
}

double *matxvec(double **A, double *x, int n){
  double *aux=(double*)malloc(n*sizeof(double));
  for(int i=0;i<n;i++){
   double suma=0; 
   for(int j=0;j<n;j++){
     suma+=A[i][j]*x[j];
   }
   aux[i]=suma;    
  }
return aux;}
//===
double *NormalizeVec(double *v, int n){
  double *aux=(double*)malloc(n*sizeof(double));
  double norma=0; 
  for(int i=0;i<n;i++) norma+=v[i]*v[i];
  for(int j=0;j<n;j++) aux[j]=v[j]/sqrt(norma);

  return aux;
}
//===
double ProdPunto(double *a, double *b, int n){
 double c=0; 
 for(int i=0;i<n;i++)
  c+=a[i]*b[i]; 
return c;  
}
//===
double kEigenValue(double **A, double *v, int n){
  double eva; 
  double *aux=(double*)malloc(n*sizeof(double));
  aux=matxvec(A,v,n);
  eva=ProdPunto(v,aux,n);
  return eva; 
}
//===
double ErrorVec(double **A, double *v, double lambda, int n){
  double mierror=0.0; //fue dejarla ir? 
  double *vaux=(double*)malloc(n*sizeof(double));
  vaux=matxvec(A,v,n);
  for(int i=0;i<n;i++) vaux[i]=vaux[i]-lambda*v[i];
  for(int k=0;k<n;k++) mierror+=vaux[k]*vaux[k];
return sqrt(mierror);}
//===
double Norma2V(double *v,int n){
  double norma=0; 
  for(int i=0;i<n;i++) norma+=v[i]*v[i];
return sqrt(norma);}
//====
double EigenValue(double **A, int n, int iter, double tol){
  double lamb,eaux; 
  double *vaux=(double*)malloc(n*sizeof(double));
  InitVec(vaux,n);
 // FILE *out; 
 // out=fopen("error1.dat","w");
  int cont;
  for(cont=0;cont<iter;cont++){
    vaux=matxvec(A,vaux,n);
    vaux=NormalizeVec(vaux,n);
    lamb=kEigenValue(A,vaux,n);
    eaux=ErrorVec(A,vaux,lamb,n);
   // fprintf(out,"%d %lf \n",cont,eaux);
    if(eaux<tol) break;
  }
 // fclose(out);
  printf("\n La matriz es de: %dx%d",n,n);
  printf("\n El valor del mayor eigenvalor es: %lf",lamb);
  printf("\n Se realizaron <%d> iteraciones",cont);
  printf("\n El error final es: %lf",eaux);
  printf("\n================================\n");
return lamb;}

double InversePower(double **A, double dlta, int n, int iter,double tol,double *err, int *Miter){
  double **aaux=createMatrix(n,n);
 for(int i=0;i<n;i++){
  for(int j=0;j<n;j++) aaux[i][j]=A[i][j]; 
 }
  double mu,rho,eaux; 
  double *vaux=(double*)malloc(n*sizeof(double));
  double *yaux=(double*)malloc(n*sizeof(double));
  double *waux=(double*)malloc(n*sizeof(double));
  double *raux=(double*)malloc(n*sizeof(double));
    for(int i=0;i<n;i++){
      aaux[i][i]=A[i][i]-dlta;
      yaux[i]=0; waux[i]=0; raux[i]=0; 
   }//}
  InitVec(vaux,n);
  int cont; 
  for(cont=0;cont<iter;cont++){
 //  printf("\nLlegas hasta aca\n");
   yaux=factLU(aaux,vaux,n); double ynorm=Norma2V(yaux,n);
   CopiaVec(waux,vaux,n); for(int k=0;k<n;k++)  waux[k]=waux[k]/ynorm; 
   CopiaVec(vaux,yaux,n); for(int k=0;k<n;k++)  vaux[k]=vaux[k]/ynorm; 
   rho=ProdPunto(vaux,waux,n);
   mu=dlta+rho;
   for(int k=0;k<n;k++) raux[k]=waux[k]-rho*vaux[k];
   eaux=Norma2V(raux,n); 
   if(eaux<tol) break; } 
  *Miter=cont;//printf("\n Se realizaron <%d> iteraciones",cont);
  *err=eaux;//printf("\n El error final es: %lf",eaux);
   freeMatrix(aaux);
   free(vaux); free(yaux); free(waux); free(raux);
return mu;}

double normaInf(double **A, int n){
  double max=0; 
  double suma=0;
  for(int j=0;j<n;j++) max+=A[0][j]; 
  for(int i=1;i<n;i++){
    for(int j=0;j<n;j++) suma+=A[i][j]; 
    if(max<suma) max=suma;
    suma=0;
  }
return max;}

void paresEigen(double **A, int n, int iter, double tol, int N){
  double eps=sqrt(DBL_EPSILON); 
  double d=normaInf(A,n); printf("Norma inf A: %lf\n",d);
  double dd=2.0*d/(double)N;
  double mu, mu0=-d*10.0;
  double err; int Miter;
  for(int i=0;i<=N;i++){
    //printf("-d+k del d: %lf\n",-d+(double)i*dd);
    mu=InversePower(A,-d+(double)i*dd,n,1000,eps,&err,&Miter);
    if(fabs(mu0-mu)>0.0001){
      mu0=mu;
      printf("\n Se realizaron <%d> iteraciones",Miter);
      printf("\n El error final es: %g",err);
      printf("\n El eigenvalor es: %lf\n",mu0);
    }
  } 

}

double encontrarMax(double **A, int n, int *mi, int *mj){
  double max=fabs(A[1][0]);
  *mi=1; *mj=0;
  for(int i=2;i<n;i++){
   for(int j=0;j<i;j++){
     if(max<fabs(A[i][j])){
       max=fabs(A[i][j]);
       *mi=i; *mj=j; }
   } }
return max;}
double sgn(double x){
if(x>0) return 1.0;
if(x<0) return -1.0;
return 0.0;
}

//======= Rotaciones de Givens ==========
double **Givens(int n, int mi, int mj, double c, double s){
  double **aux=(double**)malloc(n*sizeof(double*));
  for(int i=0;i<n;i++) aux[i]=(double*)malloc(n*sizeof(double));
  InitMat(aux,n);
  aux[mi][mi]=c; aux[mj][mj]=c;
  aux[mi][mj]=s; aux[mj][mi]=-s;
return aux;}

void mulAG(double **aux, double **A, int mi, int mj, int n, double c, double s){
 for(int i=0;i<n;i++){
  for(int j=0;j<n;j++) aux[i][j]=A[i][j]; }

 for(int i=0;i<n;i++){
   aux[i][mi]=c*A[i][mi]-s*A[i][mj];
   aux[i][mj]=c*A[i][mj]+s*A[i][mi];
 }
}

void mulGA(double **aux,double **A, int mi, int mj, int n, double c, double s){
 for(int i=0;i<n;i++){
  for(int j=0;j<n;j++) aux[i][j]=A[i][j]; }
 for(int j=0;j<n;j++){
   aux[mi][j]=c*A[mi][j]-s*A[mj][j];
   aux[mj][j]=s*A[mi][j]+c*A[mj][j];
 }
} 
//=======================================
double AVVD(double **A, double **V, double **D,int n){
  double **AV=(double**)malloc(n*sizeof(double*));  
  double **VD=(double**)malloc(n*sizeof(double*));  
  for(int i=0;i<n;i++){
    AV[i]=(double*)malloc(n*sizeof(double));
    VD[i]=(double*)malloc(n*sizeof(double));
  }
  double sumaav; double sumavd; 
  double err;
  for(int i=0;i<n;i++){
    sumaav=0; sumavd=0;
    for(int j=0;j<n;j++){
      for(int k=0;k<n;k++){
        sumaav+=A[i][k]*V[k][j];
        sumavd+=V[i][k]*D[k][j];
      }
      AV[i][j]=sumaav;
      VD[i][j]=sumavd;
    }
  } 
  for(int i=0;i<n;i++) for(int j=0;j<n;j++) AV[i][j]=AV[i][j]-VD[i][j];
  err=normaInf(AV,n);

for(int i=0;i<n;i++){free(AV[i]); free(VD[i]);}
free(AV); free(VD);
return err;}
//=======================================

void Jacobi(double **V,double **A, int n, int iter, double tol){
 double **aux=(double**)malloc(n*sizeof(double*));
 double **Aoriginal=(double**)malloc(n*sizeof(double*));
 double **Vaux=(double**)malloc(n*sizeof(double*));
 double delta,t,c,s,Amax;
 int maxi,maxj; 
  for(int i=0;i<n;i++){
  aux[i]=(double*)malloc(n*sizeof(double));
  Aoriginal[i]=(double*)malloc(n*sizeof(double));
  Vaux[i]=(double*)malloc(n*sizeof(double));}
for(int i=0;i<n;i++) for(int j=0;j<n;j++) Aoriginal[i][j]=A[i][j];
 InitMat(V,n);
 int cont;
 for(cont=0;cont<iter;cont++){
   Amax=encontrarMax(A,n,&maxi,&maxj);// printf("Par i,j: %d, %d\n",maxi,maxj);
   if(fabs(Amax)<tol) break;
   delta=(A[maxj][maxj]-A[maxi][maxi])/(2.0*A[maxi][maxj]);
   t=sgn(delta)/(fabs(delta)+sqrt(1+delta*delta));
   c=1/sqrt(1+t*t); s=c*t;
   mulAG(aux,A,maxi,maxj,n,c,s);// for(int b=0;b<n;b++)free(aux->aux[b]);
   mulGA(A,aux,maxi,maxj,n,c,s);
   mulAG(Vaux,V,maxi,maxj,n,c,s);
   for(int k=0;k<n;k++) for(int l=0;l<n;l++) V[k][l]=Vaux[k][l];
 }
 if(n<30){
 printf("Matriz de eigenvalores\n");
 printMatrix(A,n,n);
 printf("Matriz de eigenvectores\n");
 printMatrix(V,n,n);}
 else{
   writeMatrix(A,n,n,"Greateigenvalues.bin");
   writeMatrix(V,n,n,"Greateigenvectors.bin");
 }
 printf("El error es de: %g\n",AVVD(Aoriginal,V,A,n));
 
//=====Liberacion y regreso
for(int i=0;i<n;i++) {free(aux[i]); free(Vaux[i]); free(Aoriginal[i]);}
 free(aux); free(Vaux); free(Aoriginal);
}
//=====Fin funcion jacobi

