/*
* Yongxiang Huang, last modification: 10/01/2009
* yongxianghuang@gmail.com
*
*/

/* This function is to perform the second order least square fitting*/

#include <stdlib.h>
#include <stdio.h>
#include "mex.h"
#include "math.h"



/************************************************************************/
/*                                                                      */
/* MAIN FUNCTION                                                        */
/*                                                                      */
/************************************************************************/

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[]) {
  
    /* declarations */
    int i,j,Ntmp,Ny,Nx,kp,kn;
    
    double *x,*y,*mp,p0,p1,p2,tmp,A,B,C,D,E,F,G; /*time index, time series*/
    
  
    
/*     check input*/
    if (nrhs!=2)     mexErrMsgTxt("Two inputs!");

    /* get input data */
    if (mxIsEmpty(prhs[0]))mexErrMsgTxt("The data is empty!");
    if (mxIsEmpty(prhs[1]))mexErrMsgTxt("The data is empty!");
    
    x=mxGetPr(prhs[0]);
    Nx=mxGetN(prhs[0]);
    Ntmp=mxGetM(prhs[0]);
    if(Nx<Ntmp) Nx=Ntmp;
    
    
    y=mxGetPr(prhs[1]);
    Ny=mxGetN(prhs[1]);
    Ntmp=mxGetM(prhs[1]);
    if(Ny<Ntmp) Ny=Ntmp;
     
    if (Nx!=Ny) mexErrMsgTxt("The length of time index and the input time series must be the same!");
   

/*Define the output*/
    plhs[0]=mxCreateDoubleMatrix(1,3,mxREAL);
    mp=mxGetPr(plhs[0]);
    A=0.0;
    B=0.0;
    C=0.0;
    D=0.0;
    E=0.0;
    F=0.0;
    G=0.0;
    
    for(i=0;i<Ny;i++)
    {
        A=A+x[i];
        B=B+y[i];
        C=C+x[i]*x[i];
        D=D+x[i]*y[i];
        E=E+x[i]*x[i]*x[i];
        F=F+x[i]*x[i]*x[i]*x[i];
        G=G+x[i]*x[i]*y[i];
    }

     /*   for(i=0;i<7;i++)mexPrintf("temp[%d]=%f\n",i,temp[i]);*/
    tmp=F*A*A-2*A*C*E+C*C*C-Ny*F*C+Ny*E*E;
    p0=(B*E*E+C*C*G-B*C*F+A*D*F-C*D*E-A*E*G)/tmp;
    p1=(C*C*D+A*B*F-B*C*E-A*C*G-Ny*D*F+Ny*E*G)/tmp;/*second order polyfit*/
    p2=(B*C*C+A*A*G-A*B*E-A*C*D+Ny*D*E-Ny*C*G)/tmp;
    mp[0]=0.0;
    kp=0;
    kn=0;
    for(i=0;i<Ny;i++)
    {
        tmp=(y[i]-p0-p1*x[i]-p2*x[i]*x[i]);
        if(tmp>0){mp[0]+=tmp*tmp;kp+=1;}
        if(tmp<0){mp[1]+=tmp*tmp;kn+=1;}
        mp[2]+=tmp*tmp;
    }/*second order DFA*/
    mp[0]=mp[0]/kp;
    mp[1]=mp[1]/kn;
    mp[2]=mp[2]/Ny;
    mp[0]=pow(mp[0],0.5);
    mp[1]=pow(mp[1],0.5);
    mp[2]=pow(mp[2],0.5);
}
