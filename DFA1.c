/*
* Yongxiang Huang, last modification: 10/01/2009
* yongxianghuang@gmail.com
*
*/

/* This function is to perform the first order least square fitting*/

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
    
    double *x,*y,p0,p1,A,B,C,D,tmp,*mp; /*time index, time series*/
    
  
    
/*     check input*/
    if (nrhs!=2)     mexErrMsgTxt("Two parameters!");

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
    
    for(i=0;i<Ny;i++)
    {
        A=A+x[i];
        B=B+y[i];
        C=C+x[i]*x[i];
        D=D+x[i]*y[i];
       
    }

    /*First order DFA*/
    tmp=Ny*C-A*A;
    p0=(B*C-A*D)/tmp;
    p1=(Ny*D-A*B)/tmp;/*First order polyfit*/
    kp=0;
    kn=0;
    for(i=0;i<Ny;i++)
    {
        tmp=(y[i]-p0-p1*x[i]);
        if(tmp>0){mp[0]+=tmp*tmp;kp+=1;}
        if(tmp<0){mp[1]+=tmp*tmp;kn+=1;}
        mp[2]+=tmp*tmp;
    }/*First order DFA*/
    mp[0]=mp[0]/kp;
    mp[1]=mp[1]/kn;
    mp[2]=mp[2]/Ny;
    mp[0]=pow(mp[0],0.5);
    mp[1]=pow(mp[1],0.5);
    mp[2]=pow(mp[2],0.5);
}
