function [coef]=mfDFA(x,Nl,tlag,q,dq)
% [coef]=mfDFA(x,Nl,tlag,q,dq)
% This function is to perform the multifractal detrended fluctuation
% analysis
% Input
% x is the data to be analyzed
% Nl is the segment length
% tlag is the maximum time lag
% q is the statistical order
% dq is the increment order
% Output
% coef is the structure function of the detrended spectra
%      coef.p1 first order DFA
%              coef.p1(1,:) the positive part
%              coef.p1(2,:) the negative part
%              coef.p1(3,:) the original
%      coef.p2 second order DFA
%              coef.p2(1,:) the positive part
%              coef.p2(2,:) the negative part
%              coef.p2(3,:) the original
%      coef.q statistical order 
%      coef.scale the corresponding scale 
% Time series may possess its owner direction, therefore, we do not do the
% analysis by reverse the data x
% 
% Written by Yongxiang HUANG 23/10/2009
% last modified by Yongxiang HUANG 30/03/2010
% 
% See also,  waveleader, sfscaling, acfscaling, pdfscaling
% 

%   References:
%   HUANG Y., SCHMITT F.G., LU Z. LIU Y. Arbitrary order Hilbert spectral analysis 
%  for time series possessing scaling statistics: a comparison study
%  Physical Review E (submitted)

if nargin==0
    error('At least one input!');
end

Nx=length(x);

if nargin==1
    Nl=Nx;
    tlag=fix(Nx/10);
    q=[0 5];
    dq=1;
end

if nargin==2
    tlag=fix(Nx/10);
    q=[0 5];
    dq=1;
end
if nargin==3
    q=[0 5];
    dq=1;
end
if nargin==4
    dq=1;
end

if Nl>Nx
    Nl=fix(Nx/4);
end

Ns=fix(Nx/Nl);

if Ns>1
    tic;
    for i=1:Ns
        tmp=x((i-1)*Nl+1:i*Nl);
        [coef1,scale]=DFAM(tmp,tlag,q,dq);
        if i==1
            coef.p1=coef1.p1;
            coef.p2=coef1.p2;
        else
            coef.p1=coef.p1+coef1.p1;
            coef.p2=coef.p2+coef1.p2;
        end
        
        if mod(i,10)==0
            disp(['Finish the ',num2str(i),'th segment DFA analysis!']);
            disp(['Still have ',num2str(Ns-i),'th segments or ',num2str((Ns-i)/Ns*100),'%!']);
            toc
        end
    end
    coef.p1=coef.p1/Ns;
    coef.p2=coef.p2/Ns;
else
    [coef,scale]=DFAM(x,tlag,q,dq);
end
coef.Tau=scale;
coef.Q=[dq:dq:q];