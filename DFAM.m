function [df,scale]=DFAM(x,tlag,q,dq)
% [df,scale]=DFAM(data,Lwin,q,dq)
% This function is to perform the multifractal maxima detrended fluctuation
% analysis (MF-DFA), which shares the same idea with WTMM (wavelet transform modulus maxima)
% Input
% x is the data to be analyzed
% Lwin is the corresponding scale
%      if length of Lwin is one, it is start scale of the detrended windows
% q    is the highest statistical moments, default value 5
% dq   is the increment of statistical moment, default value dq=0.5
% Output
% df is the detrended spectrum
% scale is the corresponding scales (windows)

% check input
if nargin~=4
    error('Four inputs are required!');
end


    Lwin=[4 6 8 fix(10.^(1:0.1:log10(tlag)))];
    


fq=dq:dq:q;
Nq=length(fq);
scale=Lwin;
Ns=length(scale);

df1=zeros(3,Nq,Ns); % the first order DFA
df2=zeros(3,Nq,Ns); % the second order DFA


% pre-treat the data
x=x-mean(x);
y=cumsum(x); % get the cumulant of the original siginal
% y=x;
for i=1:Ns
    coef=MYDFA2(y,scale(i));
    for j=1:Nq
        if fq(j)~=0
            df1(1,j,i)=mean(coef(1,:).^fq(j)).^(1/fq(j));% first order DFA
            df1(2,j,i)=mean(coef(2,:).^fq(j)).^(1/fq(j));
            df1(3,j,i)=mean(coef(3,:).^fq(j)).^(1/fq(j));
            
            df2(1,j,i)=mean(coef(4,:).^fq(j)).^(1/fq(j));% second order DFA
            df2(2,j,i)=mean(coef(5,:).^fq(j)).^(1/fq(j));
            df2(3,j,i)=mean(coef(6,:).^fq(j)).^(1/fq(j));
        else
            df1(1,j,i)=mean(coef(1,:).^fq(j)); % first order DFA
            df1(2,j,i)=mean(coef(2,:).^fq(j));
            df1(3,j,i)=mean(coef(3,:).^fq(j));
           
            df2(1,j,i)=mean(coef(4,:).^fq(j));% second order DFA
            df2(2,j,i)=mean(coef(5,:).^fq(j));
            df2(3,j,i)=mean(coef(6,:).^fq(j));
        end
    end
    
end
df.p1=df1;
df.p2=df2;
% df.p3=df3;
% df.p4=df4;

        