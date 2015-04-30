%  fbm=fbm(H,n)
% This function is to synthesize the fractional Brownian motion by the Wood-Chan
% algorithm
% Input
% H is the given Hurst number
% n is the length of the data
% Output
% fbm is the time sequences of synthesized fractional Brownian motion
% 
%See also: wfbm
% 



% ---------------------------------------generateur------------
function fbm=fbm(H,n)
vp         =  vpropC(n, H);
m          =  length(vp);
ar        =  randn(1,m/2 + 1);      
ai        =  randn(1,m/2 + 1);            
ar(1)         =  sqrt(2) * ar(1);
ar(m/2 + 1)     =  sqrt(2) * ar(m/2 + 1);
ai(1)         =  0;
ai(m/2 + 1)     =  0;
ar        = [ar(1:m/2 + 1)  ar(m/2:-1:2)];
aic         = - ai;
ai        = [ai(1:m/2 + 1)  aic(m/2:-1:2)];
 
 
a        = ar+i*ai;
a         = sqrt(vp) .* a;
a         = real(fft(a));
a         = (1/(sqrt(2 * m))) * a;
fbm        = cumsum(a(1:n));
% ------------------------vpropC---------------------------
function    v=vpropC(n, H)
%  -------------------------------------------------------------------------
%  Input                   :    n : length of time series
%                              H : self-similarity parameter
%                         H in (0,1)     
%
%  Output                  :     eigenvalues of the Circulant Matrix of 
% covFGN
% --------------------------------------------------------------------
m=2^nextpow2(2*(n-1));        % research of the power of two greater than 2*(n-1)
u=fgncv(n, H, 1+m/2);             % u=c(0),c(1/n),...c(m/2/n)
ligneC=[u fliplr(u(2:m/2))];     % [c(0),c(1/n),...c(m/2/n)],[c(m/2-1)/n),...,c(1/n)]]
v=real(fft(ligneC));
% ------------------------cFGN----------------------
function cFGN=fgncv(n,H,M)  
%  Output         cFGN : vector of covariances of Fractional Gaussian
%     cFGN(i)=Cov(X(0),X(i/n)) = (1/2) (1/n)^2H ( |i-1|^2H -2|i|^2H 
% +|i+1|^2H )
 
k =  0:(M - 1);
H2 = 2 * H;
cFGN = (abs((k - 1)/n).^H2 - 2 * abs(k/n).^H2 + abs((k + 1)/n).^H2)/2;
