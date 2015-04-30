

function coef=MYDFA2(y,Lwin)
% to estimate the fluctuation function
if size(y,1)>size(y,2);
    y=y';
end

Ny=length(y);
Nsg=fix(Ny/Lwin); % get the number of segments
coef=zeros(6,Nsg);
tt=[1:Lwin];
for i=1:Nsg
    tmp=y((i-1)*Lwin+1:i*Lwin);
    
    coefT=DFA1(tt,tmp);% first order DFA
    coef(1,i)=coefT(1);
    coef(2,i)=coefT(2);
    coef(3,i)=coefT(3);
    
    coefT=DFA2(tt,tmp);% second order DFA
    coef(4,i)=coefT(1);
    coef(5,i)=coefT(1);
    coef(6,i)=coefT(1);
end

