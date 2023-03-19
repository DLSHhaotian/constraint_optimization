function [c,ceq,dcdx,dceqdx]=confunHimmeblau(x,p)
c=zeros(0,1);
ceq=zeros(1,1);
%c<=0  
x1=x(1,1);
x2=x(2,1);
tmp=x1+2;
%Equality constraint 
ceq=tmp^2-x2;

% computer constraint gradients
if nargout>2
    dcdx=zeros(2,0);
    dceqdx=zeros(2,1);
    %Gradient of Equality constraint
    dceqdx(1,1)=2*tmp;
    dceqdx(2,1)=-1;
end