function [ceq,dceq,d2ceq]=Hg_ceq(x)
%Syntax: [ceq,dceq,d2ceq]=Hg_ceq(x) 
%Equality constraint ceq(x1,x2)=(x1+2)^2-x2=0
%Implemention of equality constraint(ceq),gradient(dceq),Hessian(d2ceq)

x1=x(1,1);
x2=x(2,1);
tmp=x1+2;
%Equality constraint 
ceq=tmp^2-x2;
%Gradient of Equality constraint
dceq=zeros(2,1);
dceq(1,1)=2*tmp;
dceq(2,1)=-1;

%Hessian of Equality constraint
d2ceq=zeros(2,2,1);
d2ceq(1,1,1)=2;
end
