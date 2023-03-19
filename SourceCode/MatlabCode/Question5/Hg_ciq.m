function [ciq,dciq,d2ciq]=Hg_ciq(x)
%Syntax: [ciq,dciq,d2ciq]=Hg_ciq(x)
%Inequality constraint ciq(x1,x2)=-x1+3>0
%                                =x2+2>0
%Implemention of inequality constraint(ciq),gradient(dciq),Hessian(d2ciq)

x1=x(1,1);
x2=x(2,1);

%Inequality constraint 
ciq=[3-x1;2+x2];
%Gradient of inequality constraint
dciq=zeros(2,2);
dciq(1,1)=-1;
dciq(2,2)=1;
%Hessian of inequality constraint
d2ciq=zeros(2,2,2);
end
