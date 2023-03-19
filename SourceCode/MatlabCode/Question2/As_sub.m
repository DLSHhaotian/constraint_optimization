function [p,lambda]=As_sub(H,g,Ae,be)
%EQP solver for the primal active set method
ginvH=pinv(H);
[n,m]=size(Ae);
%For singular matrices caused by equality constraints
if(n>0)
    rb=Ae*ginvH*g+be;
    lambda=pinv(Ae*ginvH*Ae')*rb;
    p=ginvH*(Ae'*lambda-g);
else
    p=-ginvH*g;
    lambda=0;
end