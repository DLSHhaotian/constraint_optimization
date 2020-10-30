function [x,f,J] = SDENewtonSolver(ffun,t,dt,psi,xinit,tol,maxit,varargin)
I= eye(length(xinit));
x = xinit;
[f,J] = feval(ffun,t,x,vararginf:g);
R = x - f*dt - psi;
it = 1;
while ( (norm(R,'inf') > tol) & (it <= maxit) )
dRdx = I - J*dt;
mdx = dRdxnR;
x = x - mdx;
[f,J] = feval(ffun,t,x,vararginf:g);
R = x - f*dt - psi;
it = it+1;
end