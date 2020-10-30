function X = SDEeulerImplicitExplicit(ffun,gfun,T,x0,dW,varargin)
tol = 1.0e-8;
maxit = 10;

N = size(T,2)-1;
nx = size(x0,1);
X = zeros(nx,N+1);

X(:,1) = x0;
for k=1:N
dt = T(k+1)-T(k);
f = feval(ffun,T(k),X(:,k),vararginf:g);
g = feval(gfun,T(k),X(:,k),vararginf:g);
psi = X(:,k) + g*dW(:,k);
xinit = psi + f*dt;
[X(:,k+1),f,~] = SDENewtonSolver(...
ffun,...
T(:,k+1),dt,psi,xinit,...
tol,maxit,...
vararginf:g);
end