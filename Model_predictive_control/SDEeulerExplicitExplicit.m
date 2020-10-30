function X = SDEeulerExplicitExplicit(ffun,gfun,T,x0,W,varargin)
N = size(T,2)-1;
nx = size(x0,1);
X = zeros(nx,N+1);

X(:,1) = x0;
for k=1:N
dt = T(k+1)-T(k);
dW = W(:,k+1)-W(:,k);
f = feval(ffun,T(k),X(:,k),varargin{:});
g = feval(gfun,T(k),X(:,k),varargin{:});
X(:,k+1) = X(:,k) + f*dt + g*dW;
end