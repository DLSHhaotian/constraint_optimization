function [x,lambda]=EqualityQPSolverNullSpace(H,g,A,b)
%EqualityQPSolverNullSpace  Equality constrained convex QP

%Syntax: [x,lambda]=EqualityQPSolverNullSpace(H,g,A,b)

%dimension of A------>number of x and equality constraints
[nx,nc]=size(A);
%QR factorization
[Q,Rbar] = qr(A);
m1 = size(Rbar,2);
Q1 = Q(:,1:m1); 
Q2 = Q(:,m1+1:nx);
R = Rbar(1:m1,1:m1);

xY=R'\b;
xZ=(Q2'*H*Q2)\(-Q2'*(H*Q1*xY+g));
x=Q1*xY+Q2*xZ;
lambda=R\(Q1'*(H*x+g));
end