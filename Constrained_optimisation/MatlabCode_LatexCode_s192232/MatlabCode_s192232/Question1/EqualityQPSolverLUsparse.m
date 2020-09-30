function [x,lambda]=EqualityQPSolverLUsparse(H,g,A,b)
%EqualityQPSolverLUsparse Equality constrained convex QP

%Syntax: [x,lambda]=EqualityQPSolverLUsparse(H,g,A,b)

%dimension of A------>number of x and equality constraints
[nx,nc]=size(A);
%setup KKT system(dense)
KKT_A=[H -A;-A' zeros(nc,nc)];
KKT_b=-[g;b];
%sparse KKT system
KKT_A=sparse(KKT_A);
KKT_b=sparse(KKT_b);
% factorize and solve KKT system using the LU factorization
[L,U,p]=lu(KKT_A,'vector');
z=U\(L\KKT_b(p));
%Extract solution
x=z(1:nx);
lambda=z(nx+1:end);
end