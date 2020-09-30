function [x,lambda]=EqualityQPSolverLDLdense(H,g,A,b)
%EqualityQPSolverLDLdense  Equality constrained convex QP

%Syntax: [x,lambda]=EqualityQPSolverLDLdense(H,g,A,b)

%dimension of A------>number of x and equality constraints
[nx,nc]=size(A);
%setup KKT system(dense)
KKT_A=[H -A;-A' zeros(nc,nc)];
KKT_b=-[g;b];
% factorize and solve KKT system using the LDL factorization
z=zeros(nx+nc,1);
[L,D,p]=ldl(KKT_A,'vector');% factorization
z(p)=L'\(D\(L\KKT_b(p)));% back substitution
%Extract solution
x=z(1:nx);
lambda=z(nx+1:end);
end