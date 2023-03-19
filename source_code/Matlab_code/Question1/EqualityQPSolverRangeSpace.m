function [x,lambda]=EqualityQPSolverRangeSpace(H,g,A,b)
%EqualityQPSolverRangeSpace  Equality constrained convex QP

%Syntax: [x,lambda]=EqualityQPSolverRangeSpace(H,g,A,b)
%H is positive definite

%solve Hv=g
v=H\g;
%form HA
H_A=A'*inv(H)*A;
%solve lambda
lambda=H_A\(b+A'*v);
%solve x
x=H\(A*lambda-g);
end
