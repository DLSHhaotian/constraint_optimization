function [H,g,A,b]=u2HgAb(n,u_mean,d0)
%u2HgAb  Generate the H,g,A,b of test problem

%Syntax: [H,g,A,b]=u2HgAb(n,u_mean,d0)
%n: The number of variables
%u_mean,d0: The changed constant in test problem
H=eye(n+1);
g=-2*u_mean*ones(n+1,1);
b=zeros(n, 1);
b(1,1)=-d0;
A=[zeros(1,n-1);eye(n-1)];
A=[A zeros(n,1)]-eye(n);
A=A';
A=[A; zeros(1, n)];
A(n,1)=1;
A(n+1,n)=-1;
end