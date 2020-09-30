function [x,output]=Pri_Simple(c,A,b,base)
% Pri_Simple  Primal simplex algorithm
%          min  c'*x
%           x
%          s.t. A x  = b      
%               x >= 0   
%         base: base vector
% Syntax: [x,output]=Pri_Simple(c,A,b,base)
%         output.fval: minimum value
%         output.case: 
%                      'found' : optimal solution found
%                      'unbound' : the problem is unbounded
%         output.xarray: Iteration trajectory  
iteration_max=400;
iteration=0;
Xarray=[];
nx=size(A,2);
nc=size(A,1);
%if c>=0 the x can be calculated directly, but only if the constraint is satisfied
if c>=0
    index_c=find(c~=0,1,'last');
    test_c=inv(A(:,(nx-nc+1):nx))*b;
    if test_c>=0
        x=zeros(1,index_c);
        output.fval=0;
    else
        output.case='no optimal point';
        x=NaN;
        output.fval=NaN;
        return;
    end
end
%Initialization of nobasevector
nobase=zeros(1,1);
comp1=1:nx;
count=1;
for i=1:nx
    if(isempty(find(base==comp1(i),1)))
        nobase(count)=i;
        count=count+1;
    end
end
B=A(:,base);
x_B=inv(B)*b;
while(iteration<=iteration_max)
    B=A(:,base);
    N=A(:,nobase);
    c_B = c(base);
    c_N = c(nobase);
    x_B=inv(B)*b;
    
    lambda=inv(B)'*c_B;
    for i=1:length(nobase)
    s_N(i)=c_N(i)-N(:,i)'*lambda;
    end
    [min_q,index_q]=min(s_N);
    %optimal solution found
    if min_q>=0
        output.case='found';
        output.fval=c_B'*x_B;
        index_c=find(c~=0,1,'last');
        for i=1:index_c
            value_c=find(base==i,1);
            if isempty(value_c)
                x(i)=0;
            else
                x(i)=x_B(value_c);
            end
        end
        break;
    end
    %current value
    index_c=find(c~=0,1,'last');
        for i=1:index_c
            value_c=find(base==i,1);
            if isempty(value_c)
                x_tmp(i)=0;
            else
                x_tmp(i)=x_B(value_c);
            end
        end
    x=x_tmp;
    %fprintf("step %d,the current feasible solution is:\n",iteration)
    %fprintf('%.4f, ',x_tmp);
    %fprintf(',\nand the value of z is %.4f\n\n',c_B'*x_B);
    %solve d
    d=inv(B)*A(:,nobase(index_q));
    
    if d<=0
        x=NaN;
        output.fval=NaN;
        output.case='unbounded';
        return;
    end
    %find the leaving indices
    min_xq=inf;
    index_xq=0;
    for i=1:length(d)
        if d(i)>0
            xq=x_B(i)/d(i);
            if xq<min_xq
                min_xq=xq;
                index_xq=i;
            end
        end
    end
    %update the basevector and nobasevector
    tmp=base(index_xq);
    base(index_xq)=nobase(index_q);
    nobase(index_q)=tmp;
    %update x_B
    x_B=x_B-d*min_xq;
    iteration=iteration+1;
    Xarray=[Xarray x_tmp];
end
output.iteration=iteration;
output.xarray=Xarray;
end