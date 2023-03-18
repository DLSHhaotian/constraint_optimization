function [x,lam,exitflag,output]=Pri_AsQp(H,g,Ae,be,Ai,bi,x0)
% Pri_AsQp   Primal Active-Set Algorithm
%
%          min  0.5*H'*x*H+g'*x
%           x
%          s.t. Ae x  = be      
%               Ai x >= bi      
%
% Syntax: [x,lam,exitflag,output]=Pri_AsQp(H,g,Ae,be,Ai,bi,x0)
%         output.lam: Lagrange multiplier
%         output.x_plot: Iteration trajectory              
%Initialization
%Ax>=b;
epsilon=1.0e-9;
err=1.0e-6;
iteration=0;
x=x0;
n=length(x);
iteration_max=20;
ne=length(be);
ni=length(bi);
lam=zeros(ne+ni,1);
index=ones(ni,1);
output.x_plot=[x0'];
output.lam=[lam'];

%initialize the active set(index=1)
for(i=1:ni)
    if(Ai(i,:)*x>bi(i)+epsilon)
        index(i)=0;
    end
end
%main program
while(iteration<=iteration_max)
    
    Aee=[];
% if the start point x is on the Equality Constraint or on the edge of the
% Equality Constraint, put it in active set.
    if(ne>0)
        Aee=Ae;
    end
    for j=1:ni
        if(index(j)>0)
            Aee=[Aee;Ai(j,:)];
        end
    end
    %Solve subproblem to find p and Compute Lagrange multipliers
    gk=H*x+g;
    [m1,n1]=size(Aee);
    [p,lam]=As_sub(H,gk,Aee,zeros(m1,1));
    if(norm(p)<=err)
        lambda_min=0.0;
        if(length(lam)>ne)
            [lambda_min,jk]=min(lam(ne+1:end));
        end
        if(lambda_min>=0)
            exitflag=1;
        else
            exitflag=0;
            %remove the (Lagrange multipliers min) inequlity Constraint
            %form active set
            for(i=1:ni)
                if(index(i)&(sum(index(1:i)))==jk)
                    index(i)=0;
                    break;
                end
            end
        end
        iteration=iteration+1;
    else
        exitflag=0;
        %computer the step length
        alpha=1.0;
        tm=1.0;
        for(i=1:ni)
            if((index(i)==0)&(Ai(i,:)*p<0))
                tm1=(bi(i)-Ai(i,:)*x)/(Ai(i,:)*p);
                if(tm1<tm)
                    tm=tm1;
                    ti=i;
                end
            end
        end
        
        alpha=min(alpha,tm);
        x=x+alpha*p;
        %update the active set
        if(tm<1)
            index(ti)=1;
        end
    end
    if(exitflag==1)
        break;
    end
    iteration=iteration+1;
    output.x_plot=[output.x_plot;x'];
    if(length(lam)<(ne+ni))
        lam=[lam;zeros(ne+ni-length(lam),1)];
    end
    output.lam=[output.lam;lam'];
end
output.fval=0.5*x'*H*x+g'*x;
output.iter=iteration;
