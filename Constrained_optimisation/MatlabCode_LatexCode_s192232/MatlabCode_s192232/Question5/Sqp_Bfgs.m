function [x,output]=Sqp_Bfgs(x0,lam_eq0,lam_ineq0)
% Sqp_Bfgs  SQP algorithm with damped BFGS and line search
%
%          min  f(x)
%           x
%          s.t. ceq(x)=0  (Lagrange multiplier lam_eq)    
%               ciq(x)>=0  (Lagrange multiplier lam_ineq)  
% Syntax: [x,output]=Sqp_Bfgs(x0,lam_eq0,lam_ineq0)
%         output.fval: minimum value
%         output.dL: convergence of delta_L
%         output.ceq: convergence of ceq
%         output.Xarray=Xarray: Iteration trajectory  
epsilon=1e-5;

Xarray=[];
normdLarray=[];
normceqarray=[];
Xarray=[Xarray x0];

nx=size(x0,1);
Bk=eye(nx);
%evaluate the gradient
[~,df,~]=Hg_f(x0);
[ceq,dceq,~]=Hg_ceq(x0);
[ciq,dciq,~]=Hg_ciq(x0);
xk=x0;
lam_eqk=lam_eq0;
lam_ineqk=lam_ineq0;
iteration=0;
iteration_max=200;
while(iteration<iteration_max)
    %Solve the subproblem(QP)
    [p,lameq,lamineq]=Qpsolver_Sqp(Bk,df,dceq,ceq,dciq,ciq);
    %Take step and update x(k+1),lameq(k+1),lamineq(k+1)
    xk=xk+p;
    lam_eqk=-lameq;
    lam_ineqk=-lamineq;
    %Calculate dL(k) with x(k) and lameq(k+1),lamineq(k+1)
    Lk=df-dceq.*lam_eqk-dciq*lam_ineqk;
    %Re-evaluate
    [~,df,~]=Hg_f(xk);
    [ceq,dceq,~]=Hg_ceq(xk);
    [ciq,dciq,~]=Hg_ciq(xk);
    %Calculate L(k+1) with x(k+1) and lam(k+1)
    Lk2=df-dceq.*lam_eqk-dciq*lam_ineqk;
    %Damped BFGS update
    pk=p;
    qk=Lk2-Lk;
    judge=(pk'*qk>=0.2*pk'*(Bk*pk));
    thetak=judge.*1+(~judge).*((0.8*pk'*(Bk*pk))/(pk'*(Bk*pk)-pk'*qk));
    rk=thetak*qk+(1-thetak)*(Bk*pk);
    Bk=Bk-((Bk*pk)*(pk'*Bk))/(pk'*(Bk*pk))+(rk*rk')/(pk'*rk);
    %Check teh convergence
    if ((norm(Lk2,"inf")<epsilon)&&(norm(ceq,"inf")<epsilon)&&(max(abs(ciq)) + epsilon >= 0)&&(min(lam_ineqk) + epsilon >= 0))
        break;
    end
    iteration=iteration+1;
    Xarray=[Xarray xk];
    normdLarray=[normdLarray norm(Lk2,"inf")];
    normceqarray=[normceqarray norm(ceq,"inf")];
end
[fval,~,~]=Hg_f(xk);
x=xk;
output.fval=fval;
output.xarray=Xarray;
output.dL=normdLarray;
output.ceq=normceqarray;
output.iteration=iteration;
end