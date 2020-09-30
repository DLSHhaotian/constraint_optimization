function [x,output]=Sqp_Bfgs_Linesearch(x0,lam_eq0,lam_ineq0)
% Sqp_Bfgs_Linesearch   SQP algorithm with damped BFGS and line search
%
%          min  f(x)
%           x
%          s.t. ceq(x)=0  (Lagrange multiplier lam_eq)    
%               ciq(x)>=0  (Lagrange multiplier lam_ineq)  
% Syntax: [x,output]=Sqp_Bfgs_Linesearch(x0,lam_eq0,lam_ineq0)
%         output.fval: minimum value
%         output.dL: convergence of delta_L
%         output.ceq: convergence of ceq
%         output.Xarray=Xarray: Iteration trajectory  
epsilon=1e-6;
ro=0.99;
eta=0.55;
tau=0.99;
Xarray=[];
normdLarray=[];
normceqarray=[];
Xarray=[Xarray x0];

nx=size(x0,1);
Bk=eye(nx);

[f,df,~]=Hg_f(x0);
[ceq,dceq,~]=Hg_ceq(x0);
[ciq,dciq,~]=Hg_ciq(x0);
xk=x0;
lam_eqk=lam_eq0;
lam_ineqk=lam_ineq0;
iteration=0;
iteration_max=200;
iter_line_max=30;
while(iteration<iteration_max)
    %Solve the subproblem(QP)
    [p,lameq,lamineq]=Qpsolver_Sqp(Bk,df,dceq,ceq,dciq,ciq);
    lam_eqk_hat=-lameq;
    lam_ineqk_hat=-lamineq;
    p_lameqk=lam_eqk_hat-lam_eqk;
    p_lamineqk=lam_ineqk_hat-lam_ineqk;
    %Line search
    alpha=1;
    mu=(df'*p+0.5*p'*Bk*p)/((1-ro)*(norm(ceq,1)+norm(ciq,1)));
    [f_new,~,~]=Hg_f(xk+alpha*p);
    [ceq_new,~,~]=Hg_ceq(xk+alpha*p);
    [ciq_new,~,~]=Hg_ciq(xk+alpha*p);
    iter_line=1;
    while iter_line<iter_line_max
        phi0=f+mu*norm(ceq,1)+mu*norm(min(0,ciq),1);
        dphi0=df'*p-mu*norm(ceq,1)-mu*norm(min(0,ciq),1);
        phi_alpha=f_new+mu*norm(ceq_new,1)+mu*norm(min(0,ciq_new),1);
        if phi_alpha<=(phi0+eta*alpha*dphi0)
            break;
        else
            alpha=tau*alpha;
            [f_new,~,~]=Hg_f(xk+alpha*p);
            [ceq_new,~,~]=Hg_ceq(xk+alpha*p);
            [ciq_new,~,~]=Hg_ciq(xk+alpha*p);
            iter_line=iter_line+1;
        end
    end
    %Take step and update x(k+1),lameq(k+1),lamineq(k+1)
    xk=xk+alpha*p;
    lam_eqk=lam_eqk+alpha*p_lameqk;
    lam_ineqk=lam_ineqk+alpha*p_lamineqk;
    %Calculate dL(k) with x(k) and lameq(k+1),lamineq(k+1)
    Lk=df-dceq.*lam_eqk-dciq*lam_ineqk;
    %Re-evaluate
    [f,df,~]=Hg_f(xk);
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