function [x,output]=Sqp_Bfgs_trust(x0,lam_eq0,lam_ineq0)
% Sqp_Bfgs_trust   SQP algorithm with damped BFGS and line search
%
%          min  f(x)
%           x
%          s.t. ceq(x)=0  (Lagrange multiplier lam_eq)    
%               ciq(x)>=0  (Lagrange multiplier lam_ineq)  
% Syntax: [x,output]=Sqp_Bfgs_trust(x0,lam_eq0,lam_ineq0)
%         output.fval: minimum value
%         output.dL: convergence of delta_L
%         output.ceq: convergence of ceq
%         output.Xarray=Xarray: Iteration trajectory  
epsilon=1e-6;
ro=0.99;
gamma=0.9;
eta=0.2;
eta2=0.2;
eta3=0.5;
deltak=0.9;

deltaAarray=[];
deltaAarray=[deltaAarray deltak];
mu=30;
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
iteration_max=80;
iter_line_max=30;
while(iteration<iteration_max)
    %Solve the subproblem(QP)
    %Linearization of the elastic-mode formulation
    
    H_L1=[Bk zeros(2,2) zeros(2,2);zeros(4,6)];
    df_L1=[df;mu;mu;mu;mu];
    dciq_L1=[dciq' [1 0 0 0;0 1 0 0];zeros(4,2) eye(4);[1 0 0 0 0 0];[0 1 0 0 0 0];[-1 0 0 0 0 0];[0 -1 0 0 0 0]];
    ciq_L1=[-ciq;zeros(4,1);-deltak;-deltak;-deltak;-deltak];
    dceq_L1=[dceq' 0 0 -1 1];
    ceq_L1=[-ceq];
    [p_all,~,~,~,lamk] = quadprog(H_L1,df_L1,-dciq_L1,-ciq_L1,dceq_L1,ceq_L1);
    p=p_all(1:2);
    lamineq=lamk.ineqlin(1:2); 
    
    lameq=lamk.eqlin;
    lam_eqk_hat=-lameq;
    lam_ineqk_hat=-lamineq;
    lam_eqk=lam_eqk_hat;
    lam_ineqk=lam_ineqk_hat;
    
    %evaluate the gradient
    [f_new,~,~]=Hg_f(xk+p);
    [ceq_new,~,~]=Hg_ceq(xk+p);
    [ciq_new,~,~]=Hg_ciq(xk+p);
    
   
    mp_eq=norm((ceq+dceq'*p),1);
    mp_ineq=norm(min(0,(ciq+dciq'*p)),1);
    mp_k=mp_eq+mp_ineq;
    %Penalty update and step computation
    if mp_k==0
        mu=mu;
    else
        p_norm=norm(p,"Inf");
        mp_eq_inf=norm((ceq+dceq'*p),1);
        mp_ineq_inf=norm(min(0,(ciq+dciq'*p)),1);
        mp_k_inf=mp_eq+mp_ineq;
        if mp_k_inf==0
            mu=mu*(1+eta2);
        else
        if norm(ceq,1)+norm(min(0,ciq),1)-mp_k>=eta3*(norm(ceq,1)+norm(min(0,ciq),1)-mp_k_inf)
            mu=mu;
        else
            mu=mu*(1+eta2);
        end
        end
    end
    %calculate the ratio
    qmu_p=f+df'*p+0.5*p'*Bk*p+mu*mp_eq+mu*mp_ineq;
    qmu_zero=f+mu*norm((ceq),1)+mu*norm(min(0,ciq),1);
    predk=qmu_zero-qmu_p;
   
    phi1=f+mu*norm(ceq,1)+mu*norm(min(0,ciq),1);
    phi1_new=f_new+mu*norm(ceq_new,1)+mu*norm(min(0,ciq_new),1);
    aredk=phi1-phi1_new;
    rok=aredk/predk;
    %judge the step is accepted or rejected
    if rok>eta
        xk=xk+p;
        deltak=1.5*deltak;
    else
        deltak=gamma*norm(p,"Inf");
    end
   
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
    deltaAarray=[deltaAarray deltak];
end
[fval,~,~]=Hg_f(xk);
x=xk;
output.fval=fval;
output.xarray=Xarray;
output.dL=normdLarray;
output.ceq=normceqarray;
output.iteration=iteration;
output.delta=deltaAarray;
end