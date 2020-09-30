function [p,lameq,lamineq]=Qpsolver_Sqp(H,df,dceq,ceq,dciq,ciq)
%quadratic programming solver for subproblem in SQP method

[p,~,~,~,lamk] = quadprog(H,df,-dciq',ciq,dceq',-ceq);
lamineq=lamk.ineqlin;
lameq=lamk.eqlin;
end