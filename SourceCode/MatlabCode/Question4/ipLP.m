function [x,output]=ipLP(g,A,b)
% ipLP   Primal-dual interior-point algorithm
%
%          min  g'*x
%           x
%          s.t. A x  = b      
%               x >= 0      
%         rc = g -A*lambda - S = 0 
%         rb = Ax ? b = 0 
%         XSe=0
%         rA = Ax + b = 0 (Lagrange multiplier y)
%         rC = Cx + s + d = 0 (Lagrange multiplier z)
%         s ¡Ý 0 (Slack variables )
%         sz = 0
% Syntax: [x,output]=ipLP(g,A,b)
%         output.fval: minimum value
%         output.lam : final lambda
%         output.s : final s
%         output.z: final z
%         output.Xarray: Iteration trajectory    
iteration_max=30;

epsilon=1.0e-6;
eta=0.99;
stop_flag=0;
nx=size(g,1);%x
nc=size(b,1);%c
%starting point
e=ones(nx,1);
x_hat=A'*inv(A*A')*b;
lam_hat=inv(A*A')*A*g;
s_hat=g-A'*lam_hat;

delta_x=max(-(3/2)*min(x_hat),0);
delta_s=max(-(3/2)*min(s_hat),0);

x_hat=x_hat+delta_x*e;
s_hat=s_hat+delta_s*e;

delta_x_hat=0.5*(x_hat'*s_hat)/(e'*s_hat);
delta_s_hat=0.5*(x_hat'*s_hat)/(e'*x_hat);

x0=x_hat+delta_x_hat*e;
lam0=lam_hat;
s0=s_hat+delta_s_hat*e;

Xarray=[];
Xarray=[Xarray x0];
%Predictor-Corrector algorithm
rc=A'*lam0+s0-g;
rb=A*x0-b;

x=x0;
lam=lam0;
s=s0;
iteration=0;
while(~stop_flag&iteration<=iteration_max)
    %solving the problem to obtain delta_x_aff,delta_lam_aff and delta_s_aff
    KKT_A=[zeros(nx,nx) A' eye(nx,nx);A zeros(nc,nc) zeros(nc,nx);diag(s) zeros(nx,nc) diag(x)];
    XSe=diag(x)*diag(s)*e;
    KKT_b=[-rc;-rb;-XSe];
    [L,D,p]=ldl(KKT_A, 'lower', 'vector');
    delta_aff=zeros(size(KKT_b,1),1);
    delta_aff(p)=L'\(D\(L\KKT_b(p)));
    delta_x_aff=delta_aff(1:nx);
    delta_lam_aff=delta_aff(nx+1:nx+nc);
    delta_s_aff=delta_aff(nx+nc+1:end);
    
    %calculate the alf_pri and alf_dual
    xi_deltax=-x./delta_x_aff;
    alf_pri=min([1;xi_deltax(delta_x_aff<0)]);
    si_deltas=-s./delta_s_aff;
    alf_dual=min([1;si_deltas(delta_s_aff<0)]);
    
    %calculate mu_aff and dual gap
    mu_aff=(x+alf_pri*delta_x_aff)'*(s+alf_dual*delta_s_aff)/nx;
    mu=(x'*s)/nx;
    %calculate centering parameter
    sigma=(mu_aff/mu)^3;
    %solving the problem to obtain delta_x_step, delta_lam_step and delta_s_step
    XSe_aff=-(diag(x)*diag(s)*e)-diag(delta_x_aff)*diag(delta_s_aff)*e+sigma*mu*e;
    KKT_b=[-rc;-rb;XSe_aff];
    delta_step=zeros(size(KKT_b,1),1);
    delta_step(p)=L'\(D\(L\KKT_b(p)));
    delta_x_step=delta_step(1:nx);
    delta_lam_step=delta_step(nx+1:nx+nc);
    delta_s_step=delta_step(nx+nc+1:end);
    
    %step length
    xk_deltaxk=-x./delta_x_step;
    alf_pri_k_max=min(xk_deltaxk(delta_x_step<0));
    alf_pri_k=min([1;eta*alf_pri_k_max]);
    
    sk_deltask=-s./delta_s_step;
    alf_dual_k_max=min(sk_deltask(delta_s_step<0));
    alf_dual_k=min([1;eta*alf_dual_k_max]);
    
    x=x+alf_pri_k*delta_x_step;
    lam=lam+alf_dual_k*delta_lam_step;
    s=s+alf_dual_k*delta_s_step;
    
    rc_norm=norm(rc,2);
    ra_norm=norm(rb,2);
    dual_gap=abs(x'*s/nx);
     %stop judge
    judge = [rc_norm;ra_norm;dual_gap];
    stop_flag = (length(judge(judge < epsilon)) == 3);
    %update rc and rb
    rc=(1-alf_dual_k)*rc;
    rb=(1-alf_pri_k)*rb;
    iteration=iteration+1;
    Xarray=[Xarray x];
end
fval=g'*x;
output.fval=fval;
output.lam=lam;
output.s=s;
output.iteration=iteration;
output.xarray=Xarray;
end