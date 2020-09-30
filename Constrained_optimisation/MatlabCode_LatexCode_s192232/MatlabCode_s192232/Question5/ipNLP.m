function [x,output]=ipNLP(x0,y0,z0,s0)
% ipNLP   Primal-Dual Interior Point Algorithms for NLP
%
%          min  f(x)
%           x
%          s.t. ceq(x)=0  (Lagrange multiplier y)    
%               ciq(x)-s=0  (Lagrange multiplier z)
%               s>=0
%          rdL=df-dceq*y-dciq*z;
%          rsz=S*z-tau*e;
%          rceq=ceq;
%          rciq=ciq-s;
% Syntax: [x,output]=ipNLP(x0,y0,z0,s0)
%         output.fval: minimum value
%         output.Xarray=Xarray: Iteration trajectory  
epsilon=1.0e-5;
tau0=0.8;
eta=0.005;
ro=0.5;
Xarray=[];
Xarray=[Xarray x0];

nx=size(x0,1);
neq=size(y0,1);
niq=size(z0,1);
e=ones(niq,1);
%evaluate the gradient and hessian
[f,df,d2f]=Hg_f(x0);
[ceq,dceq,d2ceq]=Hg_ceq(x0);
[ciq,dciq,d2ciq]=Hg_ciq(x0);

xk=x0;
yk=y0;
Z0=diag(z0);
S0=diag(s0);
Zk=Z0;
Sk=S0;
tauk=tau0;
iteration=0;
iteration_out=0;
iteration_max=700;
iteration_max_out=700;
while(iteration_out<iteration_max_out)
    %convergence judge
    
    rdL=df-dceq*yk-dciq*diag(Zk);
    rsz=Sk*diag(Zk)-0*e;
    rceq=ceq;
    rciq=ciq-diag(Sk);
    Error=max([norm(rdL,1),norm(rsz,1),norm(rceq,1),norm(rciq,1)]);
    if Error<epsilon
        break;
    end
    while(iteration<iteration_max)
        %solve the subproblem
        [f,df,d2f]=Hg_f(xk);
        [ceq,dceq,d2ceq]=Hg_ceq(xk);
        [ciq,dciq,d2ciq]=Hg_ciq(xk);
        sk=diag(Sk);
        zk=diag(Zk);
        sumeq=0;
        sumiq=0;
        for i=1:length(yk)
            sumeq=sumeq+yk(i)*d2ceq(:,:,i);
        end
        for i=1:length(zk)
            sumiq=sumiq+zk(i)*d2ciq(:,:,i);
        end
        ddL=d2f-sumeq-sumiq;
        sumk=inv(Sk)*Zk;
        KKT_A=[ddL,zeros(nx,niq),dceq,dciq;zeros(niq,nx),sumk,zeros(niq,neq),-eye(niq);dceq',zeros(neq,niq),zeros(neq,neq),zeros(neq,niq);dciq',-eye(niq),zeros(niq,neq),zeros(niq,niq)];
        rdL=df-dceq*yk-dciq*zk;
        rsz=zk-tauk.*inv(Sk)*e;
        rceq=ceq;
        rciq=ciq-sk;
        KKT_b=-[rdL;rsz;rceq;rciq];
        [L,D,p]=ldl(KKT_A,'vector');
        p_all=zeros(size(KKT_b,1),1);
        p_all(p)=L'\(D\(L\KKT_b(p)));
        px=p_all(1:nx);
        ps=p_all(nx+1:nx+niq);
        py=-p_all(nx+niq+1:nx+niq+neq);
        pz=-p_all(nx+niq+neq+1:end);
        %compute the step length
        alpha_tmp_s=(eta.*sk-sk)./ps;
        alpha_s=min([1;alpha_tmp_s(ps<0)]);
        alpha_tmp_z=(eta.*zk-zk)./pz;
        alpha_z=min([1;alpha_tmp_z(pz<0)]);
        %Update x_k+1, y_k+1, s_k+1, z_k+1
        xk=xk+alpha_s*px;
        sk=sk+alpha_s*ps;
        yk=yk+alpha_z*py;
        zk=zk+alpha_z*pz;
        Sk=diag(sk);
        Zk=diag(zk);
        iteration=iteration+1;
        Xarray=[Xarray xk];
        %convergence judge
        rdL=df-dceq*yk-dciq*zk;
        rsz=Sk*zk-tauk.*e;
        rceq=ceq;
        rciq=ciq-diag(Sk);
        Error=max([norm(rdL,1),norm(rsz,1),norm(rceq,1),norm(rciq,1)]);
        if Error<tauk
            break;
        end
    end
    tauk=ro*tauk;
    iteration_out=iteration_out+1;
end
x=xk;
[fval,~,~]=Hg_f(xk);
output.xarray=Xarray;
output.iteration=iteration;
output.fval=fval;
end