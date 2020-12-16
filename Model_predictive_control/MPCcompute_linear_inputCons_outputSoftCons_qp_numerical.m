function u_new=MPCcompute_linear_inputCons_outputSoftCons_qp_numerical(R,Y,D,U,Ad,Bd,Bd_d,Css,u_min,u_max,u_delta_min,u_delta_max,z_min,z_max,N,Q_cof,u_delta_cof,eta_cof)
%%%%需要模型信息，A,B,C,D...
num_x=4;
num_u=2;
num_d=2;
num_y=4;
num_z=2;
num_eta=num_y;
num_u_delta_cons=N*num_u;
num_z_cons=N*num_y;
Phi=zeros(N*num_z,num_x);%[CA]*x0 ->>> y #####[num_z*N,num_x]
Gam=zeros(N*num_y,N*num_u+N*num_eta);%[H]*[Uk,etak] ->>> y
Gam_d=zeros(N*num_y,N*num_d);%[Hd]*Dk ->>> y

%预测值与参考值误差的权重矩阵Qz,海森矩阵Hr,一阶梯度Mr
Qz=diag(zeros(N*num_y,1));
%input rate的权重矩阵S,海森矩阵Hs,一阶梯度Mu_delta
Sz=diag(u_delta_cof*ones(num_u,1));
Hs=zeros(N*num_u+N*num_eta,N*num_u+N*num_eta);
Mu_delta=zeros(N*num_u+N*num_eta,num_u);
Mu_delta(1:num_u,:)=-Sz;
H_eta=diag([zeros(N*num_u,1);eta_cof*ones(N*num_eta,1)]);
U0=zeros(N*num_u+N*num_eta,1);
Iu=diag(ones(num_u,1));
Au_delta_cons=zeros(num_u_delta_cons,N*num_u+N*num_eta);

for i=1:N
    %fill the Phi
    Phi((i-1)*num_y+1:i*num_y,:)=Css*Ad^i;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %fill the Gam
    Gam_i=zeros(num_y,N*num_u+N*num_eta);%[num_x * num_u*N]
    for j=1:i
        Gam_i(:,(j-1)*num_u+1:j*num_u)=Css*(Ad^(i-j))*Bd;%[Hn,Hn-1,Hn-2]
    end
    Gam((i-1)*num_y+1:i*num_y,:)=Gam_i;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %fill the Gam_d
    Gam_di=zeros(num_y,N*num_d);%[num_x * num_u*N]
    for j=1:i
        Gam_di(:,(j-1)*num_d+1:j*num_d)=Css*(Ad^(i-j))*Bd_d;%[Hn,Hn-1,Hn-2]
    end
    Gam_d((i-1)*num_y+1:i*num_y,:)=Gam_di;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %fill the Qz
    qz=diag([Q_cof*ones(num_z,1);zeros(num_y-num_z,1)]);
    Qz((i-1)*num_y+1:i*num_y,(i-1)*num_y+1:i*num_y)=qz;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %fill the Hs
    if i==1
        Hs(1:num_u,1:2*num_u)=[2*Sz,-Sz];
    elseif i==N
        Hs((i-1)*num_u+1:(i)*num_u,(i-2)*num_u+1:i*num_u)=[-Sz,Sz];
    else
        Hs((i-1)*num_u+1:(i)*num_u,(i-2)*num_u+1:(i+1)*num_u)=[-Sz,2*Sz,-Sz];
    end
    %fill the U0
    U0((i-1)*num_u+1:i*num_u,:)=U;
    %fill the Au_delta_cons
    if i==1
        Au_delta_cons(1:num_u,1:num_u)=Iu;
    else
        Au_delta_cons((i-1)*num_u+1:(i)*num_u,(i-2)*num_u+1:i*num_u)=[-Iu,Iu];
    end
end

Mx0=Gam'*Qz*Phi;
Md=Gam'*Qz*Gam_d;
Mr=-Gam'*Qz;
Hr=Gam'*Qz*Gam;
Hu=Hr+Hs+H_eta;
g=Mx0*Y+Mr*R+Mu_delta*U+Md*D;
Az_cons_min=Gam;
Az_cons_min(:,N*num_u+1:end)=diag(1*ones(N*num_eta,1));
Az_cons_max=Gam;
Az_cons_max(:,N*num_u+1:end)=diag(-1*ones(N*num_eta,1));
Z_min=z_min*ones(num_z_cons,1)-Phi*Y-Gam_d*D;
Z_max_fake=5000*ones(num_z_cons,1);
Z_max=z_max*ones(num_z_cons,1)-Phi*Y-Gam_d*D;
Z_min_fake=zeros(num_z_cons,1);

U_min=[u_min*ones(N*num_u,1);zeros(N*num_eta,1)];
U_max=[u_max*ones(N*num_u,1);5000*ones(N*num_eta,1)];%eta上限给个很大的值就可以
U_delta_min=u_delta_min*ones(num_u_delta_cons,1);
U_delta_max=u_delta_max*ones(num_u_delta_cons,1);
U_delta_min(1)=u_delta_min+U(1);
U_delta_min(2)=u_delta_min+U(2);
U_delta_max(1)=u_delta_max+U(1);
U_delta_max(2)=u_delta_max+U(2);

A_ieq=[Au_delta_cons;Az_cons_min;Az_cons_max];
b_ieq_min=[U_delta_min;Z_min;Z_min_fake];
b_ieq_max=[U_delta_max;Z_max_fake;Z_max];

[u_opt,~] = qpsolver(0.5*(Hu+Hu'),g,U_min,U_max,A_ieq,b_ieq_min,b_ieq_max,U0);
u_new=u_opt(1:num_u);

end