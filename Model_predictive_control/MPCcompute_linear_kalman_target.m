function u_new=MPCcompute_linear_kalman_target(R,X0,D,U,Us,Ad,Bd,Bd_d,Czss,N,Q_cof,u_delta_cof,Wu_cof)
%%%%需要模型信息，A,B,C,D...
num_x=4;
num_u=2;
num_d=2;
num_y=4;
num_z=2;
Phi=zeros(N*num_z,num_x);%[CA]*x0 ->>> z #####[num_z*N,num_x]
Phi_d=zeros(N*num_z,num_d);
Gam=zeros(N*num_z,N*num_u);%[H]*Uk ->>> z
%Gam_d=zeros(N*num_z,N*num_d);%[Hd]*Dk ->>> z
%与输入参考值误差
Hu_ref=diag(Wu_cof*ones(N*num_u,1));
%预测值与参考值误差的权重矩阵Qz,海森矩阵Hr,一阶梯度Mr
Qz=diag(Q_cof*ones(N*num_z,1)) ;
%input rate的权重矩阵S,海森矩阵Hs,一阶梯度Mu_delta
Sz=diag(u_delta_cof*ones(num_u,1));
Hs=zeros(N*num_u,N*num_u);
Mu_delta=zeros(N*num_u,num_u);
Mu_delta(1:num_u,:)=-Sz;
U0=zeros(N*num_u,1);
for i=1:N
    %fill the Phi
    Phi((i-1)*num_z+1:i*num_z,:)=Czss*Ad^i;
    Phi_d((i-1)*num_z+1:i*num_z,:)=Czss*Ad^i*Bd_d;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %fill the Gam
    Gam_i=zeros(num_z,N*num_u);%[num_x * num_u*N]
    for j=1:i
        Gam_i(:,(j-1)*num_u+1:j*num_u)=Czss*(Ad^(i-j))*Bd;%[Hn,Hn-1,Hn-2]
    end
    Gam((i-1)*num_z+1:i*num_z,:)=Gam_i;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %fill the Gam_d
%     Gam_di=zeros(num_z,N*num_d);%[num_x * num_u*N]
%     for j=1:i
%         Gam_di(:,(j-1)*num_d+1:j*num_d)=Czss*(Ad^(i-j))*Bd_d;%[Hn,Hn-1,Hn-2]
%     end
%     Gam_d((i-1)*num_z+1:i*num_z,:)=Gam_di;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %fill the Hs
    if i==1
        Hs(1:num_u,1:2*num_u)=[2*Sz,-Sz];
    elseif i==N
        Hs(end-num_u+1:end,end-2*num_u+1:end)=[-Sz,Sz];
    else
        Hs((i-1)*num_u+1:(i)*num_u,(i-2)*num_u+1:(i+1)*num_u)=[-Sz,2*Sz,-Sz];
    end
    %fill the U0
    U0((i-1)*num_u+1:i*num_u,:)=U;
end

Mx0=Gam'*Qz*Phi;
Md=Gam'*Qz*Phi_d;
Mr=-Gam'*Qz;
Mu_ref=-Hu_ref;
Hr=Gam'*Qz*Gam;
Hu=Hr+Hs;
g=Mx0*X0+Mr*R+Mu_delta*U+Md*D+Mu_ref*Us;
[u_opt,~] = qpsolver(0.5*(Hu+Hu'),g,[],[],[],[],[],U0);
u_new=u_opt(1:num_u);

end