function u_new=MPCcompute_linear(R,Y,D,U,Ad,Bd,Bd_d,Czss,N)
%%%%需要模型信息，A,B,C,D...
num_x=4;
num_u=2;
num_d=2;
num_y=4;
num_z=2;
Phi=zeros(N*num_z,num_x);%[CA]*x0 ->>> z #####[num_z*N,num_x]
Lam=zeros(N*num_z,N*num_u);%[H]*Uk ->>> z
Lam_d=zeros(N*num_z,N*num_d);%[Hd]*Dk ->>> z

%预测值与参考值误差的权重矩阵Qz,海森矩阵Hr,一阶梯度Mr
Qz=diag(ones(N*num_z,1)) ;
%input rate的权重矩阵S,海森矩阵Hs,一阶梯度Mu_delta
Sz=diag(10*ones(num_u,1));
Hs=zeros(N*num_u,N*num_u);
Mu_delta=zeros(N*num_u,num_u);
Mu_delta(1:num_u,:)=-Sz;

for i=1:N
    %fill the Phi
    Phi((i-1)*num_z+1:i*num_z,:)=Czss*Ad^i;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %fill the Lam
    lam_i=zeros(num_z,N*num_u);%[num_x * num_u*N]
    for j=1:i
        lam_i(:,(j-1)*num_u+1:j*num_u)=Czss*(Ad^(i-j))*Bd;%[Hn,Hn-1,Hn-2]
    end
    Lam((i-1)*num_z+1:i*num_z,:)=lam_i;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %fill the Lam_d
    lam_di=zeros(num_z,N*num_d);%[num_x * num_u*N]
    for j=1:i
        lam_di(:,(j-1)*num_d+1:j*num_d)=Czss*(Ad^(i-j))*Bd_d;%[Hn,Hn-1,Hn-2]
    end
    Lam_d((i-1)*num_z+1:i*num_z,:)=lam_di;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %fill the Hs
    if i==1
        Hs(1:num_u,1:2*num_u)=[2*Sz,-Sz];
    elseif i==N
        Hs(end-num_u+1:end,end-2*num_u+1:end)=[-Sz,Sz];
    else
        Hs((i-1)*num_u+1:(i)*num_u,(i-2)*num_u+1:(i+1)*num_u)=[-Sz,2*Sz,-Sz];
    end
end

Mx0=Lam'*Qz*Phi;
Md=Lam'*Qz*Lam_d;
Mr=-Lam'*Qz;
Hr=Lam'*Qz*Lam;
Hu=Hr+Hs;
Lx0=-Hu\Mx0;
Lr=-Hu\Mr;
Lu_delta=-Hu\Mu_delta;
Ld=-Hu\Md;
Kx0=Lx0(1:num_u,:);
Kr=Lr(1:num_u,:);
Ku_delta=Lu_delta(1:num_u,:);
Kd=Ld(1:num_u,:);

u_new=Kx0*Y+Kr*R+Ku_delta*U+Kd*D;
end