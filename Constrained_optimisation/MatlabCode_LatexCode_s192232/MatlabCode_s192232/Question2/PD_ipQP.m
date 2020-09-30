function [x,output]=PD_ipQP(H,g,A,b,C,d,x0,y0,z0,s0)
% PD_ipQP   Primal-dual interior-point algorithm
%
%          min  0.5*x'*H*x+g'*x
%           x
%          s.t. A x  = b      
%               C x >= d      
%         rL = Hx + g - Ay - Cz = 0 
%         rA = Ax + b = 0 (Lagrange multiplier y)
%         rC = Cx + s + d = 0 (Lagrange multiplier z)
%         s ¡Ý 0 (Slack variables )
%         sz = 0
% Syntax: [x,output]=PD_ipQP(H,g,A,b,C,d,x0,y0,z0,s0)
%         output.fval: minimum value
%         output.y: final y
%         output.s: final s
%         output.z: final z
%         output.Xarray: Iteration trajectory  
          
iteration_max=50;
epsilon=1.0e-9;
eta=0.995;
noeq=0;
Xarray=[];
Xarray=[Xarray x0];
nh=size(H,1);%dimension of x
na=size(A,2);%A'x=b
nc=size(C,2);%C'x>=d
e=ones(nc,1);
%residual
S0=diag(s0);
Z0=diag(z0);
%if no equlity equations
if isempty(A)&isempty(b)
    noeq=1;
end

if  noeq==1
    rL=H*x0+g-C*z0;
    rA=0.*x0;
    rC=s0+d-C'*x0;
    rsz=S0*Z0*e;
else
    rL=H*x0+g-A*y0-C*z0;
    rA=b-A'*x0;
    rC=s0+d-C'*x0;
    rsz=S0*Z0*e;
end

%duality gap
mu=(z0'*s0)/nc;

stop_flag=0;
iteration=0;
x=x0;
y=y0;
S=S0;
Z=Z0;
while(~stop_flag&iteration<=iteration_max)
    
    H_bar=H+C*(inv(S)*Z)*C';
    if noeq==1
     KKT=H_bar;
    else
    KKT=[H_bar -A;-A' zeros(na,na)];
    end
    
    [L, D, p] = ldl(KKT, 'lower', 'vector');
    %Affine Direction
    rL_bar=rL-C*(inv(S)*Z)*(rC-inv(Z)*rsz);
    if noeq==1
        KKT_b=-rL_bar;
    else
    KKT_b=-[rL_bar;rA];
    end
    delta_xy=zeros(size(KKT_b,1),1);
    delta_xy(p) = L'\(D\(L\KKT_b(p)));
    delta_x=delta_xy(1:nh);
    delta_z=-(inv(S)*Z)*C'*delta_x+(inv(S)*Z)*(rC-inv(Z)*rsz);
    delta_s=-inv(Z)*rsz-inv(Z)*S*delta_z;
    
    %compute the lagest alf
    z=diag(Z);
    s=diag(S);
    delta_zs=[delta_z;delta_s];
    alf_c=-[z;s]./delta_zs;
    alf_aff=min([1;alf_c(delta_zs<0)]);
    
    %compute the affine duality gap
    mu_aff=(z+alf_aff*delta_z)'*(s+alf_aff*delta_s)/nc;
    %compute the centering parameter
    sigma=(mu_aff/mu)^3;
    
    %Affine-Centering-Correction Direction
    delta_S=diag(delta_s);
    delta_Z=diag(delta_z);
    rsz_bar=rsz+delta_S*delta_Z*e-sigma*mu*e;
    %update Affine Direction
    rL_bar=rL-C*(inv(S)*Z)*(rC-inv(Z)*rsz_bar);
    if noeq==1
        KKT_b=-rL_bar;
    else
    KKT_b=-[rL_bar;rA];
    end
    delta_xy=zeros(size(KKT_b,1),1);
    delta_xy(p) = L'\(D\(L\KKT_b(p)));
    delta_x=delta_xy(1:nh);
    delta_y=delta_xy(nh+1:end);
    delta_z=-(inv(S)*Z)*C'*delta_x+(inv(S)*Z)*(rC-inv(Z)*rsz_bar);
    delta_s=-inv(Z)*rsz_bar-inv(Z)*S*delta_z;
    %update alf
    delta_zs=[delta_z;delta_s];
    alf_c=-[z;s]./delta_zs;
    alf_aff=min([1;alf_c(delta_zs<0)]);
    alf_bar=alf_aff*eta;
    %compute new state
    x=x+alf_bar*delta_x;
    y=y+alf_bar*delta_y;
    z=z+alf_bar*delta_z;
    s=s+alf_bar*delta_s;
    Z=diag(z);
    S=diag(s);
    %compute residuals
    if noeq==1
        rL=H*x+g-C*z;
        rA=0.*x;
        rC=s+d-C'*x;
        rsz=S*Z*e;
    else
        rL=H*x+g-A*y-C*z;
        rA=b-A'*x;
        rC=s+d-C'*x;
        rsz=S*Z*e;
    end
    mu=z'*s/nc;
    %stop judge
    judge = [norm(rL,1), norm(rA,1), norm(rC,1), abs(mu)];
    stop_flag = (length(judge(judge < epsilon)) == 4);
    Xarray=[Xarray x];
    iteration=iteration+1;
end
fval=0.5*x'*H*x+g'*x;
output.fval=fval;
output.y=y;
output.s=s;
output.z=z;
output.Xarray=Xarray;
output.iteration=iteration;
end