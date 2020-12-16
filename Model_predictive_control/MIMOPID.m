function [u,I] = MIMOPID(ubar,ybar,y,yold,I,KP,KI,KD,dt,umin,umax,decop)

e = ybar-y;
e_decop=decop*e;
P = KP*e_decop;
D = -KD*(y-yold)/dt;
u = ubar + P + I + D;
flag=1;
for n=1:2
if (u(n) >= umax)
u(n) = umax;
flag=0;
elseif (u(n) <= umin)
u(n) = umin;
flag=0;
end
end
if(flag==1)
I = I + KI*e_decop*dt;
end