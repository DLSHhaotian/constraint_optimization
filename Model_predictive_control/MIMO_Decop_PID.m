function [u,I] = MIMO_Decop_PID(ubar,ybar,y,yold,I,KP,KI,KD,dt,umin,umax,gain_decop)

e = ybar-y;
p_decop=gain_decop*e;
P = KP*p_decop;
D = -KD*(y-yold)*p_decop/dt;
u = ubar + P + I + D;

if (u >= umax)
u = umax;
elseif (u <= umin)
u = umin;
else
I = I + KI*p_decop*dt;
end