function s_list=ode_traj(s0,t_list,Ts,u_angle_list,u_engin_list,u_break_list)
x=s0(1);
y=s0(2);
theta=s0(3);
v=s0(4);
L=1.516;
s_list=zeros(4,0);
for i=1:1:size(t_list,1)
    t=Ts;
    u_angle=u_angle_list(i);
    u_engin=u_engin_list(i);
    u_break=u_break_list(i);
    v_last=v;
    theta_last=theta;
    theta=theta_last+t*(v_last/L)*tan(u_angle);
    v=v_last+t*(u_engin+u_break);
    x=x+t*v_last*cos(theta_last);
    y=y+t*v_last*sin(theta_last);
    s_list=[s_list [x;y;theta;v]];
end
end