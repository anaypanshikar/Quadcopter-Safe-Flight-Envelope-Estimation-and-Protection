close all
clear all 
clc

global V theta0 Omega_f0 Omega_b0 Vx0 Vz0 T N h N_traj p_c m Iy k_r

%% Trim conditions
V = 8;
theta0 = -13.1*pi/180;
Omega_f0 = 709.5;
Omega_b0 = 860.4;
Vx0 = V*cos(theta0);        %7.792
Vz0 = V*sin(theta0);       %1.813

%% Monte carlo method variables
T = 0.15;
N = 100;
h = T/N;         %simulation time step
N_traj = 1000;
p_c = 0.1;
m = 0.375;
Iy = 1.3*10^-3;
k_r = 1;
%% simulation

for i = 1:N_traj
    [Vx_f,Vz_f,theta_f,q_f,Vel_f,gamma_f,Vx_b,Vz_b,theta_b,q_b,Vel_b,gamma_b,u1_f,u2_f] = reachability();
    envf(i,:) = [Vx_f,Vz_f,theta_f,q_f,Vel_f,gamma_f];
    envb(i,:) = [Vx_b,Vz_b,theta_b,q_b,Vel_b,gamma_b];
end    
fb = boundary(envf(:,3)*180/pi,envf(:,4)*180/pi,0.001);
bb = boundary(envb(:,3)*180/pi,envb(:,4)*180/pi,0.001);

figure()
fp = polyshape(envf(fb,3)*180/pi,envf(fb,4)*180/pi);
hold on
bp = polyshape(envb(bb,3)*180/pi,envb(bb,4)*180/pi);
polyout = intersect(fp,bp);
plot(polyout);
hold off
grid on
xlabel("theta (deg)")
ylabel("q (deg/s)")

[xv,yv] = boundary(polyout);
xq = -10;
yq = 0;
[in,on] = inpolygon(xq,yq,xv,yv);
in
on

% figure()
% fp = patch(envf(fb,3)*180/pi,envf(fb,4)*180/pi,[0.5843 0.3557 0.9982],"EdgeColor",[0.5843 0.3557 0.9982],"LineWidth",0.001);
% hold on
% bp = patch(envb(bb,3)*180/pi,envb(bb,4)*180/pi,[0.5843 0.9982 0.2557],"EdgeColor",[0.5843 0.9982 0.2557],"LineWidth",0.001);
% scatter(envf(:,3)*180/pi,envf(:,4)*180/pi,40,"b",".")
% scatter(envb(:,3)*180/pi,envb(:,4)*180/pi,40,"g",".")
% scatter(theta0*180/pi,0,100,"red","x","LineWidth",2)
% hold off
% alpha(fp,0.5)
% alpha(bp,0.5)
% grid on
% xlabel("theta (deg)")
% ylabel("q (deg/s)")
% legend("Forward Set","Backward Set",'Data forward','Data backward',"Trim","location","best")
% 
% figure()
% scatter(envf(:,5),envf(:,6)*180/pi,40,"b",".")
% hold on
% scatter(envb(:,5),envb(:,6)*180/pi,40,"g",".")
% scatter(V,0,100,"red","x","LineWidth",2)
% hold off
% grid on
% xlabel("Velocity (m/s)")
% ylabel("gamma (deg)")
% legend("Forward Set","Backward Set","Trim","location","best")
% 
% figure()
% scatter3(envf(:,5),envf(:,6)*180/pi,envf(:,3)*180/pi,40,"b",".")
% hold on
% scatter3(envb(:,5),envb(:,6)*180/pi,envb(:,3)*180/pi,40,"g",".")
% scatter(V,0,100,"red","x","LineWidth",2)
% hold off
% grid on
% zlabel("theta (deg)")
% ylabel("gamma (deg)")
% xlabel("V (m/s)")
% legend("Forward Set","Backward Set","location","best")
% 
% figure()
% scatter3(envf(:,6)*180/pi,envf(:,3)*180/pi,envf(:,4)*180/pi,40,"b",".")
% hold on
% scatter3(envb(:,6)*180/pi,envb(:,3)*180/pi,envb(:,4)*180/pi,40,"g",".")
% scatter(V,0,100,"red","x","LineWidth",2)
% hold off
% grid on
% ylabel("theta (deg)")
% xlabel("gamma (deg)")
% zlabel("q (deg/s)")
% legend("Forward Set","Backward Set","location","best")
%% functions

function [Vx_f,Vz_f,theta_f,q_f,Vel_f,gamma_f,Vx_b,Vz_b,theta_b,q_b,Vel_b,gamma_b,u1_f,u2_f] = reachability()
    global theta0 Omega_f0 Omega_b0 Vx0 Vz0 N h p_c m k_r
    
    Vx = Vx0;
    Vz = Vz0;
    theta = theta0;
    q = 0;
    u1 = 2*Omega_f0^2*k_r;
    u2 = 2*Omega_b0^2*k_r;
    u1reff = [2*0^2*k_r, 2*1200^2*k_r];
    u2reff = [2*0^2*k_r, 2*1200^2*k_r];
    u1ref = u1reff(randi(2));
    u2ref = u2reff(randi(2));
    longf = [Vx,Vz,theta,q,u1,u2];
    longb = [Vx,Vz,theta,q,u1,u2];
    Fxf = m*9.81*sin(theta0);
    Fzf = -m*9.81*cos(theta0);
    Myf = 0;
    Fxb = m*9.81*sin(theta0);
    Fzb = -m*9.81*cos(theta0);
    Myb = 0;
    
    for k = 1:N-1
        
        % forward analysis

        [k_11,k_12,k_13,k_14,k_15,k_16] = forward(longf,Fxf,Fzf,Myf,u1ref,u2ref);
        [k_21,k_22,k_23,k_24,k_25,k_26] = forward(longf+0.5*h.*[k_11,k_12,k_13,k_14,k_15,k_16],Fxf,Fzf,Myf,u1ref,u2ref);
        [k_31,k_32,k_33,k_34,k_35,k_36] = forward((longf+0.5*h.*[k_21,k_22,k_23,k_24,k_25,k_26]),Fxf,Fzf,Myf,u1ref,u2ref);
        [k_41,k_42,k_43,k_44,k_45,k_46] = forward((longf+[k_31,k_32,k_33,k_34,k_35,k_36].*h),Fxf,Fzf,Myf,u1ref,u2ref);
        longf = longf + (1/6).*([k_11,k_12,k_13,k_14,k_15,k_16]+2.*[k_21,k_22,k_23,k_24,k_25,k_26]+2.*[k_31,k_32,k_33,k_34,k_35,k_36]+[k_41,k_42,k_43,k_44,k_45,k_46]).*h;  % main equation
        
        [Cd1,Cd2,Cd3,Cd4,Cz0,Cz1,Cz2,Cm0,Cm1,Cm2] = param(longf(1),longf(2));
        [Fxf,Fzf,Myf] = FnM(Cd1,Cd2,Cd3,Cd4,Cz0,Cz1,Cz2,Cm0,Cm1,Cm2,longf(5),longf(6),longf(1),longf(2));
        
        %backward analysis

        [k_11,k_12,k_13,k_14,k_15,k_16] = backward(longb,Fxb,Fzb,Myb,u1ref,u2ref);
        [k_21,k_22,k_23,k_24,k_25,k_26] = backward(longb+0.5*h.*[k_11,k_12,k_13,k_14,k_15,k_16],Fxb,Fzb,Myb,u1ref,u2ref);
        [k_31,k_32,k_33,k_34,k_35,k_36] = backward((longb+0.5*h.*[k_21,k_22,k_23,k_24,k_25,k_26]),Fxb,Fzb,Myb,u1ref,u2ref);
        [k_41,k_42,k_43,k_44,k_45,k_46] = backward((longb+[k_31,k_32,k_33,k_34,k_35,k_36].*h),Fxb,Fzb,Myb,u1ref,u2ref);
        longb = longb + (1/6).*([k_11,k_12,k_13,k_14,k_15,k_16]+2.*[k_21,k_22,k_23,k_24,k_25,k_26]+2.*[k_31,k_32,k_33,k_34,k_35,k_36]+[k_41,k_42,k_43,k_44,k_45,k_46]).*h;  % main equation
        
        [Cd1,Cd2,Cd3,Cd4,Cz0,Cz1,Cz2,Cm0,Cm1,Cm2] = param(longb(1),longb(2));
        [Fxb,Fzb,Myb] = FnM(Cd1,Cd2,Cd3,Cd4,Cz0,Cz1,Cz2,Cm0,Cm1,Cm2,longb(5),longb(6),longb(1),longb(2));
        
        R = rand;
        if R <= (1-p_c^(1/N))
            if u1ref == 2*0^2*k_r
                u1ref = 2*1200^2*k_r;
            elseif u1ref == 2*1200^2*k_r
                u1ref = 2*0^2*k_r;
            end
        end
        R = rand;
        if R <= (1-p_c^(1/N))
            if u2ref == 2*1200^2*k_r
                u2ref = 2*0^2*k_r;
            elseif u2ref == 2*0^2*k_r
                u2ref = 2*1200^2*k_r;
            end
        end
        
    end
    
    Vx_f = longf(1);
    Vz_f = longf(2);
    theta_f = longf(3);
    q_f = longf(4);
    Vel_f = sqrt(longf(1)^2 + longf(2)^2);
    gamma_f = longf(3) - atan(longf(2)/longf(1));
    
    Vx_b = longb(1);
    Vz_b = longb(2);
    theta_b = longb(3);
    q_b = longb(4);
    Vel_b = sqrt(longb(1)^2 + longb(2)^2);
    gamma_b = longb(3) - atan(longb(2)/longb(1));
    
    u1_f = longf(5);
    u2_f = longf(6);
end

function [Cd1,Cd2,Cd3,Cd4,Cz0,Cz1,Cz2,Cm0,Cm1,Cm2] = param(Vx,Vz)
    Cd1 = -0.217;
    Cd2 = 0.0184;
    Cd3 = -9.61*10^-4;
    Cd4 = 6.17046*10^-2;
    Cz0 = 0.0298*abs(Vz)*Vz - 3.77*10^-3*Vz^3;
    Cz1 = 1.67 - 0.0858*Vx + 2.2*10^-3*Vx^2;
    Cz2 = 2.15 + 1.97*10^-2*Vx^2 + 0.0728*Vz - 6.84*10^-4*Vx^3 - 1.97*10^-4*Vx^3*Vz + 4.34*10^-3*Vx^2*Vz;
    Cm0 = 0.0103*Vx - 6.77*10^-4*Vx^2 + 8.64*10^-3*Vz + 7.17*10^-5*Vx^2*Vz + 2.63*10^-4*Vx*Vz^2;
    Cm1 = 0.152 + 1.04*10^-3*Vx^2 + 1.66*10^-3*Vx*Vz - 1.86*10^-3*Vx;
    Cm2 = -0.163 + 8.04*10^-3*Vx - 2.11*10^-4*Vx*Vz - 6.31*10^-4*Vx^2;
    
end

function [Fx,Fz,My] = FnM(Cd1,Cd2,Cd3,Cd4,Cz0,Cz1,Cz2,Cm0,Cm1,Cm2,u1,u2,Vx,Vz)
    Fx = Cd1*Vx + Cd2*Vx^2 + Cd3*Vx^3 + Cd4*Vz;
    Fz = Cz0 + Cz1*u1 + Cz2*u2;
    My = Cm0 + Cm1*u1 + Cm2*u2;

end

function [Vxdot,Vzdot,thetadot,qdot,u1dot,u2dot] = forward(long,Fx,Fz,My,u1ref,u2ref)
    global m Iy
    Vxdot = Fx*1.9*10^-6/m - long(4)*long(2) - 9.81*sin(long(3)) ;
    Vzdot = Fz*1.9*10^-6/m + long(4)*long(1) + 9.81*cos(long(3)) ;
    thetadot = long(4);
    qdot = My*1.9*10^-6*0.0875/Iy;
    %rotor dynamics
    u1dot = 30*(u1ref-long(5));
    u2dot = 30*(u2ref-long(6));     
end

function [Vxdot,Vzdot,thetadot,qdot,u1dot,u2dot] = backward(long,Fx,Fz,My,u1ref,u2ref)
    global m Iy
    Vxdot = -Fx*1.9*10^-6/m + long(4)*long(2) + 9.81*sin(long(3)); 
    Vzdot = -Fz*1.9*10^-6/m - long(4)*long(1) - 9.81*cos(long(3));
    thetadot = -long(4);
    qdot = -My*1.9*10^-6*0.0875/Iy;
    %rotor dynamics
    u1dot = 30*(u1ref-long(5));
    u2dot = 30*(u2ref-long(6));    
end
