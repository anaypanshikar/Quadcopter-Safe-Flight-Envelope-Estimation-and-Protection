close all
clear all 
clc

%%
global V theta0 Omega_f0 Omega_b0 Vx0 Vz0 T N h N_traj p_c m Iy k_r

load('qtheta_env.mat')

%% Trim conditions
V = 8;
theta0 = -13.1*pi/180;
Omega_f0 = 709.5;
Omega_b0 = 860.4;
Vx0 = V*cos(theta0);        %7.792
Vz0 = V*sin(theta0);       %1.813

%% Initial plot
figure()
plot(xv,yv);
hold on
scatter(theta0*180/pi,0,100,"red","x","LineWidth",2)
hold off
grid on
xlabel("theta (deg)")
ylabel("q (deg/s)")
legend("SFE","Trim condition")

%% Monte carlo method variables
T = 0.15;
N = 100;
h = T/N;         %simulation time step
N_traj = 1;
p_c = 0.1;
m = 0.375;
Iy = 1.3*10^-3;
k_r = 1;

%% simulation

for i = 1:N_traj
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
    Fxf = m*9.81*sin(theta0);
    Fzf = -m*9.81*cos(theta0);
    Myf = 0;
    path = [longf];
    for k = 1:N-1
        % forward analysis
        [k_11,k_12,k_13,k_14,k_15,k_16] = forward(longf,Fxf,Fzf,Myf,u1ref,u2ref);
        [k_21,k_22,k_23,k_24,k_25,k_26] = forward(longf+0.5*h.*[k_11,k_12,k_13,k_14,k_15,k_16],Fxf,Fzf,Myf,u1ref,u2ref);
        [k_31,k_32,k_33,k_34,k_35,k_36] = forward((longf+0.5*h.*[k_21,k_22,k_23,k_24,k_25,k_26]),Fxf,Fzf,Myf,u1ref,u2ref);
        [k_41,k_42,k_43,k_44,k_45,k_46] = forward((longf+[k_31,k_32,k_33,k_34,k_35,k_36].*h),Fxf,Fzf,Myf,u1ref,u2ref);
        longf = longf + (1/6).*([k_11,k_12,k_13,k_14,k_15,k_16]+2.*[k_21,k_22,k_23,k_24,k_25,k_26]+2.*[k_31,k_32,k_33,k_34,k_35,k_36]+[k_41,k_42,k_43,k_44,k_45,k_46]).*h;  % main equation
        
        %check SFE
        %append to path
        xq = longf(3)*180/pi;
        yq = longf(4)*180/pi;
        [in,on] = inpolygon(xq,yq,xv,yv);
        if in==1
            path = [path; longf];
            %next step
            [Cd1,Cd2,Cd3,Cd4,Cz0,Cz1,Cz2,Cm0,Cm1,Cm2] = param(longf(1),longf(2));
            [Fxf,Fzf,Myf] = FnM(Cd1,Cd2,Cd3,Cd4,Cz0,Cz1,Cz2,Cm0,Cm1,Cm2,longf(5),longf(6),longf(1),longf(2));
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
        else
            break
        end
           
    end
    
end    

%% Final plots
figure()
plot(xv,yv);
hold on
scatter(theta0*180/pi,0,100,"red","x","LineWidth",2)
scatter(path(:,3)*180/pi,path(:,4)*180/pi,40,"red",".")
hold off
grid on
xlabel("theta (deg)")
ylabel("q (deg/s)")
legend("SFE","Trim condition","path")

%% Functions

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