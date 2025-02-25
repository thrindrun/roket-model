clear;
clc;
%close all;

% Sayısal değerler %
g_val = 9.81; m_val = 25; rho_val = 1.225;
Ix_val = .5; Iy_val = .6; Iz_val = .7;
%phi_val = 0; theta_val = 0; psi_val = 0;
Sx_val = .3; Sy_val = .3; Sz_val = .05;
Cu_val = .1; Cv_val = .1; Cw_val = .1;
Cp_val = .01; Cq_val = .01; Cr_val = .01;

% Syms %
syms u v w p q r phi theta psi real % states
syms Fx Fy Fz L M N real % inputs
syms g m Ix Iy Iz rho real
syms Sx Sy Sz real % area
syms Cu Cv Cw Cp Cq Cr

% State vector %
x = [u; v; w; p; q; r; phi; theta; psi;];

% Input vector %
U = [Fx; Fy; Fz; L; M; N];

% Sistem denklemleri %
f = [
    Fx/m - 0.5*rho*Sx*u^2*Cu/m - q*w + r*v - g*cos(theta)*cos(psi);
    Fy/m - 0.5*rho*Sy*v^2*Cv/m - r*u + p*w - g*(sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi));
    Fz/m - 0.5*rho*Sz*w^2*Cw/m - p*w + q*u - g*(cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi));
    (L-q*r*(Iz-Iy)) / Ix - Cp*p;
    (M-r*p*(Ix-Iz)) / Iy - Cq*q;
    (N-p*q*(Iy-Ix)) / Iz - Cr*r;
    p + (q*sin(phi)+r*cos(phi)) * tan(theta);
    q*cos(phi) - r*sin(phi);
    (q*sin(phi)+r*cos(phi)) / cos(theta);
];

% Jacobian %
A_matrix = jacobian(f,x);
B_matrix = jacobian(f,U);

% Equilibrium %
x_eq = [0; 0; 0; 0; 0; 0; 0; 0; 0];
U_eq = [0; 0; 0; 0; 0; 0];

% Subs %
A_eq = subs(A_matrix, x, x_eq);
A_eq = subs(A_eq, U, U_eq);
A_eq = subs(A_eq, [m, g, rho, Ix, Iy, Iz, Cu, Cv, Cw, Cp, Cq, Cr, Sx, Sy, Sz], ...
    [m_val, g_val, rho_val, Ix_val, Iy_val, Iz_val, Cu_val, Cv_val, Cw_val, ...
    Cp_val, Cq_val, Cr_val, Sx_val, Sy_val, Sz_val]);
A = double(A_eq);

B_eq = subs(B_matrix, x, x_eq);
B_eq = subs(B_eq, U, U_eq);
B_eq = subs(B_eq, [m, g, rho, Ix, Iy, Iz, Cu, Cv, Cw, Cp, Cq, Cr, Sx, Sy, Sz], ...
    [m_val, g_val, rho_val, Ix_val, Iy_val, Iz_val, Cu_val, Cv_val, Cw_val, ...
    Cp_val, Cq_val, Cr_val, Sx_val, Sy_val, Sz_val]);
B = double(B_eq);

% Controllability (rank should be 9 since there are 9 states) %
controllability = ctrb(A,B);
rank_controllability = rank(controllability);

% LQR %
Q = diag([100,1,1,1,1,1,1,1,1]);
R = diag([.1,1,1,1,1,1]);
K = lqr(A,B,Q,R);

% Output %
C = eye(9);
D = zeros(9,6);

% Initial conditions %
x_init = [1; 1; 1; 0; 0; 0; 0; 0; 0];
t = 0:0.01:10;

% Closed loop matrices %
A_cl = A - B*K;
B_cl = B;
C_cl = C;
D_cl = D;

% State space %
sys_cl = ss(A_cl,B_cl,C_cl,D_cl);

% Simulation %
[Y, T, X] = initial(sys_cl, x_init, t);

% Results %
figure 

% Linear velocities %
subplot(3,1,1);
plot(T, X(:,1:3), 'LineWidth', 1.5);
legend('u (X ekseni hızı)', 'v (Y ekseni hızı)', 'w (Z ekseni hızı)');
xlabel('Zaman (s)');
ylabel('Hız (m/s)');
title('Roketin Gövde Eksenindeki Hızları');
grid on;

% Angular velocities %
subplot(3,1,2);
plot(T, X(:,4:6), 'LineWidth', 1.5);
legend('p (Roll hızı)', 'q (Pitch hızı)', 'r (Yaw hızı)');
xlabel('Zaman (s)');
ylabel('Açısal Hız (rad/s)');
title('Roketin Gövde Eksenindeki Açısal Hızları');
grid on;

% Euler angles %
subplot(3,1,3);
plot(T, X(:,7:9), 'LineWidth', 1.5);
legend('\phi (Roll açısı)', '\theta (Pitch açısı)', '\psi (Yaw açısı)');
xlabel('Zaman (s)');
ylabel('Açı (rad)');
title('Roketin Euler Açıları');
grid on;