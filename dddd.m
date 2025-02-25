clear;
clc;
close all;

% **Sayısal değerler**
m_val = 25;    
m_dot_val = -0.05; 
Ix_val = 0.5;   
Iy_val = 0.6;   
Iz_val = 0.7;  
cu_val = 0.1;   
cv_val = 0.1;   
cw_val = 0.1;   
cp_val = 0.01;   
cq_val = 0.01;   
cr_val = 0.01;
g_val = 9.81; % Yerçekimi ivmesi
rho_val = 1.225;
Sx_val = 0.3;
Sy_val = 0.3;
Sz_val = 0.05;
% **Sembolik değişkenleri tanımla**
syms u v w p q r phi theta psi Fx Fy Fz L M N real
syms m Ix Iy Iz g m_dot real  % Yerçekimi ve kütle değişimi tanımlandı
syms cu cv cw cp cq cr real   % Sürtünme katsayıları
syms rho Sx Sy Sz % density area

% **Durum değişkenleri**
x = [u; v; w; p; q; r; phi; theta; psi];

% **Girdi vektörü**
u_vec = [Fx; Fy; Fz; L; M; N];

% **Güncellenmiş sistem denklemleri (yerçekimi ve sürtünme dahil)**
f = [
    Fx/m - q*w + r*v - 0.5*rho*Sx*u^2*cu/m + g*sin(theta);  
    Fy/m - r*u + p*w - 0.5*rho*Sy*v^2*cv/m - g*sin(phi)*cos(theta);
    Fz/m - p*v + q*u + 0.5*rho*Sz*w^2*cw/m - g*cos(phi)*cos(theta);
    (L - q*r*(Iz - Iy))/Ix - cp*p;
    (M - r*p*(Ix - Iz))/Iy - cq*q;
    (N - p*q*(Iy - Ix))/Iz - cr*r;
    p + (q*sin(phi) + r*cos(phi)) * tan(theta);
    q*cos(phi) - r*sin(phi);
    (q*sin(phi) + r*cos(phi)) / cos(theta);
];

% **Jacobian matrislerini hesapla**
A_matrix = jacobian(f, x);
B_matrix = jacobian(f, u_vec);

% **Denge noktaları**
x_equilibrium = [0; 0; 0; 0; 0; 0; 0; 0; 0];
u_equilibrium = [0; 0; 250; 0; 0; 0];

% **A ve B matrislerini hesapla**
A_eq = subs(A_matrix, x, x_equilibrium);
A_eq = subs(A_eq, u_vec, u_equilibrium);
A_eq = subs(A_eq, [m, m_dot, Ix, Iy, Iz, g, cu, cv, cw, cp, cq, cr, rho, Sx, Sy, Sz], ...
                [m_val, m_dot_val, Ix_val, Iy_val, Iz_val, g_val, cu_val, cv_val, cw_val, cp_val, cq_val, cr_val, rho_val, Sx_val, Sy_val, Sz_val])
A = double(A_eq); % Artık tüm değişkenler sayısal, double() çağrılabilir

B_eq = subs(B_matrix, [x; u_vec], [x_equilibrium; u_equilibrium]);
B_eq = subs(B_eq, [m, m_dot, Ix, Iy, Iz, g, cu, cv, cw, cp, cq, cr], ...
                [m_val, m_dot_val, Ix_val, Iy_val, Iz_val, g_val, cu_val, cv_val, cw_val, cp_val, cq_val, cr_val]);
B = double(B_eq); % Aynı şekilde B için de double() çağrılabilir

% **Kontrol edilebilirlik matrisi**
controllability = ctrb(A, B);
rank_controllability = rank(controllability);

% **LQR ile Optimal Kontrol Kazancı K Matrisi**
Q = diag([1,1,1,1,1,1,1,1,1]);
R = diag([1,1,1,1,1,1]);
K = lqr(A,B,Q,R);

% **Çıkış ve Direkt Geçiş Matrisleri**
C = eye(9);
D = zeros(9,6);

% **Başlangıç koşulları**
x_init = [1; 1; 1; 1; 1; 1; 1; 1; 1];
t = 0:0.01:10;

% **Kapalı çevrim sistem matrisi**
A_closed = A - B*K;
B_closed = B;
C_closed = C;
D_closed = D;

% **Simülasyon için state-space modeli oluştur**
sys_cl = ss(A_closed, B_closed, C_closed, D_closed);

% **Simülasyonu çalıştır**
[Y, T, X] = initial(sys_cl, x_init, t);

% **Sonuçları Çizdir**
figure;

% **Hız Bileşenleri (u, v, w)**
subplot(3,1,1);
plot(T, X(:,1:3), 'LineWidth', 1.5);
legend('u (X ekseni hızı)', 'v (Y ekseni hızı)', 'w (Z ekseni hızı)');
xlabel('Zaman (s)');
ylabel('Hız (m/s)');
title('Roketin Gövde Eksenindeki Hızları');
grid on;

% **Açısal Hızlar (p, q, r)**
subplot(3,1,2);
plot(T, X(:,4:6), 'LineWidth', 1.5);
legend('p (Roll hızı)', 'q (Pitch hızı)', 'r (Yaw hızı)');
xlabel('Zaman (s)');
ylabel('Açısal Hız (rad/s)');
title('Roketin Gövde Eksenindeki Açısal Hızları');
grid on;

% **Euler Açıları (phi, theta, psi)**
subplot(3,1,3);
plot(T, X(:,7:9), 'LineWidth', 1.5);
legend('\phi (Roll açısı)', '\theta (Pitch açısı)', '\psi (Yaw açısı)');
xlabel('Zaman (s)');
ylabel('Açı (rad)');
title('Roketin Euler Açıları');
grid on;
