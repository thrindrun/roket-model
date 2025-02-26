clear;
clc;
close all;



% **Sayısal değerler**
m_val = 25.0;    
Ix_val = 0.5;   
Iy_val = 0.7;   
Iz_val = 0.7;
theta_val=0;
phi_val=0;
psi_val=0;

g_val = 9.81; % Yerçekimi ivmesi
% **Sembolik değişkenleri tanımla**
syms u v w p q r phi theta psi Fx Fy Fz L M N real
syms m Ix Iy Iz g real  

% **Durum değişkenleri**
x = [u; v; w; p; q; r; phi; theta; psi];

% **Girdi vektörü**
u_vec = [Fx; Fy; Fz; L; M; N];

% **Güncellenmiş sistem denklemleri (yerçekimi ve sürtünme dahil)**
f = [
    -g*cos(theta)*cos(psi)+Fx/m-0.016*u^2-q*w+r*v;  %-0.016*(u^2 + v^2 + w^2) - q*w + r*v// (-g*cos(theta)*cos(psi))
    -g*(sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi))+(Fy/m)-0.16*v^2-r*u+ p*w;
    -g*(cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi))+(Fz/m)-0.16*w^2-p*v + q*u;
    (L - q*r*(Iz - Iy))/Ix; %- cp*p
    (0.16*q^2-M)/Iy; % - cq*q
    (0.16*r^2-N)/Iz; % - cr*r
    p + (q*sin(phi) + r*cos(phi)) * tan(theta);
    q*cos(phi) - r*sin(phi);
    (q*sin(phi) + r*cos(phi)) / cos(theta);
];


% **Jacobian matrislerini hesapla**    
A_matrix = jacobian(f, x);
B_matrix = jacobian(f, u_vec);

% **Denge noktaları**
x_equilibrium = [0; 0; 0; 0; 0; 0; 0; 0; 0];
u_equilibrium = [m_val*g_val; 0; 0; 0; 0; 0];

% **A ve B matrislerini hesapla**
A_eq = subs(A_matrix, x, x_equilibrium);
A_eq = subs(A_eq, u_vec, u_equilibrium);
A_eq = subs(A_eq, [m, Ix, Iy, Iz, g,phi,theta,psi],[m_val, Ix_val, Iy_val, Iz_val, g_val,phi_val,theta_val,psi_val]);
A = double(A_eq); % Artık tüm değişkenler sayısal, double() çağrılabilir

B_eq = subs(B_matrix, [x; u_vec], [x_equilibrium; u_equilibrium]);
B_eq = subs(B_eq, [m, Ix, Iy, Iz, g,phi,theta,psi],[m_val, Ix_val, Iy_val, Iz_val, g_val,phi_val,theta_val,psi_val]);
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
x_init = [1; 1; 1; 0; 0; 0; 0; 0; 0];
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
