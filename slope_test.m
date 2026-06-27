
N = 1;  % mu_0
M = 40; % N_omega

load("A_Kn0125_eta1.mat","A_Kn0125_eta1"); load("A_Kn0125_eta2.mat","A_Kn0125_eta2");
load("A_Kn025_eta1.mat","A_Kn025_eta1"); load("A_Kn025_eta2.mat","A_Kn025_eta2");
load("A_Kn05_eta1.mat","A_Kn05_eta1"); load("A_Kn05_eta2.mat","A_Kn05_eta2");
load("A_Kn1_eta1.mat","A_Kn1_eta1"); load("A_Kn1_eta2.mat","A_Kn1_eta2");
load("A_Kn2_eta1.mat","A_Kn2_eta1"); load("A_Kn2_eta2.mat","A_Kn2_eta2");
load("A_Kn4_eta1.mat","A_Kn4_eta1"); load("A_Kn4_eta2.mat","A_Kn4_eta2");

%%
% A1 = reshape(A_Kn1,[359,N*M]);
% A1_1 = reshape(A_Kn1_1,[359,N*M]);
% A2 = reshape(A_Kn05,[359,N*M]);
% A2_1 = reshape(A_Kn05_1,[359,N*M]);
% A3 = reshape(A_Kn033,[538,N*M]);
% A3_1 = reshape(A_Kn033_1,[538,N*M]);
% A4 = reshape(A_Kn025,[717,N*M]);
% A4_1 = reshape(A_Kn025_1,[717,N*M]);
% A5 = reshape(A_Kn016,[1076,N*M]);
% A5_1 = reshape(A_Kn016_1,[1076,N*M]);

% A1 = reshape(A_Kn05,[4041,N*M]);
% A1_1 = reshape(A_Kn05_eta1,[4041,N*M]);
% A2 = reshape(A_Kn033,[6061,N*M]);
% A2_1 = reshape(A_Kn033_eta1,[6061,N*M]);
% A3 = reshape(A_Kn025,[8081,N*M]);
% A3_1 = reshape(A_Kn025_eta1,[8081,N*M]);
% A4 = reshape(A_Kn02,[10102,N*M]);
% A4_1 = reshape(A_Kn02_eta1,[10102,N*M]);

Kn_data = [1/8; 1/4; 1/2; 1; 2; 4];

dx_total = min(0.02/5, Kn_data/125);
dt_total=min(Kn_data.*dx_total/2, Kn_data.^2);

diff = zeros(6,1);

diff(1) = dt_total(1)*max(max(abs(A_Kn0125_eta1 - A_Kn0125_eta2)));
diff(2) = dt_total(2)*max(max(abs(A_Kn025_eta1 - A_Kn025_eta2)));
diff(3) = dt_total(3)*max(max(abs(A_Kn05_eta1 - A_Kn05_eta2)));
diff(4) = dt_total(4)*max(max(abs(A_Kn1_eta1 - A_Kn1_eta2)));
diff(5) = dt_total(5)*max(max(abs(A_Kn2_eta1 - A_Kn2_eta2)));
diff(6) = dt_total(6)*max(max(abs(A_Kn4_eta1 - A_Kn4_eta2)));
%diff(7) = dt_total(7)*max(max(abs(A_Kn8_eta1 - A_Kn8_eta2)));

diff_pos(1) = max(max(abs(-A_Kn0125_eta1 + A_Kn0125_eta2)));
diff_pos(2) = max(max(abs(-A_Kn025_eta1 + A_Kn025_eta2)));
diff_pos(3) = max(max(abs(-A_Kn05_eta1 + A_Kn05_eta2)));
diff_pos(4) = max(max(abs(-A_Kn1_eta1 + A_Kn1_eta2)));
diff_pos(5) = max(max(abs(-A_Kn2_eta1 + A_Kn2_eta2)));
diff_pos(6) = max(max(abs(-A_Kn4_eta1 + A_Kn4_eta2)));
%diff_pos(7) = max(max(abs(-A_Kn8+A_Kn8_eta1)));

[max1_x,max1_y] = find((-A_Kn0125_eta1 + A_Kn0125_eta2)==diff_pos(1));
[max2_x,max2_y] = find((-A_Kn025_eta1 + A_Kn025_eta2)==diff_pos(2));
[max3_x,max3_y] = find((-A_Kn05_eta1 + A_Kn05_eta2)==diff_pos(3));
[max4_x,max4_y] = find((-A_Kn1_eta1 + A_Kn1_eta2)==diff_pos(4));
[max5_x,max5_y] = find((-A_Kn2_eta1 + A_Kn2_eta2)==diff_pos(5));
[max6_x,max6_y] = find((-A_Kn4_eta1 + A_Kn4_eta2)==diff_pos(6));
%[max7_x,max7_y] = find((-A_Kn8+A_Kn8_eta1)==diff_pos(7));

[p]=polyfit(1./(Kn_data(1:end)),log(diff(1:end)),1);


figure(1)
p2 = plot((1./Kn_data),log(diff),'-o','Linewidth',2,'MarkerSize',10); hold on;
p1 = plot((1./Kn_data),p(1)*(1./Kn_data)+p(2),'--','LineWidth',1.5); 
xlabel('$1/\varepsilon$', 'Interpreter','latex')
ylabel('$log(\|\Lambda^{\varepsilon}_{\eta_1} - \Lambda^{\varepsilon}_{\eta_2}\|_{\max})$', ...
    'Interpreter','latex')
legend(p1,'fit, $y=k\frac{1}{\varepsilon} + b$','interpreter','latex')
set(gca,'FontSize',20)
hold off;

%% plot out the temperature for two etas, for all input sources
domega = 0.05; omega_max = 2; omega_min = domega; omega = [omega_min:domega:omega_max];

x=1:40; 
y_0125=1:2140; y_025 = 1:1070; y=1:535;
T_0125 = dt_total(1)*y_0125; T_025 = dt_total(2)*y_025; T_05 = dt_total(3)*y;
T_1 = dt_total(4)*y; T_2 = dt_total(5)*y; T_4 = dt_total(6)*y;

figure(20)

max_value = exp(-0.05./Kn_data)*exp(p(2));

subplot(2,5,1), mesh(omega, T_0125, dt_total(1)*A_Kn0125_eta1(y_0125,x)), view(2); 
    title('$\varepsilon$=0.125','interpreter','latex', 'FontSize', 16);  
    xlabel('$\omega$','Interpreter','latex', 'FontSize', 14); ylabel('$t$', 'Interpreter','latex', 'FontSize', 14) 
    ylim([0 T_0125(end)]); xlim([omega_min omega_max]); colorbar; %caxis([0 max_value(1)])
subplot(2,5,2), mesh(omega, T_025, dt_total(2)*A_Kn025_eta1(y_025,x)), view(2); 
    title('$\varepsilon$=0.25','interpreter','latex', 'FontSize', 16); 
    xlabel('$\omega$','Interpreter','latex', 'FontSize', 14); ylabel('$t$', 'Interpreter','latex', 'FontSize', 14) 
    ylim([0 T_025(end)]); xlim([omega_min omega_max]); colorbar; %caxis([0 max_value(2)])
subplot(2,5,3), mesh(omega, T_05, dt_total(3)*A_Kn05_eta1(y,x)), view(2); 
    title('$\varepsilon$=0.5','interpreter','latex', 'FontSize', 16); 
    xlabel('$\omega$','Interpreter','latex', 'FontSize', 14); ylabel('$t$', 'Interpreter','latex', 'FontSize', 14) 
    ylim([0 T_05(end)]); xlim([omega_min omega_max]); colorbar; caxis([0 max_value(3)])
subplot(2,5,4), mesh(omega, T_1, dt_total(4)*A_Kn1_eta1(y,x)), view(2); 
    title('$\varepsilon$=1','interpreter','latex', 'FontSize', 16); 
    xlabel('$\omega$','Interpreter','latex', 'FontSize', 14); ylabel('$t$', 'Interpreter','latex', 'FontSize', 14) 
    ylim([0 T_1(end)]); xlim([omega_min omega_max]); colorbar; caxis([0 max_value(4)])
% subplot(2,5,5), mesh(omega, T_2, dt_total(5)*A_Kn2_eta1(y,x)), view(2); 
%     title('$\varepsilon$=2','interpreter','latex', 'FontSize', 13); 
%     xlabel('$\omega$','Interpreter','latex'); ylabel('$t$', 'Interpreter','latex') 
%     colorbar; %caxis([0 max_value(5)])
subplot(2,5,5), mesh(omega, T_4, dt_total(6)*A_Kn4_eta1(y,x)), view(2); 
    title('$\varepsilon$=4','interpreter','latex', 'FontSize', 16); 
    xlabel('$\omega$','Interpreter','latex', 'FontSize', 14); ylabel('$t$', 'Interpreter','latex', 'FontSize', 14) 
    ylim([0 T_4(end)]); xlim([omega_min omega_max]); colorbar; caxis([0 max_value(6)])

subplot(2,5,6), mesh(omega,T_0125, dt_total(1)*A_Kn0125_eta2(y_0125,x)), view(2);
xlabel('$\omega$','Interpreter','latex', 'FontSize', 14); ylabel('$t$', 'Interpreter','latex', 'FontSize', 14) 
    ylim([0 T_0125(end)]); xlim([omega_min omega_max]); colorbar; %caxis([0 max_value(1)])
subplot(2,5,7), mesh(omega,T_025, dt_total(2)*A_Kn025_eta2(y_025,x)), view(2); 
xlabel('$\omega$','Interpreter','latex', 'FontSize', 14); ylabel('$t$', 'Interpreter','latex', 'FontSize', 14) 
    ylim([0 T_025(end)]); xlim([omega_min omega_max]); colorbar; %caxis([0 max_value(2)])
subplot(2,5,8), mesh(omega,T_05, dt_total(3)*A_Kn05_eta2(y,x)), view(2);
xlabel('$\omega$','Interpreter','latex', 'FontSize', 14); ylabel('$t$', 'Interpreter','latex', 'FontSize', 14) 
    ylim([0 T_05(end)]); xlim([omega_min omega_max]); colorbar; caxis([0 max_value(3)])
subplot(2,5,9), mesh(omega,T_1, dt_total(4)*A_Kn1_eta2(y,x)), view(2); 
xlabel('$\omega$','Interpreter','latex', 'FontSize', 14); ylabel('$t$', 'Interpreter','latex', 'FontSize', 14) 
    ylim([0 T_1(end)]); xlim([omega_min omega_max]); colorbar; caxis([0 max_value(4)])
% subplot(2,5,10), mesh(omega,T_2, dt_total(5)*A_Kn2_eta2(y,x)), view(2); 
% xlabel('$\omega$','Interpreter','latex'); ylabel('$t$', 'Interpreter','latex') 
%     colorbar; %caxis([0 max_value(5)])
subplot(2,5,10), mesh(omega,T_4,dt_total(6)*A_Kn4_eta2(y,x)), view(2); 
xlabel('$\omega$','Interpreter','latex', 'FontSize', 14); ylabel('$t$', 'Interpreter','latex', 'FontSize', 14) 
    ylim([0 T_4(end)]); xlim([omega_min omega_max]); colorbar; caxis([0 max_value(6)])

%%
figure(21)
max_value = exp(-0.05./Kn_data)*exp(p(2));
max_diff10 = max(max(dt_total(6)*(A_Kn4_eta1(y,x) - A_Kn4_eta2(y,x))));
min_diff10 = min(min(dt_total(6)*(A_Kn4_eta1(y,x) - A_Kn4_eta2(y,x))));
max_diff7 = max(max(dt_total(2)*(A_Kn025_eta1(y_025,x) - A_Kn025_eta2(y_025,x))));
min_diff7 = min(min(dt_total(2)*(A_Kn025_eta1(y_025,x) - A_Kn025_eta2(y_025,x))));
max_diff1 = max(max(dt_total(1)*(A_Kn0125_eta1(y_0125,x) - A_Kn0125_eta2(y_0125,x))));
min_diff1 = min(min(dt_total(1)*(A_Kn0125_eta1(y_0125,x) - A_Kn0125_eta2(y_0125,x))));
max_diff3 = max(max(dt_total(3)*(A_Kn05_eta1(y,x) - A_Kn05_eta2(y,x))));
min_diff3 = min(min(dt_total(3)*(A_Kn05_eta1(y,x) - A_Kn05_eta2(y,x))));

subplot(2,5,1), mesh(omega, T_0125, dt_total(1)*A_Kn0125_eta1(y_0125,x)), view(2); colorbar; 
    xlabel('$\omega$','Interpreter','latex'); ylabel('$t$', 'Interpreter','latex') 
    ylim([0 T_0125(end)]); xlim([omega_min omega_max]); set(gca,'fontsize',16);
    title('$\varepsilon$=0.125','interpreter','latex', 'FontSize', 18);  
subplot(2,5,2), mesh(omega, T_025, dt_total(2)*A_Kn025_eta1(y_025,x)), view(2); colorbar;
    xlabel('$\omega$','Interpreter','latex'); ylabel('$t$', 'Interpreter','latex') 
    ylim([0 T_025(end)]); xlim([omega_min omega_max]); set(gca,'fontsize',16);
    title('$\varepsilon$=0.25','interpreter','latex', 'FontSize', 18); 
subplot(2,5,3), mesh(omega, T_05, dt_total(3)*A_Kn05_eta1(y,x)), view(2); colorbar; caxis([0 max_value(3)])
    xlabel('$\omega$','Interpreter','latex'); ylabel('$t$', 'Interpreter','latex') 
    ylim([0 T_05(end)]); xlim([omega_min omega_max]);  set(gca,'fontsize',16);
    title('$\varepsilon$=0.5','interpreter','latex', 'FontSize', 18);
subplot(2,5,4), mesh(omega, T_1, dt_total(4)*A_Kn1_eta1(y,x)), view(2); colorbar; caxis([0 max_value(4)])
    xlabel('$\omega$','Interpreter','latex'); ylabel('$t$', 'Interpreter','latex') 
    ylim([0 T_1(end)]); xlim([omega_min omega_max]);  set(gca,'fontsize',16)
    title('$\varepsilon$=1','interpreter','latex', 'FontSize', 18); 
subplot(2,5,5), mesh(omega, T_4, dt_total(6)*A_Kn4_eta1(y,x)), view(2); colorbar; caxis([0 max_value(6)])
    xlabel('$\omega$','Interpreter','latex'); ylabel('$t$', 'Interpreter','latex') 
    ylim([0 T_4(end)]); xlim([omega_min omega_max]); set(gca,'fontsize',16)
    title('$\varepsilon$=4','interpreter','latex', 'FontSize', 18);

subplot(2,5,6), mesh(omega,T_0125, dt_total(1)*(A_Kn0125_eta1(y_0125,x) - A_Kn0125_eta2(y_0125,x))), view(2);
xlabel('$\omega$','Interpreter','latex'); ylabel('$t$', 'Interpreter','latex') 
    ylim([0 T_0125(end)]); xlim([omega_min omega_max]); colorbar; caxis([min_diff1 max_diff1]);
    set(gca,'fontsize',16)
subplot(2,5,7), mesh(omega,T_025, dt_total(2)*(A_Kn025_eta1(y_025,x) - A_Kn025_eta2(y_025,x))), view(2); 
xlabel('$\omega$','Interpreter','latex'); ylabel('$t$', 'Interpreter','latex') 
    ylim([0 T_025(end)]); xlim([omega_min omega_max]); colorbar; caxis([min_diff7 max_diff7])
    set(gca,'fontsize',16)
subplot(2,5,8), mesh(omega,T_05, dt_total(3)*(A_Kn05_eta1(y,x) - A_Kn05_eta2(y,x))), view(2);
xlabel('$\omega$','Interpreter','latex'); ylabel('$t$', 'Interpreter','latex') 
    ylim([0 T_05(end)]); xlim([omega_min omega_max]); colorbar; caxis([min_diff10 max_diff10])
    set(gca,'fontsize',16)
subplot(2,5,9), mesh(omega,T_1, dt_total(4)*(A_Kn1_eta1(y,x) - A_Kn1_eta2(y,x))), view(2); 
xlabel('$\omega$','Interpreter','latex'); ylabel('$t$', 'Interpreter','latex') 
    ylim([0 T_1(end)]); xlim([omega_min omega_max]); colorbar; caxis([min_diff10 max_diff10])
    set(gca,'fontsize',16)
subplot(2,5,10), mesh(omega,T_4,dt_total(6)*(A_Kn4_eta1(y,x) - A_Kn4_eta2(y,x))), view(2); 
xlabel('$\omega$','Interpreter','latex'); ylabel('$t$', 'Interpreter','latex') 
    ylim([0 T_4(end)]); xlim([omega_min omega_max]); colorbar; caxis([min_diff10 max_diff10]); 
    set(gca,'fontsize',16);

%% separate plots
figure(211)
mesh(omega, T_0125, dt_total(1)*A_Kn0125_eta1(y_0125,x)), view(2);  
ylim([0 T_0125(end)]); xlim([omega_min omega_max]); colorbar; 
set(gca,'fontsize',18);
% title('$\varepsilon$=0.125','interpreter','latex', 'FontSize', 20);  
xlabel('$\omega$','Interpreter','latex'); ylabel('$t$', 'Interpreter','latex')
pbaspect([0.5 1 1])

figure(212)
mesh(omega, T_025, dt_total(2)*A_Kn025_eta1(y_025,x)), view(2); 
xlabel('$\omega$','Interpreter','latex'); ylabel('$t$', 'Interpreter','latex') 
ylim([0 T_025(end)]); xlim([omega_min omega_max]); colorbar;
set(gca,'fontsize',18);
pbaspect([0.5 1 1])

figure(213)
mesh(omega, T_05, dt_total(3)*A_Kn05_eta1(y,x)), view(2); 
xlabel('$\omega$','Interpreter','latex'); ylabel('$t$', 'Interpreter','latex') 
ylim([0 T_05(end)]); xlim([omega_min omega_max]); colorbar;
set(gca,'fontsize',18);
pbaspect([0.5 1 1])

figure(214)
mesh(omega, T_1, dt_total(4)*A_Kn1_eta1(y,x)), view(2);  
xlabel('$\omega$','Interpreter','latex'); ylabel('$t$', 'Interpreter','latex') 
ylim([0 T_1(end)]); xlim([omega_min omega_max]); colorbar; caxis([0 max_value(4)])
set(gca,'fontsize',18);
pbaspect([0.5 1 1])

figure(215)
mesh(omega, T_4, dt_total(6)*A_Kn4_eta1(y,x)), view(2); 
xlabel('$\omega$','Interpreter','latex'); ylabel('$t$', 'Interpreter','latex') 
ylim([0 T_4(end)]); xlim([omega_min omega_max]); colorbar; caxis([0 max_value(6)])
set(gca,'fontsize',18);
pbaspect([0.5 1 1])

%%
figure(221)
mesh(omega,T_0125, dt_total(1)*(A_Kn0125_eta1(y_0125,x) - A_Kn0125_eta2(y_0125,x))), view(2); 
ylim([0 T_0125(end)]); xlim([omega_min omega_max]); colorbar; caxis([min_diff1 max_diff1]); 
set(gca,'fontsize',18); 
xlabel('$\omega$','Interpreter','latex'); ylabel('$t$', 'Interpreter','latex')
pbaspect([0.5 1 1])

figure(222)
mesh(omega,T_025, dt_total(2)*(A_Kn025_eta1(y_025,x) - A_Kn025_eta2(y_025,x))), view(2); 
xlabel('$\omega$','Interpreter','latex'); ylabel('$t$', 'Interpreter','latex') 
ylim([0 T_025(end)]); xlim([omega_min omega_max]); colorbar; caxis([min_diff7 max_diff7])
set(gca,'fontsize',18);
pbaspect([0.5 1 1])

figure(223)
mesh(omega,T_05, dt_total(3)*(A_Kn05_eta1(y,x) - A_Kn05_eta2(y,x))), view(2);
xlabel('$\omega$','Interpreter','latex'); ylabel('$t$', 'Interpreter','latex') 
ylim([0 T_05(end)]); xlim([omega_min omega_max]); colorbar; caxis([min_diff3 max_diff3])
set(gca,'fontsize',18);
pbaspect([0.5 1 1])

figure(224)
mesh(omega,T_1, dt_total(4)*(A_Kn1_eta1(y,x) - A_Kn1_eta2(y,x))), view(2); 
xlabel('$\omega$','Interpreter','latex'); ylabel('$t$', 'Interpreter','latex') 
ylim([0 T_1(end)]); xlim([omega_min omega_max]); colorbar; caxis([min_diff10 max_diff10])
set(gca,'fontsize',18);
pbaspect([0.5 1 1])

figure(225)
mesh(omega,T_4,dt_total(6)*(A_Kn4_eta1(y,x) - A_Kn4_eta2(y,x))), view(2); 
xlabel('$\omega$','Interpreter','latex'); ylabel('$t$', 'Interpreter','latex') 
ylim([0 T_4(end)]); xlim([omega_min omega_max]); c=colorbar; caxis([min_diff10 max_diff10]); 
set(gca,'fontsize',18);
pbaspect([0.5 1 1])

%% plot out the difference between measurements generated by \eta_1 and \eta_2
AA1 = reshape(A_Kn0125_eta1 - A_Kn0125_eta2, 2140, 40);
AA2 = reshape(A_Kn025_eta1 - A_Kn025_eta2, 1070, 40);
AA3 = reshape(A_Kn05_eta1 - A_Kn05_eta2, 535, 40);
AA4 = reshape(A_Kn1_eta1 - A_Kn1_eta2, 535, 40);
AA5 = reshape(A_Kn2_eta1 - A_Kn2_eta2, 535, 40);
AA6 = reshape(A_Kn4_eta1 - A_Kn4_eta2, 535, 40);

figure(30)
% max_value = exp(-1.7850./Kn_data).*5*1e-3;
subplot(1,6,1), mesh(dt_total(1)*AA1), view(2); colorbar; caxis([-diff(1) diff(1)])
subplot(1,6,2), mesh(dt_total(2)*(AA2)), view(2);  colorbar; caxis([-diff(2) diff(2)])
subplot(1,6,3), mesh(dt_total(3)*(AA3)), view(2); colorbar; caxis([-diff(3) diff(3)])
subplot(1,6,4), mesh(dt_total(4)*(AA4)), view(2); colorbar; caxis([-diff(4) diff(4)])
subplot(1,6,5), mesh(dt_total(5)*(AA5)), view(2); colorbar; caxis([-diff(5) diff(5)])
subplot(1,6,6), mesh(dt_total(6)*(AA6)), view(2); colorbar; caxis([-diff(6) diff(6)])

%%
diff_eta12_omega1 = A_Kn0125_eta1(:,1,end-11)-A_Kn0125_eta2(:,1,end-11);
tt_omega1 = dt_total(1)*(1:2140);

figure(456)
yyaxis left
plot(tt_omega1, dt_total(1)*A_Kn0125_eta1(:,1,end-11),'LineWidth',2.5); 
% hold on; plot(tt_omega1, A_Kn0125_eta2(:,1,end),'--','LineWidth',2); hold off
ylabel('$\Delta T[\eta_1]$','Interpreter','latex')

yyaxis right
plot(tt_omega1, dt_total(1)*diff_eta12_omega1,'LineWidth',2);
ylabel('$\Delta T[\eta_1]- \Delta T[\eta_2]$', 'Interpreter','latex')
xlabel('$t$', 'Interpreter','latex')
title('$\varepsilon=0.125, \omega_0=1.45$', 'Interpreter','latex')
xlim([0 0.134]); %ylim([0 0.00096])
set(gca,'FontSize',20)

%%
diff_eta12_omega1 = A_Kn025_eta1(:,1,end-11)-A_Kn025_eta2(:,1,end-11);
tt_omega1 = dt_total(2)*(1:1070);

figure(457)
yyaxis left
plot(tt_omega1, dt_total(2)*A_Kn025_eta1(:,1,end-11),'LineWidth',2.5); 
% hold on; plot(tt_omega1, A_Kn025_eta2(:,1,end),'--','LineWidth',2); hold off
ylabel('$\Delta T[\eta_1]$','Interpreter','latex')

yyaxis right
plot(tt_omega1, dt_total(2)*diff_eta12_omega1,'LineWidth',2);
ylabel('$\Delta T[\eta_1]-\Delta T[\eta_2]$', 'Interpreter','latex')
xlabel('$t$', 'Interpreter','latex')
title('$\varepsilon=0.25, \omega_0 = 1.45$', 'Interpreter','latex')
xlim([0 0.267])
set(gca,'FontSize',20)

%%
diff_eta12_omega1 = A_Kn4_eta1(:,1,end-11)-A_Kn4_eta2(:,1,end-11);
tt_omega1 = dt_total(6)*(1:535);

figure(458)
yyaxis left
plot(tt_omega1, dt_total(6)*A_Kn4_eta1(:,1,end-11),'LineWidth',2.5); 
% hold on; plot(tt_omega1, A_Kn4_eta2(:,1,end),'--','LineWidth',2); hold off
ylabel('$\Delta T[\eta_1]$','Interpreter','latex')

yyaxis right
plot(tt_omega1, dt_total(6)*diff_eta12_omega1,'LineWidth',2);
ylabel('$\Delta T[\eta_1]-\Delta T[\eta_2]$', 'Interpreter','latex')
xlabel('$t$', 'Interpreter','latex')
title('$\varepsilon=4, \omega_0 = 1.45$', 'Interpreter','latex')
xlim([0 4.28])
set(gca,'FontSize',20)

%%
diff_eta12_omega1 = A_Kn1_eta1(:,1,end-11)-A_Kn1_eta2(:,1,end-11);
tt_omega1 = dt_total(4)*(1:535);

figure(459)
yyaxis left
plot(tt_omega1, dt_total(4)*A_Kn1_eta1(:,1,end-11),'LineWidth',2.5); 
% hold on; plot(tt_omega1, A_Kn4_eta2(:,1,end),'--','LineWidth',2); hold off
ylabel('$\Delta T[\eta_1]$','Interpreter','latex')

yyaxis right
plot(tt_omega1, dt_total(4)*diff_eta12_omega1,'LineWidth',2);
ylabel('$\Delta T[\eta_1]-\Delta T[\eta_2]$', 'Interpreter','latex')
xlabel('$t$', 'Interpreter','latex')
title('$\varepsilon=1, \omega_0 = 1.45$', 'Interpreter','latex')
xlim([0 1.07])
set(gca,'FontSize',20)


%%
diff_eta12_omega1 = A_Kn2_eta1(:,1,end-11)-A_Kn2_eta2(:,1,end-11);
tt_omega1 = dt_total(5)*(1:535);

figure(459)
yyaxis left
plot(tt_omega1, dt_total(5)*A_Kn2_eta1(:,1,end-11),'LineWidth',2.5); 
% hold on; plot(tt_omega1, A_Kn4_eta2(:,1,end),'--','LineWidth',2); hold off
ylabel('$\Delta T[\eta_1]$','Interpreter','latex')

yyaxis right
plot(tt_omega1, dt_total(5)*diff_eta12_omega1,'LineWidth',2);
ylabel('$\Delta T[\eta_1]-\Delta T[\eta_2]$', 'Interpreter','latex')
xlabel('$t$', 'Interpreter','latex')
title('$\varepsilon=2, \omega_0 = 1.45$', 'Interpreter','latex')
xlim([0 2.14])
set(gca,'FontSize',20)

%%
diff_eta12_omega1 = A_Kn05_eta1(:,1,end-11)-A_Kn05_eta2(:,1,end-11);
tt_omega1 = dt_total(3)*(1:535);

figure(460)
yyaxis left
plot(tt_omega1, dt_total(3)*A_Kn05_eta1(:,1,end-11),'LineWidth',2.5); 
% hold on; plot(tt_omega1, A_Kn4_eta2(:,1,end),'--','LineWidth',2); hold off
ylabel('$\Delta T[\eta_1]$','Interpreter','latex')

yyaxis right
plot(tt_omega1, dt_total(3)*diff_eta12_omega1,'LineWidth',2);
ylabel('$\Delta T[\eta_1]-\Delta T[\eta_2]$', 'Interpreter','latex')
xlabel('$t$', 'Interpreter','latex')
title('$\varepsilon=0.5, \omega_0 = 1.45$', 'Interpreter','latex')
xlim([0 0.535])
set(gca,'FontSize',20)
