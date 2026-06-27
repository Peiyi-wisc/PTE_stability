% coupled system

Kn = 0.5;
Ls = 4;

% discretizations
dx = min(0.02/5, Kn/125); x_max = 1; x = 0:dx:x_max;  Nx = length(x); 
dy = min(0.02/5, Kn/125); y_max = Ls; y = 0:dy:y_max; Ny = length(y);
dmu = 0.01; mu = -1+dmu/2:dmu:1-dmu/2; N_mu = length(mu);
domega = 0.05; omega_max = 2; omega_min = domega; omega = omega_min:domega:omega_max;

% choose mu_0, omega_0:
mu_m = 194; ome_n = 20;

dt=min(Kn*dx/2, Kn^2); T=(2*Kn*x_max/(mu(mu_m)*omega(ome_n)))*4; 
time=0:dt:T; 
Nt = length(time);

eta1 = (tanh(10*(omega-1.5))-tanh(2*(omega-1)))/4 + 0.9;%1/2;   
zeta1= (1 - eta1);

eta_sub  = 2*eta1 - 1;
zeta_sub = 2*(1-eta1);

%% plot the (eta, zeta) pair of transducer & substrate
figure(1);
plot(omega,eta1, '-o', 'LineWidth', 2); hold on; 
plot(omega,zeta1,'--','LineWidth', 2.5); hold off
xlabel('$\omega$','Interpreter','latex');
ylabel('$\eta_t, \zeta_t$','Interpreter','latex');
legend({'$\eta_t$', '$\zeta_t$'}, 'Interpreter','latex', 'FontSize', 18)
set(gca,'Fontsize',18)
title('Transducer property', 'interpreter', 'latex', 'FontSize', 18)
%
figure(2);
p1 = plot(omega,eta_sub, '-o', 'LineWidth',2); hold on; 
p2 = plot(omega,zeta_sub,'--','LineWidth',2.5); hold off
xlabel('$\omega$','Interpreter','latex')
ylabel('$\eta_s, \zeta_s$','Interpreter','latex');
legend({'$\eta_s$', '$\zeta_s$'}, 'Interpreter','latex', 'FontSize', 18)
set(gca,'Fontsize',18)
title('Substrate property', 'interpreter', 'latex', 'FontSize',18)

%% run the coupled system (f,g)

[rho_f, rho_g] = coupled_phonon_fg(x, y, dt, T, dmu, domega, omega_min, omega_max, eta1, zeta1, eta_sub, zeta_sub, Kn, mu_m, ome_n);
% rho_f: Nx * (Nt+1); rho_g: Ny * (Nt+1)

disp('Simulation completed!');
%
figure(34)
imagesc(x, time(2:end), rho_f(:,2:end)')
set(gca, 'YDir', 'normal') 
colorbar
clim([0 30])
ylabel('$t$','Interpreter','latex')
xlabel('$x$ (Transducer)', 'Interpreter', 'latex')
set(gca,'Fontsize',18)
title('Phonon Density in the transducer ($\rho_f$)', 'Interpreter','latex','Fontsize',18)

figure(4)
imagesc(y, time(2:end), rho_g(:,2:end)')
set(gca, 'YDir', 'normal') 
colorbar
clim([0 30])
ylabel('$t$','Interpreter','latex')
xlabel('$x$ (Substrate)', 'Interpreter','latex')
set(gca,'Fontsize',18)
title('Phonon Density in the substrate ($\rho_g$)', 'Interpreter','latex','Fontsize',18)

%
rho_fg = zeros(Nx + Ny, Nt+1);
rho_fg(1:Nx,:) = rho_f; rho_fg(Nx+1:end,:) = rho_g;
xy_max = x_max + y_max;
xy = 0:dx:xy_max;

figure(5)
imagesc(xy, time(2:end), rho_fg(:,2:end)')
set(gca, 'YDir', 'normal') 
colorbar
clim([0 30])
ylabel('$t$','Interpreter','latex')
xlabel('$x$ (Transducer+Substrate)', 'Interpreter','latex')
set(gca,'Fontsize',18)
title('Phonon Density', 'Interpreter','latex','Fontsize',18)


%%

rho_f_t0 = rho_f(:,end);
figure(3)
plot(x, rho_f_t0,LineWidth=1);
xlabel('$x$','Interpreter','latex');
ylabel('$\rho$', 'Interpreter','latex')

rho_g_t0 = rho_g(:,end);
figure(444)
plot(y, rho_g_t0,LineWidth=1);
xlabel('$x$','Interpreter','latex');
ylabel('$\rho$', 'Interpreter','latex')