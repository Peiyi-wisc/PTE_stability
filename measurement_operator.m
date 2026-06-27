clc
tic

Kn_data = [1/8; 1/4; 1/2; 1; 2; 4];

Kn = Kn_data(3); 

dx = min(0.02/5, Kn/125); x_max = 0.5; x = [0:dx:x_max]; 

dmu = 0.01; mu = [-1+dmu/2:dmu:1-dmu/2]; N_mu = length(mu);

domega = 0.05; omega_max = 2; omega_min = domega; omega = [omega_min:domega:omega_max];

%M = ceil(length(mu)/50);
M = 1; mu_m = 194; %mu(3*N_mu/4-M+1)
N_omega=length(omega);
%dt=min(Kn_min*dx/2, Kn_min^2);T=(2*Kn_max*x_max/(mu(N_mu-M+1)*omega(1))); time=[0:dt:T];
dt=min(Kn*dx/2, Kn^2); T=(2*Kn*x_max/(mu(mu_m)*omega(1)))/20; time=[0:dt:T]; 

Nt = length(time);

eta1 = (tanh(10*(omega-1.5))-tanh(2*(omega-1)))/4+1/2;      % ground truth eta, eta1
eta2 = (tanh(10*(omega-1.4))-tanh(2*(omega-0.9)))/4+1/2;    % eta2

% %% 
% figure(222)
% p1 = plot(omega, eta1, '-o', 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]); hold on 
% p2 = plot(omega, eta2, '-o', 'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980]); hold off
% xlabel('$\omega$','Interpreter','latex')
% ylabel('$\eta$','Interpreter','latex')
% legend([p1 p2], {'$\eta_1$', '$\eta_2$'}, 'Interpreter','latex', 'Location', 'best')
% fontsize(gca,20,"points")

% compute temperature for 40 phi(omega)-sources 
A = zeros(Nt,M,N_omega); %% store measurement operator

for m=1:M  % choose \mu_0
    for n = 1:N_omega % choose omega(n)
        %T1=(2*Kn*x_max/(mu(3*N_mu/4-M+m)*omega(n))); t1=[0:dt:T1]; Nt1 = length(t1);
        T1 = T;
        temperature=phonon(dx, x_max, dt, T1, dmu, domega, omega_min, omega_max, eta1, mu_m, n, Kn); %% compute temperature  
        A(:,m,n) = temperature;
    end
end
toc

%%
A_Kn05_2d = reshape(A, Nt, N_omega);
figure(2205); mesh(omega,time, dt*A_Kn05_2d), view(2); colorbar; caxis([0 max_value(3)])
% save temperature for different Kn and \eta:
%A_Kn05_eta2 = A;
%save("A_Kn0125_eta2.mat","A_Kn0125_eta2")

