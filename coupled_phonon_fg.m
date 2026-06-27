function [rho_f, rho_g] = coupled_phonon_fg(x, y, dt, T, dmu, domega, omega_min, omega_max, eta, zeta, eta_sub, zeta_sub, Kn, mu_m, ome_n)
%% Compute the coupled PTE (f,g) where Ls >> Lt
%   Inputs: x, y (coordinate vectors)
%   Outputs: rho_f(t,x) - density of f, matrix of (Nx,Nt); 
%            rho_g(t,x) - density of g, matrix of (Ny,Nt)

%% Sample grids
Nx = length(x); Ny = length(y);
dx = x(2) - x(1); 
dy = y(2) - y(1);

mu = -1+dmu/2:dmu:1-dmu/2;
omega = omega_min:domega:omega_max;

[~, mumu_f, omegaomega_f] = meshgrid(x, mu, omega);
[~, mumu_g, omegaomega_g] = meshgrid(y, mu, omega);

pos_mu_f = repmat(max(mu',0), [1, Nx-1, length(omega)]);
neg_mu_f = repmat(min(mu',0), [1, Nx-1, length(omega)]);

pos_mu_g = repmat(max(mu',0), [1, Ny-1, length(omega)]);
neg_mu_g = repmat(min(mu',0), [1, Ny-1, length(omega)]);

tt = 0:dt:T; 
Nt = length(tt); 
  
normalization = sum((10*omega).^3.*exp(10*omega)./(exp(10*omega)-1).^2)*domega*dmu*sum(length(mu));

%% Initialization and boundary condition
fn   = zeros(length(mu), Nx, length(omega));
fn_1 = fn;
gn   = zeros(length(mu), Ny, length(omega));
gn_1 = gn;
col_f = fn; col_g = gn;
rho_f = zeros(Nx, Nt); rho_g = zeros(Ny, Nt);

bdy_func = @(t_s, mu_s, omega_s) (1/dmu)*(1/(4*dt))*(1/domega).*eq(mu_s, mu(mu_m)).*eq(omega_s, omega(ome_n)).*(eq(t_s, 0) + eq(t_s, dt) + eq(t_s, 2*dt) + eq(t_s, 3*dt)) + mu_s-mu_s + t_s-t_s;

% Record the initial state profiles at index j = 1
time = 0; j = 1;
bdy = bdy_func(time, max(mumu_f, 0), omegaomega_f); 
fn((mu>0), 1, :) = bdy((mu>0), 1, :);

rho_f(:, j) = sum(sum(fn, 3), 1)' * dmu * domega;
rho_g(:, j) = sum(sum(gn, 3), 1)' * dmu * domega;

%% Evolution 
while time < T
    
    % flux - f
    mu_px_f = (1/Kn)*((omegaomega_f(:,2:end,:).*fn(:,2:end,:))-(omegaomega_f(:,1:end-1,:).*fn(:,1:end-1,:)))/dx; 
    pos_flux_f = pos_mu_f.*mu_px_f;
    neg_flux_f = neg_mu_f.*mu_px_f;
    
    % flux - g
    mu_px_g = (1/Kn)*((omegaomega_g(:,2:end,:).*gn(:,2:end,:))-(omegaomega_g(:,1:end-1,:).*gn(:,1:end-1,:)))/(2*dy); 
    pos_flux_g = pos_mu_g.*mu_px_g;
    neg_flux_g = neg_mu_g.*mu_px_g;

    % BC at x=0 - f
    bdy = bdy_func(time, max(mumu_f, 0), omegaomega_f); 
    fn_1((mu>0), 1, :) = bdy((mu>0), 1, :);
    % advection for mu>0 - f
    fn_1(:,2:end,:)   = fn(:,2:end,:) - dt * pos_flux_f;
    
    % BC at x=Lt - f
    for km = 1:length(omega)
        outgoing = fn((mu>0), end, km);
        fn_1((mu<0), end, km) = eta(km)*flipud(outgoing) + zeta(km)*gn((mu<0), 1, km);
    end
    % advection for mu<0 - f
    fn_1(:,1:end-1,:) = fn_1(:,1:end-1,:) - dt * neg_flux_f; 
  
    % Collision - f
    for km = 2:Nx-1
        col_temp = omegaomega_f(:,km,:).*fn(:,km,:); col_temp = reshape(col_temp, length(mu), length(omega));
        values = sum(sum(col_temp))*dmu*domega/normalization;
        col_temp = values*ones(length(mu),1,length(omega)).*((10*omegaomega_f(:,km,:)).^3).*exp(10*omegaomega_f(:,km,:))./((exp(10*omegaomega_f(:,km,:))-1).^2);
        col_f(:,km,:) = (1/Kn^2)*( -omegaomega_f(:,km,:).*fn(:,km,:) + reshape(col_temp, [length(mu), 1, length(omega)]));
    end 
    fn_1 = fn_1 + dt * col_f;

    % BC at x=Lt - g
    for km = 1:length(omega)
        incoming = gn((mu<0), 1, km);
        gn_1((mu>0), 1, km) = eta_sub(km)*flipud(incoming) + zeta_sub(km)*fn((mu>0), end, km);
    end
    % advection for mu>0 - g
    gn_1(:,2:end,:)  = gn(:,2:end,:) - dt * pos_flux_g;

    % BC at x=Lt+Ls - g
    for km = 1:length(omega)
        gn_1((mu<0), end, km) = 0;
    end
    % advection for mu<0 - g
    gn_1(:,1:end-1,:) = gn_1(:,1:end-1,:) - dt * neg_flux_g;

    % Collision - g
    for km = 2:Ny-1
        col_temp = omegaomega_g(:,km,:).*gn(:,km,:); col_temp = reshape(col_temp, length(mu), length(omega));
        values = sum(sum(col_temp))*dmu*domega/normalization;
        col_temp = values*ones(length(mu),1,length(omega)).*((10*omegaomega_g(:,km,:)).^3).*exp(10*omegaomega_g(:,km,:))./((exp(10*omegaomega_g(:,km,:))-1).^2);
        col_g(:,km,:) = (1/Kn^2)*(1/4)*( -omegaomega_g(:,km,:).*gn(:,km,:) + reshape(col_temp, [length(mu), 1, length(omega)]));
    end 
    gn_1 = gn_1 + dt * col_g;
    
    % 5. ADVANCE TIME & RECORD OUTPUT
    time = time + dt;
    j = j + 1;
    fn = fn_1;
    gn = gn_1;
         
    rho_f(:, j) = sum(sum(fn, 3), 1)' * dmu * domega;
    rho_g(:, j) = sum(sum(gn, 3), 1)' * dmu * domega;
end
end