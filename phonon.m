function temperature=phonon(dx, x_max, dt, T, dmu, domega, omega_min, omega_max, eta_func,m,n,Kn)
% output: temperature as a function of time, Nt*1 
%% sample grids

x = [0:dx:x_max]; 
mu = [-1+dmu/2:dmu:1-dmu/2];
omega = [omega_min:domega:omega_max];
[xx,mumu,omegaomega] = meshgrid(x,mu,omega);
pos_mu = repmat(max(mu',0),[1,length(x)-1,length(omega)]);
neg_mu = repmat(min(mu',0),[1,length(x)-1,length(omega)]);
  
%eta_fun = @(omega) 0.8*sin(pi*omega/5);
eta = eta_func;
  
normalization = sum((10*omega).^3.*exp(10*omega)./(exp(10*omega)-1).^2)*domega*dmu*sum(length(mu));
%% initialization and boundary condition
f = zeros(length(mu),length(x),length(omega));
g = f;
col = f;
temperature = [];

bdy_func = @(t_s,mu_s,omega_s) (1/dmu)*(1/dt)*(1/domega).*eq(mu_s,mu(m)).*eq(omega_s,omega(n)).*eq(t_s,0) + mu_s-mu_s + t_s-t_s;
 %bdy_func = @(mu_s,omega_s) mu_s.*(1-mu_s).*exp(-omega_s);
% f(:,1,:) = bdy(:,1,:);

%% evolution
time = 0; 
%dt=min(Kn*dx/2, Kn^2);

%N_T=(T-time)/dt+1;
%P1=zeros(length(x),N_T);
j=1;
while time < T
    
bdy = bdy_func(time,max(mumu,0),omegaomega);
f((mu>0),1,:) = bdy((mu>0),1,:);
    time = time + dt;
     j=j+1;
    mu_px_f = (1/Kn)*((omegaomega(:,2:end,:).*f(:,2:end,:))-(omegaomega(:,1:end-1,:).*f(:,1:end-1,:)))/dx; %%dmu
    pos_flux = pos_mu.*mu_px_f;
    neg_flux = neg_mu.*mu_px_f;
  
    % transport
    g(:,2:end,:) = f(:,2:end,:) - dt * pos_flux;
    g(:,1:end-1,:) = g(:,1:end-1,:) - dt * neg_flux;
    % collision
    for km = 2:length(x)-1
        col_temp = omegaomega(:,km,:).*f(:,km,:); col_temp = reshape(col_temp,length(mu),length(omega));
        values = sum(sum(col_temp))*dmu*domega/normalization;
        col_temp = values*ones(length(mu),1,length(omega)).*((10*omegaomega(:,km,:)).^3).*exp(10*omegaomega(:,km,:))./((exp(10*omegaomega(:,km,:))-1).^2);
        col(:,km,:) = (1/Kn^2)*( -omegaomega(:,km,:).*f(:,km,:) + reshape(col_temp,[length(mu),1,length(omega)]));
     
%     max(P1)
  
    end
   %P1(:,j)=sum(sum(sum(col),1),3);
     
    g  = g + dt * col;
%         if time > 2
%             a=1;
%         end

    % boundaries
    bdy = bdy_func(time,max(mumu,0),omegaomega);

    g((mu>0),1,:) = bdy((mu>0),1,:);
     
    for km = 1:length(omega)
        outgoing = f((mu>0),end,km);
        g((mu<0),end,km) = eta(km)*flipud(outgoing);
 
    end
      
    f = g;
     
    temp_t =omegaomega(:,1,:).* f(:,1,:); temp_t = sum(sum(temp_t))*dmu*domega;
    temperature = [temperature;temp_t];
      
%figure(7); mesh(reshape(sum(g(:,:,:)*domega,3),[length(mu),length(x)]));title(['time = ',num2str(time)]);pause(0.001);
  
%figure(8); mesh(reshape(sum(g(:,:,:)*dmu,1),[length(x),length(omega)]));title(['time = ',num2str(time)]);pause(0.01);

end
