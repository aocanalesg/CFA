
function [bstar] = fcn_bstar_mcmc(zrs0, z0, uA0, uG0, uc0, per, perb, ...
    sigma_uA, sigma_uG, rho_A, rho_g, mu_z, zbar, ...
    Abar, phi, gbar, zrsprob_cum, beta, zeta_z)%, windex0, wprob_cum, wgrid)
% calculate conditional b*_t
% b*_t = b*(A_t, \theta_t)

% initialize the shock process
shock_uA = randn(per,1)*sqrt(sigma_uA);
shock_uG = randn(per,1)*sqrt(sigma_uG);
shock_w = rand(per,1);
shock_zrs = rand(per,1);

% initialize the vector vars
A_vec = zeros(per,1);
A_vec(1) = Abar*exp(uA0);
g_vec = zeros(per,1);
g_vec(1) = gbar*exp(uG0);
% windex_vec = zeros(per,1);
% windex_vec(1) = windex0;
zrs_vec = zeros(per,1);
zrs_vec(1) = zrs0;
z_vec = zeros(per,1);
z_vec(1) = z0;
uc_vec = zeros(per,1);
uc_vec(1) = uc0;

T_vec = zeros(per,1);
supr_vec = zeros(per,1);

% t = 2... infty
% each time period (it)   
for it = 2:per

    % productivity and govt spending  
    A = Abar^(1-rho_A)*A_vec(it-1)^rho_A*exp(shock_uA(it));
    A_vec(it) = A;
    g = gbar^(1-rho_g)*g_vec(it-1)^rho_g*exp(shock_uG(it));
    g_vec(it) = g;
    taumax = 1+phi-sqrt((1+phi)*phi*(A-g)/A);
    cmax = (A-g)*(1-taumax)/(1+phi-taumax);
    Lmax = 1-(cmax+g)/A;
    T_vec(it) = A*(1-Lmax)*taumax;
    uc_vec(it) = 1/cmax;
    
    % transfer regime
    zrs_vec(it) = min(find(zrsprob_cum(zrs_vec(it-1),:)-shock_zrs(it)>0));    
    if zrs_vec(it) == 1
        z_vec(it)= zbar+zeta_z*(A-Abar);
    else
        z_vec(it)= mu_z*z_vec(it-1)+zeta_z*(A-Abar);
    end
    z = z_vec(it);

    % political risk
    %theta = wgrid(1);
%     windex_vec(it) = min(find(wprob_cum(windex_vec(it-1),:)-shock_w(it)>0));
%     theta = wgrid(windex_vec(it));
    
    % compute discounted surplus
    supr_vec(it) = (beta)^it*uc_vec(it)*(T_vec(it) - g - z);%*theta;
end

bstar = sum(supr_vec(perb:end)*beta^(-perb)/uc_vec(perb));
   