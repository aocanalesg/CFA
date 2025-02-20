% Purpose: calibrate the model to Greek economy
% Note: mu_zy = 1.026 with p11_z = 1-1/20 or 1-1/40 --> very disbursed dist
% Last update: 09/02/2011 by Huixin Bi

%--- calibration in the paper---------
taubar = 0.32; % = tauss
govery = 0.167; % = gyss
zovery = 0.1334; % = zyss
zeta_z_log = -0.45;         % estimated elasticity (log) % = eta
gammabar = 0.42; % = gamma
rho_zy = 0.9899;            % not used 
mu_zy = 1.015;              % explosive transfers  %1.026;   % = mu
p11_z = 1-1/40;             % regime transfer parameter for rs^z % = p00
p22_z = 1-1/40;             % 1-1/20; % = p11

% wgrid = [1; 1];             % political risk regime
% p11_w = 1-1/13;             % regime transfer parameter for w
% p22_w = 1-1/13;

%----------------------------------------------------
% steady states under gamma-tax policy
Lbar = 0.75;             % Lbar is leisure, LL is labor % = 1 - nss
LL = 1-Lbar;             % in order to nail down CHI  % = nss
Abar = 1; % = Ass
beta = 0.95;    
qbar = beta; 
delta_bar = 0.03;

covery = 1 - govery;
bovery = (taubar - govery - zovery)/(1-qbar);

ybar = LL*Abar; 
zbar = ybar*zovery;
bbar= ybar*bovery;                    
cbar = ybar*covery; 
gbar = ybar*govery; 
phi = Abar*(1-taubar)*Lbar/cbar;
taumaxbar = 1+phi-sqrt((1+phi)*phi*(Abar-gbar)/Abar);
zeta_z = zeta_z_log*zbar/Abar;
mu_z = mu_zy;

%=================================================
% grid of shock u_t  (A)
nuA = 11; 
std_uA = 0.033;
sigma_uA = std_uA^2;
uA_std_max = 4;
uA_max = uA_std_max*std_uA; 
uA_min = -uA_max;
uA_step = (uA_max-uA_min)/(nuA-1);
uA_grid = uA_min:uA_step:uA_max;
pr_uA=(2*pi)^(-0.5)/std_uA*exp(-(uA_grid.^2)/(std_uA^2)/2);	
%[uAgrid,uAprob] = tauchen(nuA,0,0,sqrt(sigma_uA),uA_std_max);

% grid of shock e_t  (log g_t/g)    
nuG = 11;  
std_uG = 0.03;
sigma_uG = std_uG^2; 
uG_std_max = 4;
uG_max = uG_std_max*std_uG; 
uG_min = -uG_max;
uG_step = (uG_max-uG_min)/(nuG-1);
uG_grid = uG_min:uG_step:uG_max;
pr_uG=(2*pi)^(-0.5)/std_uG*exp(-(uG_grid.^2)/(std_uG^2)/2);
% [uGgrid,uGprob] = tauchen(nuG,0,0,sqrt(sigma_uG),uG_std_max);

% grid of A
A_std_max = 4;
nA = 5;
rho_A = 0.45;
Amin = Abar*exp(-A_std_max/sqrt(1-rho_A^2)*std_uA);  
Amax = Abar*exp(A_std_max/sqrt(1-rho_A^2)*std_uA); 
Astep = (Amax-Amin)/(nA-1);
Agrid = Amin:Astep:Amax;
	
% grid of g
g_std_max = 4;
ng = 5;
rho_g = 0.426;
gmin = gbar*exp(-g_std_max/sqrt(1-rho_g^2)*std_uG);  
gmax = gbar*exp(g_std_max/sqrt(1-rho_g^2)*std_uG); 
gstep = (gmax-gmin)/(ng-1);
ggrid = gmin:gstep:gmax;

% transfer regime shock
zrsprob = [p11_z 1-p11_z; 1-p22_z p22_z];
zrsprob_cum = cumsum(zrsprob,2);
zrsgrid = [1;2];
nzrs = length(zrsgrid);

% grid of z
nz = 21;9;
zmin  = zbar+zeta_z*(Amax-Abar);        
zmax = zbar*mu_z^(1/(1-p11_z))+zeta_z*(Amin-Abar); % zbar+zeta_z*(Amin-Abar);  %
zstep = (zmax-zmin)/(nz-1);
zgrid = zmin:zstep:zmax;

% political parameter shock 
% nw = length(wgrid);
% wprob = [p11_w 1-p11_w; 1-p22_w p22_w];
% wprob_cum = cumsum(wprob,2);


