%---------------------------------------------------------
% Purpose:  calculate b* distribution (conditional)
% Model:    z_t = z+zeta_z(A_t-A) // mu_z z_t-1+zeta_z(A_t-A)
%           stochastic \theta_t and A_t, g_t, zrs_t
% Last updated:     Aug 31/2011
%----------------------------------------------------------

clc
clear all
close all

cd('C:\Users\GaboP\OneDrive - Ministerio de Hacienda\CFA\11. L�mite de Deuda\Insumos\Bi (2012) - Sovereign default risk premia, fiscal limits and fiscal policy\MatlabCode_FiscalLimit');

compute_index = 1;
compute_vmat_index = 1;
save_v_index = 1;
save_bstar_index = 1;
per = 200; 
N = 500; 
perb = 1; 

parameter;

if compute_index == 1
%   % enable the following two lines for parallel computing
    % distcomp.feature( 'LocalUseMpiexec', false );
    % matlabpool; 
    
    bstar_mat = zeros(N,nuA,nuG,nzrs);
    
    for iuA = 1:nuA
        for iuG = 1:nuG
            for izrs = 1:nzrs
                tic
                parfor is = 1:N
                    % t = 1    
                    uA0 = uA_grid(iuA);    
                    uG0 = uG_grid(iuG);
                    zrs0 = zrsgrid(izrs);   
                    windex0 = 1;
                    z0 = zbar;

                    A0 = Abar*exp(uA0); 
                    g0 = gbar*exp(uG0);
                    taumax = 1+phi-sqrt((1+phi)*phi*(A0-g0)/A0);
                    cmax = (A0-g0)*(1-taumax)/(1+phi-taumax);
                    uc0 = 1/cmax;                  

                    bstar_mat(is,iuA,iuG,izrs) = fcn_bstar_mcmc(windex0, zrs0, z0, uA0, uG0, uc0, per, perb, ...
                            sigma_uA, sigma_uG, rho_A, rho_g, mu_z, zbar, ...
                            Abar, phi, gbar, zrsprob_cum, wprob_cum, beta, zeta_z, wgrid);
                end
                toc 
            end
        end
    end    
%     matlabpool close
end

compute_v
plot_v
save_v



