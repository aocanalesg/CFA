!****************************************************************************
!
!  Modules: bstar_sim
!
!****************************************************************************

MODULE mod_bstar_sim

!use imsl_libraries
use mod_global_constant
use mod_parameter_fl
use mod_interp_locate
use rnun_int
use rnnoa_int

IMPLICIT NONE

!****************************************************************************
CONTAINS


!****************************************************************************
!
!  FUNCTION: fcn_nonlinear
!           to solve pnstar_rule and cstar_rule
!
!****************************************************************************

function fcn_nonlinear(vector_laffer,index_laffer)

	implicit none
	
	! specification
	real(dp) :: fcn_nonlinear(nfcn_star)
	real(dp), intent(in) :: vector_laffer(nfcn_star)
	integer, intent(in) :: index_laffer(nstate_star)
	
	real(dp) :: pn_laffer,ln_laffer,kn_laffer,g_laffer_p,kn_t1_laffer,kt_t1_laffer,tau_laffer_p,tot_laffer_p,a_laffer_p,an_laffer,at_laffer,g_laffer,tau_laffer,tot_laffer,z_laffer
	real(dp) :: s_laffer,px_laffer,pg_laffer,tauk_laffer,taul_laffer,yn_laffer,wn_laffer,rkn_laffer,coeff_1,coeff_2,coeff,lt_laffer,yt_laffer,wt_laffer,rkt_laffer,l_laffer,w_laffer
	real(dp) :: lambda_laffer,c_laffer,y_laffer,dn_laffer,in_laffer,qn_laffer,aa,bb,cc,it_laffer,kt_laffer,qt_laffer,i_laffer
	real(dp) :: ap_laffer_p,anp_laffer,atp_laffer,totp_laffer_p,totp_laffer,gp_laffer_p,gp_laffer,taup_laffer_p,taup_laffer,zp_laffer,En_rhs_nug(nug),Et_rhs_nug(nug),En_rhs_nua(nua),Et_rhs_nua(nua),En_rhs_nutau(ntauv),Et_rhs_nutau(ntauv),En_rhs_nutot(nutot),Et_rhs_nutot(nutot)
	real(dp) :: pnp_laffer,lnp_laffer,knp_laffer,sp_laffer,pxp_laffer,pgp_laffer,taukp_laffer,taulp_laffer,ynp_laffer,wnp_laffer
	real(dp) :: coeff_1p,coeff_2p,coeffp,ltp_laffer,ytp_laffer,wtp_laffer,lp_laffer,wp_laffer
	real(dp) :: lambdap_laffer,cp_laffer,rknp_laffer,inp_laffer,qnp_laffer,rktp_laffer,dnp_laffer,aap,bbp,ccp,itp_laffer,qtp_laffer
	real(dp) :: EV_n, EV_t	
	integer  :: ia_laffer,ig_laffer,itau_laffer,itot_laffer,ikn_laffer,ikt_laffer,iuap_laffer,iugp_laffer,iutaup_laffer,itaup_laffer,iutotp_laffer
	real(dp) :: const_1,const_2,const_1p,const_2p
        
    ! inputs
    pn_laffer  = vector_laffer(1)         ! guess for the pn and c rule     
    ln_laffer  = vector_laffer(2)
    kn_laffer = vector_laffer(3)
               
    ia_laffer = index_laffer(1)
    ig_laffer = index_laffer(2)
    itau_laffer = index_laffer(3)
    itot_laffer = index_laffer(4)
    ikn_laffer  = index_laffer(5)
    ikt_laffer  = index_laffer(6)
    
    a_laffer_p = av_grid(ia_laffer)        !percentage change
    g_laffer_p = gv_grid(ig_laffer)        !percentage change
    tau_laffer_p = tauv_grid(itau_laffer)  
    tot_laffer_p = totv_grid(itot_laffer)  !percentage change
    kn_t1_laffer = knv_grid(ikn_laffer)
    kt_t1_laffer = ktv_grid(ikt_laffer)

    z_laffer = z_bar
    
    an_laffer = an_bar*exp(a_laffer_p) 
    at_laffer = at_bar*exp(a_laffer_p) 
    g_laffer = g_bar*exp(g_laffer_p)
    tau_laffer = tau_laffer_p !tau_laffer_bar*exp(tau_laffer_p)
    tot_laffer = tot_bar*exp(tot_laffer_p)  
    
    tauk_laffer = tau_laffer
    taul_laffer = tauk_laffer

               
    !! For test
    !an_laffer = an_bar
    !at_laffer = at_bar
    !g_laffer = g_bar 
    !tau_laffer = tau_bar
    !tot_laffer = tot_bar
    !kn_t1_laffer = kn_bar
    !kt_t1_laffer = kt_bar   
    !tauk_laffer = tau_laffer
    !taul_laffer = tauk_laffer

    
    ! compute pg_t, px_t, s_t
    s_laffer = ((1.0_dp-varphi*pn_laffer**(1.0_dp-chi_c))/(1.0_dp -varphi))**(1.0_dp/(1.0_dp-chi_c))
    px_laffer = s_laffer*tot_laffer
    !pg_laffer = (varphi_g*pn_laffer**(1.0_dp-chi_c)+(1.0_dp-varphi_g)*s_laffer**(1.0_dp-chi_c))**(1.0_dp/(1.0_dp-chi_c))
  
    !! tax rate
    !tauk_laffer = tau_laffer_bar
    !taul_laffer = tauk_laffer           ! tax on labor income

    ! wage, labor and output
    yn_laffer = an_laffer*(kn_t1_laffer**(1.0_dp-alpha_n))*(ln_laffer**alpha_n)
    wn_laffer = alpha_n*pn_laffer*yn_laffer/ln_laffer 
    rkn_laffer = (1.0_dp-alpha_n)*pn_laffer*yn_laffer/kn_t1_laffer

    coeff_1 = pn_laffer*alpha_n*an_laffer*(kn_t1_laffer)**(1.0_dp-alpha_n)*(varphi_l)**(1.0_dp/chi_l)
    coeff_2 = px_laffer*alpha_t*at_laffer*(kt_t1_laffer)**(1.0_dp-alpha_t)*(1.0_dp-varphi_l)**(1.0_dp/chi_l)
    coeff = (coeff_1/coeff_2*ln_laffer**(alpha_n-alpha_t))**(1.0_dp/(alpha_t-1.0_dp-1.0_dp/chi_l))
    lt_laffer = ln_laffer*coeff
    yt_laffer = at_laffer*kt_t1_laffer**(1.0_dp-alpha_t)*lt_laffer**alpha_t
    wt_laffer = px_laffer*alpha_t*yt_laffer/lt_laffer
    rkt_laffer = (1.0_dp-alpha_t)*px_laffer*yt_laffer/kt_t1_laffer
    
    l_laffer = varphi_l**(-1.0_dp/chi_l)*ln_laffer**((1.0_dp+chi_l)/chi_l) + &
            (1.0_dp-varphi_l)**(-1.0_dp/chi_l)*lt_laffer**((1.0_dp+chi_l)/chi_l)
    l_laffer = l_laffer**(chi_l/(1.0_dp+chi_l))
    w_laffer = varphi_l*wn_laffer**(1.0_dp+chi_l)+(1.0_dp-varphi_l)*wt_laffer**(1.0_dp+chi_l)
    w_laffer = w_laffer**(1.0_dp/(1.0_dp+chi_l))
      
    ! solve for consumption    
!    lambda_laffer = phi/((1.0_dp-taul_laffer)*w_laffer*(1.0_dp-l_laffer))
    if (index_logl .eq. 1) then
        lambda_laffer = phi/(1.0_dp-l_laffer)/((1.0_dp-taul_laffer)*w_laffer)
    elseif (index_logl .eq. 2) then    
    lambda_laffer = phi*l_laffer**(sig_l)/((1.0_dp-taul_laffer)*w_laffer)
    end if
    c_laffer = lambda_laffer**(-1.0_dp/sig_c)
    !const_1 = g_laffer**(1.0_dp-1.0_dp/v)
    !const_2 = lambda_laffer
    !call solve_c_bisect(c_laffer,0.0_dp,10.0_dp,-10000.0_dp,1000.0_dp,const_1,const_2)
    !call solve_c_bisect(c_laffer,0.5_dp,2.0_dp,-10000.0_dp,1000.0_dp,const_1,const_2)
  
    ! solve for capital
    y_laffer = pn_laffer*yn_laffer+px_laffer*yt_laffer    
    dn_laffer =  yn_laffer*((pn_laffer)**chi_c)
    in_laffer = kn_laffer - (1-delta_kn)*kn_t1_laffer
    qn_laffer = 1.0_dp+kappa_kn*(in_laffer/kn_t1_laffer-delta_kn)
    
    ! compute investment with adjustment cost (appendix A.45)    
    aa = kappa_kt/(2.0_dp*kt_t1_laffer) 
    bb = 1.0_dp-delta_kt*kappa_kt
    cc = kappa_kt*(delta_kt**2.0_dp)*kt_t1_laffer/2.0_dp+c_laffer+in_laffer+kappa_kn/2.0_dp*(in_laffer/kn_t1_laffer-delta_kn)**2.0_dp*kn_t1_laffer +g_laffer-dn_laffer/varphi
    it_laffer = (-bb+sqrt(bb**2.0_dp-4.0_dp*aa*cc))/(2.0_dp*aa)
    kt_laffer = (1.0_dp-delta_kt)*kt_t1_laffer+it_laffer    
    qt_laffer = 1.0_dp+kappa_kt*(it_laffer/kt_t1_laffer-delta_kt)
    
    i_laffer = in_laffer+it_laffer

    ! compute expectation in t+1
    !------------------------------------------------------------------------------------------------------------------    
    zp_laffer = z_bar             	
    
    En_rhs_nua = 0.0_dp
    Et_rhs_nua = 0.0_dp
    ! shock on a_t+1
	do iuap_laffer = 1,nua
	ap_laffer_p = rho_a*a_laffer_p+ua_grid(iuap_laffer)
    !ap_laffer_p =0.0_dp
	anp_laffer = an_bar*exp(ap_laffer_p)
	atp_laffer = at_bar*exp(ap_laffer_p)

       En_rhs_nug = 0.0_dp
       Et_rhs_nug = 0.0_dp
       ! shock on g_t+1
       do iugp_laffer = 1,nug
          gp_laffer_p = rho_g*g_laffer_p+eta_g*log(y_laffer/y_bar)+ug_grid(iugp_laffer)
          gp_laffer = g_bar*exp(gp_laffer_p)
          !gp_laffer_p = 0.0_dp
          !gp_laffer = g_bar*exp(gp_laffer_p)
          
             En_rhs_nutau = 0.0_dp
             Et_rhs_nutau = 0.0_dp
             ! shock on tau_t+1
             !do iutaup_laffer = 1,nutau
             !taup_laffer_p = rho_tau*tau_laffer_p+utau_grid(iutaup_laffer)
             !taup_laffer_p = 0.0_dp
             !taup_laffer = tau_laffer_bar*exp(taup_laffer_p) 
             !taukp_laffer = taup_laffer
             !taulp_laffer = taukp_laffer
             do itaup_laffer = 1,ntauv
                 taup_laffer = tauv_grid(itaup_laffer)
                 taukp_laffer = taup_laffer
                 taulp_laffer = taukp_laffer
             
          En_rhs_nutot = 0.0_dp
          Et_rhs_nutot = 0.0_dp
          ! shock on tot_t+1
          do iutotp_laffer = 1,nutot
          totp_laffer_p = rho_tot*tot_laffer_p+utot_grid(iutotp_laffer)
          totp_laffer = tot_bar*exp(totp_laffer_p)
          !totp_laffer_p = 0.0_dp
          !totp_laffer = tot_bar*exp(totp_laffer_p)
            
            pnp_laffer = interp6(av_grid,nav,gv_grid,ngv,tauv_grid,ntauv,totv_grid,ntotv,knv_grid,nknv,ktv_grid,nktv,ap_laffer_p,gp_laffer_p,taup_laffer,totp_laffer_p,kn_laffer,kt_laffer,pnstar_rule)
            lnp_laffer = interp6(av_grid,nav,gv_grid,ngv,tauv_grid,ntauv,totv_grid,ntotv,knv_grid,nknv,ktv_grid,nktv,ap_laffer_p,gp_laffer_p,taup_laffer,totp_laffer_p,kn_laffer,kt_laffer,lnstar_rule)                                                 
            knp_laffer = interp6(av_grid,nav,gv_grid,ngv,tauv_grid,ntauv,totv_grid,ntotv,knv_grid,nknv,ktv_grid,nktv,ap_laffer_p,gp_laffer_p,taup_laffer,totp_laffer_p,kn_laffer,kt_laffer,knstar_rule)    
                        
            sp_laffer = ((1.0_dp-varphi*pnp_laffer**(1.0_dp-chi_c))/(1.0_dp -varphi))**(1.0_dp/(1.0_dp-chi_c))                
            pxp_laffer = sp_laffer*totp_laffer
            
            ynp_laffer = anp_laffer*(kn_laffer**(1.0_dp-alpha_n))*(lnp_laffer**alpha_n)
            wnp_laffer = alpha_n*pnp_laffer*ynp_laffer/lnp_laffer 

            coeff_1p = pnp_laffer*alpha_n*anp_laffer*(kn_laffer)**(1.0_dp-alpha_n)*(varphi_l)**(1.0_dp/chi_l)
            coeff_2p = pxp_laffer*alpha_t*atp_laffer*(kt_laffer)**(1.0_dp-alpha_t)*(1.0_dp-varphi_l)**(1.0_dp/chi_l)
            coeffp = (coeff_1p/coeff_2p*lnp_laffer**(alpha_n-alpha_t))**(1.0_dp/(alpha_t-1.0_dp-1.0_dp/chi_l))
            ltp_laffer = lnp_laffer*coeffp
            ytp_laffer = atp_laffer*kt_laffer**(1.0_dp-alpha_t)*ltp_laffer**alpha_t
            wtp_laffer = pxp_laffer*alpha_t*ytp_laffer/ltp_laffer

            lp_laffer = varphi_l**(-1.0_dp/chi_l)*lnp_laffer**((1.0_dp+chi_l)/chi_l) + &
                        (1.0_dp-varphi_l)**(-1.0_dp/chi_l)*ltp_laffer**((1.0_dp+chi_l)/chi_l)
            lp_laffer = lp_laffer**(chi_l/(1.0_dp+chi_l))
            wp_laffer = varphi_l*wnp_laffer**(1.0_dp+chi_l)+(1.0_dp-varphi_l)*wtp_laffer**(1.0_dp+chi_l)
            wp_laffer = wp_laffer**(1.0_dp/(1.0_dp+chi_l))                
                
            ! solve for consumption
            if (index_logl .eq. 1) then
            lambdap_laffer = phi/(1.0_dp-lp_laffer)/((1.0_dp-taulp_laffer)*wp_laffer)
            elseif (index_logl .eq. 2) then    
            lambdap_laffer = phi*lp_laffer**(sig_l)/((1.0_dp-taulp_laffer)*wp_laffer)
            end if
            cp_laffer = lambdap_laffer**(-1.0_dp/sig_c)
            !const_1p = gp_laffer**(1.0_dp-1.0_dp/v)
            !const_2p = lambdap_laffer           
            !call solve_c_bisect(cp_laffer,0.0_dp,10.0_dp,-10000.0_dp,1000.0_dp,const_1p,const_2p)
                
            rknp_laffer = (1.0_dp-alpha_n)*pnp_laffer*ynp_laffer/kn_laffer
            inp_laffer = knp_laffer-(1.0_dp-delta_kn)*kn_laffer
            qnp_laffer = 1.0_dp+kappa_kn*(inp_laffer/kn_laffer-delta_kn)
            rktp_laffer = (1.0_dp-alpha_t)*pxp_laffer*ytp_laffer/kt_laffer
            dnp_laffer =  ynp_laffer*((pnp_laffer)**chi_c)
                            
            ! compute investment with adjustment cost (appendix A.45)
            aap = kappa_kt/(2.0_dp*kt_laffer) 
            bbp = 1.0_dp-delta_kt*kappa_kt
            ccp = kappa_kt*(delta_kt**2.0_dp)*kt_laffer/2.0_dp+cp_laffer+inp_laffer+kappa_kn/2.0_dp*(inp_laffer/kn_laffer-delta_kn)**2.0_dp*kn_laffer +gp_laffer-dnp_laffer/varphi
            itp_laffer = (-bbp+sqrt(bbp**2.0_dp-4.0_dp*aap*ccp))/(2.0_dp*aap) 
            qtp_laffer = 1.0_dp+kappa_kt*(itp_laffer/kt_laffer-delta_kt)
              
            !-------------------------------------------------------------------------------------------------------------- 
            En_rhs_nutot(iutotp_laffer) = beta*lambdap_laffer/lambda_laffer*((1.0_dp-taukp_laffer)*rknp_laffer &
                -kappa_kn/2.0_dp*((inp_laffer/kn_laffer-delta_kn)**2)+kappa_kn*(inp_laffer/kn_laffer-delta_kn)*(inp_laffer/kn_laffer) &
                +qnp_laffer*(1.0_dp-delta_kn))*pr_ua(iuap_laffer)*pr_ug(iugp_laffer)*pr_utot(iutotp_laffer)*pr_tau(itaup_laffer)
                
            Et_rhs_nutot(iutotp_laffer) = beta*lambdap_laffer/lambda_laffer*((1.0_dp-taukp_laffer)*rktp_laffer &
                -kappa_kt/2.0_dp*((itp_laffer/kt_laffer-delta_kt)**2)+kappa_kt*(itp_laffer/kt_laffer-delta_kt)*(itp_laffer/kt_laffer) &
                +qtp_laffer*(1.0_dp-delta_kt))*pr_ua(iuap_laffer)*pr_ug(iugp_laffer)*pr_utot(iutotp_laffer)*pr_tau(itaup_laffer)
             end do ! for iutaup
          
         En_rhs_nutau(itaup_laffer) = (2.0_dp*sum(En_rhs_nutot)-En_rhs_nutot(1)-En_rhs_nutot(nutot))*utot_step/2.0_dp
         Et_rhs_nutau(itaup_laffer) = (2.0_dp*sum(Et_rhs_nutot)-Et_rhs_nutot(1)-Et_rhs_nutot(nutot))*utot_step/2.0_dp
         end do ! for iutaup               
         !
         !   
         En_rhs_nug(iugp_laffer) = sum(En_rhs_nutau)
         Et_rhs_nug(iugp_laffer) = sum(Et_rhs_nutau)
         end do ! for iugp   
                
        En_rhs_nua(iuap_laffer) = (2.0_dp*sum(En_rhs_nug)-En_rhs_nug(1)-En_rhs_nug(nug))*ug_step/2.0_dp
        Et_rhs_nua(iuap_laffer) = (2.0_dp*sum(Et_rhs_nug)-Et_rhs_nug(1)-Et_rhs_nug(nug))*ug_step/2.0_dp


        end do ! for iuap
    
	EV_n = (2.0_dp*sum(En_rhs_nua)-En_rhs_nua(1)-En_rhs_nua(nua))*ua_step/2.0_dp
	EV_t = (2.0_dp*sum(Et_rhs_nua)-Et_rhs_nua(1)-Et_rhs_nua(nua))*ua_step/2.0_dp 
	
    fcn_nonlinear(1) = taul_laffer*(w_laffer*l_laffer)+tauk_laffer*(rkn_laffer*kn_t1_laffer+rkt_laffer*kt_t1_laffer)-z_laffer+c_laffer+i_laffer+kappa_kn/2.0_dp*(in_laffer/kn_t1_laffer-delta_kn)**2.0_dp*kn_t1_laffer &                        
                       &+kappa_kt/2.0_dp*(it_laffer/kt_t1_laffer-delta_kt)**2.0_dp*kt_t1_laffer-y_laffer
    fcn_nonlinear(2) = qn_laffer-EV_n
    fcn_nonlinear(3) = qt_laffer-EV_t
    
    !print *, fcn_nonlinear
    !pause
    
 
end function fcn_nonlinear


!****************************************************************************
!
!  FUNCTION: compute_Tmax
!       compute Tmax and cmax using pnstar_rule, lnstar_rule, and knstar_rule
!
!****************************************************************************
subroutine compute_Tmax(vector_out_Tmax, vector_in_Tmax,index_in_Tmax)

    use mod_interp_locate
	implicit none
	
	! specification
	real(dp), intent(out):: vector_out_Tmax(nrule_Tmax)
	real(dp), intent(in) :: vector_in_Tmax(nfcn_star)
	integer, intent(in) :: index_in_Tmax(nstate_star)

    real(dp) :: pn_Tmax,ln_Tmax,kn_Tmax,kn_t1_Tmax,kt_t1_Tmax,g_Tmax,tot_Tmax,a_Tmax,g_Tmax_p,tau_Tmax_p,tau_Tmax,tot_Tmax_p,a_Tmax_p,an_Tmax,at_Tmax,z_Tmax
    real(dp) :: s_Tmax,px_Tmax,pg_Tmax,taul_Tmax,tauk_Tmax,yn_Tmax,wn_Tmax,rkn_Tmax,coeff_1_Tmax,coeff_2_Tmax,coeff_Tmax,lt_Tmax,yt_Tmax,wt_Tmax,rkt_Tmax,l_Tmax,w_Tmax
    real(dp) :: lambda_Tmax,c_Tmax,y_Tmax,dn_Tmax,in_Tmax,qn_Tmax
    real(dp) :: aa_Tmax,bb_Tmax,cc_Tmax,it_Tmax,kt_Tmax,qt_Tmax,i_Tmax
    real(dp) :: const_1_Tmax,const_2_Tmax
    integer  :: ia_Tmax,ig_Tmax,itau_Tmax,itot_Tmax,ikn_Tmax,ikt_Tmax
        
    ! inputs
    pn_Tmax  = vector_in_Tmax(1)         ! guess for the pn and ln rule     
    ln_Tmax  = vector_in_Tmax(2)
    kn_Tmax  = vector_in_Tmax(3)
    
    ia_Tmax  = index_in_Tmax(1)
    ig_Tmax  = index_in_Tmax(2)
    itau_Tmax = index_in_Tmax(3)
    itot_Tmax  = index_in_Tmax(4)
    ikn_Tmax = index_in_Tmax(5)
    ikt_Tmax = index_in_Tmax(6)
    
    a_Tmax_p = av_grid(ia_Tmax)
    g_Tmax_p = gv_grid(ig_Tmax)
    tau_Tmax_p = tauv_grid(itau_Tmax)
    tot_Tmax_p = totv_grid(itot_Tmax)    
    kn_t1_Tmax = knv_grid(ikn_Tmax)
    kt_t1_Tmax = ktv_grid(ikt_Tmax)

    z_Tmax = z_bar 
       
    an_Tmax = an_bar*exp(a_Tmax_p)  
    at_Tmax = at_bar*exp(a_Tmax_p)  
    g_Tmax = g_bar*exp(g_Tmax_p)
    tau_Tmax = tau_Tmax_p !tau_laffer_bar*exp(tau_Tmax_p)    
    tot_Tmax = tot_bar*exp(tot_Tmax_p)         
    
    ! compute pg_t, px_t, s_t
    s_Tmax = ((1.0_dp-varphi*pn_Tmax**(1.0_dp-chi_c))/(1.0_dp -varphi))**(1.0_dp/(1.0_dp-chi_c))
    px_Tmax = s_Tmax*tot_Tmax
  
    ! tax rate
    !tauk_Tmax = tau_laffer_bar
    !taul_Tmax = tauk_Tmax           ! tax on labor income
    tauk_Tmax = tau_Tmax
    taul_Tmax = tauk_Tmax

    ! wage, labor and output
    yn_Tmax = an_Tmax*(kn_t1_Tmax**(1.0_dp-alpha_n))*(ln_Tmax**alpha_n)
    wn_Tmax = alpha_n*pn_Tmax*yn_Tmax/ln_Tmax 
    rkn_Tmax = (1.0_dp-alpha_n)*pn_Tmax*yn_Tmax/kn_t1_Tmax

    coeff_1_Tmax = pn_Tmax*alpha_n*an_Tmax*(kn_t1_Tmax)**(1.0_dp-alpha_n)*(varphi_l)**(1.0_dp/chi_l)
    coeff_2_Tmax = px_Tmax*alpha_t*at_Tmax*(kt_t1_Tmax)**(1.0_dp-alpha_t)*(1.0_dp-varphi_l)**(1.0_dp/chi_l)
    coeff_Tmax = (coeff_1_Tmax/coeff_2_Tmax*ln_Tmax**(alpha_n-alpha_t))**(1.0_dp/(alpha_t-1.0_dp-1.0_dp/chi_l))
    lt_Tmax = ln_Tmax*coeff_Tmax
    yt_Tmax = at_Tmax*kt_t1_Tmax**(1.0_dp-alpha_t)*lt_Tmax**alpha_t
    wt_Tmax = px_Tmax*alpha_t*yt_Tmax/lt_Tmax
    rkt_Tmax = (1.0_dp-alpha_t)*px_Tmax*yt_Tmax/kt_t1_Tmax
    
    l_Tmax = varphi_l**(-1.0_dp/chi_l)*ln_Tmax**((1.0_dp+chi_l)/chi_l) + &
            (1.0_dp-varphi_l)**(-1.0_dp/chi_l)*lt_Tmax**((1.0_dp+chi_l)/chi_l)
    l_Tmax = l_Tmax**(chi_l/(1.0_dp+chi_l))
    w_Tmax = varphi_l*wn_Tmax**(1.0_dp+chi_l)+(1.0_dp-varphi_l)*wt_Tmax**(1.0_dp+chi_l)
    w_Tmax = w_Tmax**(1.0_dp/(1.0_dp+chi_l))
   
    ! solve for consumption    
    if (index_logl .eq. 1) then
    lambda_Tmax = phi/((1.0_dp-taul_Tmax)*w_Tmax*(1.0_dp-l_Tmax))
    elseif (index_logl .eq. 2) then    
    lambda_Tmax = phi*l_Tmax**(sig_l)/((1.0_dp-taul_Tmax)*w_Tmax)
    end if
    c_Tmax = lambda_Tmax**(-1.0_dp/sig_c)
    !const_1_Tmax = g_Tmax**(1.0_dp-1.0_dp/v)
    !const_2_Tmax = lambda_Tmax
    !call solve_c_bisect(c_Tmax,0.0_dp,10.0_dp,-10000.0_dp,1000.0_dp,const_1_Tmax,const_2_Tmax)
    !call solve_c_bisect(c_Tmax,0.5_dp,2.0_dp,-10000.0_dp,1000.0_dp,const_1_Tmax,const_2_Tmax)
    
    ! solve for capital
    y_Tmax = pn_Tmax*yn_Tmax+px_Tmax*yt_Tmax    
    dn_Tmax =  yn_Tmax*((pn_Tmax)**chi_c)
    in_Tmax = kn_Tmax - (1-delta_kn)*kn_t1_Tmax
    qn_Tmax = 1.0_dp+kappa_kn*(in_Tmax/kn_t1_Tmax-delta_kn)
    
    ! compute investment with adjustment cost (appendix A.45)    
    aa_Tmax = kappa_kt/(2.0_dp*kt_t1_Tmax) 
    bb_Tmax = 1.0_dp-delta_kt*kappa_kt
    cc_Tmax = kappa_kt*(delta_kt**2.0_dp)*kt_t1_Tmax/2.0_dp+c_Tmax+in_Tmax+kappa_kn/2.0_dp*(in_Tmax/kn_t1_Tmax-delta_kn)**2.0_dp*kn_t1_Tmax +g_Tmax-dn_Tmax/varphi
    it_Tmax = (-bb_Tmax+sqrt(bb_Tmax**2.0_dp-4.0_dp*aa_Tmax*cc_Tmax))/(2.0_dp*aa_Tmax)
    kt_Tmax = (1.0_dp-delta_kt)*kt_t1_Tmax+it_Tmax    
    qt_Tmax = 1.0_dp+kappa_kt*(it_Tmax/kt_t1_Tmax-delta_kt)
    
    i_Tmax = in_Tmax+it_Tmax
    
    vector_out_Tmax(1) = taul_Tmax*(w_Tmax*l_Tmax)+tauk_Tmax*(rkn_Tmax*kn_t1_Tmax+rkt_Tmax*kt_t1_Tmax)
    vector_out_Tmax(2) = s_Tmax  
    vector_out_Tmax(3) = kt_Tmax
    vector_out_Tmax(4) = y_Tmax
  
end subroutine compute_Tmax

!****************************************************************************
!
!  FUNCTION: compute_bstar
!       compute the conditional distribution
!
!****************************************************************************

subroutine compute_dist(bstar_mat_dist,sstar_mat_dist)

	implicit none
	real(dp),intent(out)::bstar_mat_dist(nbstar_sim),sstar_mat_dist(nbstar_sim,per_bstar_sim)!,ngv,nknv,nktv
	
	real(dp):: a_dist,g_dist,z_dist,tau_dist,tot_dist,s_dist,y_dist,pg_dist,lambda_dist,r_dist,uc_dist,T_bar,T_dist,kn_dist_t1,kn_dist,kt_dist_t1,kt_dist,supr_dist	
	real(dp), dimension(per_bstar_sim):: shock_ua_dist,shock_ug_dist,shock_utot_dist,shock_utau_dist,tau_random
	real(dp), dimension(per_bstar_sim):: svec_dist,avec_dist,gvec_dist,totvec_dist,tauvec_dist,lambdavec_dist,rvec_dist
	real(dp), dimension(per_bstar_sim):: knvec_dist,ktvec_dist,yvec_dist,pgvec_dist,suprvec_dist,suprvec_z_dist
	integer :: is_dist,it_dist,ia,ig,itau,itot,ikn,ikt
	
	!$OMP single     
    !do ig = 1,ngv 
    !do ikn = 1,nknv
    !do ikt = 1,nktv 
	    !print *, ig,ikn,ikt
	    call itime(time_now)
        write (*,10), time_now
        10 format ( ' time ', i2.2, ':', i2.2, ':', i2.2 )

    call load_data(pr_tau, cdf_tau, tau_random)

    ! at each state, do simulations of # nbstar_sim 
        
    do is_dist = 1,nbstar_sim	    
	    ! generate shocks
        call d_rnnoa(shock_ua_dist(:))
        shock_ua_dist = std_ua*shock_ua_dist  	    
        call d_rnnoa(shock_ug_dist(:))
        shock_ug_dist = std_ug*shock_ug_dist  
        !call d_rnnoa(shock_utau_dist(:))
        !shock_utau_dist = std_utau*shock_utau_dist                
        call d_rnnoa(shock_utot_dist(:))
        shock_utot_dist = std_utot*shock_utot_dist  
	    
        ! period t=1
        it_dist = 1	
        avec_dist(it_dist) = 0.0_dp
        gvec_dist(it_dist) =  0.0_dp !gv_grid(ig)
        !tauvec_dist(it_dist) = 0.0_dp
        totvec_dist(it_dist) = 0.0_dp 
        knvec_dist(it_dist) = kn_bar  !knv_grid(ikn)
        ktvec_dist(it_dist) = kt_bar !ktv_grid(ikt)
        yvec_dist(it_dist) = y_bar
                
        T_dist =             interp6(av_grid,nav,gv_grid,ngv,tauv_grid,ntauv,totv_grid,ntotv,knv_grid,nknv,ktv_grid,nktv,avec_dist(it_dist),gvec_dist(it_dist),tau_random(it_dist),totvec_dist(it_dist),knvec_dist(it_dist),ktvec_dist(it_dist),Tstar_rule)
        lambda_dist =        interp6(av_grid,nav,gv_grid,ngv,tauv_grid,ntauv,totv_grid,ntotv,knv_grid,nknv,ktv_grid,nktv,avec_dist(it_dist),gvec_dist(it_dist),tau_random(it_dist),totvec_dist(it_dist),knvec_dist(it_dist),ktvec_dist(it_dist),lambdastar_rule)       
        svec_dist(it_dist) = interp6(av_grid,nav,gv_grid,ngv,tauv_grid,ntauv,totv_grid,ntotv,knv_grid,nknv,ktv_grid,nktv,avec_dist(it_dist),gvec_dist(it_dist),tau_random(it_dist),totvec_dist(it_dist),knvec_dist(it_dist),ktvec_dist(it_dist),sstar_rule)
        !sstar_mat_dist(is_dist,1)=svec_dist(it_dist)
        
        lambdavec_dist(it_dist)=lambda_dist
        r_dist=1.0_dp/beta
        rvec_dist(it_dist)=r_dist

        !uc_dist = lambda_dist
        !ucvec_dist(it_dist) = uc_dist
        
        !suprvec_dist(it_dist)=(beta**it_dist)*uc_dist*(T_dist-g_bar*exp(gvec_dist(it_dist))-z_bar)
        !suprvec_dist(it_dist)=(beta**it_dist)*(T_dist-pg_dist*g_bar*exp(gvec_dist(it_dist))-z_bar)
        !suprvec_dist(it_dist)=(beta**it_dist)*(T_dist-pg_dist*g_bar*exp(gvec_dist(it_dist))-z_bar)/svec_dist(it_dist)
        suprvec_dist(it_dist)=(beta**it_dist)*(T_dist-g_bar*exp(gvec_dist(it_dist))-z_bar)/svec_dist(it_dist)  
        
        
        ! from t = 2 on             
        do it_dist = 2,per_bstar_sim            
         
            ! exog states
            a_dist = rho_a*avec_dist(it_dist-1)+shock_ua_dist(it_dist)
            g_dist = rho_g*gvec_dist(it_dist-1)+eta_g*log(yvec_dist(it_dist-1)/y_bar)+shock_ug_dist(it_dist)
            !tau_dist = rho_tau*tauvec_dist(it_dist-1)+shock_utau_dist(it_dist)
            tot_dist = rho_tot*totvec_dist(it_dist-1)+shock_utot_dist(it_dist)
	        z_dist = z_bar
	        
            ! interpolate Tstar_rule and cstar_rule 
            kn_dist_t1 = knvec_dist(it_dist-1)
            kt_dist_t1 = ktvec_dist(it_dist-1)
            kn_dist = interp6(av_grid,nav,gv_grid,ngv,tauv_grid,ntauv,totv_grid,ntotv,knv_grid,nknv,ktv_grid,nktv,a_dist,g_dist,tau_random(it_dist),tot_dist,kn_dist_t1,kt_dist_t1,knstar_rule)
            kt_dist = interp6(av_grid,nav,gv_grid,ngv,tauv_grid,ntauv,totv_grid,ntotv,knv_grid,nknv,ktv_grid,nktv,a_dist,g_dist,tau_random(it_dist),tot_dist,kn_dist_t1,kt_dist_t1,ktstar_rule)
            T_dist =  interp6(av_grid,nav,gv_grid,ngv,tauv_grid,ntauv,totv_grid,ntotv,knv_grid,nknv,ktv_grid,nktv,a_dist,g_dist,tau_random(it_dist),tot_dist,kn_dist_t1,kt_dist_t1,Tstar_rule)
            s_dist =  interp6(av_grid,nav,gv_grid,ngv,tauv_grid,ntauv,totv_grid,ntotv,knv_grid,nknv,ktv_grid,nktv,a_dist,g_dist,tau_random(it_dist),tot_dist,kn_dist_t1,kt_dist_t1,sstar_rule)      
            y_dist =  interp6(av_grid,nav,gv_grid,ngv,tauv_grid,ntauv,totv_grid,ntotv,knv_grid,nknv,ktv_grid,nktv,a_dist,g_dist,tau_random(it_dist),tot_dist,kn_dist_t1,kt_dist_t1,ystar_rule)  

            !uc_dist = lambda_dist
            !supr_dist=(beta**it_dist)*uc_dist*(T_dist-g_bar*exp(g_dist)-z_dist)
            !supr_dist=(beta**it_dist)*(T_dist-pg_dist*g_bar*exp(g_dist)-z_dist)            
            !supr_dist=(beta**it_dist)*(T_dist-pg_dist*g_bar*exp(g_dist)-z_dist)/s_dist
            !r_dist=1.0_dp/beta*lambdavec_dist(it_dist-1)/lambda_dist         !endogenous interest rate
            r_dist=1.0_dp/beta                                                !exogenous interest rate
            supr_dist=(beta**it_dist)*(T_dist-g_bar*exp(g_dist)-z_dist)/s_dist
            
	        
            ! save the paths
            knvec_dist(it_dist) = kn_dist
            ktvec_dist(it_dist) = kt_dist
            yvec_dist(it_dist) = y_dist
            avec_dist(it_dist) = a_dist
            gvec_dist(it_dist) = g_dist
            lambdavec_dist(it_dist)=lambda_dist
            rvec_dist(it_dist)=r_dist
            tauvec_dist(it_dist) = tau_dist            
            totvec_dist(it_dist) = tot_dist
            svec_dist(it_dist) = s_dist  	
            suprvec_dist(it_dist) = supr_dist
            
            !sstar_mat_dist(is_dist,it_dist)=svec_dist(it_dist)
            !print *, sstar_mat_dist(is_dist,it_dist)
            
        end do	
               
        !bstar_mat_dist(is_dist,ig,ikn,ikt)=sum(suprvec_dist(perb_bstar_sim:per_bstar_sim)*(beta**(-perb_bstar_sim))/ucvec_dist(perb_bstar_sim))*beta_pol
        !bstar_mat_dist(is_dist,ig,ikn,ikt)=sum(suprvec_dist(perb_bstar_sim:per_bstar_sim)*(beta**(-perb_bstar_sim)))*beta_pol
        bstar_mat_dist(is_dist)=sum(suprvec_dist(perb_bstar_sim:per_bstar_sim)*(beta**(-perb_bstar_sim)))*beta_pol
        sstar_mat_dist(is_dist,:)=svec_dist(:)
        !
        !print *, is_dist,sstar_mat_dist(is_dist,:)
        !pause
        
    end do 
    
    !end do
    !end do
    !end do

    !$OMP end single  
   
end subroutine compute_dist


END MODULE mod_bstar_sim

