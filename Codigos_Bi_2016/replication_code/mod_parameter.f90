!****************************************************************************
!
!  Modules: assign parameter values: Argentina
!
!****************************************************************************

MODULE mod_parameter_fl

use mod_global_constant
use mod_interp_locate

use neqnf_int 
use umach_int

IMPLICIT NONE

!-------------------------------------------------------------------------------
! Model Parameters
!-------------------------------------------------------------------------------
integer :: openstatus
integer, parameter  :: index_load    = 1                          ! = 1: rule; =2: ss; =3: ini =4:
integer, parameter  :: index_logl    = 2                          ! = 1, logl; =2, crra L

integer, parameter  :: nav          = 5  
integer, parameter  :: ngv          = 5 
integer, parameter  :: ntauv        = 11
integer, parameter  :: ntotv        = 7   
integer, parameter	:: nknv 	    = 5                        
integer, parameter	:: nktv 	    = 5                        

integer, parameter  :: nav_ini      = 5  
integer, parameter  :: ngv_ini      = 5 
integer, parameter  :: ntauv_ini    = 11
integer, parameter  :: ntotv_ini    = 5   
integer, parameter	:: nknv_ini 	= 5                        
integer, parameter	:: nktv_ini 	= 5                        

!--------------------------------------------------------------------------	
! Argentina
!--------------------------------------------------------------------------	
! deep parameters 
real(dp), parameter :: beta_pol     = 1.00_dp               ! political discount factor       
real(dp), parameter :: beta 	    = 0.981_dp	            ! discount rate
real(dp), parameter :: q_bar 	    = beta	   
real(dp), parameter :: delta_kn     = 0.03_dp               ! depreciation rate (no-tradable)
real(dp), parameter :: delta_kt     = 0.03_dp               ! depreciation rate (tradable)
real(dp), parameter :: chi_c        = 0.76_dp               ! elasticity of substitution between tradable and non-tradable
real(dp), parameter :: chi_l        = 0.69_dp !0.681_dp                ! elasticity between the two labors
real(dp), parameter :: varphi       = 0.53_dp               ! home goods ratio in c and i
real(dp), parameter :: alpha_n      = 0.47_dp                ! labor share in the nontradable sector
real(dp), parameter :: alpha_t      = 0.47_dp               ! labor share in the tradable sector
real(dp), parameter :: kappa_kn     = 1.66_dp !1.632_dp                ! capital ajustment cost (non-tradable)
real(dp), parameter :: kappa_kt     = 1.66_dp !1.632_dp                ! capital ajustment cost (tradable)
real(dp), parameter :: sig_c        = 2.68_dp               ! inverse of elasticity of substitution in consumption
real(dp), parameter :: sig_l        = 3.54_dp                ! inverse of frisch elasticity of labor supply

real(dp), parameter :: L_bar        = 0.26_dp                ! steady state labor
real(dp), parameter :: tau_bar      = 0.194_dp  !0.227_dp
real(dp), parameter :: govery       = 0.121_dp  !0.1476_dp
real(dp), parameter :: sb           = 0.448_dp*4.0_dp        ! quarterly debt to GDP ratio

real(dp), parameter :: an_bar       = 1.0_dp
real(dp), parameter :: s_bar        = 1.0_dp                ! ss real exchange rate
real(dp), parameter :: pn_bar       = 1.0_dp
real(dp), parameter :: tot_bar      = 1.0_dp 

real(dp), parameter :: tau_laffer_bar = tau_bar !0.207_dp       ! maximum tax rate

!! estimated parameter
!real(dp), parameter :: eta_g        = 0.24_dp  ! procyclical response of g
!real(dp), parameter :: rho_g        = 0.39_dp  
!real(dp), parameter :: rho_tau      = 0.80_dp 
!real(dp), parameter :: rho_a 	    = 0.85_dp 
!real(dp), parameter :: rho_tot      = 0.78_dp  
!real(dp), parameter :: std_ug 	    = 2.42_dp*0.01_dp 
!real(dp), parameter :: std_utau 	= 5.55_dp*0.01_dp
!real(dp), parameter :: std_ua 	    = 2.54_dp*0.01_dp       
!real(dp), parameter :: std_utot 	= 3.44_dp*0.01_dp     

real(dp), parameter :: eta_g        = 0.12_dp  ! procyclical response of g
real(dp), parameter :: rho_g        = 0.47_dp  
real(dp), parameter :: rho_tau      = 0.80_dp 
real(dp), parameter :: rho_a 	    = 0.84_dp 
real(dp), parameter :: rho_tot      = 0.78_dp  
real(dp), parameter :: std_ug 	    = 2.41_dp*0.01_dp 
real(dp), parameter :: std_utau 	= 5.55_dp*0.01_dp
real(dp), parameter :: std_ua 	    = 2.54_dp*0.01_dp       
real(dp), parameter :: std_utot 	= 3.44_dp*0.01_dp     

!-------------------------------------------------------------------------------
! shock (ua,ug,uz,utot)
!-------------------------------------------------------------------------------
! grid parameter (ua)
real(dp), parameter	:: ua_std_max    = 6.0_dp                 ! max std (uagrid)	
integer, parameter	:: nua           = 13 
real(dp), parameter :: ua_max        = ua_std_max*std_ua
real(dp), parameter :: ua_min        = -ua_max
real(dp), parameter :: ua_step	     = (ua_max-ua_min)/(nua-1.0_dp)

! grid parameter (ug)	
real(dp), parameter :: ug_std_max    = 6.0_dp !6.0_dp         ! max std (uggrid)		
integer, parameter  :: nug           = 13  	                		
real(dp), parameter :: ug_max        = ug_std_max*std_ug
real(dp), parameter :: ug_min        = -ug_max
real(dp), parameter :: ug_step	     = (ug_max-ug_min)/(nug-1.0_dp)	

!! grid parameter (utau)
!real(dp), parameter  :: utau_std_max  = 6.0_dp      	
!integer, parameter  :: nutau         = 13 	         		
!real(dp), parameter :: utau_max      = utau_std_max*std_utau
!real(dp), parameter :: utau_min      = -utau_max
!real(dp), parameter :: utau_step	 = (utau_max-utau_min)/(nutau-1.0_dp)

! grid parameter (utot)
real(dp), parameter  :: utot_std_max = 6.0_dp      	
integer, parameter  :: nutot         = 13 	         		
real(dp), parameter :: utot_max      = utot_std_max*std_utot
real(dp), parameter :: utot_min      = -utot_max
real(dp), parameter :: utot_step	 = (utot_max-utot_min)/(nutot-1.0_dp)

!-------------------------------------------------------------------------------
! state variables (a_t,g_t,z_t,tot_t,kn_{t-1},kt_{t-1})
! (1-delta_t)b_{t-1}* if nonlinear model
!-------------------------------------------------------------------------------
! grid parameter (a_t)
real(dp), parameter :: a_std_max    = 4.0_dp

! grid parameter (g_t)
real(dp), parameter :: g_std_max    = 4.0_dp  !7.0_dp

!! grid parameter (tau_t)
!real(dp), parameter :: tau_std_max  = 4.0_dp
!integer, parameter  :: ntauv        = 5

! grid parameter (tot_t)
real(dp), parameter :: tot_std_max  = 6.5_dp

! grid parameter (kn_{t-1})         
real(dp), parameter :: kn_margin_l  = 0.3_dp                   ! deviation from the ss when construct gird
real(dp), parameter :: kn_margin_h  = 0.3_dp                   ! deviation from the ss when construct gird

! grid parameter (kt_{t-1})         
real(dp), parameter :: kt_margin_l  = 0.3_dp                   ! deviation from the ss when construct gird
real(dp), parameter :: kt_margin_h  = 0.3_dp                   ! deviation from the ss when construct gird

!-------------------------------------------------------------------------------
! parameter to compute fiscal limits
!-------------------------------------------------------------------------------

integer, parameter  :: per_bstar_sim    = 1000                  ! simulation periods
integer, parameter  :: perb_bstar_sim   = 1                     ! starting period
integer, parameter  :: nbstar_sim       = 10000                ! repeated simulations

!-------------------------------------------------------------------------------
! rules for computingv the fiscal limits
integer, parameter  :: nfcn_star = 3                     ! number of equations in fcn_nonlinear
integer, parameter  :: nstate_star = 6                   ! number of state var	
integer, parameter  :: nrule_star = 3                    ! number of rules 
integer, parameter  :: nrule_Tmax = 4                    ! number of rules in Tmax
real(dp), dimension(nav,ngv,ntauv,ntotv,nknv,nktv):: pnstar_rule,lnstar_rule,knstar_rule
real(dp), dimension(nav,ntauv,ntotv,nknv,nktv):: pnstar_rule_nog,lnstar_rule_nog,knstar_rule_nog
real(dp), dimension(nav_ini,ngv_ini,ntauv_ini,ntotv_ini,nknv_ini,nktv_ini) :: pnstar_rule_ini,lnstar_rule_ini,knstar_rule_ini
real(dp), dimension(nav,ngv,ntauv,ntotv,nknv,nktv):: pnstar_new,lnstar_new,knstar_new
real(dp), dimension(nav,ngv,ntauv,ntotv,nknv,nktv):: Tstar_rule,lambdastar_rule,ktstar_rule,ystar_rule,sstar_rule,pgstar_rule
real(dp), dimension(nav,ngv,ntauv,ntotv,nknv,nktv,nrule_Tmax):: rules_star
real(dp) :: bstar_mat(nbstar_sim),sstar_mat(nbstar_sim,per_bstar_sim)!,ngv,nknv,nktv

!-------------------------------------------------------------------------------
! paramters for shock simulation (determinancy)
!-------------------------------------------------------------------------------

integer,  parameter :: T_shock        = 3080 
real(dp):: bsim_shock(T_shock),gsim_shock(T_shock),shock_ug(T_shock)
real(dp):: tausim_shock(T_shock),ssim_shock(T_shock),ynsim_shock(T_shock),lnsim_shock(T_shock),wnsim_shock(T_shock),ytsim_shock(T_shock),ltsim_shock(T_shock),wtsim_shock(T_shock),qsim_shock(T_shock),rsim_shock(T_shock)
integer:: tsim_shock
integer:: index_g

real(dp), parameter :: tol_bsim_indet   = 0.1_dp
integer,  parameter :: T_noshock        = 4000 
real(dp):: bsim_noshock(T_noshock),gsim_noshock(T_noshock),asim_noshock(T_noshock),tausim_noshock(T_noshock),knsim_noshock(T_noshock),ktsim_noshock(T_noshock)
real(dp):: diff_bsim_noshock
integer:: index_noeqm, tsim_noshock


!! other parameters
!!------------------------------------------------------------------ 	
integer     :: ib,ia,ig,itau,itot,ikn,ikt,index_compute	    
integer     :: iter, index_t_star(nstate_star),rc_nsolve
real(dp)    :: vector_t_star(nfcn_star),vector_out_t_star(nrule_Tmax)
real(dp)    :: pnstar_max_diff,lnstar_max_diff,knstar_max_diff
real(dp), dimension(nav,ngv,ntauv,ntotv,nknv,nktv):: pnstar_diff_iter,lnstar_diff_iter,knstar_diff_iter
integer     :: time_now(3)                          ! measure the timingv

!-------------------------------------------------------------------------------
! paramters for solution precision
!-------------------------------------------------------------------------------

real(dp), parameter :: tol_conv     = (0.1_dp)**6                ! conv crit for the whole loop (all state var)
real(dp), parameter :: tol_nsolve   = (0.1_dp)**7                ! con crit for nonlinear solver
integer, parameter  :: nmax_nsolve  = 10                         ! max number of iterations in nonlinear solver


! Construct Grid and Prob matrix
!-------------------------------------------------------------------------------
real(dp) :: ua_grid(nua),ug_grid(nug),utot_grid(nutot) 
real(dp) :: av_grid(nav),gv_grid(ngv),tauv_grid(ntauv),totv_grid(ntotv),knv_grid(nknv),ktv_grid(nktv)
real(dp) :: av_grid_ini(nav_ini),gv_grid_ini(ngv_ini),tauv_grid_ini(ntauv_ini),totv_grid_ini(ntotv_ini),knv_grid_ini(nknv_ini),ktv_grid_ini(nktv_ini)
real(dp) :: pr_ua(nua), pr_ug(nug),pr_tau(ntauv),pr_utot(nutot),cdf_tau(ntauv),tau_random(per_bstar_sim)

contains

!****************************************************************************

!subroutine: solve the nonlinear system in steady state

!****************************************************************************

subroutine solve_ss(x, fnorm)

    integer n
    parameter (n=2)
    
    integer k, nout
    real(dp) :: fnorm, x(n), xguess(n)

    data xguess/1.0_dp, 0.3_dp/   !guess for total output in units
   
    call umach (2, nout)
    call neqnf (fcn_ss, x, xguess=xguess, fnorm=fnorm)
     
end subroutine solve_ss   

!****************************************************************************

!subroutine: optimization condition in steady state

!****************************************************************************

subroutine fcn_ss(x, f, n) 

    implicit none
    
    !specification
    integer n
    real(dp) :: x(n), f(n)
    real(dp) :: y_ss,b_ss,bstar_ss,g_ss,ci_ss,dn_ss,yn_r_ss,yn_n_ss,yt_n_ss,yt_r_ss
    real(dp) :: varphi_l_ss,rkn_ss,rkt_ss,kn_ss,kt_ss,ln_ss,wn_ss,w_ss,wt_ss,lt_ss
    real(dp) :: rt_kt_ss  
   
    ! guess y and yn_real
    y_ss = x(1)
    yn_r_ss = x(2)

    ! steady state from calibrationss
    b_ss = sb*y_ss
    bstar_ss = b_ss/s_bar
    g_ss = govery*y_ss
    !g_r_ss = g_n_ss/pg_bar
    ci_ss = s_bar*(q_bar*bstar_ss-bstar_ss)-g_ss+y_ss
    dn_ss = varphi*(ci_ss+g_ss)       
    
    yn_n_ss = yn_r_ss*pn_bar
    yt_n_ss = y_ss - yn_n_ss
    yt_r_ss = yt_n_ss/(s_bar*tot_bar)
    varphi_l_ss = alpha_n*yn_n_ss/(alpha_n*yn_n_ss+alpha_t*yt_n_ss)         ! assign the value (calibrate) 

    ! capital    
    rkn_ss = (1.0_dp/beta-1.0_dp+delta_kn)/(1.0_dp-tau_bar)
    rkt_ss = (1.0_dp/beta-1.0_dp+delta_kt)/(1.0_dp-tau_bar)
    kn_ss = (1.0_dp-alpha_n)*yn_n_ss/rkn_ss
    kt_ss = (1.0_dp-alpha_t)*yt_n_ss/rkt_ss
                                               
    ! labor and wage
    ln_ss = (yn_r_ss/an_bar/(kn_ss**(1.0_dp-alpha_n)))**(1.0_dp/alpha_n)
    wn_ss = alpha_n*yn_n_ss/ln_ss
    w_ss = (ln_ss/L_bar/varphi_l_ss)**(-1.0_dp/chi_l)*wn_ss
    wt_ss = (alpha_t*yt_n_ss/(1.0_dp-varphi_l_ss)*(w_ss**chi_l)/L_bar)**(1.0_dp/(chi_l+1.0_dp))   
    lt_ss = alpha_t*yt_n_ss/wt_ss
       
    ! equation to solve y_ss
    F(1) = yn_r_ss - pn_bar**(-chi_c)*dn_ss !-yn_r_ss+an_bar*kn_ss**(1-alpha_n)*ln_ss**(alpha_n)
    F(2) = -w_ss + (varphi_l_ss*wn_ss**(1.0_dp+chi_l)+(1.0_dp-varphi_l_ss)*wt_ss**(1.0_dp+chi_l))**(1/(1.0_dp+chi_l))
end subroutine fcn_ss   
 
!****************************************************************

!subroutine read steady state value

!****************************************************************
 
 subroutine read_ss(y_bar,b_bar,bstar_bar,g_bar,dn_bar,yn_r_bar,yt_r_bar, &
            & rkn_bar,rkt_bar,kn_bar,kt_bar,in_bar,it_bar,i_bar,c_bar,ln_bar,lt_bar,wn_bar,wt_bar,w_bar, &
            & varphi_l,at_bar,z_bar,zovery,tb_bar,syn,lambda_bar,phi)
 
    implicit none 
    
    real(dp), intent(out):: y_bar,b_bar,bstar_bar,g_bar,dn_bar,yn_r_bar,yt_r_bar
    real(dp), intent(out):: rkn_bar,rkt_bar,kn_bar,kt_bar,in_bar,it_bar,i_bar,c_bar,ln_bar,lt_bar,wn_bar,wt_bar,w_bar
    real(dp), intent(out):: varphi_l,at_bar,z_bar,zovery,tb_bar,syn,lambda_bar,phi
    
    integer, parameter :: n = 2
    real(dp):: fnorm, x_input(n)    
    real(dp):: ci_bar,yn_n_bar,yt_n_bar,cin_bar
    
    ! solve for y_bar and yn_r_bar
    call solve_ss(x_input, fnorm)   !solve y and yn    
    y_bar = x_input(1) 
    yn_r_bar = x_input(2)
    
    ! steady state from calibrationss
    b_bar = sb*y_bar
    bstar_bar = b_bar/s_bar
    g_bar = govery*y_bar                ! g_bar is g_n_bar
    !g_r_bar = g_bar/pg_bar

    ci_bar = s_bar*(q_bar*bstar_bar-bstar_bar)-g_bar+y_bar
    dn_bar = varphi*(ci_bar+g_bar)
    yn_n_bar = yn_r_bar*pn_bar
    yt_n_bar = y_bar - yn_n_bar
    yt_r_bar = yt_n_bar/(s_bar*tot_bar)
    varphi_l = alpha_n*yn_n_bar/(alpha_n*yn_n_bar+alpha_t*yt_n_bar)         ! assign the value (calibrate) 

    ! capital    
    rkn_bar = (1.0_dp/beta-1.0_dp+delta_kn)/(1.0_dp-tau_bar)
    rkt_bar = (1.0_dp/beta-1.0_dp+delta_kt)/(1.0_dp-tau_bar)
    kn_bar = (1.0_dp-alpha_n)*yn_n_bar/rkn_bar
    kt_bar = (1.0_dp-alpha_t)*yt_n_bar/rkt_bar
    in_bar = delta_kn*kn_bar
    it_bar = delta_kt*kt_bar
    i_bar = in_bar+it_bar
            
    ! consumption
    c_bar = ci_bar-i_bar
        
    ! labor and wage
    ln_bar = (yn_r_bar/an_bar/(kn_bar**(1.0_dp-alpha_n)))**(1.0_dp/alpha_n)
    wn_bar = alpha_n*yn_n_bar/ln_bar
    w_bar = (ln_bar/L_bar/varphi_l)**(-1.0_dp/chi_l)*wn_bar
    wt_bar = (alpha_t*yt_n_bar/(1.0_dp-varphi_l)*(w_bar**chi_l)/L_bar)**(1.0_dp/(chi_l+1.0_dp))   
    lt_bar = alpha_t*yt_n_bar/wt_bar
    
    ! tradable sector
    at_bar = yt_r_bar/(kt_bar**(1.0_dp-alpha_t)*lt_bar**alpha_t)
    
    ! other variables
    z_bar = tau_bar*y_bar+s_bar*(beta-1.0_dp)*bstar_bar-g_bar   
    zovery = z_bar/y_bar
    tb_bar=y_bar-c_bar-i_bar-g_bar   
    syn = yn_n_bar/y_bar    
    
    !ctemp_bar =  c_bar**((v-1.0_dp)/v)
    !gtemp_bar =  g_r_bar**((v-1.0_dp)/v)
    !ctil_bar =  (omega*ctemp_bar+(1.0_dp-omega)*gtemp_bar)**(v/(v-1))    ! composed consumption basket
    lambda_bar = c_bar**(-sig_c) !omega*c_bar**(-1.0_dp/v)*ctil_bar**(1.0_dp/v-sig_c)
    if (index_logl .eq. 1) then
       phi = lambda_bar*(1.0_dp-tau_bar)*w_bar*(1.0_dp-L_bar)
    elseif (index_logl .eq. 2) then
       phi = lambda_bar*(1.0_dp-tau_bar)*w_bar*L_bar**(-sig_l)
    end if

end subroutine read_ss


!************************************************************************
!
!  SUBROUTINE: construct grid
!
!************************************************************************
subroutine construct_grid(tauv_grid,av_grid,gv_grid,totv_grid,knv_grid,ktv_grid,ua_grid,ug_grid,utot_grid,pr_ua,pr_ug,pr_utot,pnstar_rule,lnstar_rule,knstar_rule)

	implicit none
    
    real(dp), intent(out) :: tauv_grid(ntauv),av_grid(nav),gv_grid(ngv),totv_grid(ntotv),knv_grid(nknv),ktv_grid(nktv)
    real(dp), intent(out) :: ua_grid(nua),ug_grid(nug),utot_grid(nutot)
    real(dp), intent(out) :: pr_ug(nug),pr_ua(nua),pr_utot(nutot)
    real(dp), dimension(nav,ngv,ntauv,ntotv,nknv,nktv), intent(out) :: pnstar_rule,lnstar_rule,knstar_rule
	
	real(dp) :: x1, x2, sqrt_2
	integer  :: i, j, k	
    real(dp) :: a_min, a_max, a_step, g_min, g_max, g_step, tau_min, tau_max, tau_step
    real(dp) :: tot_min, tot_max, tot_step, kn_min, kn_max, kn_step, kt_min, kt_max, kt_step

!-------------------------------------------------------------------------------
! state variable range: a,g,z,tot,kn,kt
!-------------------------------------------------------------------------------
    open(unit = 8, file = "./dat_files/tauv_grid.dat", status = "old", action = "read", iostat = openstatus)
    read(8,*) tauv_grid
    close(8)	
    
    a_min       = -a_std_max/dsqrt(1.0_dp-rho_a**2.0_dp)*std_ua  
    a_max       = a_std_max/dsqrt(1.0_dp-rho_a**2.0_dp)*std_ua
    a_step      = (a_max-a_min)/(nav-1.0_dp)
	
    g_min       = -g_std_max/dsqrt(1.0_dp-rho_g**2.0_dp)*std_ug
    g_max       = g_std_max/dsqrt(1.0_dp-rho_g**2.0_dp)*std_ug
    g_step      = (g_max-g_min)/(ngv-1.0_dp)

    !tau_min     = -tau_std_max/dsqrt(1.0_dp-rho_tau**2.0_dp)*std_utau  
    !tau_max     = tau_std_max/dsqrt(1.0_dp-rho_tau**2.0_dp)*std_utau
    !tau_step    = (tau_max-tau_min)/(ntauv-1.0_dp)   
    
    tot_min     = -tot_std_max/dsqrt(1.0_dp-rho_tot**2.0_dp)*std_utot  
    tot_max     = tot_std_max/dsqrt(1.0_dp-rho_tot**2.0_dp)*std_utot
    tot_step    = (tot_max-tot_min)/(ntotv-1.0_dp)
        
    kn_min      = kn_bar*(1.0_dp-kn_margin_l)
    kn_max      = kn_bar*(1.0_dp+kn_margin_h)
    kn_step     = (kn_max-kn_min)/(nknv-1.0_dp)
    
    kt_min      = kt_bar*(1.0_dp-kt_margin_l)
    kt_max      = kt_bar*(1.0_dp+kt_margin_h)
    kt_step     = (kt_max-kt_min)/(nktv-1.0_dp)
    
!-------------------------------------------------------------------------------
! construct state varialbe grids: a,g,z,tot,kn,kt,ua,ug,uz,utot
!-------------------------------------------------------------------------------
	do i = 1,nav
		av_grid(i) = a_min+(i-1)*a_step
	end do
	do i = 1,ngv
		gv_grid(i) = g_min+(i-1)*g_step
    end do
 !   do i = 1,ntauv
	!	tauv_grid(i) = tau_min+(i-1)*tau_step
	!end do	    
	do i = 1,ntotv
		totv_grid(i) = tot_min+(i-1)*tot_step
	end do	
	do i = 1,nknv
		knv_grid(i) = kn_min+(i-1)*kn_step
	end do	
	do i = 1,nktv
		ktv_grid(i) = kt_min+(i-1)*kt_step
	end do	
	
	do i = 1,nua
		ua_grid(i) = ua_min+(i-1)*ua_step
	end do	
	do i = 1,nug
		ug_grid(i) = ug_min+(i-1)*ug_step
    end do
 !   do i = 1,nutau
	!	utau_grid(i) = utau_min+(i-1)*utau_step
	!end do    
	do i = 1,nutot
		utot_grid(i) = utot_min+(i-1)*utot_step
	end do

!-------------------------------------------------------------------------------
! compute transition prob matrix: ua,ug,uz,utot
!-------------------------------------------------------------------------------
	do i = 1,nua
	    pr_ua(i)=(2.0_dp*pi)**(-0.5_dp)/std_ua*exp(-(ua_grid(i)**2.0_dp)/(std_ua**2.0_dp)/2.0_dp)
	end do
	do i = 1,nug
	    pr_ug(i)=(2.0_dp*pi)**(-0.5_dp)/std_ug*exp(-(ug_grid(i)**2.0_dp)/(std_ug**2.0_dp)/2.0_dp)
    end do
 !   do i = 1,nutau
	!    pr_utau(i)=(2.0_dp*pi)**(-0.5_dp)/std_utau*exp(-(utau_grid(i)**2.0_dp)/(std_utau**2.0_dp)/2.0_dp)
	!end do    
	do i = 1,nutot
	    pr_utot(i)=(2.0_dp*pi)**(-0.5_dp)/std_utot*exp(-(utot_grid(i)**2.0_dp)/(std_utot**2.0_dp)/2.0_dp)
	end do

!-------------------------------------------------------------------------------
! save grids
!-------------------------------------------------------------------------------
	open(unit = 1, file = "./dat_files/av_grid.dat", status = "replace", action = "write", iostat = openstatus)
        write(1,'(F19.16)') av_grid
      	close(1)   
	open(unit = 1, file = "./dat_files/gv_grid.dat", status = "replace", action = "write", iostat = openstatus)
        write(1,'(F19.16)') gv_grid
      	close(1)
    !open(unit = 1, file = "./dat_files/tauv_grid.dat", status = "replace", action = "write", iostat = openstatus)
    !    write(1,'(F19.16)') tauv_grid
    !  	close(1)        
    open(unit = 1, file = "./dat_files/totv_grid.dat", status = "replace", action = "write", iostat = openstatus)
        write(1,'(F19.16)') totv_grid
      	close(1)
	open(unit = 1, file = "./dat_files/knv_grid.dat", status = "replace", action = "write", iostat = openstatus)
        write(1,'(F19.16)') knv_grid
      	close(1)
	open(unit = 1, file = "./dat_files/ktv_grid.dat", status = "replace", action = "write", iostat = openstatus)
        write(1,'(F19.16)') ktv_grid
      	close(1)
	
	if (index_load==1) then	
		open(unit = 8, file = "./dat_files/pnstar_rule.dat", status = "old", action = "read", iostat = openstatus)
		read(8,*) pnstar_rule
		close(8)
		open(unit = 8, file = "./dat_files/lnstar_rule.dat", status = "old", action = "read", iostat = openstatus)
		read(8,*) lnstar_rule
		close(8)
		open(unit = 8, file = "./dat_files/knstar_rule.dat", status = "old", action = "read", iostat = openstatus)
		read(8,*) knstar_rule
		close(8)
	elseif (index_load==2) then                  
		pnstar_rule = pn_bar
        lnstar_rule = ln_bar
        knstar_rule = kn_bar
	elseif (index_load == 3) then 
		! read in initial files
		open(unit = 8, file = "./dat_files/ini/pnstar_rule.dat", status = "old", action = "read", iostat = openstatus)
		read(8,*) pnstar_rule_ini
		close(8)
		open(unit = 8, file = "./dat_files/ini/lnstar_rule.dat", status = "old", action = "read", iostat = openstatus)
		read(8,*) lnstar_rule_ini
		close(8)
		open(unit = 8, file = "./dat_files/ini/knstar_rule.dat", status = "old", action = "read", iostat = openstatus)
		read(8,*) knstar_rule_ini
		close(8)
		open(unit = 8, file = "./dat_files/ini/av_grid.dat", status = "old", action = "read", iostat = openstatus)
		read(8,*) av_grid_ini
		close(8)
		open(unit = 8, file = "./dat_files/ini/gv_grid.dat", status = "old", action = "read", iostat = openstatus)
		read(8,*) gv_grid_ini
		close(8)    
		open(unit = 8, file = "./dat_files/ini/tauv_grid.dat", status = "old", action = "read", iostat = openstatus)
		read(8,*) tauv_grid_ini
		close(8)
		open(unit = 8, file = "./dat_files/ini/totv_grid.dat", status = "old", action = "read", iostat = openstatus)
		read(8,*) totv_grid_ini
		close(8)        
		open(unit = 8, file = "./dat_files/ini/knv_grid.dat", status = "old", action = "read", iostat = openstatus)
		read(8,*) knv_grid_ini
		close(8)
		open(unit = 8, file = "./dat_files/ini/ktv_grid.dat", status = "old", action = "read", iostat = openstatus)
		read(8,*) ktv_grid_ini
		close(8)
        do ia = 1,nav
        do ig = 1,ngv
		do itau = 1,ntauv
        do itot = 1,ntotv
		do ikn = 1,nknv
		do ikt = 1,nktv
			pnstar_rule(ia,ig,itau,itot,ikn,ikt) = interp6(av_grid_ini,nav_ini,gv_grid_ini,ngv_ini,tauv_grid_ini,ntauv,totv_grid_ini,ntotv_ini,knv_grid_ini,nknv_ini,ktv_grid_ini,nktv_ini,&
					& av_grid(ia),gv_grid(ig),tauv_grid(itau),totv_grid(itot),knv_grid(ikn),ktv_grid(ikt),pnstar_rule_ini)
			lnstar_rule(ia,ig,itau,itot,ikn,ikt) = interp6(av_grid_ini,nav_ini,gv_grid_ini,ngv_ini,tauv_grid_ini,ntauv,totv_grid_ini,ntotv_ini,knv_grid_ini,nknv_ini,ktv_grid_ini,nktv_ini,&
					& av_grid(ia),gv_grid(ig),tauv_grid(itau),totv_grid(itot),knv_grid(ikn),ktv_grid(ikt),lnstar_rule_ini)
			knstar_rule(ia,ig,itau,itot,ikn,ikt) = interp6(av_grid_ini,nav_ini,gv_grid_ini,ngv_ini,tauv_grid_ini,ntauv,totv_grid_ini,ntotv_ini,knv_grid_ini,nknv_ini,ktv_grid_ini,nktv_ini,&
					& av_grid(ia),gv_grid(ig),tauv_grid(itau),totv_grid(itot),knv_grid(ikn),ktv_grid(ikt),knstar_rule_ini)        
		end do
        end do
        end do
        end do
        end do
        end do
        
    elseif (index_load == 4) then    
        ! read in initial files
        open(unit = 8, file = "./dat_files/pnstar_rule_nog.dat", status = "old", action = "read", iostat = openstatus)
        read(8,*) pnstar_rule_nog
        close(8)
        open(unit = 8, file = "./dat_files/lnstar_rule_nog.dat", status = "old", action = "read", iostat = openstatus)
        read(8,*) lnstar_rule_nog
        close(8)
        open(unit = 8, file = "./dat_files/knstar_rule_nog.dat", status = "old", action = "read", iostat = openstatus)
        read(8,*) knstar_rule_nog
        close(8)
        
        do ig = 1,ngv
            pnstar_rule(:,ig,:,:,:,:) = pnstar_rule_nog(:,:,:,:,:)
            lnstar_rule(:,ig,:,:,:,:) = lnstar_rule_nog(:,:,:,:,:)
            knstar_rule(:,ig,:,:,:,:) = knstar_rule_nog(:,:,:,:,:)
        end do

	end if  
		

            end subroutine construct_grid

!************************************************************************
!
!  SUBROUTINE: load data
!
!************************************************************************

subroutine load_data(pr_tau, cdf_tau, tau_random)
    !use mod_interp_locate
    implicit none
    real(dp), intent(out) :: pr_tau(ntauv), cdf_tau(ntauv), tau_random(per_bstar_sim)

    !! save v_grid.dat 
    !open(unit = 8, file = "./dat_files/tauv_grid.dat", status = "old", action = "read", iostat = openstatus)
    !read(8,*) tauv_grid
    !close(8)	
    ! save pr_v
    open(unit = 8, file = "./dat_files/pr_tau.dat", status = "old", action = "read", iostat = openstatus)
    read(8,*) pr_tau
    close(8)
    ! save cdf_v
    open(unit = 8, file = "./dat_files/cdf_tau.dat", status = "old", action = "read", iostat = openstatus)
    read(8,*) cdf_tau
    close(8)
    open(unit = 8, file = "./dat_files/tau_random.dat", status = "old", action = "read", iostat = openstatus)
    read(8,*) tau_random
    close(8)
    
end subroutine load_data    
!!************************************************************************
!!
!!  SUBROUTINE: load steady state
!!
!!************************************************************************
!
!subroutine load_ss(y_bar,b_bar,bstar_bar,g_bar,dn_bar,yn_r_bar,yt_r_bar, &
!            & rkn_bar,rkt_bar,kn_bar,kt_bar,in_bar,it_bar,i_bar,c_bar,ln_bar,lt_bar,wn_bar,wt_bar,w_bar, &
!            & varphi_l,at_bar,z_bar,zovery,tb_bar,syn,lambda_bar,phi)
!
!    implicit none
!
!    real(dp) :: steady_state_mat(28)
!    real(dp), intent(out) :: y_bar,b_bar,bstar_bar,g_bar,dn_bar,yn_r_bar,yt_r_bar,rkn_bar,rkt_bar,kn_bar,kt_bar,in_bar,it_bar,i_bar,c_bar,ln_bar,lt_bar,wn_bar,wt_bar,w_bar
!    real(dp), intent(out) :: varphi_l,at_bar,z_bar,zovery,tb_bar,syn,lambda_bar,phi
!
!    open(unit = 8, file = "./dat_files/steady_state_mat.dat", status = "old", action = "read", iostat = openstatus)
!    read(8,*) steady_state_mat
!    close(8)
!    
!    y_bar = steady_state_mat(1)
!    b_bar = steady_state_mat(2)
!    bstar_bar = steady_state_mat(3)
!    g_bar = steady_state_mat(4)
!    dn_bar = steady_state_mat(5)
!    yn_r_bar = steady_state_mat(6)
!    yt_r_bar = steady_state_mat(7)
!    rkn_bar = steady_state_mat(8)
!    rkt_bar = steady_state_mat(9)
!    kn_bar = steady_state_mat(10)
!    kt_bar = steady_state_mat(11)
!    in_bar = steady_state_mat(12)
!    it_bar = steady_state_mat(13)
!    i_bar = steady_state_mat(14)
!    c_bar = steady_state_mat(15)
!    ln_bar = steady_state_mat(16)
!    lt_bar = steady_state_mat(17)
!    wn_bar = steady_state_mat(18)
!    wt_bar = steady_state_mat(19)
!    w_bar = steady_state_mat(20)
!    varphi_l = steady_state_mat(21)
!    at_bar = steady_state_mat(22)
!    z_bar = steady_state_mat(23)
!    zovery = steady_state_mat(24)
!    tb_bar = steady_state_mat(25)
!    syn = steady_state_mat(26)
!    lambda_bar = steady_state_mat(27)  
!    phi = steady_state_mat(28) 
!    
!    
!    end subroutine load_ss

END MODULE mod_parameter_fl

