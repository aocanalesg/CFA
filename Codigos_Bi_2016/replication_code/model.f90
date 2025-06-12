!****************************************************************************
!  PROGRAM: model.f90
!  
! 07/03/2012
! need to adjust bgrid when change the fiscal limit distribution (vgrid)
! increase ui_std_max from 4 to 5.0_dp (to increase precision)
!
! 07/03/12
! change fcn_eulor and solve_model for unconditional fiscal limit
! need to lower haircut_bar (b/c steeper distribution)
! raise tol doesn't change precision much
! increasing grid points or reducing variance can increase precision 

! 07/21/12
! add capital, add two state variables: kn and kt, one choice variable:investment
! at this stage, assume capital adjustment cost is 0, v=0
! in this code: state variable: b, an, kn, kt, #grids:7,7,7,7
! add eulor equation to fcn_eulor
!
! 07/24/12
! add capital, add two state variables: kn_{t-1} and kt_{t-1}, two choice variable: kn_t, kt_t
! can only achieve convergence when std(an) is small and grids for kn and kt are very narrow, otherwise, randomize activated

! 08/08/12
! change model to FDI, state variable (b_{t-1},an_t,g_t,kn_{t-1}), choice variable (pn_t,ln_t,b_t)
! 08/14/12
! add state variable (g_t)

! 09/4/12
! add capital adjustment cost

! 10/8/12
! change shocks to g and tot
! 
! 10/26/12
! change the fiscal limit in terms of local currency
! see mod_fcn_eulor and mod_solve_model
!
! b_t-1 * s_t-1 >=< bstar_t
!
! 11/14/12
! based on the model of Nov_14 version: g in utility, procyclical g and z, different calibration
! all shocks are consistent with Argentina estimation, except sigma_tot
! new (ad-hoc) fiscal limit distribution: (0.45 - 0.7) of ybar
!
! 11/19/12
! use long-term gamma_tau
!
! 11 /21/12
! change v = 0.4, rho_g from 0.32 to 0.35, sigma_a from 0.88 to 0.83 (unused) --- re-estimation results
!
! 11/28/12
! change g_std_max and ng, ug_std_max and nuG to make ggrid finer
! 
! 12/2/12
! no tot grid
!
! 3/12/13
! change to CRRA utility
!
! 3/13/13
! add sector specific capital, no default
!
! 3/25/13
! fiscal limit, a_g_tot shock
!****************************************************************************
program model

use mod_global_constant
use mod_nonlinear_solver
use mod_parameter_fl
use mod_bstar_sim
use mod_interp_locate

!** include 'link_fnl_shared.h' 

use neqnf_int 
use umach_int

implicit none 

real(dp) :: steady_state_mat(28)
integer, parameter :: index_laffer = 1              ! solve for pnstar, lnstar, knstar rules
integer, parameter :: index_compute_Tmax = 1        ! solve for Tmax rule
integer, parameter :: index_bstar = 1               ! whether or not compute bstar 


!****************************************************************************
! compute steady state
!****************************************************************************
    
call read_ss(y_bar,b_bar,bstar_bar,g_bar,dn_bar,yn_r_bar,yt_r_bar, &
            & rkn_bar,rkt_bar,kn_bar,kt_bar,in_bar,it_bar,i_bar,c_bar,ln_bar,lt_bar,wn_bar,wt_bar,w_bar, &
            & varphi_l,at_bar,z_bar,zovery,tb_bar,syn,lambda_bar,phi)
 
!print *, y_bar,b_bar,bstar_bar,g_bar,dn_bar,yn_r_bar,yt_r_bar
!print *, rkn_bar,rkt_bar,kn_bar,kt_bar,in_bar,it_bar,i_bar,c_bar,ln_bar,lt_bar,wn_bar,wt_bar,w_bar
!print *, varphi_l,at_bar,z_bar,zovery,tb_bar,syn,lambda_bar,L_bar,phi

steady_state_mat = (/y_bar,b_bar,bstar_bar,g_bar,dn_bar,yn_r_bar,yt_r_bar, &
            & rkn_bar,rkt_bar,kn_bar,kt_bar,in_bar,it_bar,i_bar,c_bar,ln_bar,lt_bar,wn_bar,wt_bar,w_bar, &
            & varphi_l,at_bar,z_bar,zovery,tb_bar,syn,lambda_bar,phi/)
open(unit = 1, file = "./dat_files/steady_state_mat.dat", status = "replace", action = "write", iostat = openstatus)
    write(1,'(F24.16)') steady_state_mat
    close(1)
    
!call load_ss(y_bar,b_bar,bstar_bar,g_bar,g_r_bar,dn_bar,yn_r_bar,yt_r_bar, &
!            & rkn_bar,rkt_bar,kn_bar,kt_bar,in_bar,it_bar,i_bar,c_bar,ln_bar,lt_bar,wn_bar,wt_bar,w_bar, &
!            & varphi_l,at_bar,z_bar,zovery,tb_bar,syn,ctil_bar,lambda_bar,phi)

print *, y_bar
print *, kn_bar
print *, kt_bar
print *, bstar_bar/y_bar/4.0_dp

!*****************************************************************************
! initialization for grids
call construct_grid(tauv_grid,av_grid,gv_grid,totv_grid,knv_grid,ktv_grid,ua_grid,ug_grid,utot_grid,pr_ua,pr_ug,pr_utot,pnstar_rule,lnstar_rule,knstar_rule)

call load_data(pr_tau, cdf_tau, tau_random)


! use nonlinear solver to solve the model
if (index_laffer==1) then
    ! initialization
    pnstar_max_diff = 100.0_dp
    lnstar_max_diff = 100.0_dp
    knstar_max_diff = 100.0_dp  

    pnstar_new = pnstar_rule
    lnstar_new = lnstar_rule
    knstar_new = knstar_rule

    iter = 1
    print *, '-----------------------------------------' 
    print *, 'start solving pnstar_rule, lnstar_rule, knstar_rule'
    
    ! to compute the decision rule at tau_max
    do
        ! exit if it converges
	    if ((pnstar_max_diff < tol_conv).AND.(lnstar_max_diff<tol_conv).AND.(knstar_max_diff<tol_conv)) exit	    	        
	    
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(pnstar_new,lnstar_new,knstar_new,pnstar_rule,lnstar_rule,knstar_rule,pnstar_diff_iter,lnstar_diff_iter,knstar_diff_iter)
        !$OMP DO
        do ia = 1,nav 
        do ig = 1,ngv 
        do itau = 1,ntauv
        do itot = 1,ntotv   
        do ikn = 1,nknv 
        do ikt = 1,nktv
		    ! call nonlinear solver				
		    index_t_star = (/ia,ig,itau,itot,ikn,ikt/)
		    vector_t_star = (/pnstar_rule(ia,ig,itau,itot,ikn,ikt),lnstar_rule(ia,ig,itau,itot,ikn,ikt),knstar_rule(ia,ig,itau,itot,ikn,ikt)/)			
		    call csolve(vector_t_star,index_t_star, nfcn_star, tol_nsolve,nmax_nsolve,rc_nsolve,nstate_star)
		    if (rc_nsolve>1) then
		        print *, '|---- rc_nsolve ----|', rc_nsolve
		        write (*,*), '|---- ia,ig,itau,itot,ikn,ikt ----|', ia,ig,itau,itot,ikn,ikt
		    end if	
		    pnstar_new(ia,ig,itau,itot,ikn,ikt) = vector_t_star(1)	
		    lnstar_new(ia,ig,itau,itot,ikn,ikt) = vector_t_star(2) 
		    knstar_new(ia,ig,itau,itot,ikn,ikt) = vector_t_star(3)  

            ! update the rules
		    pnstar_diff_iter(ia,ig,itau,itot,ikn,ikt) = dabs(pnstar_new(ia,ig,itau,itot,ikn,ikt)-pnstar_rule(ia,ig,itau,itot,ikn,ikt))
		    lnstar_diff_iter(ia,ig,itau,itot,ikn,ikt) = dabs(lnstar_new(ia,ig,itau,itot,ikn,ikt)-lnstar_rule(ia,ig,itau,itot,ikn,ikt))
		    knstar_diff_iter(ia,ig,itau,itot,ikn,ikt) = dabs(knstar_new(ia,ig,itau,itot,ikn,ikt)-knstar_rule(ia,ig,itau,itot,ikn,ikt))
            
	    end do
	    end do
	    end do
	    end do
        end do
        end do
	    !$OMP END DO
	    !$OMP END PARALLEL
	    
	        pnstar_max_diff = maxval(pnstar_diff_iter)
            lnstar_max_diff = maxval(lnstar_diff_iter)
            knstar_max_diff = maxval(knstar_diff_iter)


	    pnstar_rule = pnstar_new
	    lnstar_rule = lnstar_new
	    knstar_rule = knstar_new
	    write(*,*), iter,pnstar_max_diff,lnstar_max_diff,knstar_max_diff
	    iter = iter +1
	    
	    ! elapsed time
	    call itime(time_now)
        write (*,10), time_now
        10 format ( ' time ', i2.2, ':', i2.2, ':', i2.2 )        	
    
    ! save the rules
    open(unit = 1, file = "./dat_files/pnstar_rule.dat", status = "replace", action = "write", iostat = openstatus)
        write(1,'(F24.16)') pnstar_rule
      	close(1)
    open(unit = 1, file = "./dat_files/lnstar_rule.dat", status = "replace", action = "write", iostat = openstatus)
        write(1,'(F24.16)') lnstar_rule
      	close(1)
    open(unit = 1, file = "./dat_files/knstar_rule.dat", status = "replace", action = "write", iostat = openstatus)
        write(1,'(F24.16)') knstar_rule
      	close(1)
        
    end do
    
else    
    !load results 
    print *, 'load pnstar_rule...'
    open(unit = 8, file = "./dat_files/pnstar_rule.dat", status = "old", action = "read", iostat = openstatus)
    read(8,*) pnstar_rule
    close(8)    
    print *, 'load lnstar_rule...'
    open(unit = 8, file = "./dat_files/lnstar_rule.dat", status = "old", action = "read", iostat = openstatus)
    read(8,*) lnstar_rule
    close(8)
    print *, 'load knstar_rule...'
    open(unit = 8, file = "./dat_files/knstar_rule.dat", status = "old", action = "read", iostat = openstatus)
    read(8,*) knstar_rule
    close(8)     
end if 


! to compute T_max for the solved decision rules
if (index_compute_Tmax == 1) then
    print *, '-----------------------------------------' 
    print *, 'start computing Tmax ...'

	    open(unit = 8, file = "./dat_files/pnstar_rule.dat", status = "old", action = "read", iostat = openstatus)
	    read(8,*) pnstar_rule
	    close(8)
	    open(unit = 8, file = "./dat_files/lnstar_rule.dat", status = "old", action = "read", iostat = openstatus)
	    read(8,*) lnstar_rule
	    close(8)
	    open(unit = 8, file = "./dat_files/knstar_rule.dat", status = "old", action = "read", iostat = openstatus)
	    read(8,*) knstar_rule
	    close(8)
    
    do ia = 1,nav
    do ig = 1,ngv 
    do itau = 1,ntauv
    do itot = 1,ntotv 
    do ikn = 1,nknv	
    do ikt = 1,nktv
	    ! call nonlinear solver				
	    index_t_star = (/ia,ig,itau,itot,ikn,ikt/)
	    vector_t_star = (/pnstar_rule(ia,ig,itau,itot,ikn,ikt),lnstar_rule(ia,ig,itau,itot,ikn,ikt),knstar_rule(ia,ig,itau,itot,ikn,ikt)/)	
	    call compute_Tmax(vector_out_t_star, vector_t_star,index_t_star)	
	    rules_star(ia,ig,itau,itot,ikn,ikt,1) = vector_out_t_star(1)		
        rules_star(ia,ig,itau,itot,ikn,ikt,2) = vector_out_t_star(2)  
        rules_star(ia,ig,itau,itot,ikn,ikt,3) = vector_out_t_star(3)  
        rules_star(ia,ig,itau,itot,ikn,ikt,4) = vector_out_t_star(4)
    end do
    end do
    end do  
    end do
    end do 
    end do
    Tstar_rule = rules_star(:,:,:,:,:,:,1)
    sstar_rule = rules_star(:,:,:,:,:,:,2)
    ktstar_rule = rules_star(:,:,:,:,:,:,3)
    ystar_rule = rules_star(:,:,:,:,:,:,4)
    
    open(unit = 1, file = "./dat_files/rules_star.dat", status = "replace", action = "write", iostat = openstatus)
        write(1,'(F24.16)') rules_star
      	close(1)
      	
end if
         
! compute the fiscal limit distribution
if (index_bstar==1) then
    print *, '-----------------------------------------' 
    print *, 'start simulating bstar_dist'    
    call compute_dist(bstar_mat,sstar_mat)   
    open(unit = 1, file = "./dat_files/bstar_mat.dat", status = "replace", action = "write", iostat = openstatus)
        write(1,'(F24.16)') bstar_mat
      	close(1)
    !open(unit = 1, file = "./dat_files/sstar_mat.dat", status = "replace", action = "write", iostat = openstatus)
    !    write(1,'(F24.16)') sstar_mat
    !  	close(1)
        
end if


end program model

