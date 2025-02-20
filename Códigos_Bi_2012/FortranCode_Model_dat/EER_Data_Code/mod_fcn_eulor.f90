!****************************************************************************
!
!  Modules: fcn_eulor
!
!****************************************************************************

MODULE mod_fcn_eulor

use mod_global_constant
use mod_parameter
use mod_interp_locate
IMPLICIT NONE

!****************************************************************************
CONTAINS

!****************************************************************************
!
!  FUNCTION: fcn_eulor
!       solve the nonlinear system (1-equation)
!
!****************************************************************************

function fcn_eulor(vector_eulor,index_eulor)
    use mod_interp_locate
	implicit none
	
	! specification
	real(dp) :: fcn_eulor(nfcn)
	real(dp), intent(in) :: vector_eulor(nfcn)
	integer, intent(in) :: index_eulor(nstate)
	
	integer :: ib_eulor, iA_eulor, ig_eulor, iz1_eulor, izrs_eulor, izrsp_eulor, iuAp_eulor
	integer :: iuGp_eulor, ivp_eulor, idelp_eulor
	real(dp) :: z1_eulor, z_eulor, g_eulor, A_eulor, b_d_eulor, gama_eulor, bp_eulor, tauL_eulor
	real(dp) :: c_eulor, q_eulor, lhs_eulor, gama_p_eulor, bp_d_eulor, gp_eulor, Ap_eulor, vp_eulor, delp_eulor
	real(dp) :: tauLp_eulor, cp_eulor, E_rhs_eulor
	real(dp) :: E_rhs_nzrs(nzrs), E_rhs_nuA(nuA), E_rhs_nuG(nuG), E_rhs_nv(nv),cp_delp_ndelta_eulor(ndelta)
	real(dp) :: prob_d_zrsp1_eulor,prob_d_zrsp2_eulor,tauLp_d_eulor,cp_d_eulor,cp_delp_eulor        
	
	! read in 
	bp_eulor = vector_eulor(1)          ! b_t
	ib_eulor = index_eulor(1)           !index_eulor(1)
	iA_eulor = index_eulor(2)           !index_eulor(2)
	ig_eulor = index_eulor(3)
	iz1_eulor = index_eulor(4)
	izrs_eulor = index_eulor(5)          

	! exog and endog states
	z1_eulor = z_grid(iz1_eulor)
	A_eulor = A_grid(iA_eulor)
	b_d_eulor = b_grid(ib_eulor)        ! b^d_t (post-default govt liability at t)
	gama_eulor = gama_bar
	
	! allow middel steps to facilitate convergence
	if (load_index==2) then
		g_eulor = g_bar  
	else
	    g_eulor = g_grid(ig_eulor)
	end if
	
	if (izrs_eulor==1) then
	    z_eulor = z_bar+zeta_z*(A_eulor-A_bar)
	else
	    z_eulor = mu_z*z1_eulor+zeta_z*(A_eulor-A_bar) 
	end if	

	! compute y,c,T at time t
	tauL_eulor = tauL_bar+gama_eulor*(b_d_eulor-b_bar)
	c_eulor = (1.0_dp-tauL_eulor)*(A_eulor-g_eulor)/(1.0_dp-tauL_eulor+chi_N)
	q_eulor = (g_eulor+b_d_eulor+z_eulor-tauL_eulor*(c_eulor+g_eulor))/bp_eulor
	lhs_eulor = q_eulor

	
	! Numerical integration
	gama_p_eulor = gama_bar	
	
	E_rhs_nuA = 0.0_dp
    do iuAp_eulor = 1,nuA
    Ap_eulor = (A_bar**(1.0_dp-rho_A))*(A_eulor**rho_A)*exp(uA_grid(iuAp_eulor))
    
    E_rhs_nuG = 0.0_dp
    do iuGp_eulor = 1,nuG
    if (load_index==2) then
		gp_eulor = g_bar  
	else
	    gp_eulor = (g_bar**(1.0_dp-rho_g))*(g_eulor**rho_g)*exp(uG_grid(iuGp_eulor))
	end if       
    
    ! compute using the cdf function
    ! default probability         
    prob_d_zrsp1_eulor = interp3(v_grid,nv,vA_grid,nvA,vg_grid,nvg,bp_eulor,Ap_eulor,gp_eulor,cdf_v(:,:,:,1))
    prob_d_zrsp2_eulor = interp3(v_grid,nv,vA_grid,nvA,vg_grid,nvg,bp_eulor,Ap_eulor,gp_eulor,cdf_v(:,:,:,2))
    
    ! using delta_grid and cdf_delta (change the decision rule at high level of debt)
    cp_delp_ndelta_eulor = 0.0_dp
    do idelp_eulor = 1,ndelta
        delp_eulor = delta_grid(idelp_eulor)
        bp_d_eulor = bp_eulor*(1.0_dp-delp_eulor)
        tauLp_d_eulor = tauL_bar+gama_p_eulor*(bp_d_eulor-b_bar)
        cp_d_eulor = (1.0_dp-tauLp_d_eulor)/(1.0_dp-tauLp_d_eulor+chi_N)*(Ap_eulor-gp_eulor)
        cp_delp_ndelta_eulor(idelp_eulor)=(1.0_dp-delp_eulor)/cp_d_eulor*pr_delta(idelp_eulor)
    end do
    cp_delp_eulor = sum(cp_delp_ndelta_eulor)        
    
    ! if not default next period
    tauLp_eulor = tauL_bar+gama_p_eulor*(bp_eulor-b_bar)
    cp_eulor = (1.0_dp-tauLp_eulor)/(1.0_dp-tauLp_eulor+chi_N)*(Ap_eulor-gp_eulor)
        
    ! expectation
    E_rhs_nuG(iuGp_eulor) = beta*c_eulor*(cp_delp_eulor*(prob_d_zrsp1_eulor*pr_zrs(izrs_eulor,1) & 
                        & +prob_d_zrsp2_eulor*pr_zrs(izrs_eulor,2)) &
                        & + 1.0_dp/cp_eulor*((1.0_dp-prob_d_zrsp1_eulor)*pr_zrs(izrs_eulor,1) & 
                        & +(1.0_dp-prob_d_zrsp2_eulor)*pr_zrs(izrs_eulor,2)))*pr_uA(iuAp_eulor)*pr_uG(iuGp_eulor)                         

    end do  ! for iuGp
    
    E_rhs_nuA(iuAp_eulor) = (2.0_dp*sum(E_rhs_nuG)-E_rhs_nuG(1)-E_rhs_nuG(nuG))*uG_step/2.0_dp
    end do  ! for iuAp
	
	E_rhs_eulor = (2.0_dp*sum(E_rhs_nuA)-E_rhs_nuA(1)-E_rhs_nuA(nuA))*uA_step/2.0_dp	
	
	! 1-equation
	fcn_eulor(1) = lhs_eulor - E_rhs_eulor
	
    
end function fcn_eulor


END MODULE mod_fcn_eulor

