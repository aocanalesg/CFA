!****************************************************************************
!
!  Modules: assign parameter values
!
!****************************************************************************

MODULE mod_parameter

use mod_global_constant
use mod_interp_locate
IMPLICIT NONE

!-------------------------------------------------------------------------------
! Model Parameters
!-------------------------------------------------------------------------------

!------------------------------------------------------------------ 	
! other parameter
real(dp), parameter :: tol_conv     = (0.1_dp)**6               ! conv crit for the whole loop (all state var)
real(dp), parameter :: tol_nsolve   = (0.1_dp)**6               ! con crit for nonlinear solver
integer, parameter  :: nv           = 160  
integer, parameter	:: nb 	        = 150                       ! grid points
integer, parameter  :: nb_ini       = 120
integer, parameter  :: load_index   = 3                         ! 3: load n dim; 
                                                                ! 4: load n dim, and interpolate over bgrid;

!--------------------------------------------------------------------------	
! deep parameters
real(dp), parameter :: delta_bar    = 0.0947_dp		        ! default rate
real(dp), parameter :: gama_bar	    = 0.37_dp               !0.42_dp
real(dp), parameter :: mu_z         = 1.015                 ! transfer growth
real(dp), parameter :: zeta_z_log   = -0.45_dp
real(dp), parameter :: beta 	    = 0.95_dp	            ! discount rate
real(dp), parameter :: pi           = 3.14159265358979_dp   ! define pi

! steady-state parameters
real(dp), parameter :: N_bar	    = 0.25_dp	    ! steady-state labor
real(dp), parameter :: tauL_bar 	= 0.32_dp	    ! ss tax rate
real(dp), parameter :: govery	    = 0.167_dp	    ! ss g/y
real(dp), parameter :: zovery	    = 0.1334_dp	    ! ss z/y
real(dp), parameter :: A_bar	    = 1.0_dp	    ! ss productivity

real(dp), parameter :: covery 	    = 1.0_dp-govery	! ss c/y
real(dp), parameter :: y_bar	    = N_bar*A_bar	! ss output
real(dp), parameter :: c_bar	    = y_bar*covery	! ss consumption
real(dp), parameter :: g_bar	    = y_bar*govery
real(dp), parameter :: z_bar 	    = y_bar*zovery
real(dp), parameter :: q_bar	    = beta	        ! ss bond price
real(dp), parameter :: bovery 	    = (tauL_bar-govery-zovery)/(1.0_dp-q_bar)
real(dp), parameter :: b_bar	    = bovery*y_bar
real(dp), parameter :: T_bar	    = tauL_bar*y_bar
real(dp), parameter :: chi_N	    = (1-tauL_bar)*A_bar*(1-N_bar)/c_bar	! utility coeff
real(dp), parameter :: zeta_z       = zeta_z_log*z_bar/A_bar

! grid parameter (A_t)
real(dp), parameter :: std_uA 	    = 0.033_dp	    ! std for uA_t
integer, parameter	:: uA_std_max   = 4
integer, parameter	:: nuA          = 11
real(dp), parameter :: uA_max       = uA_std_max*std_uA
real(dp), parameter :: uA_min       = -uA_max
real(dp), parameter :: uA_step	    = (uA_max-uA_min)/(nuA-1)

integer, parameter  :: A_std_max    = 4
integer, parameter  :: nA           = 5
real(dp), parameter :: rho_A 	    = 0.45_dp	! persistent for productivity shock
real(dp), parameter :: A_min        = A_bar*exp(-A_std_max/dsqrt(1.0_dp-rho_A**2)*std_uA)  
real(dp), parameter :: A_max        = A_bar*exp(A_std_max/dsqrt(1.0_dp-rho_A**2)*std_uA)
real(dp), parameter :: A_step       = (A_max-A_min)/(nA-1)

! grid parameter (g_t)	
real(dp), parameter :: std_uG 	    = 0.03_dp	    ! std for uG_t	
integer, parameter  :: uG_std_max   = 4	            ! max std (uAgrid uGgrid)		
integer, parameter  :: nuG          = 11	        ! grid points 		
real(dp), parameter :: uG_max       = uG_std_max*std_uG
real(dp), parameter :: uG_min       = -uG_max
real(dp), parameter :: uG_step	    = (uG_max-uG_min)/(nuG-1)	

integer, parameter  :: g_std_max    = 4
integer, parameter  :: ng           = 5
real(dp), parameter :: rho_g        = 0.426_dp      ! iid normal shock to g
real(dp), parameter :: g_min        = g_bar*exp(-g_std_max/dsqrt(1.0_dp-rho_g**2)*std_uG) ! 
real(dp), parameter :: g_max        = g_bar*exp(g_std_max/dsqrt(1.0_dp-rho_g**2)*std_uG)
real(dp), parameter :: g_step       = (g_max-g_min)/(ng-1)	

! grid parameter (zrs_t)
integer, parameter	:: nzrs	            = 2
integer, parameter  :: zrs_grid(nzrs)   = (/1,2/)   
real(dp), parameter :: p11_zrs          = 0.975_dp  !0.95_dp   
real(dp), parameter :: p22_zrs          = 0.975_dp  !0.95_dp   
real(dp), parameter :: pr_zrs(nzrs,nzrs) = reshape((/p11_zrs, 1.0_dp-p22_zrs, 1.0_dp-p11_zrs, p22_zrs/),(/nzrs,nzrs/))    ! 40 years
!!! reshape (a, b, c, d) to 2*2 --> (a,c ; b, d), NOT (a, b; c,d)

! grid parameter (z_t)
integer, parameter  :: nz           = 21 !9
real(dp), parameter :: z_min        = z_bar+zeta_z*(A_max-A_bar)        ! countercyclical
real(dp), parameter :: z_max        = z_bar*mu_z**(1.0_dp/(1.0_dp-p11_zrs))+zeta_z*(A_min-A_bar)    ! z_bar+zeta_z*(A_min-A_bar) !   
real(dp), parameter :: z_step       = (z_max-z_min)/(nz-1)
	
! grid parameter (v_t) 	
integer, parameter  :: nvA          = 11
integer, parameter  :: nvg          = 11        ! depend on the fiscal limit computation (matlab file)
real(dp) :: pr_v(nv,nvA,nvg,nzrs), cdf_v(nv,nvA,nvg,nzrs)
real(dp) :: v_step, pr_v_vec(nv*nvA*nvg*nzrs), cdf_v_vec(nv*nvA*nvg*nzrs)
real(dp) :: vA_grid(nvA), vg_grid(nvg), v_grid(nv)

! grid parameter (delta_t)
integer, parameter :: ndelta    = 26
real(dp) :: pr_delta(ndelta), delta_grid(ndelta), cdf_delta(ndelta)

! grid parameter (b_t)	
!	real(dp), parameter :: b_min    = 0.05_dp
!	real(dp), parameter :: b_max    = 0.75_dp	    ! bgrid upper and lower limit

! initialization (rule)
real(dp), dimension(nb,nA,ng,nz,nzrs) :: b_rule
real(dp), dimension(nb,nA,ng,nz,nzrs) :: b_rule_new

! other paramters
integer, parameter  :: nfcn = 1                     ! number of equations in fcn_eulor
integer, parameter  :: nstate = 5                   ! number of state var	
integer, parameter  :: nmax_nsolve = 10             ! max number of iterations in nonlinear solver

integer     :: openstatus, ib, iA,ig,iz,izrs,index_compute	    
integer     :: iter, index_t(nstate), rc_nsolve
real(dp)    :: b_max_diff,b_diff_iter
real(dp)    :: vector_t(nfcn)
integer     :: time_now(3)                          ! measure the timing

! parameter for loading initial guess (with n-1 dimension)
integer, parameter  :: ngrid_ini_short = nb_ini*nA*nz*nzrs
real(dp)            :: b_rule_ini_mat_short(ngrid_ini_short)
real(dp), dimension(nb_ini,nA,nz,nzrs) :: b_rule_ini_short

! parameter for loading initial guess (with n dimension and shorter grid)
integer, parameter  :: ngrid_ini = nb_ini*nA*nz*ng*nzrs
real(dp)            :: b_rule_ini_mat(ngrid_ini)
real(dp), dimension(nb_ini,nA,ng,nz,nzrs) :: b_rule_ini	

! parameter for loading initial guess (with n dimension)
integer, parameter  :: ngrid = nb*nA*nz*ng*nzrs
real(dp)            :: b_rule_ini_mat_nd(ngrid)

!-------------------------------------------------------------------------------
! Construct Grid and Prob matrix
!-------------------------------------------------------------------------------
real(dp) :: uG_grid(nuG), uA_grid(nuA), b_grid(nb), g_grid(ng), A_grid(nA),z_grid(nz)
real(dp) :: pr_uG(nuG), pr_uA(nuA), b_grid_ini(nb_ini)


contains

!************************************************************************
!
!  SUBROUTINE: load data
!
!************************************************************************

subroutine load_data(b_rule, b_grid, v_grid, vA_grid, vg_grid, pr_v_vec, cdf_v_vec, delta_grid, pr_delta, cdf_delta)
    use mod_interp_locate
    implicit none
    real(dp), intent(out) :: b_grid(nb), v_grid(nv), vA_grid(nvA), vg_grid(nvg), pr_v_vec(nv*nuA*nuG*nzrs)
    real(dp), intent(out) :: cdf_v_vec(nv*nuA*nuG*nzrs), b_rule(nb,nA,ng,nz,nzrs)
    real(dp), intent(out) :: delta_grid(ndelta), pr_delta(ndelta), cdf_delta(ndelta)

    ! save b_grid.dat b_grid -ascii -double (in Matlab)
    open(unit = 8, file = "b_grid.dat", status = "old", action = "read", iostat = openstatus)
    read(8,*) b_grid
    close(8)		
    ! save v_grid.dat 
    open(unit = 8, file = "v_grid.dat", status = "old", action = "read", iostat = openstatus)
    read(8,*) v_grid
    close(8)	
    ! save vA_grid.dat 
    open(unit = 8, file = "vA_grid.dat", status = "old", action = "read", iostat = openstatus)
    read(8,*) vA_grid
    close(8)
    ! save vg_grid.dat 
    open(unit = 8, file = "vg_grid.dat", status = "old", action = "read", iostat = openstatus)
    read(8,*) vg_grid
    close(8)
    ! save pr_v
    open(unit = 8, file = "pr_v_vec.dat", status = "old", action = "read", iostat = openstatus)
    read(8,*) pr_v_vec
    close(8)
    ! save cdf_v
    open(unit = 8, file = "cdf_v_vec.dat", status = "old", action = "read", iostat = openstatus)
    read(8,*) cdf_v_vec
    close(8)
    ! save v_grid.dat 
    open(unit = 8, file = "delta_grid.dat", status = "old", action = "read", iostat = openstatus)
    read(8,*) delta_grid
    close(8)	
    ! save pr_v
    open(unit = 8, file = "pr_delta.dat", status = "old", action = "read", iostat = openstatus)
    read(8,*) pr_delta
    close(8)
    ! save cdf_v
    open(unit = 8, file = "cdf_delta.dat", status = "old", action = "read", iostat = openstatus)
    read(8,*) cdf_delta
    close(8)

    !-------------------------------------------------------------------------------
    ! Initialization (rules): specify or load
    !-------------------------------------------------------------------------------

    if (load_index == 1) then
	! initialization following bgrid
		do ib = 1,nb
		    b_rule(ib,:,:,:,:) = b_grid(ib)
	    end do	    
	    
	elseif (load_index == 2) then
	! load ini mat and expand to 5 dimensions
	    open(unit = 8, file = "b_rule_ini_mat.dat", status = "old", action = "read", iostat = openstatus)
	    read(8,*) b_rule_ini_mat_short
	    close(8)
	    b_rule_ini_short = reshape(b_rule_ini_mat_short,(/nb_ini,nA,nz,nzrs/))	
	    ! generate the initial guess of b_rule 	    
	    do ib = 1,nb
	    do iA = 1,nA
	    do ig = 1,ng
	    do iz = 1,nz
	    do izrs = 1,nzrs
	        b_rule(ib,iA,ig,iz,izrs) = b_rule_ini_short(ib,iA,iz,izrs)
	    end do
	    end do 
	    end do 
	    end do
	    end do   
	
	elseif (load_index == 3) then
	! load ini mat
	    open(unit = 8, file = "b_rule_ini_mat.dat", status = "old", action = "read", iostat = openstatus)
	    read(8,*) b_rule_ini_mat_nd
	    close(8)
	    b_rule = reshape(b_rule_ini_mat_nd,(/nb,nA,ng,nz,nzrs/))	
	        
	elseif (load_index == 4) then    
	! load ini mat and interp
	    open(unit = 8, file = "b_grid_ini.dat", status = "old", action = "read", iostat = openstatus)
	    read(8,*) b_grid_ini
	    close(8)
	    open(unit = 8, file = "b_rule_ini_mat.dat", status = "old", action = "read", iostat = openstatus)
	    read(8,*) b_rule_ini_mat
	    close(8)
	    b_rule_ini = reshape(b_rule_ini_mat,(/nb_ini,nA,ng,nz,nzrs/))    	
	    ! generate the initial guess of b_rule
	    do ib = 1,nb
	    do iA = 1,nA
	    do ig = 1,ng
	    do iz = 1,nz
	    do izrs = 1,nzrs
	        b_rule(ib,iA,ig,iz,izrs) = interp1(b_grid_ini,nb_ini,b_grid(ib),b_rule_ini(:,iA,ig,iz,izrs))
	    end do
	    end do
	    end do 
	    end do 
	    end do
	end if  

end subroutine load_data


!************************************************************************
!
!  SUBROUTINE: construct grid
!
!************************************************************************
subroutine construct_grid(g_grid_local,A_grid_local,z_grid_local,uG_grid_local, uA_grid_local,pr_uG_local,pr_uA_local)

	implicit none
  
    real(dp), intent(out) :: pr_uG_local(nuG), pr_uA_local(nuA), g_grid_local(ng),A_grid_local(nA),z_grid_local(nz)	
	real(dp), intent(out) :: uG_grid_local(nuG), uA_grid_local(nuA)
	
	real(dp) :: x1_local, x2_local, sqrt_2
	integer :: i_local, j_local, k_local
		
	! g grid, uA grid, uG grid, v grid, z grid, uZ grid
	do i_local = 1,ng
		g_grid_local(i_local) = g_min+(i_local-1)*g_step
	end do	
	do i_local = 1,nuG
		uG_grid_local(i_local) = uG_min+(i_local-1)*uG_step
	end do	
	do i_local = 1,nA
		A_grid_local(i_local) = A_min+(i_local-1)*A_step
	end do
	do i_local = 1,nuA
		uA_grid_local(i_local) = uA_min+(i_local-1)*uA_step
	end do	
	do i_local = 1,nz
		z_grid_local(i_local) = z_min+(i_local-1)*z_step
	end do
	
	! compute transition prob matrix
	do i_local = 1,nuG
	    pr_uG_local(i_local)=(2.0_dp*pi)**(-0.5_dp)/std_uG*exp(-(uG_grid_local(i_local)**2.0_dp)/(std_uG**2.0_dp)/2.0_dp)
	end do	
	do i_local = 1,nuA
	    pr_uA_local(i_local)=(2.0_dp*pi)**(-0.5_dp)/std_uA*exp(-(uA_grid_local(i_local)**2.0_dp)/(std_uA**2.0_dp)/2.0_dp)
	end do
	

end subroutine construct_grid


END MODULE mod_parameter

