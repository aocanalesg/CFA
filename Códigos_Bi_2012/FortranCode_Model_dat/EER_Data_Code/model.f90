!****************************************************************************************************
!
!  PROGRAM: model.f90
!  
!  Fortran code for "Sovereign Default Risk Premia, Fiscal Limits, and Fiscal Policy" (Huixin Bi): 
!       published at European Economic Review, 56(3), April 2012: 389-410.
!  
!****************************************************************************************************
!include 'mkl_vsl.fi'
program model

use mod_global_constant
use mod_nonlinear_solver        ! nonlinear solver from Chris Sims' csolve
use mod_parameter
use mod_fcn_eulor
use mod_interp_locate           ! interpolate

implicit none

! construct the grids 
call construct_grid(g_grid,A_grid,z_grid,uG_grid, uA_grid,pr_uG, pr_uA)

! load the data
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+++ Note: need to specify the initial guesses for the decision in order to run the code properly ++++
call load_data(b_rule, b_grid, v_grid, vA_grid, vg_grid, pr_v_vec, cdf_v_vec, delta_grid, pr_delta, cdf_delta)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


! process the input matrix
v_step = v_grid(2)-v_grid(1)
pr_v = reshape(pr_v_vec,(/nv,nuA,nuG,nzrs/))
cdf_v = reshape(cdf_v_vec,(/nv,nuA,nuG,nzrs/))

!-------------------------------------------------------------------------------
! Solve the model
!-------------------------------------------------------------------------------
print *, "Computing the solution..."
print *, "iter max_diff_b "
iter = 2
b_max_diff = 100_dp
b_rule_new = b_rule

index_compute = 2              ! whether or not compute    
if (index_compute == 2) then

    ! elapsed time
    print *, '-----------------------------------------'
	call itime(time_now)
    write (*,12), time_now
    12 format ( ' time ', i2.2, ':', i2.2, ':', i2.2 )
        
    !-----------------------------------------------------------------------------------
    ! iterate until it converges
	do
	    ! exit if it converges
		if (b_max_diff < tol_conv) exit
		
		! initialization
		iter = iter +1
		b_max_diff = -1.0_dp

        ! loop over state vars
        do ib = 1,nb
		do iA = 1,nA
		do ig = 1,ng
		do iz = 1,nz
		do izrs = 1,nzrs
			
			! call nonlinear solver				
			index_t = (/ib,iA,ig,iz,izrs/)
			vector_t = (/b_rule(ib,iA,ig,iz,izrs)/)				
			call csolve(vector_t,index_t, nfcn, tol_nsolve,nmax_nsolve,rc_nsolve,nstate)	
			if (rc_nsolve>1) then
			    print *, '|---- rc_nsolve ----|', rc_nsolve
			    write (*,*), '|---- ib,iA,ig,iz,izrs ----|', ib,iA,ig,iz,izrs
			end if	
			b_rule_new(ib,iA,ig,iz,izrs) = vector_t(1)	    

            ! update the rules
			b_diff_iter = dabs(b_rule_new(ib,iA,ig,iz,izrs)-b_rule(ib,iA,ig,iz,izrs))
			if (b_diff_iter>b_max_diff) then
				b_max_diff = b_diff_iter
			end if
			
		end do
		end do
		end do
		end do
		end do
	
		b_rule = b_rule_new
		write(*,*), iter,b_max_diff		
			
		! elapsed time
		call itime(time_now)
        write (*,10), time_now
        10 format ( ' time ', i2.2, ':', i2.2, ':', i2.2 )
        print *, '-----------------------------------------'            

	end do
    !-----------------------------------------------------------------------------------      
    ! save final results
    open(unit = 1, file = "b_rule_mat.dat", status = "replace", action = "write", iostat = openstatus)
        write(1,'(F24.16)') b_rule
      	close(1)

end if



end program model
