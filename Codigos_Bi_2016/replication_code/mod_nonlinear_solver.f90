!****************************************************************************
!
!  Modules: csolve
!
!****************************************************************************
include 'mkl_vsl.fi'

MODULE mod_nonlinear_solver

use mod_global_constant
use mod_bstar_sim          ! if compute the fiscal limit  
!use mod_fcn_eulor         

IMPLICIT NONE

CONTAINS

!****************************************************************************
!
!  subroutine: csolve 
!
!****************************************************************************

subroutine csolve(x_csolve,index_csolve, nx_csolve, tol_csolve,nmax_csolve,rc_csolve,nstate_csolve)

    USE MKL_VSL_TYPE
    USE MKL_VSL

    implicit none        
    ! parameter specification        
    real(dp), parameter :: delta_csolve = 0.000001_dp        ! differencing interval for numerical gradient
    real(dp), parameter :: alpha_csolve = 0.001_dp                    ! tolerance on rate of descent 
    real(dp), parameter :: nm_mean = 0.0_dp, nm_std = 1.0_dp          ! random number generator
       
    ! input and output
    integer, intent(in)     :: nstate_csolve
    integer, intent(in)     :: nx_csolve
    integer, intent(in)     :: index_csolve(nstate_csolve)
    real(dp), intent(inout) :: x_csolve(nx_csolve)
    real(dp), intent(in)    :: tol_csolve                ! convergence crit
    integer, intent(in)     :: nmax_csolve               ! max # of iterations
    integer, intent(out)    :: rc_csolve          
    ! var specification
    real(dp),dimension(nx_csolve,nx_csolve) :: tvec_csolve, grad_csolve
    real(dp),dimension(nx_csolve) :: f0_csolve,dx0_csolve,fmin_csolve,dx_csolve,xmin_csolve
    real(dp),dimension(nx_csolve) :: xph_csolve,h_csolve,xj_csolve,f_csolve
    real(dp)    :: af0_csolve,af00_csolve,xnorm_csolve, af_csolve,d_af0_csolve,up_af0_csolve,up_inv_af0_csolve
    real(dp)    :: lambda_csolve,lambdamin_csolve,afmin_csolve,dxsize_csolve,factor_csolve,factor_06_csolve
    real(dp)    :: up_abs_lambda_csolve, up_delta_csolve
    integer     :: ix_csolve,done_csolve,itct_csolve,itctmod_csolve,randomize_csolve,shrink_csolve,subdone_csolve
    ! lapack routine parameters
    real(dp)    :: b_csolve(nx_csolve,1)
    integer     :: infoinv_csolve, lda_csolve, ldb_csolve, nrhs_csolve, ipiv_csolve(nx_csolve)
    ! random number generator parameters
    real(dp), dimension(nx_csolve) :: randnormx_csolve
    integer :: errcode,brng,method,seed
    type (VSL_STREAM_STATE) :: stream
    
    
    ! step matrix (eye matrix)
    tvec_csolve = 0.0_dp
    do ix_csolve = 1,nx_csolve
        tvec_csolve(ix_csolve,ix_csolve) = delta_csolve
    end do        
    ! initialization
    f0_csolve = fcn_nonlinear(x_csolve,index_csolve)    
    af0_csolve = sum(dabs(f0_csolve))
    af00_csolve = af0_csolve
    itct_csolve = 0
    factor_06_csolve = 0.6_dp      ! Powell's hybrid method
    
    !-------------------------------------------------------------------
    ! computation loop
    done_csolve = 0
    do while (done_csolve==0)
        ! whether randomize_csolve
        itctmod_csolve = mod(itct_csolve,2)
        if ((itct_csolve>3) .AND. (af00_csolve-af0_csolve<tol_csolve*max(1.0_dp,af0_csolve)) .AND. (itctmod_csolve==1)) then
            randomize_csolve =1

        else
            !----------------compute the jacobian matrix grad------------------------
            ! use numerical recipe fdjac. Forward-difference approximation to Jacobian
            xj_csolve=x_csolve
            h_csolve=delta_csolve  !delta_csolve*dabs(x_csolve)
            where (h_csolve == 0.0_dp) h_csolve=delta_csolve
            xph_csolve=x_csolve+h_csolve
            !h_csolve=xph_csolve-x_csolve
            do ix_csolve=1,nx_csolve
                xj_csolve(ix_csolve)=xph_csolve(ix_csolve)
                grad_csolve(:,ix_csolve)=(fcn_nonlinear(xj_csolve,index_csolve)-f0_csolve(:))/h_csolve(ix_csolve)                
                xj_csolve(ix_csolve)=x_csolve(ix_csolve)
            end do
            
            ! update dx
            nrhs_csolve=1
            lda_csolve=nx_csolve
            ldb_csolve=nx_csolve
            b_csolve(:,1)=-f0_csolve
            call DGESV(nx_csolve, nrhs_csolve, grad_csolve, lda_csolve, ipiv_csolve, b_csolve, ldb_csolve, infoinv_csolve)
            !print *, 'infoinv_csolve', infoinv_csolve
            dx0_csolve = b_csolve(:,1)   
            
            randomize_csolve = 0

        end if   
        
        !-----update dx0 through randomization
        if (randomize_csolve == 1) then                
            ! calculate norm(x_csolve)
            xnorm_csolve = dsqrt(sum(x_csolve*x_csolve))
            print *, '|---- randomization is activated ---|'
!            pause
            ! update dx0
            brng= VSL_BRNG_MT19937               ! brng =VSL_BRNG_SOBOL //VSL_BRNG_MRG32K3A    
            method= VSL_RNG_METHOD_GAUSSIAN_ICDF!VSL_METHOD_SGAUSSIAN_ICDF
            seed= 10100             
            errcode=vslnewstream(stream, brng, seed)            ! Initialize 
            errcode=vdrnggaussian(method, stream, nx_csolve, randnormx_csolve, nm_mean, nm_std)   ! call rng      
            errcode=vsldeletestream(stream)                     ! Deinitialize 
            
            do ix_csolve = 1,nx_csolve
                dx0_csolve(ix_csolve) = xnorm_csolve/randnormx_csolve(ix_csolve)
            end do
        end if

        !-------------------------------------------------------------------
        ! Powell's hybrid method
        lambda_csolve = 1.0_dp
        lambdamin_csolve = 1.0_dp
        fmin_csolve = f0_csolve
        xmin_csolve = x_csolve
        afmin_csolve = af0_csolve
        dxsize_csolve =dsqrt(sum(dx0_csolve*dx0_csolve))
        factor_csolve = 0.6_dp
        shrink_csolve = 1
        subdone_csolve = 0
        do while (subdone_csolve==0)
            dx_csolve = lambda_csolve*dx0_csolve
            f_csolve = fcn_nonlinear(x_csolve+dx_csolve,index_csolve)
            af_csolve = sum(dabs(f_csolve))
            
            if (af_csolve<afmin_csolve) then
                afmin_csolve = af_csolve
                fmin_csolve = f_csolve
                lambdamin_csolve = lambda_csolve
                xmin_csolve = x_csolve + dx_csolve
            end if
            
            d_af0_csolve = af0_csolve-af_csolve
            up_af0_csolve = alpha_csolve*lambda_csolve*af0_csolve
            up_inv_af0_csolve = (1-alpha_csolve)*lambda_csolve*af0_csolve
           
             !==========================================================================
             if (((lambda_csolve>0.0_dp) .AND. (d_af0_csolve<up_af0_csolve)) .OR. ((lambda_csolve<0.0_dp) .AND. (d_af0_csolve<0.0_dp))) then
!            if (((lambda_csolve>0.0_dp) .AND. (af0_csolve-af_csolve<alpha_csolve*lambda_csolve*af0_csolve)) &
!            .OR. ((lambda_csolve<0.0_dp) .AND. (af0_csolve-af_csolve<0.0_dp))) then
                
                !---------------------------------------------------------------------------
                if (shrink_csolve == 0) then
                    factor_csolve = factor_csolve**(0.6_dp)
                    shrink_csolve = 1
                end if
                
                up_abs_lambda_csolve = dabs(lambda_csolve*(1-factor_csolve))*dxsize_csolve
                up_delta_csolve = delta_csolve*0.1_dp
                
                if (up_abs_lambda_csolve>up_delta_csolve) then
                !if ((dabs(lambda_csolve*(1-factor_csolve))*dxsize_csolve)>delta_csolve*0.1_dp) then
                    lambda_csolve = lambda_csolve*factor_csolve
                elseif ((lambda_csolve>0.0_dp).AND.(factor_csolve==factor_06_csolve)) then
                    lambda_csolve = -0.3_dp
                else
                    subdone_csolve = 1
                    if (lambda_csolve>0.0_dp) then
                        if (factor_csolve==factor_06_csolve) then
                            rc_csolve =2
                        else 
                            rc_csolve = 1
                        end if
                    else
                        rc_csolve = 3
                    end if
                end if
             
             !==========================================================================
             elseif ((lambda_csolve>0.0_dp).AND.(d_af0_csolve>up_inv_af0_csolve)) then
                !print *, 'step 2'   
            !elseif ((lambda_csolve>0.0_dp).AND.((af_csolve-af0_csolve)>(1-alpha_csolve)*lambda_csolve*af0_csolve)) then
                if (shrink_csolve == 1) then
                    factor_csolve = factor_csolve**0.6_dp
                    shrink_csolve = 0
                end if
                lambda_csolve = lambda_csolve/factor_csolve
            
            !==========================================================================
            else    
                !print *, 'step 3'            
                subdone_csolve = 1
                rc_csolve = 0
            end if
        end do
       
        itct_csolve = itct_csolve +1
        x_csolve = xmin_csolve
        f0_csolve = fmin_csolve
        af00_csolve = af0_csolve
        af0_csolve = afmin_csolve
        if (itct_csolve>=nmax_csolve) then
            done_csolve = 1
            rc_csolve = 4
        elseif (af0_csolve<tol_csolve) then
            done_csolve = 1
            rc_csolve =0
        end if
    
    end do

end subroutine csolve

END MODULE mod_nonlinear_solver

