!****************************************************************************
!
!  Modules: interp1, interp5, locate
!
!****************************************************************************

MODULE mod_interp_locate

use mod_global_constant
IMPLICIT NONE

CONTAINS

!************************************************************************
!
!  FUNCTION: locate (Written by Troy Davig)
!
!************************************************************************

function locate(xi,nx,x_grid)

	implicit none

	integer, intent(in) :: nx
    real(dp), intent(in) :: xi
	real(dp), intent(in) :: x_grid(nx)

	integer :: i_low, i_mid, i_up, locate
	logical :: increasing

    increasing = (x_grid(nx) >= x_grid(1))
	i_low = 0
	i_up = nx+1

	do
    	if (i_up-i_low <= 1) exit
		i_mid = (i_up + i_low)/2
		if (increasing .eqv. (xi > x_grid(i_mid))) then
          	i_low = i_mid
		else
		  	i_up = i_mid
		end if
    end do

	if (xi <= x_grid(1)) then
	    locate = 1
	else if (xi >= x_grid(nx)) then
	    locate = nx-1 	 
	else
	    locate = i_low
	end if

end function locate



!************************************************************************
!
!  FUNCTION: interp1 and 6 (written by Troy Davig)
!
!	- 1/6-D linear interpolation/extrapolation
!
!************************************************************************

function interp1(x1,nx1,x1i,z)

	implicit none
	integer, intent(in) :: nx1
    real(dp), intent(in), dimension(nx1) :: x1
	real(dp), intent(in), dimension(nx1) :: z
	
	real(dp) :: interp1, w1(2),x1i
	integer :: m1, loc1

    loc1 = locate(x1i,nx1,x1)
	w1(2) = (x1i-x1(loc1)) / (x1(loc1+1)-x1(loc1))
    w1(1) = 1-w1(2)
	interp1 = 0.0_dp

	do m1 = 0, 1	  
	    interp1 = interp1 + w1(m1+1)*z(loc1+m1)		            
	end do

end function interp1

!------------------------------------------------------------------------------------------
function interp2(x1,nx1,x2,nx2,x1i,x2i,z)

	implicit none
	integer, intent(in) :: nx1, nx2
    real(dp), intent(in), dimension(nx1) :: x1
	real(dp), intent(in), dimension(nx2) :: x2
	real(dp), intent(in), dimension(nx1,nx2) :: z
	real(dp) :: interp2, w1(2), w2(2),x1i,x2i

	integer :: m1, m2, loc1, loc2

    loc1 = locate(x1i,nx1,x1)
	loc2 = locate(x2i,nx2,x2)

	w1(2) = (x1i-x1(loc1)) / (x1(loc1+1)-x1(loc1))
    w1(1) = 1-w1(2)

    w2(2) = (x2i-x2(loc2)) / (x2(loc2+1)-x2(loc2))
    w2(1) = 1-w2(2)	
    
	interp2 = 0.0_dp

	do m1 = 0, 1
	do m2 = 0, 1
	    interp2 = interp2 + w1(m1+1)*w2(m2+1)*z(loc1+m1,loc2+m2)
	end do
	end do	

end function interp2

!------------------------------------------------------------------------------------------
function interp3(x1,nx1,x2,nx2,x3,nx3,x1i,x2i,x3i,z)

	implicit none
	integer, intent(in) :: nx1, nx2, nx3
    real(dp), intent(in), dimension(nx1) :: x1
	real(dp), intent(in), dimension(nx2) :: x2
	real(dp), intent(in), dimension(nx3) :: x3
	real(dp), intent(in), dimension(nx1,nx2,nx3) :: z
	real(dp) :: interp3, w1(2), w2(2), w3(2), x1i,x2i,x3i

	integer :: m1, m2, m3, loc1, loc2, loc3

    loc1 = locate(x1i,nx1,x1)
	loc2 = locate(x2i,nx2,x2)
	loc3 = locate(x3i,nx3,x3)

	w1(2) = (x1i-x1(loc1)) / (x1(loc1+1)-x1(loc1))
    w1(1) = 1-w1(2)

    w2(2) = (x2i-x2(loc2)) / (x2(loc2+1)-x2(loc2))
    w2(1) = 1-w2(2)
	
	w3(2) = (x3i-x3(loc3)) / (x3(loc3+1)-x3(loc3))
    w3(1) = 1-w3(2)		
    
	interp3 = 0.0_dp

	do m1 = 0, 1
	do m2 = 0, 1
    do m3 = 0, 1
	    interp3 = interp3 + w1(m1+1)*w2(m2+1)*w3(m3+1)*z(loc1+m1,loc2+m2,loc3+m3)
	end do	
	end do
	end do	

end function interp3

!------------------------------------------------------------------------------------

function interp4(x1,nx1,x2,nx2,x3,nx3,x4,nx4,x1i,x2i,x3i,x4i,z)

	implicit none
	integer, intent(in) :: nx1, nx2, nx3, nx4
    real(dp), intent(in), dimension(nx1) :: x1
	real(dp), intent(in), dimension(nx2) :: x2
	real(dp), intent(in), dimension(nx3) :: x3
	real(dp), intent(in), dimension(nx4) :: x4
	real(dp), intent(in), dimension(nx1,nx2,nx3,nx4) :: z
	real(dp) :: interp4, w1(2), w2(2), w3(2), w4(2), x1i,x2i,x3i,x4i

	integer :: m1, m2, m3, m4, loc1, loc2, loc3, loc4

    loc1 = locate(x1i,nx1,x1)
	loc2 = locate(x2i,nx2,x2)
	loc3 = locate(x3i,nx3,x3)
	loc4 = locate(x4i,nx4,x4)

	w1(2) = (x1i-x1(loc1)) / (x1(loc1+1)-x1(loc1))
    w1(1) = 1-w1(2)

    w2(2) = (x2i-x2(loc2)) / (x2(loc2+1)-x2(loc2))
    w2(1) = 1-w2(2)
	
	w3(2) = (x3i-x3(loc3)) / (x3(loc3+1)-x3(loc3))
    w3(1) = 1-w3(2)	
    
    w4(2) = (x4i-x4(loc4)) / (x4(loc4+1)-x4(loc4))
    w4(1) = 1-w4(2)
    	
    
	interp4 = 0.0_dp

	do m1 = 0, 1
	do m2 = 0, 1
    do m3 = 0, 1
	do m4 = 0, 1
	    interp4 = interp4 + w1(m1+1)*w2(m2+1)*w3(m3+1)*w4(m4+1)*z(loc1+m1,loc2+m2,loc3+m3,loc4+m4)
	end do
	end do	
	end do
	end do

end function interp4

!--------------------------------------------------------------------------------------

function interp5(x1,nx1,x2,nx2,x3,nx3,x4,nx4,x5,nx5,x1i,x2i,x3i,x4i,x5i,z)

	implicit none
	integer, intent(in) :: nx1, nx2, nx3, nx4, nx5
    real(dp), intent(in), dimension(nx1) :: x1
	real(dp), intent(in), dimension(nx2) :: x2
	real(dp), intent(in), dimension(nx3) :: x3
	real(dp), intent(in), dimension(nx4) :: x4
	real(dp), intent(in), dimension(nx5) :: x5
	real(dp), intent(in), dimension(nx1,nx2,nx3,nx4,nx5) :: z
	real(dp) :: interp5, w1(2), w2(2), w3(2), w4(2), w5(2), x1i,x2i,x3i,x4i,x5i

	integer :: m1, m2, m3, m4, m5, loc1, loc2, loc3, loc4, loc5

    loc1 = locate(x1i,nx1,x1)
	loc2 = locate(x2i,nx2,x2)
	loc3 = locate(x3i,nx3,x3)
	loc4 = locate(x4i,nx4,x4)
	loc5 = locate(x5i,nx5,x5)

	w1(2) = (x1i-x1(loc1)) / (x1(loc1+1)-x1(loc1))
    w1(1) = 1-w1(2)

    w2(2) = (x2i-x2(loc2)) / (x2(loc2+1)-x2(loc2))
    w2(1) = 1-w2(2)
	
	w3(2) = (x3i-x3(loc3)) / (x3(loc3+1)-x3(loc3))
    w3(1) = 1-w3(2)	
    
    w4(2) = (x4i-x4(loc4)) / (x4(loc4+1)-x4(loc4))
    w4(1) = 1-w4(2)
    
    w5(2) = (x5i-x5(loc5)) / (x5(loc5+1)-x5(loc5))
    w5(1) = 1-w5(2)		
    
	interp5 = 0.0_dp

	do m1 = 0, 1
	do m2 = 0, 1
    do m3 = 0, 1
	do m4 = 0, 1
	do m5 = 0, 1
	    interp5 = interp5 + w1(m1+1)*w2(m2+1)*w3(m3+1)*w4(m4+1)*w5(m5+1)*z(loc1+m1,loc2+m2,loc3+m3,loc4+m4,loc5+m5)
	end do	
	end do
	end do	
	end do
	end do

end function interp5

!------------------------------------------------------------------ 	

function interp6(x1,nx1,x2,nx2,x3,nx3,x4,nx4,x5,nx5,x6,nx6,x1i,x2i,x3i,x4i,x5i,x6i,z)

	implicit none
	integer, intent(in) :: nx1, nx2, nx3, nx4, nx5, nx6
    real(dp), intent(in), dimension(nx1) :: x1
	real(dp), intent(in), dimension(nx2) :: x2
	real(dp), intent(in), dimension(nx3) :: x3
	real(dp), intent(in), dimension(nx4) :: x4
	real(dp), intent(in), dimension(nx5) :: x5
	real(dp), intent(in), dimension(nx6) :: x6
	real(dp), intent(in), dimension(nx1,nx2,nx3,nx4,nx5,nx6) :: z
	real(dp) :: interp6, w1(2), w2(2), w3(2), w4(2), w5(2), w6(2), x1i,x2i,x3i,x4i,x5i,x6i

	integer :: m1, m2, m3, m4, m5, m6, loc1, loc2, loc3, loc4, loc5, loc6

    loc1 = locate(x1i,nx1,x1)
	loc2 = locate(x2i,nx2,x2)
	loc3 = locate(x3i,nx3,x3)
	loc4 = locate(x4i,nx4,x4)
	loc5 = locate(x5i,nx5,x5)
	loc6 = locate(x6i,nx6,x6)

	w1(2) = (x1i-x1(loc1)) / (x1(loc1+1)-x1(loc1))
    w1(1) = 1-w1(2)

    w2(2) = (x2i-x2(loc2)) / (x2(loc2+1)-x2(loc2))
    w2(1) = 1-w2(2)
	
	w3(2) = (x3i-x3(loc3)) / (x3(loc3+1)-x3(loc3))
    w3(1) = 1-w3(2)	
    
    w4(2) = (x4i-x4(loc4)) / (x4(loc4+1)-x4(loc4))
    w4(1) = 1-w4(2)
    
    w5(2) = (x5i-x5(loc5)) / (x5(loc5+1)-x5(loc5))
    w5(1) = 1-w5(2)		
    
    w6(2) = (x6i-x6(loc6)) / (x6(loc6+1)-x6(loc6))
    w6(1) = 1-w6(2)		
    
	interp6 = 0.0_dp

	do m1 = 0, 1
	do m2 = 0, 1
    do m3 = 0, 1
	do m4 = 0, 1
	do m5 = 0, 1
	do m6 = 0, 1
	    interp6 = interp6 + w1(m1+1)*w2(m2+1)*w3(m3+1)*w4(m4+1)*w5(m5+1)*w6(m6+1)*z(loc1+m1,loc2+m2,loc3+m3,loc4+m4,loc5+m5,loc6+m6)
	end do	
	end do
	end do	
	end do
	end do
	end do

end function interp6

!------------------------------------------------------------------ 	

function interp7(x1,nx1,x2,nx2,x3,nx3,x4,nx4,x5,nx5,x6,nx6,x7,nx7,x1i,x2i,x3i,x4i,x5i,x6i,x7i,z)

	implicit none
	integer, intent(in) :: nx1, nx2, nx3, nx4, nx5, nx6,nx7
    real(dp), intent(in), dimension(nx1) :: x1
	real(dp), intent(in), dimension(nx2) :: x2
	real(dp), intent(in), dimension(nx3) :: x3
	real(dp), intent(in), dimension(nx4) :: x4
	real(dp), intent(in), dimension(nx5) :: x5
	real(dp), intent(in), dimension(nx6) :: x6
	real(dp), intent(in), dimension(nx7) :: x7
	real(dp), intent(in), dimension(nx1,nx2,nx3,nx4,nx5,nx6,nx7) :: z
	real(dp) :: interp7, w1(2), w2(2), w3(2), w4(2), w5(2), w6(2), w7(2),x1i,x2i,x3i,x4i,x5i,x6i,x7i

	integer :: m1, m2, m3, m4, m5, m6, m7, loc1, loc2, loc3, loc4, loc5, loc6, loc7

    loc1 = locate(x1i,nx1,x1)
	loc2 = locate(x2i,nx2,x2)
	loc3 = locate(x3i,nx3,x3)
	loc4 = locate(x4i,nx4,x4)
	loc5 = locate(x5i,nx5,x5)
	loc6 = locate(x6i,nx6,x6)
	loc7 = locate(x7i,nx7,x7)

	w1(2) = (x1i-x1(loc1)) / (x1(loc1+1)-x1(loc1))
    w1(1) = 1-w1(2)

    w2(2) = (x2i-x2(loc2)) / (x2(loc2+1)-x2(loc2))
    w2(1) = 1-w2(2)
	
	w3(2) = (x3i-x3(loc3)) / (x3(loc3+1)-x3(loc3))
    w3(1) = 1-w3(2)	
    
    w4(2) = (x4i-x4(loc4)) / (x4(loc4+1)-x4(loc4))
    w4(1) = 1-w4(2)
    
    w5(2) = (x5i-x5(loc5)) / (x5(loc5+1)-x5(loc5))
    w5(1) = 1-w5(2)		
    
    w6(2) = (x6i-x6(loc6)) / (x6(loc6+1)-x6(loc6))
    w6(1) = 1-w6(2)		

    w7(2) = (x7i-x7(loc7)) / (x7(loc7+1)-x7(loc7))
    w7(1) = 1-w7(2)		
   
	interp7 = 0.0_dp

	do m1 = 0, 1
	do m2 = 0, 1
    do m3 = 0, 1
	do m4 = 0, 1
	do m5 = 0, 1
	do m6 = 0, 1
	do m7 = 0, 1
	    interp7 = interp7 + w1(m1+1)*w2(m2+1)*w3(m3+1)*w4(m4+1)*w5(m5+1)* &
                   w6(m6+1)*w7(m7+1)*z(loc1+m1,loc2+m2,loc3+m3,loc4+m4, &
                   loc5+m5,loc6+m6,loc7+m7)
	end do	
	end do
	end do	
	end do
	end do
	end do
	end do

end function interp7
!------------------------------------------------------------------ 	


END MODULE mod_interp_locate

