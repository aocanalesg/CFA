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


END MODULE mod_interp_locate

