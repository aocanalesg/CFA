!****************************************************************************
!
!  Modules: define global constants
!
!****************************************************************************

MODULE mod_global_constant

IMPLICIT NONE
integer, parameter :: dp            = kind(1.0d0)
real(dp),parameter :: sqrt_2        = dsqrt(2.0_dp) 
real(dp), parameter :: pi           = 3.14159265358979_dp   ! define pi


real(dp):: y_bar,b_bar,bstar_bar,g_bar,dn_bar,yn_r_bar,yt_r_bar
real(dp):: rkn_bar,rkt_bar,kn_bar,kt_bar,in_bar,it_bar,i_bar,c_bar,ln_bar,lt_bar,wn_bar,wt_bar,w_bar
real(dp):: varphi_l,at_bar,z_bar,zovery,tb_bar,syn,lambda_bar,phi
real(dp):: ci_bar,yn_n_bar,yt_n_bar,cin_bar

END MODULE mod_global_constant

