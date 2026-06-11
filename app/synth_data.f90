program benchmark
use kinds_mod
use rndgen_mod
use empirical_pl_mod
implicit none
integer(i4) , parameter :: N = 50000 , n_tail(5) = [ 5*9000 , 5*8000 , 5*7000 , 5*6000 , 5*5000 ]
real(dp) , parameter :: alpha = 2._dp , p_tail(5) = [ 0.9_dp , 0.8_dp , 0.7_dp , 0.6_dp , 0.5_dp ]
real(dp) , allocatable :: r_data( : , : )
real(dp) :: x_min , mle_x_min , mle_alpha , lambda
type(empirical_pl) :: pl
type(rndgen) :: gen_1 , gen_2
integer(i4) :: seed = 294727492 , time_values(8)
integer(i4) :: i , j , case , io 
integer(i4) , allocatable :: arr( : )
  


end program benchmark