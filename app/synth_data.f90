program synth_data
use kinds_mod
use rndgen_mod
use power_law_mle_fit_mod
implicit none
integer(i4) , parameter :: N = 50000 , n_tail(5) = [ 5*9000 , 5*8000 , 5*7000 , 5*6000 , 5*5000 ]
real(dp) , parameter :: alpha = 2._dp , p_tail(5) = [ 0.9_dp , 0.8_dp , 0.7_dp , 0.6_dp , 0.5_dp ]
real(dp) , allocatable :: r_data( : , : )
real(dp) :: x_min , mle_x_min , mle_alpha , lambda
type(empirical_pl) :: pl
type(rndgen) :: gen_1 , gen_2
integer(i4) :: seed = 294727492 , time_values(8)
integer(i4) :: i , j , case
  
call gen_1%init(seed)
call date_and_time(values=time_values)
call gen_2%init( time_values(8) + 1000*time_values(7) + 60000*time_values(6) + 3600000*time_values(5))

allocate( r_data( 1:size(n_tail) , 1:N ) )

do j=1,size(n_tail)
    x_min = real(N-n_tail(j),dp)
    do i=1,N
    if ( i <= n_tail(j) ) then
        r_data( j , i ) = x_min*(1.d0-gen_2%rnd())**(-1.d0/(alpha-1.d0))
    else
        r_data( j , i ) = gen_2%real( 0._dp , x_min )
    endif
    enddo
    call pl%init( r_data( j , : ) )
    call pl%fast_fit( r_data( j , : ) , xmin=mle_x_min , alpha=mle_alpha )
    print*, ""
    print*, "Ajuste MLE obteve:"
    print*, "x_min = ", mle_x_min
    print*, "alpha = ", mle_alpha
    print*, "Valores reais:"
    print*, "x_min = ", x_min
    print*, "alpha = ", alpha
    print*, ""
enddo



do j=1,size(n_tail)
x_min = real(N-n_tail(j),dp)
lambda = 0.0_dp ; case = 0
do
    case = case + 1
    lambda = lambda + 0.01_dp
    if ( lambda > 0.1_dp ) exit
    call pl%fast_fit( r_data(j,:) , xmin=mle_x_min , alpha=mle_alpha , lambda_in=lambda )
    print*, "Valor de lambda:" , lambda
    print*, "Ajuste MLE+tail obteve :"
    print*, "x_min = ", mle_x_min
    print*, "alpha = ", mle_alpha
    call pl%wfast_fit( r_data(j,:) , xmin=mle_x_min , alpha=mle_alpha , lambda_in=lambda )
    print*, "Valor de lambda:" , lambda
    print*, "Ajuste weighted MLE+tail obteve :"
    print*, "x_min = ", mle_x_min
    print*, "alpha = ", mle_alpha
    print*, "Valores reais:"
    print*, "x_min = ", x_min
    print*, "alpha = ", alpha
    print*, ""
enddo
enddo

end program synth_data