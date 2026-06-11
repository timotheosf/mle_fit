program benchmark
use kinds_mod
use rndgen_mod
use empirical_pl_mod
implicit none
integer(i4) , parameter :: N(3) = [ 1000 , 10000 , 50000 ]
real(dp) , parameter :: alpha = 2._dp
real(dp) , allocatable :: r_data( : , : )
real(dp) :: x_min , mle_x_min , mle_alpha , lambda , time_fit , time_pvalued
type(clock_time) :: clock_control
type(empirical_pl) :: pl_1 , pl_2
type(rndgen) :: gen , gen
integer(i4) :: seed = 294727492 , time_values(8)
integer(i4) :: i , j , case , io , n_tail(5) , samples , time_OPM_unit , fit_unit
integer(i4) , allocatable :: arr( : )

open(newunit=time_OPM_unit, file="benchmark_OPM.dat", status="replace", action="write")
open(newunit=fit_unit     , file="benchmark_fit.dat", status="replace", action="write")

call date_and_time(values=time_values) !> Get a random seed based on time
call gen%init( time_values(8) + 1000*time_values(7) + 60000*time_values(6) + 3600000*time_values(5) )

write( time_OPM_unit , * ) "N case time_fit time_pvalued"
write( fit_unit , * ) "N case x_min alpha p_value"

do samples = 1 , size(N)
    n_tail = [ (N(samples)*9)/10 , (N(samples)*8)/10 , (N(samples)*7)/10 , (N(samples)*6)/10 , (N(samples)*5)/10 ]
    
    allocate( r_data(size(n_tail),N(samples)) )
    do i = 1 , size(n_tail)
        x_min = real(N(samples) - n_tail(i),dp)
        do j=1,N(samples)
            if ( j <= n_tail(i) ) then
                r_data(i,j) = x_min*(1.d0-gen%rnd())**(-1.d0/(alpha-1.d0))
            else
                r_data(i,j) = gen%real(0.001_dp,x_min)
            endif
        enddo
        !> Time control for fit 
        call clock_control%start()
        call pl_1%init( r_data( i , : ) )
        call pl_1%fast_fit( )
        call clock_control%stop()
        time_fit = clock_control%elapsed
        call clock_control%start()
        call pl_1%p_value( N_samples=2500 )
        call clock_control%stop()
        time_pvalued = clock_control%elapsed
        write( time_OPM_unit , * ) N(samples) , "pl_1" , time_fit , time_pvalued
        print*, N(samples) , "pl_1" , time_fit , time_pvalued
        call clock_control%start()
        call pl_2%init( r_data( i , : ) )
        call pl_2%fast_fit( )
        call clock_control%stop()
        time_fit = clock_control%elapsed
        call clock_control%start()
        call pl_2%p_value( N_samples=2500 )
        call clock_control%stop()
        time_pvalued = clock_control%elapsed
        write( time_OPM_unit , * ) N(samples) , "pl_2" , time_fit , time_pvalued
        print*, N(samples) , "pl_2" , time_fit , time_pvalued

        write( fit_unit , * ) N(samples) , "real" , x_min , alpha , ""
        write( fit_unit , * ) N(samples) , "pl_1" , pl_1%x_min , pl_1%alpha , pl_1%goodness_of_fit
        write( fit_unit , * ) N(samples) , "pl_2" , pl_2%x_min , pl_2%alpha , pl_2%goodness_of_fit

    enddo
    deallocate( r_data )
enddo



end program benchmark