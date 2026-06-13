program check
use mle_kinds_mod
use empirical_pl_mod
implicit none
type(clock_time) :: clock
type(empirical_pl) :: pl_1 , pl_2
integer(i4) :: file_freq_mb , file_rank_mb
real(dp) :: time_fit , time_pvalue

open( newunit=file_freq_mb , file='test/mb_freq_data.dat' , action='read' , status='old' )
open( newunit=file_rank_mb , file='test/mb_rank_data.dat' , action='read' , status='old' )

call clock%start()
call pl_1%data%from_file( file_freq_mb )
call pl_1%fast_fit()
call clock%stop()
time_fit = clock%elapsed
call clock%start()
call pl_1%p_value( N_samples=100 )
call clock%stop()
time_pvalue = clock%elapsed
call pl_1%report()
print*, ""
print*, "Consumed time for fitting data:", time_fit
print*, "Consumed time for p-valuing data:", time_pvalue

call clock%start()
call pl_2%data%from_file( file_rank_mb )
call pl_2%fast_fit()
call clock%stop()
time_fit = clock%elapsed
call clock%start()
call pl_2%p_value( N_samples=100 )
call clock%stop()
time_pvalue = clock%elapsed
call pl_2%report()
print*, ""
print*, "Consumed time for fitting data:", time_fit
print*, "Consumed time for p-valuing data:", time_pvalue

end program check