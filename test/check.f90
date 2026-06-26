program check
use mle_kinds_mod
use power_law_mod   
use mle_fit_mod     
implicit none
type(mle) :: engine_1, engine_2
type(power_law) :: pl_1, pl_2
integer(i4) :: file_freq_mb, file_rank_mb

open( newunit=file_freq_mb , file='test/mb_freq_data.dat' , action='read' , status='old' )
open( newunit=file_rank_mb , file='test/mb_rank_data.dat' , action='read' , status='old' )

!> 1. Wake-up distributions
call pl_1%start_pl()
call pl_2%start_pl()


call engine_1%data%from_file( file_freq_mb )
call engine_1%greed_fit( dist=pl_1 , look_whole=.true. )
call engine_1%report()
call pl_1%start_pl()
call engine_1%fast_fit( dist=pl_1 )
call engine_1%p_value( N_samples=2500 )
call engine_1%report()


close(file_freq_mb)
close(file_rank_mb)

end program check