program check
use mle_kinds_mod
use power_law_mod   !> O módulo da Power Law
use mle_fit_mod     !> O módulo do Motor Central
implicit none

type(clock_time) :: clock

!> Agora instanciamos os motores e as distribuições separadamente
type(mle) :: engine_1, engine_2
type(power_law) :: pl_1, pl_2

integer(i4) :: file_freq_mb, file_rank_mb
real(dp) :: time_fit, time_pvalue

open( newunit=file_freq_mb , file='test/mb_freq_data.dat' , action='read' , status='old' )
open( newunit=file_rank_mb , file='test/mb_rank_data.dat' , action='read' , status='old' )

!> 1. Inicializar (Acordar) as distribuições
call pl_1%start_pl()
call pl_2%start_pl()

! ==========================================
! AJUSTE DO DATASET 1 (freq_mb)
! ==========================================
call clock%start()
call engine_1%data%from_file( file_freq_mb )

!> Passamos a distribuição usando a keyword 'distribuiton' pois 'r_data' era opcional antes dela
call engine_1%fast_fit( distribuiton=pl_1 )

call clock%stop()
time_fit = clock%elapsed

call clock%start()
call engine_1%p_value( N_samples=100 )
call clock%stop()
time_pvalue = clock%elapsed

call engine_1%report()
print*, ""
print*, "Consumed time for fitting data (external clock):", time_fit
print*, "Consumed time for p-valuing data (external clock):", time_pvalue

! ==========================================
! AJUSTE DO DATASET 2 (rank_mb)
! ==========================================
call clock%start()
call engine_2%data%from_file( file_rank_mb )

call engine_2%fast_fit( distribuiton=pl_2 )

call clock%stop()
time_fit = clock%elapsed

call clock%start()
call engine_2%p_value( N_samples=100 )
call clock%stop()
time_pvalue = clock%elapsed

call engine_2%report()
print*, ""
print*, "Consumed time for fitting data (external clock):", time_fit
print*, "Consumed time for p-valuing data (external clock):", time_pvalue

!> Boa prática de fechar os arquivos
close(file_freq_mb)
close(file_rank_mb)

end program check