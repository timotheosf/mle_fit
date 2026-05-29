program check
use mle_kinds_mod
use pl_mod
use rndgen_mod
use mle_fit_mod
implicit none

real(kind=dp) , parameter :: x_min=1.d0 , alpha=2.d0
type(power_law) :: pl
type(empirical_pl) :: pl_fit1 , pl_fit2
type(rndgen) :: generator
integer(kind=i4) :: seed = 294727492
integer(kind=i4) , allocatable :: dados_inteiros(:)
real(kind=dp) , allocatable :: dados_reais(:)

call generator%init(seed)

!print *, "Put some tests in here!"

!> Test 1 - Power Law type is working?
pl%x_min = x_min ; pl%alpha=alpha
print*, "Power law x_min=", pl%x_min
print*, "Power law alpha=", pl%alpha
print*, ""
print*, "Power law PDF at x=x_min is P(x)=" , pl%pdf( x_min )
print*, "Power law CCDF at x=x_min is P(x)=" , pl%ccdf( x_min )
print*, ""
print*, "Power law real generated x_r=", pl%r_rnd( generator )
print*, "Power law int generated x_r=", pl%i_rnd( generator )

print*, ""
print*, "Generic KS test com valores reais:"
print*, "Com declaração generic KS=", pl%ks_stats([1.1_dp, 2.3_dp, 3.5_dp, 4.8_dp, 5.9_dp])
print*, "Com declaração explícita KS=", pl%ks_stats_real([1.1_dp, 2.3_dp, 3.5_dp, 4.8_dp, 5.9_dp])
print*, "Generic KS test com valores inteiros:"
print*, "Com declaração generic KS=", pl%ks_stats([1, 2, 4, 8, 16])
print*, "Com declaração explícita KS=", pl%ks_stats_int([1, 2, 4, 8, 16])
print*, ""

dados_inteiros = [ &
        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &
        1, 1, 1, 1, 1, 1, 2, 2, 2, 2, &
        2, 2, 2, 2, 3, 3, 3, 3, 4, 4, &
        4, 5, 5, 6, 7, 8, 9, 10, 12, 14, &
        18, 23, 31, 45, 72, 130, 290, 850, 4200, 15000 &
    ]


dados_reais = [ &
        1.01_dp, 1.03_dp, 1.05_dp, 1.08_dp, 1.11_dp, 1.13_dp, 1.15_dp, 1.18_dp, 1.22_dp, 1.25_dp, &
        1.29_dp, 1.33_dp, 1.37_dp, 1.41_dp, 1.45_dp, 1.50_dp, 1.55_dp, 1.62_dp, 1.68_dp, 1.75_dp, &
        1.83_dp, 1.90_dp, 1.98_dp, 2.07_dp, 2.15_dp, 2.24_dp, 2.36_dp, 2.48_dp, 2.61_dp, 2.76_dp, &
        2.92_dp, 3.10_dp, 3.31_dp, 3.55_dp, 3.82_dp, 4.14_dp, 4.52_dp, 4.98_dp, 5.54_dp, 6.25_dp, &
        7.15_dp, 8.35_dp, 10.02_dp, 12.45_dp, 16.32_dp, 23.15_dp, 36.42_dp, 68.51_dp, 185.34_dp, 850.12_dp &
    ]

call pl_fit1%fit(dados_inteiros)
call pl_fit2%fit(dados_reais)

print*, ""
print*, "Teste com o ajuste inteiro:"
print*, "x_min" , pl_fit1%x_min
print*, "alpha" , pl_fit1%alpha
print*, "ks_stats" , pl_fit1%ks_stats( dados_inteiros )
print*, "Possíveis mínimos de KS:" , pl_fit1%ks_statistics_arr
print*, "Teste com o ajuste real:"
print*, "x_min" , pl_fit2%x_min
print*, "alpha" , pl_fit2%alpha
print*, "ks_stats" , pl_fit2%ks_stats( dados_reais )
print*, "Possíveis mínimos de KS:" , pl_fit2%ks_statistics_arr


end program check