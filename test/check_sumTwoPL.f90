program check_convolutions
use mle_kinds_mod
use mle_fit_mod     
use rndgen_mod
implicit none
integer(i4), parameter :: rnd_num = 100000
type(rndgen) :: gen_1 , gen_2
real(dp) :: z_r(1:rnd_num) , x_r(1:rnd_num) , y_r(1:rnd_num)
integer :: j
type(mle_t) :: engine
type(power_law_t) :: pl

call gen_1%init( 12940871 )
call gen_2%init( 80976324 )

do j = 1, rnd_num
        x_r(j) = gen_1%rnd()
        y_r(j) = gen_2%rnd()
end do

x_r = 5._dp * (1.d0 - x_r)**(-1.d0 / (2.2 - 1.d0))
y_r = 5._dp * (1.d0 - y_r)**(-1.d0 / (4.5 - 1.d0))
z_r = x_r + y_r

call pl%start_pl()
call engine%greed_fit( r_data=z_r , dist=pl , look_whole=.true. )
call engine%report()

end program check_convolutions