module pl_mod
  use mle_kinds_mod
  use rndgen_mod
  !use mle_fit_mod
  implicit none

  private

    type :: power_law
    real(kind=dp) :: x_min , alpha
     
    contains
        procedure :: pdf => power_law_pdf
        procedure :: cdf => power_law_cdf
        procedure :: ccdf => power_law_complementary_cdf
        procedure :: i_rnd => int_power_law_random_number_generator
        procedure :: r_rnd => real_power_law_random_number_generator
        procedure :: ks_stats_int => empirical_kolmogorov_stats_measure_int
        procedure :: ks_stats_real => empirical_kolmogorov_stats_measure_real
        !procedure :: mle => mle_adjust

        generic :: ks_stats => ks_stats_int, ks_stats_real

    end type

    public :: power_law
    
contains

function power_law_cdf( this , x ) result( cdf )
    class(power_law) , intent(in) :: this
    real(kind=dp) , intent(in) :: x
    real(kind=dp) cdf
    if ( x < this%x_min) then
        cdf = 0.0d0
    else
        cdf = 1.0d0 - (this%x_min / x)**(this%alpha - 1.0d0)
    end if
end function

function power_law_complementary_cdf( this , x ) result( ccdf )
    class(power_law) , intent(in) :: this
    real(kind=dp) , intent(in) :: x
    real(kind=dp) ccdf
    if ( x < this%x_min) then
        ccdf = 1.0d0
    else
        ccdf = (this%x_min / x)**(this%alpha - 1.0d0)
    end if
end function

function power_law_pdf( this , x ) result( pdf )
    class(power_law) , intent(in) :: this
    real(kind=dp) , intent(in) :: x
    real(kind=dp) pdf
    if ( x < this%x_min) then
        pdf = 0.0d0
    else
        pdf = ( this%alpha - 1.0d0 ) * this%x_min**( this%alpha - 1.0d0 ) * x**( -this%alpha )
    end if
end function

function real_power_law_random_number_generator( this , rnd_gen ) result( pl_random_number )
    class(power_law) , intent(in) :: this
    class(rndgen) , intent(in) :: rnd_gen
    real(kind=dp) pl_random_number 
    real(kind=dp) :: r_r
        !> r_r \in [0,1)
        r_r = rnd_gen%rnd()
        !> CDF invertion , x_r in [x_min,\infty) real PL 
        pl_random_number = this%x_min*(1.d0-r_r)**(-1.d0/(this%alpha-1.d0))
endfunction

function int_power_law_random_number_generator( this , rnd_gen ) result( pl_random_number )
    class(power_law) , intent(in) :: this
    class(rndgen) , intent(in) :: rnd_gen
    integer(kind=i4) pl_random_number 
    real(kind=dp) :: r_r
        !> r_r \in [0,1)
        r_r = rnd_gen%rnd()
        !> CDF invertion , x_r in [x_min,\infty) int PL 
        pl_random_number = floor( (this%x_min-0.5)*(1.d0-r_r)**(-1.d0/(this%alpha-1.d0))+0.5 )
endfunction

function empirical_kolmogorov_stats_measure_real( this , data_array ) result( ks_statistics )
    class(power_law) , intent(in) :: this
    real(kind=dp) , intent(in) :: data_array(:)
    real(kind=dp) KS_plus , KS_minus , ks_statistics
    integer(kind=i4) i , data_length
    KS_plus = 0.d0 ; KS_minus = 0.d0
    data_length = size(data_array)
    do i = 1,data_length
        KS_plus = max( KS_plus , (real(i,dp)/real(data_length,dp) - this%cdf(data_array(i))))
        KS_minus = max( KS_minus , (this%cdf(data_array(i))-(real(i,dp)-1.d0)/real(data_length,dp)))
    enddo
    ks_statistics = max(KS_plus,KS_minus)
endfunction

function empirical_kolmogorov_stats_measure_int( this , data_array ) result( ks_statistics )
    class(power_law) , intent(in) :: this
    integer(kind=i4) , intent(in) :: data_array(:)
    real(kind=dp) KS_plus , KS_minus , ks_statistics
    integer(kind=i4) i , data_length
    KS_plus = 0.d0 ; KS_minus = 0.d0
    data_length = size(data_array)
    do i = 1,data_length
        KS_plus = max( KS_plus , (real(i,dp)/real(data_length,dp) - this%cdf(real(data_array(i),dp))))
        KS_minus = max( KS_minus , (this%cdf(real(data_array(i),dp))-(real(i,dp)-1.d0)/real(data_length,dp)))
    enddo
    ks_statistics = max(KS_plus,KS_minus)
endfunction

end module pl_mod
