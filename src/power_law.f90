module pl_mod
  use kinds_mod
  use rndgen_mod
  !use mle_fit_mod
  implicit none

  private

    type :: power_law
    !> Distribution parameters
    real(kind=dp) :: x_min , alpha
    !> Intrinsic random number generator
    type(rndgen) :: rnd_num_gen
    contains
        !> Core functions in dp
        procedure, private :: core_pdf_dp
        procedure, private :: core_cdf_dp
        !> Gerneric interface for pdf
        procedure :: pdf_sp; procedure :: pdf_dp
        procedure :: pdf_i1; procedure :: pdf_i2
        procedure :: pdf_i4; procedure :: pdf_i8
        generic :: pdf => pdf_sp, pdf_dp, pdf_i1, pdf_i2, pdf_i4, pdf_i8
        !> Gerneric interface for cdf
        procedure :: cdf_sp; procedure :: cdf_dp
        procedure :: cdf_i1; procedure :: cdf_i2
        procedure :: cdf_i4; procedure :: cdf_i8
        generic :: cdf => cdf_sp, cdf_dp, cdf_i1, cdf_i2, cdf_i4, cdf_i8
        !> Gerneric interface for ccdf
        procedure :: ccdf_sp; procedure :: ccdf_dp
        procedure :: ccdf_i1; procedure :: ccdf_i2
        procedure :: ccdf_i4; procedure :: ccdf_i8
        generic :: ccdf => ccdf_sp, ccdf_dp, ccdf_i1, ccdf_i2, ccdf_i4, ccdf_i8
        
        !> PL random number generator
        procedure :: real_rnd => single_real_pl_rnd_number
        procedure :: int_rnd => single_int_pl_rnd_number
        procedure :: real_array_rnd => real_array_pl_rnd_number
        procedure :: int_array_rnd => int_array_pl_rnd_number


        procedure :: ks_stats_int => empirical_kolmogorov_stats_measure_int
        procedure :: ks_stats_real => empirical_kolmogorov_stats_measure_real
        procedure :: ad_stats => empirical_adw_kolmogorov_stats_measure_int !empirical_anderson_darling_weighted_kolmogorov_stats_measure_int !Anderson–Darling stats
        !procedure :: ad_stats_real => empirical_anderson_darling_weighted_kolmogorov_stats_measure_real !Anderson–Darling stats
        procedure :: tail_1_stats => empirical_t2w_kolmogorov_stats_measure_int !empirical_tail_1_weighted_kolmogorov_stats_measure_int ! weight in tail: 1/F(x)
        !procedure :: tail_1_stats_real => empirical_tail_1_weighted_kolmogorov_stats_measure_real ! weight in tail: 1/F(x)
        procedure :: tail_2_stats => empirical_t2w_kolmogorov_stats_measure_int !empirical_tail_2_weighted_kolmogorov_stats_measure_int ! weight in tail: 1/sqrt(F(x))
        !procedure :: tail_2_stats_real => empirical_tail_2_weighted_kolmogorov_stats_measure_real ! weight in tail: 1/sqrt(F(x))
        
        !procedure :: mle => mle_adjust

        generic :: ks_stats => ks_stats_int, ks_stats_real
        !generic :: ad_stats => ad_stats_int, ad_stats_real  
        !generic :: tail_1_stats => tail_1_stats_int, tail_1_stats_real  
        !generic :: tail_2_stats => tail_2_stats_int, tail_2_stats_real  

    end type

    public :: power_law
      
contains
  
elemental function core_pdf_dp( this , x ) result( res )
    class(power_law) , intent(in) :: this
    real(dp) , intent(in) :: x
    real(dp) :: res
    if ( x < this%x_min ) then
        res = 0._dp
    else
        res = ( this%alpha - 1.0_dp ) * this%x_min**( this%alpha - 1.0_dp ) * x**( -this%alpha )
    endif
end function

elemental function core_cdf_dp( this , x ) result( res )
    class(power_law) , intent(in) :: this
    real(dp) , intent(in) :: x
    real(dp) :: res
    if ( x < this%x_min ) then
        res = 0._dp
    else
        res = 1.0_dp - (this%x_min / x )**(this%alpha - 1.0_dp)
    endif
end function

elemental function pdf_dp( this, x ) result(res)
    class(power_law), intent(in) :: this; real(dp), intent(in) :: x; real(dp) :: res
    res = this%core_pdf_dp(x)
end function
elemental function pdf_sp( this, x ) result(res)
    class(power_law), intent(in) :: this; real(sp), intent(in) :: x; real(dp) :: res
    res = this%core_pdf_dp(real(x, dp))
end function
elemental function pdf_i1( this, x ) result(res)
    class(power_law), intent(in) :: this; integer(i1), intent(in) :: x; real(dp) :: res
    res = this%core_pdf_dp(real(x, dp))
end function
elemental function pdf_i2( this, x ) result(res)
    class(power_law), intent(in) :: this; integer(i2), intent(in) :: x; real(dp) :: res
    res = this%core_pdf_dp(real(x, dp))
end function
elemental function pdf_i4( this, x ) result(res)
    class(power_law), intent(in) :: this; integer(i4), intent(in) :: x; real(dp) :: res
    res = this%core_pdf_dp(real(x, dp))
end function
elemental function pdf_i8( this, x ) result(res)
    class(power_law), intent(in) :: this; integer(i8), intent(in) :: x; real(dp) :: res
    res = this%core_pdf_dp(real(x, dp))
end function

elemental function cdf_dp( this, x ) result(res)
    class(power_law), intent(in) :: this; real(dp), intent(in) :: x; real(dp) :: res
    res = this%core_cdf_dp(x)
end function
elemental function cdf_sp( this, x ) result(res)
    class(power_law), intent(in) :: this; real(sp), intent(in) :: x; real(dp) :: res
    res = this%core_cdf_dp(real(x, dp))
end function
elemental function cdf_i1( this, x ) result(res)
    class(power_law), intent(in) :: this; integer(i1), intent(in) :: x; real(dp) :: res
    res = this%core_cdf_dp(real(x, dp))
end function
elemental function cdf_i2( this, x ) result(res)
    class(power_law), intent(in) :: this; integer(i2), intent(in) :: x; real(dp) :: res
    res = this%core_cdf_dp(real(x, dp))
end function
elemental function cdf_i4( this, x ) result(res)
    class(power_law), intent(in) :: this; integer(i4), intent(in) :: x; real(dp) :: res
    res = this%core_cdf_dp(real(x, dp))
end function
elemental function cdf_i8( this, x ) result(res)
    class(power_law), intent(in) :: this; integer(i8), intent(in) :: x; real(dp) :: res
    res = this%core_cdf_dp(real(x, dp))
end function

elemental function ccdf_dp( this, x ) result(res)
    class(power_law), intent(in) :: this; real(dp), intent(in) :: x; real(dp) :: res
    res = 1.0_dp - this%cdf(x)
end function
elemental function ccdf_sp( this, x ) result(res)
    class(power_law), intent(in) :: this; real(sp), intent(in) :: x; real(dp) :: res
    res = 1.0_dp - this%cdf(x)
end function
elemental function ccdf_i1( this, x ) result(res)
    class(power_law), intent(in) :: this; integer(i1), intent(in) :: x; real(dp) :: res
    res = 1.0_dp - this%cdf(x)
end function
elemental function ccdf_i2( this, x ) result(res)
    class(power_law), intent(in) :: this; integer(i2), intent(in) :: x; real(dp) :: res
    res = 1.0_dp - this%cdf(x)
end function
elemental function ccdf_i4( this, x ) result(res)
    class(power_law), intent(in) :: this; integer(i4), intent(in) :: x; real(dp) :: res
    res = 1.0_dp - this%cdf(x)
end function
elemental function ccdf_i8( this, x ) result(res)
    class(power_law), intent(in) :: this; integer(i8), intent(in) :: x; real(dp) :: res
    res = 1.0_dp - this%cdf(x)
end function

function single_real_pl_rnd_number( this ) result( pl_real_number )
    class(power_law) , intent(in) :: this
    real(dp) :: pl_real_number
        pl_real_number = this%x_min*(1.d0-this%rnd_num_gen%rnd())**(-1.d0/(this%alpha-1.d0))
end function

function real_array_pl_rnd_number( this , array_size ) result( pl_real_number )
    class(power_law) , intent(in) :: this
    integer(i4) , intent(in) :: array_size
    real(dp) , allocatable :: pl_real_number( : )
        allocate( pl_real_number(array_size) )
        pl_real_number = this%x_min*(1.d0-this%rnd_num_gen%rnd_array( array_size ))**(-1.d0/(this%alpha-1.d0))
end function

function single_int_pl_rnd_number( this ) result( pl_int_number )
    class(power_law) , intent(in) :: this
    integer(i4) :: pl_int_number
        pl_int_number = floor( (this%x_min-0.5_dp)*(1.d0-this%rnd_num_gen%rnd())**(-1.d0/(this%alpha-1.d0)) + 0.5_dp )
end function

function int_array_pl_rnd_number( this , array_size ) result( pl_int_number )
    class(power_law) , intent(in) :: this
    integer(i4) , intent(in) :: array_size
    integer(i4) , allocatable :: pl_int_number( : )
        allocate( pl_int_number(array_size) )
        pl_int_number = floor( (this%x_min-0.5_dp)*(1.d0-this%rnd_num_gen%rnd_array( array_size ))**(-1.d0/(this%alpha-1.d0)) + 0.5_dp )
end function



function core_ks_dp( this , r_data ) result( ks_statistics )
    class(power_law) , intent(in) :: this
    class(random_data) , intent(in) :: r_data
    integer(i4) :: i
    real(dp) :: KS_plus , KS_minus , ks_statistics
    real(dp) :: real_data_length , current_cdf
    if ( .not. r_data%sorted_data ) then
        print*, "Error: data ins't sorted"
        return
    endif

     select type( a => r_data%arr )
        
        type is ( real(dp) )
            do i = 1 , r_data%length
                !> Agora o compilador sabe que 'a' é real(dp) e aceita chamar o CDF!
                current_cdf = this%cdf( a(i) ) 
                KS_plus  = max( KS_plus , (real(i,dp)/real_data_length - current_cdf) )
                KS_minus = max( KS_minus , (current_cdf - (real(i,dp)-1.0_dp)/real_data_length) )
            enddo
            
        type is ( real(sp) )
            do i = 1 , r_data%length
                current_cdf = this%cdf( a(i) )
                KS_plus  = max( KS_plus , (real(i,dp)/real_data_length - current_cdf) )
                KS_minus = max( KS_minus , (current_cdf - (real(i,dp)-1.0_dp)/real_data_length) )
            enddo

        type is ( integer(i4) )
            do i = 1 , r_data%length
                current_cdf = this%cdf( a(i) )
                KS_plus  = max( KS_plus , (real(i,dp)/real_data_length - current_cdf) )
                KS_minus = max( KS_minus , (current_cdf - (real(i,dp)-1.0_dp)/real_data_length) )
            enddo

        type is ( integer(i8) )
            do i = 1 , r_data%length
                current_cdf = this%cdf( a(i) )
                KS_plus  = max( KS_plus , (real(i,dp)/real_data_length - current_cdf) )
                KS_minus = max( KS_minus , (current_cdf - (real(i,dp)-1.0_dp)/real_data_length) )
            enddo

        class default
            print *, "Error: KS Test not supported for this data type."
            return
    end select

    ks_statistics = max(KS_plus, KS_minus)
end function

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

function empirical_adw_kolmogorov_stats_measure_int( this , data_array) result( w_statistics )
    class(power_law) , intent(in) :: this
    integer(kind=i4) , intent(in) :: data_array(:)
    real(kind=dp) KS_plus , KS_minus , w_statistics , w
    integer(kind=i4) i , data_length
    KS_plus = 0.d0 ; KS_minus = 0.d0
    data_length = size(data_array)
    do i = 1,data_length
        w = sqrt( ( this%cdf(data_array(i)) * ( 1._dp - this%cdf(data_array(i))) + 1.d-20 ) )
        KS_plus = max( KS_plus , (real(i,dp)/real(data_length,dp) - this%cdf(data_array(i)))/w )
        KS_minus = max( KS_minus , (this%cdf(real(data_array(i),dp))-(real(i,dp)-1.d0)/real(data_length,dp))/w )
    enddo
   w_statistics = max(KS_plus,KS_minus)
endfunction

function empirical_t1w_kolmogorov_stats_measure_int( this , data_array) result( w_statistics )
    class(power_law) , intent(in) :: this
    integer(kind=i4) , intent(in) :: data_array(:)
    real(kind=dp) KS_plus , KS_minus , w_statistics , w
    integer(kind=i4) i , data_length
    KS_plus = 0.d0 ; KS_minus = 0.d0
    data_length = size(data_array)
    do i = 1,data_length
        w = ( this%cdf(data_array(i)) + 1.d-20 )
        KS_plus = max( KS_plus , (real(i,dp)/real(data_length,dp) - this%cdf(data_array(i)))/w )
        KS_minus = max( KS_minus , (this%cdf(real(data_array(i),dp))-(real(i,dp)-1.d0)/real(data_length,dp))/w )
    enddo
    w_statistics = max(KS_plus,KS_minus)
endfunction

function empirical_t2w_kolmogorov_stats_measure_int( this , data_array) result( w_statistics )
    class(power_law) , intent(in) :: this
    integer(kind=i4) , intent(in) :: data_array(:)
    real(kind=dp) KS_plus , KS_minus , w_statistics , w
    integer(kind=i4) i , data_length
    KS_plus = 0.d0 ; KS_minus = 0.d0
    data_length = size(data_array)
    do i = 1,data_length
        w = sqrt( ( this%cdf(data_array(i)) + 1.d-20 ) )
        KS_plus = max( KS_plus , (real(i,dp)/real(data_length,dp) - this%cdf(data_array(i)))/w )
        KS_minus = max( KS_minus , (this%cdf(real(data_array(i),dp))-(real(i,dp)-1.d0)/real(data_length,dp))/w )
    enddo
    w_statistics = max(KS_plus,KS_minus)
endfunction

end module pl_mod