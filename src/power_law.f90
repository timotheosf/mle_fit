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

        procedure :: ks_stats => empirical_ks_stats
        procedure :: ad_stats => empirical_ks_ad_weighted_stats
        procedure :: t1_stats => empirical_ks_tail1_weighted_stats
        procedure :: t2_stats => empirical_ks_tail2_weighted_stats

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

function empirical_ks_stats( this , r_data ) result( ks_statistics )
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

function empirical_ks_ad_weighted_stats( this , r_data ) result( w_statistics )
    class(power_law) , intent(in) :: this
    class(random_data) , intent(in) :: r_data
    integer(i4) :: i
    real(dp) :: KS_plus , KS_minus , w_statistics , w
    real(dp) :: real_data_length , current_cdf
    if ( .not. r_data%sorted_data ) then
        print*, "Error: data ins't sorted"
        return
    endif

     select type( a => r_data%arr )
        
        type is ( real(dp) )
            do i = 1 , r_data%length
                current_cdf = this%cdf( a(i) ) 
                w = sqrt( ( current_cdf * ( 1._dp - current_cdf) + epsilon(w) ) )
                KS_plus = max( KS_plus , (real(i,dp)/real_data_length - current_cdf)/w )
                KS_minus = max( KS_minus , (current_cdf-(real(i,dp)-1.d0)/real_data_length)/w )
            enddo
            
        type is ( real(sp) )
            do i = 1 , r_data%length
                current_cdf = this%cdf( a(i) ) 
                w = sqrt( ( current_cdf * ( 1._dp - current_cdf) + epsilon(w) ) )
                KS_plus = max( KS_plus , (real(i,dp)/real_data_length - current_cdf)/w )
                KS_minus = max( KS_minus , (current_cdf-(real(i,dp)-1.d0)/real_data_length)/w )
            enddo

        type is ( integer(i4) )
            do i = 1 , r_data%length
                current_cdf = this%cdf( a(i) ) 
                w = sqrt( ( current_cdf * ( 1._dp - current_cdf) + epsilon(w) ) )
                KS_plus = max( KS_plus , (real(i,dp)/real_data_length - current_cdf)/w )
                KS_minus = max( KS_minus , (current_cdf-(real(i,dp)-1.d0)/real_data_length)/w )
            enddo

        type is ( integer(i8) )
            do i = 1 , r_data%length
                current_cdf = this%cdf( a(i) ) 
                w = sqrt( ( current_cdf * ( 1._dp - current_cdf) + epsilon(w) ) )
                KS_plus = max( KS_plus , (real(i,dp)/real_data_length - current_cdf)/w )
                KS_minus = max( KS_minus , (current_cdf-(real(i,dp)-1.d0)/real_data_length)/w )
            enddo

        class default
            print *, "Error: Stats Test not supported for this data type."
            return
    end select
    w_statistics = max(KS_plus, KS_minus)
end function

function empirical_ks_tail1_weighted_stats( this , r_data ) result( w_statistics )
    class(power_law) , intent(in) :: this
    class(random_data) , intent(in) :: r_data
    integer(i4) :: i
    real(dp) :: KS_plus , KS_minus , w_statistics , w
    real(dp) :: real_data_length , current_cdf
    if ( .not. r_data%sorted_data ) then
        print*, "Error: data ins't sorted"
        return
    endif

     select type( a => r_data%arr )
        
        type is ( real(dp) )
            do i = 1 , r_data%length
                current_cdf = this%cdf( a(i) ) 
                w = sqrt( ( current_cdf + epsilon(w) ) )
                KS_plus = max( KS_plus , (real(i,dp)/real_data_length - current_cdf)/w )
                KS_minus = max( KS_minus , (current_cdf-(real(i,dp)-1.d0)/real_data_length)/w )
            enddo
            
        type is ( real(sp) )
            do i = 1 , r_data%length
                current_cdf = this%cdf( a(i) ) 
                w = sqrt( ( current_cdf + epsilon(w) ) )
                KS_plus = max( KS_plus , (real(i,dp)/real_data_length - current_cdf)/w )
                KS_minus = max( KS_minus , (current_cdf-(real(i,dp)-1.d0)/real_data_length)/w )
            enddo

        type is ( integer(i4) )
            do i = 1 , r_data%length
                current_cdf = this%cdf( a(i) ) 
                w = sqrt( ( current_cdf + epsilon(w) ) )
                KS_plus = max( KS_plus , (real(i,dp)/real_data_length - current_cdf)/w )
                KS_minus = max( KS_minus , (current_cdf-(real(i,dp)-1.d0)/real_data_length)/w )
            enddo

        type is ( integer(i8) )
            do i = 1 , r_data%length
                current_cdf = this%cdf( a(i) ) 
                w = sqrt( ( current_cdf + epsilon(w) ) )
                KS_plus = max( KS_plus , (real(i,dp)/real_data_length - current_cdf)/w )
                KS_minus = max( KS_minus , (current_cdf-(real(i,dp)-1.d0)/real_data_length)/w )
            enddo

        class default
            print *, "Error: Stats Test not supported for this data type."
            return
    end select
    w_statistics = max(KS_plus, KS_minus)
end function

function empirical_ks_tail2_weighted_stats( this , r_data ) result( w_statistics )
    class(power_law) , intent(in) :: this
    class(random_data) , intent(in) :: r_data
    integer(i4) :: i
    real(dp) :: KS_plus , KS_minus , w_statistics , w
    real(dp) :: real_data_length , current_cdf
    if ( .not. r_data%sorted_data ) then
        print*, "Error: data ins't sorted"
        return
    endif

     select type( a => r_data%arr )
        
        type is ( real(dp) )
            do i = 1 , r_data%length
                current_cdf = this%cdf( a(i) ) 
                w = ( current_cdf + epsilon(w) )
                KS_plus = max( KS_plus , (real(i,dp)/real_data_length - current_cdf)/w )
                KS_minus = max( KS_minus , (current_cdf-(real(i,dp)-1.d0)/real_data_length)/w )
            enddo
            
        type is ( real(sp) )
            do i = 1 , r_data%length
                current_cdf = this%cdf( a(i) ) 
                w = ( current_cdf + epsilon(w) )
                KS_plus = max( KS_plus , (real(i,dp)/real_data_length - current_cdf)/w )
                KS_minus = max( KS_minus , (current_cdf-(real(i,dp)-1.d0)/real_data_length)/w )
            enddo

        type is ( integer(i4) )
            do i = 1 , r_data%length
                current_cdf = this%cdf( a(i) ) 
                w = ( current_cdf + epsilon(w) )
                KS_plus = max( KS_plus , (real(i,dp)/real_data_length - current_cdf)/w )
                KS_minus = max( KS_minus , (current_cdf-(real(i,dp)-1.d0)/real_data_length)/w )
            enddo

        type is ( integer(i8) )
            do i = 1 , r_data%length
                current_cdf = this%cdf( a(i) ) 
                w = ( current_cdf + epsilon(w) )
                KS_plus = max( KS_plus , (real(i,dp)/real_data_length - current_cdf)/w )
                KS_minus = max( KS_minus , (current_cdf-(real(i,dp)-1.d0)/real_data_length)/w )
            enddo

        class default
            print *, "Error: Stats Test not supported for this data type."
            return
    end select
    w_statistics = max(KS_plus, KS_minus)
end function

end module pl_mod