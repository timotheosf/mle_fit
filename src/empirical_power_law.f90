module power_law_mle_fit_mod
    use kinds_mod
    use pl_mod
    implicit none
    
    private
    
    type, extends(power_law) :: empirical_pl
        real(dp) , allocatable :: std_alpha
        real(dp) , allocatable :: x_min_arr( : )
        real(dp) , allocatable :: alpha_arr( : )
        real(dp) , allocatable :: std_alpha_arr( : )
        real(dp) , allocatable :: ks_statistics_arr( : )
        real(dp) , allocatable :: ad_statistics_arr( : )
        real(dp) , allocatable :: t1_statistics_arr( : )
        real(dp) , allocatable :: t2_statistics_arr( : )
    
        type(random_data) :: data

    contains
        
        procedure :: init => start_adjust_parameters
        !procedure :: fast_fit => fast_best_parameters
        !procedure :: mle_alpha => mle_alpha_estimation
        !procedure :: fit => best_parameters

    end type

public :: empirical_pl
    
contains

function mle_alpha_estimation( this ) result( best_alpha )
    class(empirical_pl) , intent(in) :: this
    integer(i4) :: i
    real(dp) :: best_alpha

    if ( .not. this%data%sorted_data ) then
        print*, "Error: data ins't sorted"
        return
    endif

    best_alpha = 0._dp
    select type( a => this%data%arr )    
        type is ( real(dp) )
            do i = 1 , this%data%length
                best_alpha = best_alpha + log( a(i) / this%x_min )
            enddo
            
        type is ( real(sp) )
            do i = 1 , this%data%length
                best_alpha = best_alpha + log( a(i) / this%x_min )
            enddo

        type is ( integer(i4) )
            do i = 1 , this%data%length
                best_alpha = best_alpha + log( real(a(i),dp) / this%x_min )
            enddo

        type is ( integer(i8) )
            do i = 1 , this%data%length
                best_alpha = best_alpha + log( real(a(i),dp) / this%x_min )
            enddo

        class default
            print *, "Error: Stats Test not supported for this data type."
            return
    end select
    best_alpha = 1._dp - this%data%real_length/best_alpha

end function

subroutine start_adjust_parameters( this , r_data , pre_ordering )    
    class(empirical_pl) , intent(inout) :: this
    class(*) , intent(in) :: r_data(:)
    logical, intent(in), optional :: pre_ordering
    !> Receive generic data
    call this%data%receive_data( r_data )
    !> Sort data
    if (present(pre_ordering)) then
        call this%data%sort_data( pre_ordering )
    else
        call this%data%sort_data( )
    endif
end subroutine

end module power_law_mle_fit_mod