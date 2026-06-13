module mle_fit_mod
  use mle_kinds_mod
  use rndgen_mod
  use pl_mod
  implicit none

  
  private
type, extends(power_law) :: empirical_pl
        real(kind=dp) , allocatable :: x_min_arr( : )
        real(kind=dp) , allocatable :: alpha_arr( : )
        real(kind=dp) , allocatable :: std_alpha_arr( : )
        real(kind=dp) , allocatable :: ks_statistics_arr( : )
    contains
        procedure :: fit_int => fitting_int_power_law
        procedure :: fit_real => fitting_real_power_law

        generic :: fit => fit_int , fit_real 

end type

public :: empirical_pl
    
contains

subroutine power_law_real_alpha_mle_adjust( data_array , x_min , alpha )
    real(kind=dp) , intent(in) :: data_array(:) , x_min
    real(kind=dp) , intent(out) :: alpha
    integer(i4) :: i , data_length
    real(kind=dp) :: sum
    sum = 0.d0 ; data_length = size(data_array)
    do i=1,data_length
        sum = sum + log(data_array(i)/x_min)
    enddo
    alpha = 1.d0 + real(data_length,kind=dp)/sum
endsubroutine

subroutine power_law_int_alpha_mle_adjust( data_array , x_min , alpha )
    integer(kind=i4) , intent(in) :: data_array(:)
    real(kind=dp) , intent(in) :: x_min
    real(kind=dp) , intent(out) :: alpha
    integer(i4) :: i , data_length
    real(kind=dp) :: sum
    sum = 0.d0 ; data_length = size(data_array)
    do i=1,data_length
        sum = sum + log(real(data_array(i),dp)/(x_min-0.5))
    enddo
    alpha = 1.d0 + real(data_length,kind=dp)/sum
endsubroutine

subroutine mle_fitting_int( data_array, array_x_min, array_alpha, &
                            array_std_alpha, array_ks_statistics, &
                            best_x_min, best_alpha )
    
    integer(kind=i4) , intent(in) :: data_array(:)
    real(kind=dp) , allocatable , intent(out) :: array_x_min(:) , array_alpha(:) , array_std_alpha(:) , array_ks_statistics(:)
    real(kind=dp) , intent(out) :: best_x_min , best_alpha
    
    type(power_law) :: pl
    integer(kind=i4) , allocatable :: tail_data(:)
    integer(kind=i4) :: i , data_length
    real(kind=dp) :: ks_min , x_min , alpha
    ks_min = huge(ks_min)
    data_length = size( data_array )
    do i = 1 , data_length
        allocate(tail_data(i:data_length))
        tail_data = data_array(i:data_length)
        
        x_min=real(data_array(i),dp)
        call power_law_int_alpha_mle_adjust( tail_data , x_min , alpha )

        pl%x_min = x_min ; pl%alpha = alpha
        if ( ks_min > pl%ks_stats( tail_data ) ) then
            ks_min = pl%ks_stats( tail_data )
            best_x_min = x_min
            best_alpha = alpha
            if (allocated(array_alpha)) then
                array_x_min = [array_x_min, x_min] 
                array_alpha = [array_alpha, alpha]
                array_std_alpha = [ array_std_alpha , (alpha-1.d0)/real(size(tail_data))]
                array_ks_statistics = [array_ks_statistics, ks_min]
            else
                array_x_min = [x_min] 
                array_alpha = [alpha] 
                array_std_alpha = [(alpha-1.d0)/real(size(tail_data))]
                array_ks_statistics = [ks_min]
            endif
        endif
        deallocate( tail_data )
    enddo
endsubroutine

subroutine mle_fitting_real( data_array, array_x_min, array_alpha, &
                            array_std_alpha, array_ks_statistics, &
                            best_x_min, best_alpha )
    
    real(kind=dp) , intent(in) :: data_array(:)
    real(kind=dp) , allocatable , intent(out) :: array_x_min(:) , array_alpha(:) , array_std_alpha(:) , array_ks_statistics(:)
    real(kind=dp) , intent(out) :: best_x_min , best_alpha
    
    type(power_law) :: pl
    real(kind=dp) , allocatable :: tail_data(:)
    integer(kind=i4) :: i , data_length
    real(kind=dp) :: ks_min , x_min , alpha
    ks_min = huge(ks_min)
    data_length = size( data_array )
    do i = 1 , data_length
        allocate(tail_data(i:data_length))
        tail_data = data_array(i:data_length)
        
        x_min=data_array(i)
        call power_law_real_alpha_mle_adjust( tail_data , x_min , alpha )

        pl%x_min = x_min ; pl%alpha = alpha
        if ( ks_min > pl%ks_stats( tail_data ) ) then
            ks_min = pl%ks_stats( tail_data )
            best_x_min = x_min
            best_alpha = alpha
            if (allocated(array_alpha)) then
                array_x_min = [array_x_min, x_min] 
                array_alpha = [array_alpha, alpha]
                array_std_alpha = [ array_std_alpha , (alpha-1.d0)/real(size(tail_data))]
                array_ks_statistics = [array_ks_statistics, ks_min]
            else
                array_x_min = [x_min] 
                array_alpha = [alpha] 
                array_std_alpha = [(alpha-1.d0)/real(size(tail_data))]
                array_ks_statistics = [ks_min]
            endif
        endif
        deallocate( tail_data )
    enddo
endsubroutine

subroutine fitting_int_power_law( this , data_array)
    class(empirical_pl) :: this
    integer(kind=i4) , intent(in) :: data_array(:)
    real(kind=dp) , allocatable :: array_x_min(:) , array_alpha(:) , array_std_alpha(:) , array_ks_statistics(:)
    real(kind=dp) :: best_x_min , best_alpha
    
    call mle_fitting_int(data_array,array_x_min,array_alpha,array_std_alpha,array_ks_statistics,best_x_min,best_alpha)
        this%x_min_arr = array_x_min
        this%alpha_arr = array_alpha
        this%std_alpha_arr = array_std_alpha
        this%ks_statistics_arr = array_ks_statistics
        this%x_min = best_x_min
        this%alpha = best_alpha
endsubroutine

subroutine fitting_real_power_law( this , data_array)
    class(empirical_pl) :: this
    real(kind=dp) , intent(in) :: data_array(:)
    real(kind=dp) , allocatable :: array_x_min(:) , array_alpha(:) , array_std_alpha(:) , array_ks_statistics(:)
    real(kind=dp) :: best_x_min , best_alpha
    
    call mle_fitting_real(data_array,array_x_min,array_alpha,array_std_alpha,array_ks_statistics,best_x_min,best_alpha)
        this%x_min_arr = array_x_min
        this%alpha_arr = array_alpha
        this%std_alpha_arr = array_std_alpha
        this%ks_statistics_arr = array_ks_statistics
        this%x_min = best_x_min
        this%alpha = best_alpha
endsubroutine


end module mle_fit_mod