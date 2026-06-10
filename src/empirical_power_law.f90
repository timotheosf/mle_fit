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
        procedure :: fast_fit => fast_find_best_parameters
        procedure :: wfast_fit => weight_fast_find_best_parameters
        !procedure :: complete_fit => fast_find_best_parameters

    end type

public :: empirical_pl
    
contains

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

subroutine fast_find_best_parameters( this , r_data , xmin , alpha , std_alpha , ks , lambda_in )
    class(empirical_pl) , intent(inout) :: this
    class(*) , intent(in) , optional :: r_data(:)
    real(dp) , intent(in) , optional :: lambda_in
    real(dp) , intent(out) :: alpha , xmin
    real(dp) , intent(out) , optional :: ks , std_alpha
    real(dp) , parameter :: eps = epsilon(1.0_dp)
    real(dp) , allocatable :: sum_log_x(:) , ks_plus_arr(:) , ks_minus_arr(:) , current_cdf(:) , seq(:)
    real(dp) :: ks_statistics , current_ks, candidate_xmin , candidate_alpha, log_sum , N_tail , candidate_std_alpha , offset , lambda
    integer(i4) :: i , N , n_tail_int

    !> Initializing: receive and sort the generic data
    if (present(r_data)) then
        call this%init( r_data )
    else if (.not. allocated(this%data%arr)) then
        error stop "Error: data not present in Empirical PL class"
    endif
    N = this%data%len
    if (this%data%data_is_discrete) then
        offset = 0.5_dp
    else
        offset = 0.0_dp
    endif
    if (present(lambda_in)) then
        lambda = lambda_in
    else
        lambda = 0.0_dp
    endif
    

    allocate( sum_log_x(N) , ks_plus_arr(N) , ks_minus_arr(N) , seq(N) , current_cdf(N) )
    sum_log_x(N) = log(this%data%arr(N)) !> Pre-calculation tail logarithm data
    seq(N) = real(N,dp) !> Pre-calculation seq
    do i = N-1 , 1 , -1
        !> O(N) cumulative summa for alpha computation in O(1) (avoiging O(N^2) calculation of alpha)
        sum_log_x(i) = sum_log_x(i+1) + log(this%data%arr(i))
        !
        seq(i) = real(i,dp)
    enddo
    
    ks_statistics = huge(1.0_dp) !> Starting ks stats

    !> Main loop
    do i = 1 , N-1
        if ( i > 1 ) then !> Avoid repeated x_min candidates
            if ( abs(candidate_xmin - this%data%arr(i)) < eps ) cycle
        endif
        !> New candidates values and tail length
        candidate_xmin = this%data%arr(i)
        n_tail_int = N - i + 1            !> Integer aux for array slicing
        N_tail = real(n_tail_int, dp)     !> Real aux for calculations

        !> O(1) calculation of alpha
        log_sum = sum_log_x(i) - N_tail*log( candidate_xmin-offset )
        candidate_alpha = 1.0_dp + N_tail/log_sum
         
        !> O(N_tail) calculation of ks_stats
        call this%update_internals( candidate_xmin , candidate_alpha )
        current_cdf( 1:n_tail_int ) = this%cdf_dp( this%data%arr( i:N ) )
        ks_plus_arr( 1:n_tail_int ) = (seq( 1:n_tail_int ) / N_tail) - current_cdf( 1:n_tail_int )
        ks_minus_arr( 1:n_tail_int ) = current_cdf( 1:n_tail_int ) - ((seq( 1:n_tail_int ) - 1.0_dp) / N_tail)
        current_ks = max( maxval(ks_minus_arr( 1:n_tail_int )), maxval(ks_plus_arr( 1:n_tail_int )) )
        current_ks = current_ks - lambda*((N_tail/real(N))**2)
        if ( current_ks <= ks_statistics ) then
            alpha = candidate_alpha     !> Update alpha value
            xmin = candidate_xmin       !> Update x_min
            candidate_std_alpha = (candidate_alpha-1.0_dp)/sqrt( N_tail ) !> Update std alpha
            ks_statistics = current_ks  !> Update ks_stats
        endif
    enddo
    call this%update_internals( xmin , alpha )
    if (present(ks)) ks = ks_statistics
    if (present(std_alpha)) std_alpha = candidate_std_alpha

    deallocate( sum_log_x, ks_plus_arr, ks_minus_arr, current_cdf, seq )
end subroutine

subroutine weight_fast_find_best_parameters( this , r_data , xmin , alpha , std_alpha , ks , lambda_in )
    class(empirical_pl) , intent(inout) :: this
    class(*) , intent(in) , optional :: r_data(:)
    real(dp) , intent(in) , optional :: lambda_in
    real(dp) , intent(out) :: alpha , xmin
    real(dp) , intent(out) , optional :: ks , std_alpha
    real(dp) , parameter :: eps = epsilon(1.0_dp)
    real(dp) , allocatable :: sum_log_x(:) , ks_plus_arr(:) , ks_minus_arr(:) , current_cdf(:) , seq(:) , w(:)
    real(dp) :: ks_statistics , current_ks, candidate_xmin , candidate_alpha, log_sum , N_tail , candidate_std_alpha , offset , lambda
    integer(i4) :: i , N , n_tail_int

    !> Initializing: receive and sort the generic data
    if (present(r_data)) then
        call this%init( r_data )
    else if (.not. allocated(this%data%arr)) then
        error stop "Error: data not present in Empirical PL class"
    endif
    N = this%data%len
    if (this%data%data_is_discrete) then
        offset = 0.5_dp
    else
        offset = 0.0_dp
    endif
    if (present(lambda_in)) then
        lambda = lambda_in
    else
        lambda = 0.0_dp
    endif
    

    allocate( sum_log_x(N) , ks_plus_arr(N) , ks_minus_arr(N) , seq(N) , current_cdf(N) , w(N) )
    sum_log_x(N) = log(this%data%arr(N)) !> Pre-calculation tail logarithm data
    seq(N) = real(N,dp) !> Pre-calculation seq
    do i = N-1 , 1 , -1
        !> O(N) cumulative summa for alpha computation in O(1) (avoiging O(N^2) calculation of alpha)
        sum_log_x(i) = sum_log_x(i+1) + log(this%data%arr(i))
        !
        seq(i) = real(i,dp)
    enddo
    
    ks_statistics = huge(1.0_dp) !> Starting ks stats

    !> Main loop
    do i = 1 , N-1
        if ( i > 1 ) then !> Avoid repeated x_min candidates
            if ( abs(candidate_xmin - this%data%arr(i)) < eps ) cycle
        endif
        !> New candidates values and tail length
        candidate_xmin = this%data%arr(i)
        n_tail_int = N - i + 1            !> Integer aux for array slicing
        N_tail = real(n_tail_int, dp)     !> Real aux for calculations

        !> O(1) calculation of alpha
        log_sum = sum_log_x(i) - N_tail*log( candidate_xmin-offset )
        candidate_alpha = 1.0_dp + N_tail/log_sum
         
        !> O(N_tail) calculation of ks_stats
        call this%update_internals( candidate_xmin , candidate_alpha )
        current_cdf( 1:n_tail_int ) = this%cdf_dp( this%data%arr( i:N ) )
        w( 1:n_tail_int ) = 1._dp/sqrt( (current_cdf(1:n_tail_int)*(1._dp-current_cdf(1:n_tail_int))+eps) )
        ks_plus_arr( 1:n_tail_int ) = ((seq( 1:n_tail_int ) / N_tail) - current_cdf( 1:n_tail_int ))*w( 1:n_tail_int )
        ks_minus_arr( 1:n_tail_int ) = (current_cdf( 1:n_tail_int ) - ((seq( 1:n_tail_int ) - 1.0_dp) / N_tail))*w( 1:n_tail_int )
        current_ks = max( maxval(ks_minus_arr( 1:n_tail_int )), maxval(ks_plus_arr( 1:n_tail_int )) )
        current_ks = current_ks - lambda*((N_tail/real(N))**2)
        if ( current_ks <= ks_statistics ) then
            alpha = candidate_alpha     !> Update alpha value
            xmin = candidate_xmin       !> Update x_min
            candidate_std_alpha = (candidate_alpha-1.0_dp)/sqrt( N_tail ) !> Update std alpha
            ks_statistics = current_ks  !> Update ks_stats
        endif
    enddo
    call this%update_internals( xmin , alpha )
    if (present(ks)) ks = ks_statistics
    if (present(std_alpha)) std_alpha = candidate_std_alpha

    deallocate( sum_log_x, ks_plus_arr, ks_minus_arr, current_cdf, seq )
end subroutine


end module power_law_mle_fit_mod