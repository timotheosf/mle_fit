module empirical_pl_mod
    use mle_kinds_mod
    use pl_mod
    use rndgen_mod
    implicit none
    
    private
    
    type, extends(power_law) :: empirical_pl    
        type(random_data) :: data

        real(dp) :: std_alpha                   !> Standard deviation of alpha
        real(dp) :: stats                       !> Kolmogorov-Smirnov/Anderson-Dalirng Statistics
        integer(i4) :: n_tail                   !> Empirical distribuiton tail length 
        real(dp) :: goodness_of_fit             !> P-value

        real(dp) , allocatable :: x_min_arr( : )
        real(dp) , allocatable :: alpha_arr( : )
        real(dp) , allocatable :: std_alpha_arr( : )
        real(dp) , allocatable :: stats_arr( : )

        !> Private variables
        logical , private  :: weighted_adjust   !> Flag to save the weighted used in the adjust
        logical , private  :: was_fitted        !> Control flag: was this fitted?
        logical , private  :: was_pvalued       !> Control flag: was this p_valued?
        real(dp) , private :: lambda_used       !> Internal variable

    contains

        procedure :: init => start_adjust_parameters
        procedure :: fast_fit => find_best_parameters
        procedure :: p_value => p_value_test
        procedure :: report => print_report

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
    this%was_pvalued = .FALSE. ; this%was_fitted = .FALSE.
end subroutine

subroutine find_best_parameters( this , r_data , xmin , alpha , std_alpha , ks , lambda_in , use_weight , synth_data_treat_as_discrete , track_history )
    !> This is a long subroutine for fitting power laws parameters in empirical data
    !   It uses a vectorized operation as a fast method to find the stats,
    !   and simplify the alpha_exponent calculation by a O(1) (effective O(N) complexity)
    class(empirical_pl) , intent(inout) :: this
    class(*) , intent(in) , optional :: r_data(:)
    logical  , intent(in) , optional :: use_weight , synth_data_treat_as_discrete , track_history
    real(dp) , intent(in) , optional :: lambda_in
    real(dp) , intent(out) , optional :: alpha , xmin , ks , std_alpha
    real(dp) , parameter :: eps = epsilon(1.0_dp)
    real(dp) , allocatable :: mle_x_min_arr( : ) , mle_alpha_arr( : ) , mle_std_alpha_arr( : ) , mle_stats_arr( : )
    real(dp) , allocatable :: sum_log_x(:) , ks_plus_arr(:) , ks_minus_arr(:) , current_cdf(:) , seq(:) , w(:) , log_x(:)
    real(dp) :: log_xmin, alpha_minus_1
    real(dp) :: ks_statistics , current_ks, candidate_xmin , candidate_alpha, log_sum , N_tail , candidate_std_alpha , offset , lambda , mle_xmin , mle_alpha
    integer(i4) :: i , N , n_tail_int , tail_len , non_zero_idx
    logical  :: apply_weight , save_history
    !--- Initializing ---!
    ! 1. Receive and sort the generic data
    if (present(r_data)) then
        call this%init( r_data )
    else if (.not. allocated(this%data%arr)) then
        error stop "Error: data not present in Empirical PL class"
    endif
    N = this%data%len
    ! 2. Defines the offset based in the data original type
    if (this%data%data_is_discrete) then
        offset = 0.5_dp
    else
        offset = 0.0_dp
    endif
    if (present(synth_data_treat_as_discrete)) then
        if (synth_data_treat_as_discrete) offset = 0.5_dp
    endif
    ! 3. Applies the lambda_in penality for the adjust
    if (present(lambda_in)) then
        lambda = lambda_in
    else
        lambda = 0.0_dp
    endif
    ! 4. Defines the weight, if it's the case
    if (present(use_weight)) then
        apply_weight = use_weight
    else
        apply_weight = .TRUE. 
    endif
    if (apply_weight) then
        allocate(w(N))
    endif
    ! 5. Allocationg variables
    allocate( sum_log_x(N), ks_plus_arr(N), ks_minus_arr(N), seq(N), current_cdf(N), log_x(N) )
    log_x(1:N) = log(this%data%arr(1:N))    !> Optimized CDF evaluation
    sum_log_x(N) = log(this%data%arr(N))    !> Pre-calculation tail logarithm data
    seq(N) = real(N,dp) !> Pre-calculation seq
    do i = N-1 , 1 , -1
        !> O(N) cumulative summa for alpha computation in O(1) (avoiging O(N^2) calculation of alpha)
        sum_log_x(i) = sum_log_x(i+1) + log(this%data%arr(i))
        !> Determining the seq variable
        seq(i) = real(i,dp)
    enddo
    ks_statistics = huge(1.0_dp) !> Starting ks stats
    ! 6. Tracking minima history
    if (present(track_history)) then
        save_history = track_history
    else
        save_history = .FALSE.
    endif
    if (save_history) then
        non_zero_idx = 0
        allocate( mle_x_min_arr(N) , mle_alpha_arr(N) , mle_std_alpha_arr(N) , mle_stats_arr(N) )
    endif
    !--- Main Loop ---!
    mle_main_loop: do i = 1 , N-1
        if ( i > 1 ) then !> Avoid repeated x_min candidates
            if ( abs(candidate_xmin - this%data%arr(i)) < eps ) cycle mle_main_loop
        endif
        !> New candidates values and tail length
        candidate_xmin = this%data%arr(i)
        n_tail_int = N - i + 1            !> Integer aux for array slicing
        N_tail = real(n_tail_int, dp)     !> Real aux for calculations

        !> O(1) calculation of alpha
        log_sum = sum_log_x(i) - N_tail*log( candidate_xmin-offset )
        candidate_alpha = 1.0_dp + N_tail/log_sum
         
        !> Fully vectorized O(N_tail) calculation of ks_stats
        log_xmin = log(candidate_xmin)
        alpha_minus_1 = candidate_alpha - 1.0_dp
        current_cdf( 1:n_tail_int ) = 1.0_dp - exp( alpha_minus_1 * (log_xmin - log_x(i:N)) ) !> Optimized form for CDF evaluation
        if ( apply_weight ) then
            w( 1:n_tail_int ) = 1._dp/sqrt( (current_cdf(1:n_tail_int)*(1._dp-current_cdf(1:n_tail_int))+eps) ) !> Anderson-Darling weight
            ks_plus_arr( 1:n_tail_int ) = ((seq( 1:n_tail_int ) / N_tail) - current_cdf( 1:n_tail_int ))*w( 1:n_tail_int )
            ks_minus_arr( 1:n_tail_int ) = (current_cdf( 1:n_tail_int ) - ((seq( 1:n_tail_int ) - 1.0_dp) / N_tail))*w( 1:n_tail_int )
        else
            ks_plus_arr( 1:n_tail_int ) = (seq( 1:n_tail_int ) / N_tail) - current_cdf( 1:n_tail_int )
            ks_minus_arr( 1:n_tail_int ) = current_cdf( 1:n_tail_int ) - ((seq( 1:n_tail_int ) - 1.0_dp) / N_tail)
        endif
        !> The current stats is update by this functional
        current_ks = max( maxval(ks_minus_arr( 1:n_tail_int )), maxval(ks_plus_arr( 1:n_tail_int )) ) - lambda*((N_tail/real(N))**2)
        if ( current_ks <= ks_statistics ) then
            tail_len = n_tail_int           !> Update the tail_len
            mle_alpha = candidate_alpha     !> Update alpha value
            mle_xmin = candidate_xmin       !> Update x_min
            candidate_std_alpha = (candidate_alpha-1.0_dp)/sqrt( N_tail ) !> Update std alpha
            ks_statistics = current_ks  !> Update ks_stats
            if ( save_history ) then
                non_zero_idx = non_zero_idx + 1
                mle_x_min_arr(non_zero_idx) = mle_xmin
                mle_alpha_arr(non_zero_idx) = mle_alpha
                mle_std_alpha_arr(non_zero_idx) = candidate_std_alpha
                mle_stats_arr(non_zero_idx) = current_ks
            endif
        endif
    enddo mle_main_loop
    !> Update the empirical PL
    call this%update_internals( mle_xmin , mle_alpha )
    this%stats = ks_statistics ; this%n_tail = tail_len ; this%std_alpha = candidate_std_alpha           
    this%weighted_adjust = apply_weight ; this%lambda_used = lambda
    this%was_pvalued = .FALSE. ; this%was_fitted = .TRUE.
    !> In the case if one uses external variables
    if (present(xmin)) xmin = mle_xmin
    if (present(alpha)) alpha = mle_alpha
    if (present(ks)) ks = ks_statistics
    if (present(std_alpha)) std_alpha = candidate_std_alpha
    !> If one wants to track history
    if (save_history) then
        if (allocated(this%x_min_arr)) deallocate(this%x_min_arr, this%alpha_arr, this%std_alpha_arr, this%stats_arr)
        allocate( this%x_min_arr(non_zero_idx) , this%alpha_arr(non_zero_idx) , this%std_alpha_arr(non_zero_idx) , this%stats_arr(non_zero_idx) )
        this%x_min_arr(1:non_zero_idx) = mle_x_min_arr( 1:non_zero_idx )
        this%alpha_arr(1:non_zero_idx) = mle_alpha_arr( 1:non_zero_idx )
        this%std_alpha_arr(1:non_zero_idx) = mle_std_alpha_arr( 1:non_zero_idx )
        this%stats_arr(1:non_zero_idx) = mle_stats_arr( 1:non_zero_idx )
        deallocate( mle_x_min_arr , mle_alpha_arr , mle_std_alpha_arr , mle_stats_arr )
    endif
    !> Deallocate all the internal arrays
    deallocate( sum_log_x, ks_plus_arr, ks_minus_arr, current_cdf, seq , log_x )
    if (apply_weight) deallocate(w)
end subroutine

subroutine p_value_test( this , N_samples , p_value )
    !$ use omp_lib    !> Includes parallel processing
    class(empirical_pl) , intent(inout) :: this
    integer(i4) , intent(in) , optional :: N_samples !> Number of samples is optional; standard = 1000
    real(dp) , intent(out) , optional :: p_value
    type(empirical_pl) :: synth_pl  !> PL fitted in the synthetic data
    type(rndgen) :: thread_gen_1 , thread_gen_2      !> Thread independent random number generator for OpenMP
    real(dp) :: ks , p_tail , synth_ks , real_ks , lambda
    integer(i4) :: i , N_trials , N , time_values(8) , hits , j , base_seed
    integer(i4) :: thread_id  , max_noise_idx , synth_head_size
    real(dp) , allocatable :: random_chooses( : ) , synth_data( : )
    !--- Initializing ---!
    ! 1. Samples used
    if (present(N_samples)) then
        N_trials = N_samples
    else
        N_trials = 1000
    endif

    real_ks = this%stats 
    N = this%data%len
    p_tail = real(this%n_tail,dp)/real(N ,dp)
    hits = 0
    max_noise_idx = max(1_i4, N - this%n_tail) !> A secure index to avoid creating 0-dimentional arrays
    thread_id = 0 !> Secure default value 
    call date_and_time(values=time_values) !> Get a random seed based on time
    base_seed = time_values(8) + 1000*time_values(7) + 60000*time_values(6) + 3600000*time_values(5)

    !--- Parallel proceding ---!
    !$omp parallel private(thread_id, j, i, random_chooses, synth_data, synth_pl, synth_ks, thread_gen_1 , thread_gen_2, synth_head_size)
    !> This clones the variables in each thread
    !$ thread_id = omp_get_thread_num() !> Gets the thread ID used in OpenMP
    call thread_gen_1%init( base_seed + thread_id * 1999  ) !> Initializes each generator by a seed deppending on the thread_id
    call thread_gen_2%init( base_seed + thread_id * 3999 + 104729  ) !> Initializes each generator by a seed deppending on the thread_id
    allocate( random_chooses(N) , synth_data(N) )

    !$omp do reduction(+:hits) !> This creates privates hits variables and safelly summation the result
    !--- Main loop ---!
    sampling_loop: do j = 1 , N_trials
        synth_head_size = 0
        do i=1,N
            random_chooses(i) = thread_gen_1%rnd()  !> Random N-arr    
            synth_data(i) = thread_gen_2%rnd() !> Random N-arr
            if ( random_chooses(i) > p_tail ) synth_head_size = synth_head_size + 1
        enddo
        if (this%data%data_is_discrete) then
            synth_data( synth_head_size+1 : N ) = floor((this%x_min-0.5_dp)*(1.d0-synth_data( synth_head_size+1 : N ))**(-1.d0/(this%alpha-1.d0))+0.5_dp)
        else
            synth_data( synth_head_size+1 : N ) = this%x_min*(1.d0-synth_data( synth_head_size+1 : N ))**(-1.d0/(this%alpha-1.d0))
        endif
        synth_data( 1:synth_head_size ) = this%data%arr( floor(real((max_noise_idx),dp)*synth_data( 1:synth_head_size ),kind=i4) + 1 )
        
        !> Each synthetic data is fitted independented
        call synth_pl%fast_fit( synth_data , ks=synth_ks , lambda_in=this%lambda_used , use_weight=this%weighted_adjust , synth_data_treat_as_discrete=this%data%data_is_discrete )
        if ( real_ks <= synth_ks ) hits = hits + 1
    enddo sampling_loop
    !$omp end do
    !$omp end parallel
    !--- End parallel proceding ---!
    !> Returns the p-value
    if (present(p_value)) p_value = real(hits,dp)/real(N_trials,dp)
    this%goodness_of_fit = real(hits,dp)/real(N_trials,dp) !> Saves the p_value in the internal variable
    this%was_pvalued = .TRUE.
end subroutine

subroutine print_report( this )
    class(empirical_pl) , intent(in) :: this
    if (.not. this%was_fitted ) then
        error stop "Empirical Power-Law is not fitted"
    endif
    print*, "============================"
    print*, "  --Empirical PL Fitted--   "
    print*, "   x_min=", this%x_min
    print*, "   alpha=", this%alpha
    print*, "   std_alpha=", this%std_alpha
    if (this%was_pvalued) print*, "   p_value=", this%goodness_of_fit 
    print*, "============================"
end subroutine

end module empirical_pl_mod