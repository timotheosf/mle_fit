module empirical_pl_mod
    use mle_kinds_mod
    use pl_mod
    use rndgen_mod
    implicit none
    
    private
    
    type, extends(power_law) :: empirical_pl    
        type(random_data) :: data

        real(dp) :: std_alpha                   !> Standard deviation of alpha
        real(dp) :: stats                       !> Used statistics
        real(dp) :: ks                          !> Kolmogorov-Smirnov
        integer(i4) :: n_tail                   !> Empirical distribuiton tail length 
        real(dp) :: goodness_of_fit             !> P-value
        real(dp) :: mle_time                    !> Time costed for fitting
        real(dp) :: hypothesis_time             !> Time costed for hypothesis testing
        real(dp) :: p_value_eps                 !> P-value error

        real(dp) , allocatable , private :: x_min_arr( : )
        real(dp) , allocatable , private :: alpha_arr( : )
        real(dp) , allocatable , private :: std_alpha_arr( : )
        real(dp) , allocatable , private :: stats_arr( : )
        integer(i4) , allocatable , private :: n_tail_arr( : )

        !> Private variables
        logical , private  :: weighted_adjust   !> Flag to save the weighted used in the adjust
        logical , private  :: was_fitted        !> Control flag: was this fitted?
        logical , private  :: was_pvalued       !> Control flag: was this p_valued?
        real(dp) , private :: lambda_used       !> Internal variable
        type(clock_time) , private :: internal_clock !> Internal clock for benchmark
        real(dp) , allocatable , private :: wrk_sort_buffer(:)

    contains

        procedure :: fast_fit => fast_find_best_parameters
        procedure :: lamb_fit => find_best_parameters_with_cost_functional
        procedure :: greed_fit => find_greed_parameters_at_all_cost
        procedure :: p_value => p_value_test
        procedure :: report => print_report

        !> Private procedures
        procedure , private :: init => start_adjust_parameters
        procedure , private :: core_fit => internal_engine_to_find_best_parameters
        procedure , private :: null_hypothesis_test

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
    if (allocated(this%wrk_sort_buffer)) then
        if (present(pre_ordering)) then
            call this%data%sort_data( pre_ordering, work_buffer=this%wrk_sort_buffer )
        else
            call this%data%sort_data( work_buffer=this%wrk_sort_buffer )
        endif
    else
        if (present(pre_ordering)) then
            call this%data%sort_data( pre_ordering )
        else
            call this%data%sort_data( )
        endif
    endif
    !> Initial flags
    this%was_pvalued = .FALSE.
    this%was_fitted = .FALSE.
end subroutine

subroutine fast_find_best_parameters( this , r_data , xmin , alpha , std_alpha , ks , use_weight )
    class(empirical_pl) , intent(inout) :: this
    class(*) , intent(in) , optional :: r_data(:)
    logical  , intent(in) , optional :: use_weight
    real(dp) , intent(out) , optional :: alpha , xmin , ks , std_alpha
    call this%core_fit( r_data=r_data , xmin=xmin , alpha=alpha , std_alpha=std_alpha , ks=ks , use_weight=use_weight )
end subroutine fast_find_best_parameters

subroutine find_best_parameters_with_cost_functional( this , r_data , xmin , alpha , std_alpha , ks , lambda_in , use_weight )
    class(empirical_pl) , intent(inout) :: this
    class(*) , intent(in) , optional :: r_data(:)
    logical  , intent(in) , optional :: use_weight
    real(dp) , intent(in) , optional :: lambda_in
    real(dp) , intent(out) , optional :: alpha , xmin , ks , std_alpha
    real(dp) :: lambda_passed
    lambda_passed = 0.15_dp
    if (present(lambda_in)) lambda_passed = lambda_in
    call this%core_fit( r_data=r_data , xmin=xmin , alpha=alpha , std_alpha=std_alpha , ks=ks , lambda_in=lambda_passed , use_weight=use_weight )
end subroutine find_best_parameters_with_cost_functional

subroutine find_greed_parameters_at_all_cost( this , r_data , greed_xmin , greed_alpha , greed_std_alpha , greed_ks , use_weight )
    class(empirical_pl) , intent(inout) :: this
    class(*) , intent(in) , optional :: r_data(:)
    logical  , intent(in) , optional :: use_weight
    real(dp) , intent(out) , optional :: greed_xmin , greed_alpha , greed_std_alpha , greed_ks
    type(empirical_pl) :: pLaw
    integer(i4) :: idx

    call this%core_fit( r_data=r_data , track_history=.TRUE. , use_weight=use_weight , xmin=greed_xmin , alpha=greed_alpha , std_alpha=greed_std_alpha , ks=greed_ks )
    
    if (size( this%x_min_arr )==1) return
    
    pLaw%data = this%data
    pLaw%was_pvalued = .FALSE. 
    pLaw%was_fitted = .TRUE.
    pLaw%weighted_adjust = this%weighted_adjust 
    pLaw%lambda_used = 0.0_dp
    find_greed_parameteres: do idx = 1 , size( this%x_min_arr )

        call pLaw%update_internals( this%x_min_arr(idx) , this%alpha_arr(idx) )
        pLaw%stats = this%stats_arr(idx) 
        pLaw%std_alpha = this%std_alpha_arr(idx)
        pLaw%n_tail = this%n_tail_arr(idx)
        
        !> Two factor filter pt. 1 -> fast p_value (error=0.3)
        call pLaw%p_value( N_samples=100 )
        if ( pLaw%goodness_of_fit >= 0.07_dp ) then
            call pLaw%p_value( N_samples=2500 ) !> Two factor filter pt. 2 -> p_value (error=0.1)
            if ( pLaw%goodness_of_fit >= 0.1_dp ) then !> pLaw is a powerlaw, in fact
                if (present(greed_xmin)) greed_xmin = pLaw%x_min
                if (present(greed_alpha)) greed_alpha = pLaw%alpha
                if (present(greed_std_alpha)) greed_std_alpha = pLaw%std_alpha
                if (present(greed_ks)) greed_ks = pLaw%ks
                
                call this%update_internals( pLaw%x_min , pLaw%alpha )
                this%std_alpha = pLaw%std_alpha
                this%stats = pLaw%stats
                this%goodness_of_fit = pLaw%goodness_of_fit
                this%n_tail = pLaw%n_tail
                this%was_pvalued = .TRUE.

                call this%internal_clock%stop()
                this%mle_time = this%internal_clock%elapsed

                exit find_greed_parameteres
            endif
        endif
    enddo find_greed_parameteres
    call this%internal_clock%stop()
    this%mle_time = this%internal_clock%elapsed
    call this%p_value( N_samples=2500 )
end subroutine find_greed_parameters_at_all_cost


subroutine p_value_test( this , N_samples , p_value )
    class(empirical_pl) , intent(inout) :: this
    integer(i4) , intent(in) , optional :: N_samples
    real(dp) , intent(out) , optional :: p_value
    call this%null_hypothesis_test( N_samples , p_value=p_value )
end subroutine p_value_test


subroutine internal_engine_to_find_best_parameters( this , r_data , xmin , alpha , std_alpha , ks , lambda_in , use_weight , synth_data_treat_as_discrete , track_history )
    !> This is a long subroutine for fitting power laws parameters in empirical data
    !   It uses a vectorized operation as a fast method to find the stats,
    !   and simplify the alpha_exponent calculation by a O(1) (effective O(N) complexity)
    class(empirical_pl) , intent(inout) :: this
    class(*) , intent(in) , optional :: r_data(:)
    logical  , intent(in) , optional :: use_weight , synth_data_treat_as_discrete , track_history
    real(dp) , intent(in) , optional :: lambda_in
    real(dp) , intent(out) , optional :: alpha , xmin , ks , std_alpha
    real(dp) , parameter :: eps = epsilon(1.0_dp)
    real(dp) :: bin_tolerance
    real(dp) , allocatable :: mle_x_min_arr( : ) , mle_alpha_arr( : ) , mle_std_alpha_arr( : ) , mle_stats_arr( : )
    integer(i4) , allocatable :: mle_n_tail_arr( : )
    real(dp) , allocatable :: sum_log_x(:) , ks_plus_arr(:) , ks_minus_arr(:) , current_cdf(:) , seq(:) , w(:) , log_x(:)
    real(dp) :: log_xmin, alpha_minus_1 , prev_ks, prev_xmin, prev_alpha, prev_std_alpha
    real(dp) :: ks_statistics , current_ks, candidate_xmin , candidate_alpha, log_sum , N_tail , candidate_std_alpha , offset , lambda , mle_xmin , mle_alpha
    integer(i4) :: i , N , n_tail_int , tail_len , non_zero_idx , prev_tail_len
    logical  :: apply_weight , save_history , is_decreasing
    real(dp) , allocatable :: first_idx(:), last_idx(:)
    real(dp) :: S_right, S_left
    real(dp) , allocatable :: first_idx(:)
    !--- Initializing ---!
    ! 0. Start clock benchmark
    call this%internal_clock%start()
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
        bin_tolerance = 0.5_dp
    else
        offset = 0.0_dp
        bin_tolerance = (this%data%arr(N)-this%data%arr(1))/N
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
        apply_weight = .FALSE. 
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
    allocate(first_idx(N), last_idx(N))

! Preenche last_idx (indo de trás pra frente)
    last_idx(N) = real(N, dp)
    do i = N-1, 1, -1
    if (this%data%arr(i) == this%data%arr(i+1)) then
        last_idx(i) = last_idx(i+1)
    else
        last_idx(i) = real(i, dp)
    end if
    end do

    ! Preenche first_idx (indo de frente pra trás)
    first_idx(1) = 1.0_dp
    do i = 2, N
    if (this%data%arr(i) == this%data%arr(i-1)) then
        first_idx(i) = first_idx(i-1)
    else
        first_idx(i) = real(i, dp)
    end if
    end do

    allocate(first_idx(N))
    first_idx(1) = 1.0_dp
    do i = 2, N
        if (this%data%arr(i) == this%data%arr(i-1)) then
            first_idx(i) = first_idx(i-1)
        else
            first_idx(i) = real(i, dp)
        end if
    end do

    ks_statistics = huge(1.0_dp) !> Starting ks stats
    ! 6. Tracking minima history
    if (present(track_history)) then
        save_history = track_history
    else
        save_history = .FALSE.
    endif
    if (save_history) then
        non_zero_idx = 0
        allocate( mle_x_min_arr(N) , mle_alpha_arr(N) , mle_std_alpha_arr(N) , mle_stats_arr(N) , mle_n_tail_arr(N) )
    endif
    !--- Main Loop ---!
    mle_main_loop: do i = 1 , N-1
        if ( i > 1 ) then !> Avoid repeated x_min candidates
            if ( abs(candidate_xmin - this%data%arr(i)) < bin_tolerance ) cycle mle_main_loop
        endif
        !> New candidates values and tail length
        candidate_xmin = this%data%arr(i)
        n_tail_int = N - i + 1            !> Integer aux for array slicing
        N_tail = real(n_tail_int, dp)     !> Real aux for calculations

        !> O(1) calculation of alpha
        log_sum = sum_log_x(i) - N_tail*log( candidate_xmin-offset )
        candidate_alpha = 1.0_dp + N_tail/log_sum

        if ( (candidate_alpha - 1.0_dp) / sqrt( N_tail ) >= 0.1_dp ) then !> Python powerlaw package schizoid flag
            exit mle_main_loop
        endif
         
        !> Fully vectorized O(N_tail) calculation of ks_stats
        if ( this%data%data_is_discrete ) then
            current_cdf( 1:n_tail_int ) = 1.0_dp - &
                    zeta_function(candidate_alpha,this%data%arr(i:N))/zeta_function(candidate_alpha,candidate_xmin)
        else
            log_xmin = log(candidate_xmin)
            alpha_minus_1 = candidate_alpha - 1.0_dp
            current_cdf( 1:n_tail_int ) = 1.0_dp - exp( alpha_minus_1 * (log_xmin - log_x(i:N)) ) !> Optimized form for CDF evaluation
        endif
        if ( apply_weight ) then
            w( 1:n_tail_int ) = 1._dp/sqrt( (current_cdf(1:n_tail_int)*(1._dp-current_cdf(1:n_tail_int))+eps) ) !> Anderson-Darling weight
            ks_plus_arr( 1:n_tail_int ) = ((seq( 1:n_tail_int ) / N_tail) - current_cdf( 1:n_tail_int ))*w( 1:n_tail_int )
            ks_minus_arr( 1:n_tail_int ) = (current_cdf( 1:n_tail_int ) - ((seq( 1:n_tail_int ) - 1.0_dp) / N_tail))*w( 1:n_tail_int )
        else
            ks_plus_arr( 1:n_tail_int )  = current_cdf( 1:n_tail_int ) - (first_idx(i:N) - real(i, dp))/N_tail
            ks_minus_arr( 1:n_tail_int ) = (first_idx(i:N) - real(i, dp))/N_tail - current_cdf( 1:n_tail_int )
            
        endif
        !> The current stats is update by this functional
        current_ks = max( maxval(ks_minus_arr( 1:n_tail_int )), maxval(ks_plus_arr( 1:n_tail_int )) ) - lambda*((N_tail/real(N)))
        if ( current_ks < ks_statistics ) then
            tail_len = n_tail_int           !> Update the tail_len
            mle_alpha = candidate_alpha     !> Update alpha value
            mle_xmin = candidate_xmin       !> Update x_min
            candidate_std_alpha = (candidate_alpha-1.0_dp)/sqrt( N_tail ) !> Update std alpha
            ks_statistics = current_ks  !> Update ks_stats
            if ( save_history ) then
                if ( current_ks < prev_ks ) then
                    is_decreasing = .TRUE.
                else if ( current_ks >= prev_ks .and. is_decreasing ) then
                    non_zero_idx = non_zero_idx + 1
                    mle_x_min_arr(non_zero_idx) = prev_xmin
                    mle_alpha_arr(non_zero_idx) = prev_alpha
                    mle_std_alpha_arr(non_zero_idx) = prev_std_alpha
                    mle_stats_arr(non_zero_idx) = prev_ks
                    mle_n_tail_arr(non_zero_idx) = prev_tail_len
                    is_decreasing = .FALSE.
                endif
                prev_ks = current_ks
                prev_xmin = candidate_xmin
                prev_alpha = candidate_alpha
                prev_std_alpha = (candidate_alpha-1.0_dp)/sqrt( N_tail )
                prev_tail_len = n_tail_int
        endif
        endif
    enddo mle_main_loop
    !> Update the empirical PL
    call this%update_internals( mle_xmin , mle_alpha )
    this%stats = ks_statistics ; this%n_tail = tail_len ; this%std_alpha = candidate_std_alpha           
    this%weighted_adjust = apply_weight ; this%lambda_used = lambda ; this%ks = ks_statistics + lambda*((tail_len/real(N)))
    this%was_pvalued = .FALSE. ; this%was_fitted = .TRUE.
    !> In the case if one uses external variables
    if (present(xmin)) xmin = mle_xmin
    if (present(alpha)) alpha = mle_alpha
    if (present(ks)) ks = ks_statistics
    if (present(std_alpha)) std_alpha = candidate_std_alpha
    !> If one wants to track history
    if (save_history) then
        if (allocated(this%x_min_arr)) deallocate( this%x_min_arr , this%alpha_arr , this%std_alpha_arr , this%stats_arr , this%n_tail_arr )
        allocate( this%x_min_arr(non_zero_idx) , this%alpha_arr(non_zero_idx) , this%std_alpha_arr(non_zero_idx) , this%stats_arr(non_zero_idx) , this%n_tail_arr(non_zero_idx) )
        this%x_min_arr(1:non_zero_idx) = mle_x_min_arr( 1:non_zero_idx )
        this%alpha_arr(1:non_zero_idx) = mle_alpha_arr( 1:non_zero_idx )
        this%std_alpha_arr(1:non_zero_idx) = mle_std_alpha_arr( 1:non_zero_idx )
        this%stats_arr(1:non_zero_idx) = mle_stats_arr( 1:non_zero_idx )
        this%n_tail_arr(1:non_zero_idx) = mle_n_tail_arr(1:non_zero_idx)
        deallocate( mle_x_min_arr , mle_alpha_arr , mle_std_alpha_arr , mle_stats_arr , mle_n_tail_arr )
    endif
    !> Deallocate all the internal arrays
    deallocate( sum_log_x, ks_plus_arr, ks_minus_arr, current_cdf, seq , log_x )
    if (apply_weight) deallocate(w)
    ! End clock benchmark
    call this%internal_clock%stop()
    this%mle_time = this%internal_clock%elapsed
end subroutine internal_engine_to_find_best_parameters

subroutine null_hypothesis_test( this , N_samples , track_penalities , p_value )
    !$ use omp_lib    !> Includes parallel processing
    class(empirical_pl) , intent(inout) :: this
    integer(i4) , intent(in) , optional :: N_samples !> Number of samples is optional; standard = 1000
    logical , intent(in) , optional :: track_penalities
    real(dp) , intent(out) , optional :: p_value
    type(empirical_pl) :: synth_pl  !> PL fitted in the synthetic data
    type(rndgen) :: thread_gen_1 , thread_gen_2      !> Thread independent random number generator for OpenMP
    real(dp) :: p_tail , synth_ks , real_ks
    integer(i4) :: i , N_trials , N , time_values(8) , hits , j , base_seed
    integer(i4) :: thread_id  , max_noise_idx , synth_head_size
    real(dp) , allocatable :: random_chooses( : ) , synth_data( : )
    !--- Initializing ---!
    ! 0. Start time benchmark
    call this%internal_clock%start()
    ! 1. Samples used
    if (present(N_samples)) then
        N_trials = N_samples
    else
        N_trials = 1000
    endif
    real_ks = this%ks
    if (present(track_penalities)) then
        if (track_penalities) real_ks = this%stats
    endif
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
    allocate( random_chooses(N) , synth_data(N) , synth_pl%wrk_sort_buffer(N) )

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
        if (present(track_penalities)) then
            if (track_penalities) then
                call synth_pl%core_fit( r_data=synth_data , ks=synth_ks , lambda_in=this%lambda_used , use_weight=this%weighted_adjust , synth_data_treat_as_discrete=this%data%data_is_discrete )
            else
                call synth_pl%core_fit( r_data=synth_data , ks=synth_ks , use_weight=this%weighted_adjust , synth_data_treat_as_discrete=this%data%data_is_discrete )
            endif
        else
            call synth_pl%core_fit( r_data=synth_data , ks=synth_ks , use_weight=this%weighted_adjust , synth_data_treat_as_discrete=this%data%data_is_discrete )
        endif
            
        if ( real_ks <= synth_ks ) hits = hits + 1
    enddo sampling_loop
    !$omp end do
    !$omp end parallel
    !--- End parallel proceding ---!
    !> Returns the p-value
    if (present(p_value)) p_value = real(hits,dp)/real(N_trials,dp)
    this%goodness_of_fit = real(hits,dp)/real(N_trials,dp) !> Saves the p_value in the internal variable
    this%p_value_eps = 1.0_dp / (2.0_dp * sqrt(real(N_trials, dp))) !> Saves the p_value precision in the internal variable
    this%was_pvalued = .TRUE.
    ! End clock benchmark
    call this%internal_clock%stop()
    this%hypothesis_time = this%internal_clock%elapsed
end subroutine

subroutine print_report( this )
    class(empirical_pl) , intent(in) :: this
    logical :: openMP
    if (.not. this%was_fitted ) then
        error stop "Empirical Power-Law is not fitted"
    endif
    openMP = .FALSE.
    !$ openMP = .TRUE.
    print*, ""
    print '(A)', " ========================================"
    print '(A)', "         --Empirical PL Fitted--         "
    print '(A)', " ----------------------------------------"
    if (this%data%data_is_discrete) then
        print '("  ", A18, " = ", I12)', "x_min", int(this%x_min)
    else
        print '("  ", A18, " = ", F12.4)', "x_min", this%x_min
    endif
    print '("  ", A18, " = ", F12.4)', "alpha", this%alpha
    print '("  ", A18, " = ", F12.4)', "std_alpha", this%std_alpha
    if (this%was_pvalued) print '("  ", A18, " = ", F12.4)', "p_value", this%goodness_of_fit 
    if (this%was_pvalued) print '("  ", A18, " = ", F12.4)', "p_value error", this%p_value_eps
    print '("  ", A18, " = ", I12)',   "Data length", this%data%len
    print '("  ", A18, " = ", I12)',   "Tail length", this%n_tail
    print '("  ", A18, " = ", F12.5)', "time for fit (s)", this%mle_time
    if (this%was_pvalued) then
        print '("  ", A18, " = ", F12.5)', "time p_value (s)", this%hypothesis_time
        if (openMP) then
            print '(A)', " ----------------------------------------"
            print '(A)', "  * Parallel execution (OpenMP enabled)  "
        endif
    endif
    print '(A)', " ========================================"
    print*, ""
end subroutine

end module empirical_pl_mod