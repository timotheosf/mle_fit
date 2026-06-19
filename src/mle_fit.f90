module mle_fit_mod
use mle_kinds_mod
use rndgen_mod
use dist_interface_mod
implicit none

private
public ::  mle

    type :: mle
    
        !> This method needs a data array and a allocated distribution
        type(random_data) :: data
        class(empirical_distribution), allocatable :: dist

        integer(i4) :: n_tail                   !> Empirical distribuiton tail length 
        real(dp) :: goodness_of_fit             !> P-value
        real(dp) :: mle_time                    !> Time costed for fitting
        real(dp) :: hypothesis_time             !> Time costed for hypothesis testing
        real(dp) :: p_value_eps                 !> P-value error

        !> Private control variables
        logical, private  :: weighted_adjust                   !> Flag to save the weighted used in the adjust
        logical, private  :: was_fitted                        !> Control flag: was this fitted?
        logical, private  :: was_pvalued                       !> Control flag: was this p_valued?
        real(dp), private :: lambda_used                       !> Lambda used in lamb_fit method
        type(clock_time), private :: internal_clock            !> Internal clock for benchmark
        real(dp), allocatable, private :: wrk_sort_buffer(:)   !> Internal buffer for optimize radix_sort in p_value subroutine

        !> Needed internal variables for greed searching
        real(dp), allocatable, private :: x_min_arr( : )
        real(dp), allocatable, private :: theta_arr( :, : )
        real(dp), allocatable, private :: std_theta_arr( :, : )
        real(dp), allocatable, private :: stats_arr( : )
        integer(i4), allocatable, private :: n_tail_arr( : )

    contains

        !> Starting the mle class
        procedure, private :: bind_dist
        procedure, private :: init => start_adjust_parameters
        !> Private fitting routines
        procedure, private :: core_fit => internal_core_fit

        !> Public
        procedure :: fast_fit => fast_find_best_parameters
        procedure :: lamb_fit => lambda_find_best_parameters
        procedure :: report => print_report
        procedure :: p_value => null_hypothesis_test 
    
    end type mle

contains

subroutine fast_find_best_parameters( this, r_data, distribuiton, xmin, theta, std_theta, ks, use_weight , fixed_xmin )
    class(mle), intent(inout) :: this
    class(empirical_distribution), intent(inout) :: distribuiton
    class(*), intent(in), optional :: r_data(:)
    logical , intent(in), optional :: use_weight
    real(dp), intent(out), optional :: theta(:), xmin, ks, std_theta(:)
    real(dp), intent(in), optional :: fixed_xmin
    call this%bind_dist( distribuiton )
    call this%core_fit( r_data=r_data, xmin=xmin, theta=theta, std_theta=std_theta, ks=ks, use_weight=use_weight , x_min_in=fixed_xmin )
end subroutine fast_find_best_parameters

subroutine lambda_find_best_parameters( this, r_data, distribuiton, xmin, theta, std_theta, ks, use_weight, lambda_in )
    class(mle), intent(inout) :: this
    class(empirical_distribution), intent(inout) :: distribuiton
    class(*), intent(in), optional :: r_data(:)
    logical , intent(in), optional :: use_weight
    real(dp), intent(out), optional :: theta(:), xmin, ks, std_theta(:)
    real(dp), intent(in), optional :: lambda_in
    real(dp) :: lambda
    lambda = 0.15_dp
    if (present(lambda_in)) lambda=lambda_in
    call this%bind_dist( distribuiton )
    call this%core_fit( r_data=r_data, xmin=xmin, theta=theta, std_theta=std_theta, ks=ks, use_weight=use_weight, lambda_in=lambda )
end subroutine lambda_find_best_parameters

subroutine bind_dist( this, distribuiton )
    class(mle), intent(inout) :: this
    class(empirical_distribution), intent(in) :: distribuiton
    if ( allocated(this%dist) ) deallocate( this%dist )
    allocate( this%dist, source=distribuiton ) !> Dynamic binding the distribution
end subroutine bind_dist

subroutine start_adjust_parameters( this, r_data, pre_ordering )    
    class(mle), intent(inout) :: this
    class(*), intent(in) :: r_data(:)
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

subroutine internal_core_fit( this, r_data, xmin, theta, std_theta, ks, lambda_in, use_weight, track_history , x_min_in )
    class(mle), intent(inout) :: this
    class(*), intent(in), optional :: r_data(:)
    logical , intent(in), optional :: use_weight, track_history
    real(dp), intent(in), optional :: lambda_in , x_min_in
    real(dp), intent(out), optional :: xmin, ks 
    real(dp), intent(out), optional :: theta(:), std_theta(:)

    real(dp), parameter :: eps = epsilon(1.0_dp)
    real(dp) :: bin_tolerance
    
    !> History tracking arrays
    real(dp), allocatable :: mle_x_min_arr(:), mle_stats_arr(:)
    real(dp), allocatable :: mle_theta_arr(:,:), mle_std_theta_arr(:,:)
    integer(i4), allocatable :: mle_n_tail_arr(:)
    
    !> Calculation arrays
    real(dp), allocatable :: ks_plus_arr(:), ks_minus_arr(:), current_cdf(:), seq(:), w(:)
    
    !> Working variables
    real(dp) :: prev_ks, prev_xmin, ks_statistics, current_ks, candidate_xmin
    real(dp) :: N_tail, lambda 
    integer(i4) :: i, N, n_tail_int, tail_len, non_zero_idx, prev_tail_len, n_p , x_min_pos 
    logical  :: apply_weight, save_history, is_decreasing

    !> Dynamic parameter arrays
    real(dp) :: best_xmin
    real(dp), allocatable :: cand_theta(:), cand_std(:)
    real(dp), allocatable :: best_theta(:), best_std(:)
    real(dp), allocatable :: prev_theta(:), prev_std(:)

    !--- Initializing ---!
    call this%internal_clock%start()

    if (present(r_data)) then
        call this%init( r_data )
    else if (.not. allocated(this%data%arr)) then
        error stop "Error: data not present in MLE class"
    endif
    
    if (.not. allocated(this%dist)) error stop "Error: Distribution not binded to MLE class"
    
    N = this%data%len
    n_p = this%dist%num_params

    !> Defines the offset based in the data original type
    if (this%data%data_is_discrete) then
        bin_tolerance = 0.5_dp
    else
        bin_tolerance = (this%data%arr(N)-this%data%arr(1))/N
    endif

    !> Applies the lambda_in penalty for the adjust
    if (present(lambda_in)) then
        lambda = lambda_in
    else
        lambda = 0.0_dp
    endif

    !> Defines the weight flags
    if (present(use_weight)) then
        apply_weight = use_weight
    else
        apply_weight = .FALSE. 
    endif
    if (apply_weight) allocate(w(N))

    !> Allocate dynamic parameter arrays based on distribution type
    allocate(cand_theta(n_p), cand_std(n_p))
    allocate(best_theta(n_p), best_std(n_p))
    allocate(prev_theta(n_p), prev_std(n_p))

    !> Allocating common calculation variables
    allocate( ks_plus_arr(N), ks_minus_arr(N), seq(N), current_cdf(N) )
    seq(N) = real(N,dp)
    do i = N-1, 1, -1
        seq(i) = real(i,dp)
    enddo
    
    ks_statistics = huge(1.0_dp) 

    !> Tracking minima history
    if (present(track_history)) then
        save_history = track_history
    else
        save_history = .FALSE.
    endif

    if (save_history) then
        non_zero_idx = 0
        allocate( mle_x_min_arr(N), mle_stats_arr(N), mle_n_tail_arr(N) )
        allocate( mle_theta_arr(n_p, N), mle_std_theta_arr(n_p, N) )
    endif

    !> Call the pre-computation routine of the specific distribution
    call this%dist%pre_compute( this%data )

    present_xmin_IF: if (present(x_min_in)) then
        candidate_xmin = x_min_in
        x_min_pos = 0
        find_x_min_pos_loop: do i = 1, N-1
        if ( abs(candidate_xmin - this%data%arr(i)) < bin_tolerance ) then
            x_min_pos = i
            exit find_x_min_pos_loop
        else
            cycle find_x_min_pos_loop
        endif
        enddo find_x_min_pos_loop
        if ( x_min_pos == 0 ) error stop "Error: passed x_min not found in data"
        i = x_min_pos
        n_tail_int = N - i + 1            
        N_tail = real(n_tail_int, dp)     

        !> DELEGATE TO DISTRIBUTION: Evaluate current tail and get parameters + CDF
        call this%dist%evaluate_tail( i, N, candidate_xmin, this%data, &
                                      current_cdf(1:n_tail_int),       &
                                      cand_theta, cand_std )

        !> Calculate KS Statistics
        if ( apply_weight ) then
            w( 1:n_tail_int ) = 1.0_dp / sqrt( (current_cdf(1:n_tail_int) * (1.0_dp - current_cdf(1:n_tail_int)) + eps) ) 
        else
            w( 1:n_tail_int ) = 1.0_dp
        endif
        if ( this%data%data_is_discrete ) then
            ks_plus_arr( 1:n_tail_int ) = abs( (seq( 1:n_tail_int ) / N_tail) - current_cdf( 1:n_tail_int ) ) * w( 1:n_tail_int )
            current_ks = maxval( ks_plus_arr( 1:n_tail_int ) )
        else
            ks_plus_arr( 1:n_tail_int )  = abs( (seq( 1:n_tail_int ) / N_tail) - current_cdf( 1:n_tail_int ) ) * w( 1:n_tail_int )
            ks_minus_arr( 1:n_tail_int ) = abs( current_cdf( 1:n_tail_int ) - ((seq( 1:n_tail_int ) - 1.0_dp) / N_tail) ) * w( 1:n_tail_int )
            current_ks = max( maxval(ks_plus_arr( 1:n_tail_int )), maxval(ks_minus_arr( 1:n_tail_int )) )
        endif

        !> The current stats is updated by this functional
        current_ks = current_ks - lambda*((N_tail/real(N)))
        tail_len = n_tail_int           
        best_xmin = candidate_xmin      
        best_theta = cand_theta
        best_std = cand_std
        ks_statistics = current_ks
    else
    !--- Main Loop ---!
    mle_main_loop: do i = 1, N-1
        if ( i > 1 ) then 
            !> Avoid repeated x_min candidates
            if ( abs(candidate_xmin - this%data%arr(i)) < bin_tolerance ) cycle mle_main_loop
        endif

        candidate_xmin = this%data%arr(i)
        n_tail_int = N - i + 1            
        N_tail = real(n_tail_int, dp)     

        !> DELEGATE TO DISTRIBUTION: Evaluate current tail and get parameters + CDF
        call this%dist%evaluate_tail( i, N, candidate_xmin, this%data, &
                                      current_cdf(1:n_tail_int),       &
                                      cand_theta, cand_std )

        !> Calculate KS Statistics
        if ( apply_weight ) then
            w( 1:n_tail_int ) = 1.0_dp / sqrt( (current_cdf(1:n_tail_int) * (1.0_dp - current_cdf(1:n_tail_int)) + eps) ) 
        else
            w( 1:n_tail_int ) = 1.0_dp
        endif
        if ( this%data%data_is_discrete ) then
            ks_plus_arr( 1:n_tail_int ) = abs( (seq( 1:n_tail_int ) / N_tail) - current_cdf( 1:n_tail_int ) ) * w( 1:n_tail_int )
            current_ks = maxval( ks_plus_arr( 1:n_tail_int ) )
        else
            ks_plus_arr( 1:n_tail_int )  = abs( (seq( 1:n_tail_int ) / N_tail) - current_cdf( 1:n_tail_int ) ) * w( 1:n_tail_int )
            ks_minus_arr( 1:n_tail_int ) = abs( current_cdf( 1:n_tail_int ) - ((seq( 1:n_tail_int ) - 1.0_dp) / N_tail) ) * w( 1:n_tail_int )
            current_ks = max( maxval(ks_plus_arr( 1:n_tail_int )), maxval(ks_minus_arr( 1:n_tail_int )) )
        endif

        !> The current stats is updated by this functional
        current_ks = current_ks - lambda*((N_tail/real(N)))
        
        if ( current_ks <= ks_statistics ) then
            tail_len = n_tail_int           
            best_xmin = candidate_xmin      
            best_theta = cand_theta
            best_std = cand_std
            ks_statistics = current_ks  
            
            if ( save_history ) then
                if ( current_ks < prev_ks ) then
                    is_decreasing = .TRUE.
                else if ( current_ks >= prev_ks .and. is_decreasing ) then
                    non_zero_idx = non_zero_idx + 1
                    mle_x_min_arr(non_zero_idx) = prev_xmin
                    mle_theta_arr(:, non_zero_idx) = prev_theta
                    mle_std_theta_arr(:, non_zero_idx) = prev_std
                    mle_stats_arr(non_zero_idx) = prev_ks
                    mle_n_tail_arr(non_zero_idx) = prev_tail_len
                    is_decreasing = .FALSE.
                endif
                prev_ks = current_ks
                prev_xmin = candidate_xmin
                prev_theta = cand_theta
                prev_std = cand_std
                prev_tail_len = n_tail_int
            endif
        endif
    enddo mle_main_loop
    endif present_xmin_IF

    !> Inject best parameters into the distribution
    call this%dist%save_params( best_xmin, (ks_statistics + lambda*(tail_len/real(N))), tail_len, best_theta, best_std )
    
    !> Update MLE internal flags
    this%n_tail = tail_len
    this%weighted_adjust = apply_weight 
    this%lambda_used = lambda 
    this%was_pvalued = .FALSE. 
    this%was_fitted = .TRUE.

    !> Return external variables if requested
    if (present(xmin)) xmin = best_xmin
    if (present(theta)) theta = best_theta
    if (present(std_theta)) std_theta = best_std
    if (present(ks)) ks = ks_statistics

    !> Track history
    if (save_history) then
        if (allocated(this%x_min_arr)) deallocate( this%x_min_arr, this%theta_arr, this%std_theta_arr, this%stats_arr, this%n_tail_arr )
        allocate( this%x_min_arr(non_zero_idx), this%stats_arr(non_zero_idx), this%n_tail_arr(non_zero_idx) )
        allocate( this%theta_arr(n_p, non_zero_idx), this%std_theta_arr(n_p, non_zero_idx) )
        
        this%x_min_arr(1:non_zero_idx) = mle_x_min_arr( 1:non_zero_idx )
        this%stats_arr(1:non_zero_idx) = mle_stats_arr( 1:non_zero_idx )
        this%n_tail_arr(1:non_zero_idx) = mle_n_tail_arr(1:non_zero_idx)
        this%theta_arr(:, 1:non_zero_idx) = mle_theta_arr(:, 1:non_zero_idx )
        this%std_theta_arr(:, 1:non_zero_idx) = mle_std_theta_arr(:, 1:non_zero_idx )
    endif

    call this%internal_clock%stop()
    this%mle_time = this%internal_clock%elapsed
end subroutine internal_core_fit


subroutine null_hypothesis_test( this, N_samples, track_penalities, p_value )
    !$ use omp_lib
    class(mle), intent(inout) :: this
    integer(i4), intent(in), optional :: N_samples 
    logical, intent(in), optional :: track_penalities
    real(dp), intent(out), optional :: p_value
    
    type(rndgen) :: thread_gen_1, thread_gen_2      
    real(dp) :: p_tail, synth_ks, real_ks
    integer(i4) :: i, N_trials, N, time_values(8), hits, j, base_seed
    integer(i4) :: thread_id , max_noise_idx, synth_head_size
    
    !> Thread specific arrays and instances
    real(dp), allocatable :: random_chooses( : ), synth_data_arr( : )
    type(mle) :: synth_mle

    !--- Initializing ---!
    call this%internal_clock%start()

    if (present(N_samples)) then
        N_trials = N_samples
    else
        N_trials = 1000
    endif

    real_ks = this%dist%ks
    if (present(track_penalities)) then
        !> Re-apply the penalization if tracking penalities
        if (track_penalities) real_ks = real_ks - this%lambda_used*((real(this%n_tail,dp)/real(this%data%len,dp)))
    endif

    N = this%data%len
    p_tail = real(this%n_tail,dp)/real(N ,dp)
    hits = 0
    max_noise_idx = max(1_i4, N - this%n_tail) 
    
    thread_id = 0  
    call date_and_time(values=time_values) 
    base_seed = time_values(8) + 1000*time_values(7) + 60000*time_values(6) + 3600000*time_values(5)

    !--- Parallel processing ---!
    !$omp parallel private(thread_id, j, i, random_chooses, synth_data_arr, synth_mle, synth_ks, thread_gen_1, thread_gen_2, synth_head_size)
    !$ thread_id = omp_get_thread_num() 
    
    !> Independent RNGs per thread
    call thread_gen_1%init( base_seed + thread_id * 1999  ) 
    call thread_gen_2%init( base_seed + thread_id * 3999 + 104729  ) 
    
    !> Deep copy of the polymorphic distribution to the thread-local MLE manager
    allocate(synth_mle%dist, source=this%dist)
    call synth_mle%dist%gen%init( base_seed + thread_id * 7777 )
    allocate(random_chooses(N), synth_data_arr(N))
    !> Ensure sorting buffer is ready for synth_mle
    allocate(synth_mle%wrk_sort_buffer(N)) 

    !$omp do reduction(+:hits) 
    sampling_loop: do j = 1, N_trials
        synth_head_size = 0
        
        !> Step 1: Bootstrap the empirical head
        do i=1,N
            random_chooses(i) = thread_gen_1%rnd()  
            if ( random_chooses(i) > p_tail ) synth_head_size = synth_head_size + 1
        enddo
        
        do i=1, synth_head_size
            synth_data_arr(i) = this%data%arr( floor(real((max_noise_idx),dp)*thread_gen_2%rnd(),kind=i4) + 1 )
        enddo

        !> Step 2: Generate the synthetic tail using the polymorphic distribution method!
        if (synth_head_size < N) then
            call synth_mle%dist%generate_rnd_array( synth_data_arr(synth_head_size+1 : N) )
        endif
        
        !> Step 3: Evaluate synthetic data
        call synth_mle%init( synth_data_arr ) 

        if (present(track_penalities)) then
            if (track_penalities) then
                call synth_mle%core_fit( ks=synth_ks, lambda_in=this%lambda_used, use_weight=this%weighted_adjust )
            else
                call synth_mle%core_fit( ks=synth_ks, use_weight=this%weighted_adjust )
            endif
        else
            call synth_mle%core_fit( ks=synth_ks, use_weight=this%weighted_adjust )
        endif
            
        if ( real_ks <= synth_ks ) hits = hits + 1
    enddo sampling_loop
    !$omp end do
    
    !> Deallocate thread-local polymorphic variable to avoid memory leaks
    if (allocated(synth_mle%dist)) deallocate(synth_mle%dist)
    
    !$omp end parallel
    !--- End parallel processing ---!

    !> Returns the p-value
    if (present(p_value)) p_value = real(hits,dp)/real(N_trials,dp)
    this%goodness_of_fit = real(hits,dp)/real(N_trials,dp) 
    this%p_value_eps = 1.0_dp / (2.0_dp * sqrt(real(N_trials, dp))) 
    this%was_pvalued = .TRUE.
    
    call this%internal_clock%stop()
    this%hypothesis_time = this%internal_clock%elapsed
end subroutine null_hypothesis_test


subroutine print_report(this)
    class(mle), intent(in) :: this
    character(len=20), allocatable :: p_names(:)
    integer(i4) :: k

    if (.not. this%was_fitted) then
        print *, "Error: MLE model is not fitted yet."
        return
    end if

    print *, ""
    print '(A)', " ========================================"
    print '(A)', "         -- Empirical MLE Fitted --      "
    print '(A)', " ----------------------------------------"

    !> Converts x_min to int if the distribution is discrete
    if (this%data%data_is_discrete) then
        print '("  ", A18, " = ", I12)', "x_min", int(this%dist%x_min)
    else
        print '("  ", A18, " = ", F12.4)', "x_min", this%dist%x_min
    end if

    !> Get the parameters names of the distribution
    p_names = this%dist%get_param_names()
    
    !> Print each param name and corresponding MLE value
    do k = 1, this%dist%num_params
        print '("  ", A18, " = ", F12.4)', trim(p_names(k)), this%dist%theta(k)
        print '("  ", A18, " = ", F12.4)', "std_"//trim(p_names(k)), this%dist%std_theta(k)
    end do

    !> Print goodness of fit if the null hypothesis was tested
    if (this%was_pvalued) then
        print '("  ", A18, " = ", F12.4)', "p_value", this%goodness_of_fit 
        print '("  ", A18, " = ", F12.4)', "p_value error", this%p_value_eps
    end if

    !> Print benchmark informations
    print '("  ", A18, " = ", I12)',   "Data length", this%data%len
    print '("  ", A18, " = ", I12)',   "Tail length", this%dist%n_tail
    print '("  ", A18, " = ", F12.5)', "time for fit (s)", this%mle_time

    if (this%was_pvalued) then
        print '("  ", A18, " = ", F12.5)', "time p_value (s)", this%hypothesis_time
    !> Print if OMP parallel execution was used in p_value test
    !$        print '(A)', " ----------------------------------------"
    !$        print '(A)', "  * Parallel execution (OpenMP enabled)  "
    end if
    
    print '(A)', " ========================================"
    print *, ""
end subroutine print_report

end module mle_fit_mod