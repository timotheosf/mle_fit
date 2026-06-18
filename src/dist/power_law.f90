module power_law_mod
use mle_kinds_mod
use rndgen_mod
use dist_interface_mod
implicit none

private
public :: power_law

    type, extends(empirical_distribution) :: power_law
        
        !> Auxiliary work variables
        real(dp), allocatable, private :: log_x(:)
        real(dp), allocatable, private :: sum_log_x(:)
        real(dp), private :: offset

    
    contains

        !> Get parameters by their names
        procedure :: alpha
        procedure :: std_alpha

        !> PDF, CDF and CCDF calculations
        procedure, private :: core_pdf , core_ccdf , &
                              real_pdf , int_pdf , &
                              real_ccdf , int_ccdf , &
                              real_cdf , int_cdf
        generic :: pdf => real_pdf , int_pdf
        generic :: cdf => real_cdf , int_cdf
        generic :: ccdf => real_ccdf , int_ccdf

        !> Allocating theta list
        procedure :: start_pl => wake_up_power_law

        !> Abstract interface procedures for power_law distribution
        procedure :: pre_compute => pl_pre_compute
        procedure :: evaluate_tail => pl_evaluate_tail
        procedure :: get_param_names => pl_get_param_names
        procedure :: generate_rnd_array => pl_generate_rnd_array

        
        
    end type
      
contains

subroutine wake_up_power_law(this, seed)
    class(power_law), intent(inout) :: this
    integer(i4), intent(in), optional :: seed

    if (present(seed)) then
        call this%base_init(n_params = 1_i4, seed = seed)
    else
        call this%base_init(n_params = 1_i4)
    endif
end subroutine wake_up_power_law

subroutine pl_pre_compute(this, r_data)
    class(power_law), intent(inout) :: this
    type(random_data), intent(in) :: r_data
    integer(i4) :: N, i

    N = r_data%len
    
    !> Offset definition for MLE in discrete/continuos data
    if (r_data%data_is_discrete) then
        this%offset = 0.5_dp
    else
        this%offset = 0.0_dp
    endif

    !> Allocate auxiliary work arrays 
    if (allocated(this%sum_log_x)) deallocate(this%sum_log_x, this%log_x)
    allocate(this%sum_log_x(N), this%log_x(N))

    !> O(N) cumulative sum for log_x
    this%log_x(1:N) = log(r_data%arr(1:N))
    this%sum_log_x(N) = this%log_x(N)
    do i = N-1, 1, -1
        this%sum_log_x(i) = this%sum_log_x(i+1) + this%log_x(i)
    enddo
end subroutine

subroutine pl_evaluate_tail(this, i, N, x_min_candidate, r_data, cdf_out, theta_out, std_theta_out)
    class(power_law), intent(inout) :: this
    integer(i4), intent(in) :: i, N
    real(dp), intent(in) :: x_min_candidate
    type(random_data), intent(in) :: r_data
    real(dp), intent(out) :: cdf_out(:)
    real(dp), intent(out) :: theta_out(:)
    real(dp), intent(out) :: std_theta_out(:)

    integer(i4) :: n_tail_int
    real(dp) :: N_tail, log_sum, candidate_alpha, log_xmin, alpha_minus_1 , zeta_denom

    n_tail_int = N - i + 1
    N_tail = real(n_tail_int, dp)

    !> Parameters calculation in O(1) time
    log_sum = this%sum_log_x(i) - N_tail * log(x_min_candidate - this%offset)
    candidate_alpha = 1.0_dp + N_tail / log_sum
    theta_out(1) = candidate_alpha
    std_theta_out(1) = (candidate_alpha - 1.0_dp) / sqrt(N_tail)

    !> Full vectorized cdf calulation 
    if (r_data%data_is_discrete) then
        
        zeta_denom = zeta_function(candidate_alpha, x_min_candidate)
        cdf_out(1) = 1.0_dp - zeta_function(candidate_alpha, r_data%arr(i)) / zeta_denom
        do j = 2, n_tail_int
            if (r_data%arr(i+j-1) == r_data%arr(i+j-2)) then
                cdf_out(j) = cdf_out(j-1)
            else
                cdf_out(j) = 1.0_dp - zeta_function(candidate_alpha, r_data%arr(i+j-1)) / zeta_denom
            endif
        enddo
        !cdf_out(1:n_tail_int) = 1.0_dp - &
        !    zeta_function(candidate_alpha, r_data%arr(i:N)) / zeta_function(candidate_alpha, x_min_candidate)
    else
        log_xmin = log(x_min_candidate)
        alpha_minus_1 = candidate_alpha - 1.0_dp
        cdf_out(1:n_tail_int) = 1.0_dp - exp(alpha_minus_1 * (log_xmin - this%log_x(i:N)))
    endif
end subroutine

subroutine pl_generate_rnd_array(this, rnd_array)
    class(power_law), intent(in) :: this
    real(dp), intent(inout) :: rnd_array(:)
    
    integer(i4) :: arr_size, j
    real(dp) :: current_alpha, current_xmin

    arr_size = size(rnd_array)

    current_alpha = this%theta(1) 
    current_xmin  = this%x_min

    do j = 1, arr_size
        rnd_array(j) = this%gen%rnd()
    end do

    !> Inverting CDF transformation
    if (this%offset > 0.0_dp) then
        !> Discrete case
        rnd_array = floor((current_xmin - 0.5_dp) * (1.d0 - rnd_array)**(-1.d0 / (current_alpha - 1.d0)) + 0.5_dp)
    else
        !> Continuos case
        rnd_array = current_xmin * (1.d0 - rnd_array)**(-1.d0 / (current_alpha - 1.d0))
    end if
end subroutine

function pl_get_param_names(this) result(names)
    class(power_law), intent(in) :: this
    character(len=20), allocatable :: names(:)
    allocate(names(1))
    names(1) = "alpha"
end function

function alpha( this ) result(res)
    class(power_law), intent(in) :: this
    real(dp) :: res
    res = this%theta(1)
end function

function std_alpha( this ) result(res)
    class(power_law), intent(in) :: this
    real(dp) :: res
    res = this%std_theta(1)
end function

elemental function core_pdf( this , x , discrete ) result(res)
    class(power_law), intent(in) :: this
    real(dp), intent(in) :: x
    logical, intent(in), optional :: discrete
    logical :: x_is_discrete
    real(dp) :: res

    if (x < this%x_min) then
        res = 0.0_dp
        return
    endif

    x_is_discrete=.FALSE.
    if (present(discrete)) x_is_discrete=discrete
    if (x_is_discrete) then
        res = x**( -this%theta(1) )/zeta_function(this%theta(1), this%x_min)
    else
        res = (this%theta(1)-1._dp) * (this%x_min**(this%theta(1)-1._dp) )* x**( -this%theta(1) )
    endif
end function core_pdf

elemental function core_ccdf( this , x , discrete ) result(res)
    class(power_law), intent(in) :: this
    real(dp), intent(in) :: x
    logical, intent(in), optional :: discrete
    logical :: x_is_discrete
    real(dp) :: res

    if (x < this%x_min) then
        res = 1.0_dp
        return
    endif

    x_is_discrete=.FALSE.
    if (present(discrete)) x_is_discrete=discrete
    if (x_is_discrete) then
        res = zeta_function(this%theta(1), x)/zeta_function(this%theta(1), this%x_min)
    else
        res = ( x/this%x_min) ** (-this%theta(1)+1._dp)
    endif
end function core_ccdf

elemental function real_pdf( this , x ) result(res)
    class(power_law), intent(in) :: this
    real(dp), intent(in) :: x
    real(dp) :: res
    res = this%core_pdf( x )
end function real_pdf

elemental function int_pdf( this , x ) result(res)
    class(power_law), intent(in) :: this
    integer(i4), intent(in) :: x
    real(dp) :: res
    res = this%core_pdf( real(x,dp) , discrete=.true. )
end function int_pdf

elemental function real_ccdf( this , x ) result(res)
    class(power_law), intent(in) :: this
    real(dp), intent(in) :: x
    real(dp) :: res
    res = this%core_ccdf( x )
end function real_ccdf

elemental function int_ccdf( this , x ) result(res)
    class(power_law), intent(in) :: this
    integer(i4), intent(in) :: x
    real(dp) :: res
    res = this%core_ccdf( real(x,dp) , discrete=.true. )
end function int_ccdf

elemental function real_cdf( this , x ) result(res)
    class(power_law), intent(in) :: this
    real(dp), intent(in) :: x
    real(dp) :: res
    res = 1.0_dp - this%ccdf(x)
end function real_cdf

elemental function int_cdf( this , x ) result(res)
    class(power_law), intent(in) :: this
    integer(i4), intent(in) :: x
    real(dp) :: res
    res = 1.0_dp - this%ccdf(x)
end function int_cdf



end module power_law_mod