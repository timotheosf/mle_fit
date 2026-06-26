module power_law_mod
use mle_kinds_mod
use rndgen_mod
use dist_interface_mod
implicit none

private
public :: power_law

    type, extends(distribution) :: power_law
        
        !> Auxiliary work variables
        real(dp), allocatable, private :: log_x(:)
        real(dp), allocatable, private :: sum_log_x(:)
        real(dp), private :: offset
    
    contains

        !> Get parameters by their names
        procedure :: alpha
        procedure :: std_alpha

        !> Allocating theta list
        procedure :: start_pl => wake_up_power_law

        !> Static binding procedures for power_law distribution
        procedure :: pre_compute => pl_pre_compute
        procedure :: evaluate_tail => pl_evaluate_tail
        procedure :: get_param_names => pl_get_param_names
        procedure :: generate_rnd_array => pl_generate_rnd_array
        procedure :: core_pdf => pl_core_pdf
        procedure :: core_cdf => pl_core_cdf

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
    type(empirical_data), intent(in) :: r_data
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
    type(empirical_data), intent(in) :: r_data
    real(dp), intent(out) :: cdf_out(:)
    real(dp), intent(out) :: theta_out(:)
    real(dp), intent(out) :: std_theta_out(:)

    integer(i4) :: n_tail_int , j , k_int
    real(dp) :: N_tail, log_sum, candidate_alpha, log_xmin, alpha_minus_1 , zeta_denom , last_x, last_cdf, gap , current_x

    n_tail_int = N - i + 1
    N_tail = real(n_tail_int, dp)

    !> Parameters calculation in O(1) time
    log_sum = this%sum_log_x(i) - N_tail * log(x_min_candidate - this%offset)
    candidate_alpha = 1.0_dp + N_tail / log_sum
    theta_out(1) = candidate_alpha
    std_theta_out(1) = (candidate_alpha - 1.0_dp) / sqrt(N_tail)

    !> Full vectorized cdf calulation 
    if (r_data%data_is_discrete) then
        zeta_denom = 1._dp/zeta_function(candidate_alpha, x_min_candidate)
        last_x = r_data%arr(i)
        last_cdf = 1.0_dp - (zeta_function(candidate_alpha, last_x + 1.0_dp) * zeta_denom)
        cdf_out(1) = last_cdf
        do j = i + 1, N
            current_x = r_data%arr(j)
            if (current_x > last_x) then
                gap = current_x - last_x
                if (gap <= 5.0_dp) then
                    do k_int = int(last_x) + 1, int(current_x)
                        last_cdf = last_cdf + (real(k_int, dp)**(-candidate_alpha)) * zeta_denom
                    end do
                else
                    last_cdf = 1.0_dp - (zeta_function(candidate_alpha, current_x + 1.0_dp) * zeta_denom)
                endif
                last_x = current_x
            endif
            cdf_out(j - i + 1) = last_cdf
        enddo
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

elemental function pl_core_pdf(this, x, discrete) result(res)
    class(power_law), intent(in) :: this
    real(dp), intent(in) :: x
    logical, intent(in)  :: discrete
    real(dp) :: res
    if (x < this%x_min) then
        res = 0.0_dp
        return
    endif
    if (discrete) then
        res = x**( -this%theta(1) ) / zeta_function(this%theta(1), this%x_min)
    else
        res = (this%theta(1) - 1._dp) * (this%x_min**(this%theta(1) - 1._dp)) * x**( -this%theta(1) )
    endif
end function pl_core_pdf

elemental function pl_core_cdf(this, x, discrete) result(res)
    class(power_law), intent(in) :: this
    real(dp), intent(in) :: x
    logical, intent(in)  :: discrete
    real(dp) :: res
    if (x < this%x_min) then
        res = 0.0_dp
        return
    endif
    if (discrete) then
        res = 1.0_dp - (zeta_function(this%theta(1), x) / zeta_function(this%theta(1), this%x_min))
    else
        res = 1.0_dp - ((x / this%x_min)**(-this%theta(1) + 1.0_dp))
    endif
end function pl_core_cdf

end module power_law_mod