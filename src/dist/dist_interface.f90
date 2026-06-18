module dist_interface_mod
use mle_kinds_mod
use rndgen_mod
implicit none
    
private
public :: empirical_distribution

    !> Abstract interface for any distribution
    type, abstract :: empirical_distribution
        
        !> Number of parameters
        integer(i4) :: num_params
        !> Intrinsec random number generator
        type(rndgen) :: gen
        
        !> Best founded values in fit
        integer(i4) :: n_tail                   !> Fitted tail size
        real(dp)    :: x_min                    !> Best x_min found
        real(dp)    :: ks                       !> ks stats in fit
        real(dp)    :: stats                    !> If use Anderson-Darling or other weighted stats
        real(dp), allocatable :: theta(:)       !> Maximum likelihook estimated parameters (theta is standard notation)
        real(dp), allocatable :: std_theta(:)   !> Std for MLE parameters

    contains
        
        !> Deffered procedures in empirical_distribution
        !> 1. For optimaze MLE routines
        procedure(pre_compute_intf), deferred :: pre_compute
        !> 2. Evaluate distribution cdf and parameters for a given x_min in data
        procedure(evaluate_tail_intf), deferred :: evaluate_tail
        !> 3. Distribution parameters names
        procedure(get_names_intf), deferred :: get_param_names
        !> 4. To gen a distribution-distributed array
        procedure(gen_rnd_array_intf), deferred :: generate_rnd_array

        !> Basic procedures
        procedure :: base_init => initialize_base_distribution
        procedure :: save_params => save_current_parameters
        
    end type empirical_distribution


    abstract interface

        subroutine pre_compute_intf(this, r_data)
            import :: empirical_distribution, random_data
            class(empirical_distribution), intent(inout) :: this
            type(random_data), intent(in) :: r_data
        end subroutine

        subroutine evaluate_tail_intf(this, i, N, x_min_candidate, r_data, cdf_out, theta_out, std_theta_out)
            import :: empirical_distribution, random_data, i4, dp
            class(empirical_distribution), intent(inout) :: this
            integer(i4), intent(in) :: i                !> Loop index
            integer(i4), intent(in) :: N                !> Total data size
            real(dp), intent(in) :: x_min_candidate     !> Candidate for x_min
            type(random_data), intent(in) :: r_data     !> Empirical data
            real(dp), intent(out) :: cdf_out(:)         !> CDF for x>=x_min
            real(dp), intent(out) :: theta_out(:)       !> MLE fitted parameters
            real(dp), intent(out) :: std_theta_out(:)   !> Std MLE parameters
        end subroutine

        function get_names_intf(this) result(names)
            import :: empirical_distribution
            class(empirical_distribution), intent(in) :: this
            character(len=20), allocatable :: names(:) 
        end function

        subroutine gen_rnd_array_intf(this, rnd_array)
            import :: empirical_distribution, i4, dp
            class(empirical_distribution), intent(in) :: this
            real(dp) , intent(inout) :: rnd_array(:)
        end subroutine

    end interface

contains

subroutine initialize_base_distribution(this, n_params, seed)
    class(empirical_distribution), intent(inout) :: this
    integer(i4), intent(in) :: n_params
    integer(i4), intent(in), optional :: seed
    type(clock_time) :: clock

    this%num_params = n_params

    if (present(seed)) then 
        call this%gen%init( seed )
    else
        call this%gen%init( clock%get_now() )
    endif
    
    if (allocated(this%theta)) deallocate(this%theta, this%std_theta)
    allocate( this%theta(this%num_params), this%std_theta(this%num_params) )
end subroutine

subroutine save_current_parameters(this, x_min_in, ks_in, n_tail_in, theta_in, std_theta_in)
    class(empirical_distribution), intent(inout) :: this
    real(dp), intent(in) ::  x_min_in, ks_in, theta_in(:), std_theta_in(:)
    integer(i4), intent(in) :: n_tail_in

        this%x_min  = x_min_in
        this%ks     = ks_in
        this%n_tail = n_tail_in

        if (.not. allocated(this%theta) ) allocate( this%theta(this%num_params) , this%std_theta(this%num_params) )
        if ( size(this%theta)/=this%num_params ) then
            deallocate( this%theta , this%std_theta )
            allocate( this%theta(this%num_params) , this%std_theta(this%num_params) )
        endif
        if ( size(theta_in)/=this%num_params ) error stop "Passed parameter array does not match in size"
        
        this%theta     = theta_in
        this%std_theta = std_theta_in
end subroutine
            

end module dist_interface_mod