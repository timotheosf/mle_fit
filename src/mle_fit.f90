module mle_fit_mod
!> Import types and interfaces from every module
use mle_kinds_mod
use dist_interface_mod
use run_mle_mod, only: mle_t
use power_law_mod, only: power_law_t, power_law
implicit none

    !> Library's facade
    private 
    public :: fit
    public :: mle_t
    public :: power_law_t, power_law

    interface fit
        module procedure quick_fit_powerlaw
    end interface fit

contains

subroutine private_quick_fit_routine( raw_data, fit_dist , x_min_out, theta_out, std_theta_out, run_pvalue, p_value_out, p_value_samples )
    class(*), intent(in) :: raw_data(:)
    class(distribution), intent(in) :: fit_dist
    real(dp), intent(out), optional :: x_min_out
    real(dp), intent(out), optional :: theta_out(:), std_theta_out(:)
    logical, intent(in), optional :: run_pvalue
    real(dp), intent(out), optional :: p_value_out
    integer(i4), intent(in), optional :: p_value_samples
    
    type(mle_t) :: engine
    class(distribution), allocatable :: mle_dist
    logical :: do_pvalue
    integer(i4) :: samples
    
    do_pvalue = .false.
    if (present(run_pvalue)) do_pvalue = run_pvalue
    samples = 1000
    if (present(p_value_samples)) samples = p_value_samples
    allocate( mle_dist , source=fit_dist )
    call engine%fast_fit( r_data = raw_data , dist = mle_dist )
    if (do_pvalue) call engine%p_value(N_samples=samples)

    if (present(x_min_out)) x_min_out = fit_dist%x_min
    if (present(theta_out)) theta_out(1:fit_dist%num_params) = fit_dist%theta
    if (present(std_theta_out)) std_theta_out(1:fit_dist%num_params) = fit_dist%std_theta
    if (present(p_value_out) .and. do_pvalue) p_value_out = engine%goodness_of_fit
    
    call engine%report()
end subroutine private_quick_fit_routine

subroutine quick_fit_powerlaw( r_data, dist, x_min , alpha, std_alpha, run_pvalue, p_value)
    class(*), intent(inout) :: r_data(:)
    class(power_law_t), intent(in) :: dist
    real(dp), intent(out), optional :: x_min, alpha, std_alpha
    logical, intent(in), optional :: run_pvalue
    real(dp), intent(out), optional :: p_value
    real(dp) :: temp_theta(1), temp_std(1)
    call private_quick_fit_routine( raw_data=r_data , fit_dist=dist , x_min_out=x_min, theta_out=temp_theta, std_theta_out=temp_std, run_pvalue=run_pvalue, p_value_out=p_value )
    if (present(alpha)) alpha = temp_theta(1)
    if (present(std_alpha)) std_alpha = temp_std(1)
end subroutine quick_fit_powerlaw


end module mle_fit_mod