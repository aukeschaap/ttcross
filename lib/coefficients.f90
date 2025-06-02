module coefficients_mod
    use s_vector_mod
    use funcs
    use constants
    implicit none

    double precision, allocatable :: coeff_mu(:)
    double precision, allocatable :: coeff_sigma(:,:)
    double precision :: lower_bound, upper_bound

contains

    subroutine init_coefficients(n_dimensions, mean, cov, lower, upper)
        integer, intent(in) :: n_dimensions
        double precision, intent(in) :: mean(n_dimensions)
        double precision, intent(in) :: cov(n_dimensions, n_dimensions)
        double precision, intent(in) :: lower, upper

        if (.not. allocated(coeff_mu)) then
            allocate(coeff_mu(n_dimensions))
        end if
        if (.not. allocated(coeff_sigma)) then
            allocate(coeff_sigma(n_dimensions, n_dimensions))
        end if

        coeff_mu(:) = mean(:)
        coeff_sigma(:,:) = cov(:,:)
        lower_bound = lower
        upper_bound = upper
    end subroutine init_coefficients


    double precision function calc_coefficient(n_dimensions, ind, mode_sizes) result(f)
        implicit none
        integer, intent(in) :: n_dimensions
        integer, intent(in) :: ind(n_dimensions)
        integer, intent(in) :: mode_sizes(n_dimensions)

        integer :: i, j, n_s
        double precision :: factor, t(n_dimensions)
        double precision :: k(n_dimensions), bounds
        complex*16 :: phi_val, exp_term, im
        double precision :: real_sum, one_over_bounds

        one_over_bounds = 1 / (upper_bound - lower_bound)
        im = (0.d0, 1.d0)

        factor = 2.d0 * one_over_bounds**n_dimensions

        ! Sum over s-vectors
        n_s = size(s_vectors, 2)
        real_sum = 0.d0

        do i = 1, n_s
            do j = 1, n_dimensions
                t(j) = (((pi * dble(s_vectors(j, i))) * dble(ind(j) - 1)) * one_over_bounds)
            end do

            phi_val = gaussian_chf_nd(n_dimensions, t, coeff_mu, coeff_sigma)
            exp_term = exp(-im * lower_bound * sum(t))
            real_sum = real_sum + real(exp_term * phi_val)
        end do

        f = factor * real_sum
    end function calc_coefficient

end module coefficients_mod
