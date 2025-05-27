module cos_approx_mod
    implicit none
    private
    public :: cos_approximate, cos_approximate_array

contains

    ! -----------------------------------------------------------------
    ! Approximate a function at a single point using its characteristic
    ! function.
    !
    ! Computes the following approximation:
    !   pdf_val = sum_{k=0}^{n_terms-1} coeff_k * cos(omega_k * (x - lower_bound))
    ! where:
    !   coeff_k = 2 / (upper_bound - lower_bound) * Re(phis[k] * exp(-i * omega_k * lower_bound))
    !
    ! Arguments:
    !   x: point at which to evaluate the function.
    !   phis: array of characteristic function values. Should be of size
    !       >= `n_terms`.
    !   lower_bound: lower bound of the interval.
    !   upper_bound: upper bound of the interval.
    !   n_terms: number of terms to use in the approximation. If not
    !       provided, defaults to the size of `phis`.
    !
    ! Returns:
    !   pdf_val: the approximate value of the function at `x`.
    ! -----------------------------------------------------------------
    function cos_approximate(x, phis, lower_bound, upper_bound, n_terms) result(pdf_val)
        implicit none
        double precision, intent(in) :: x
        double complex, intent(in) :: phis(:)
        double precision, intent(in) :: lower_bound, upper_bound
        integer, intent(in), optional :: n_terms
        double precision :: pdf_val

        integer :: k, n
        double precision :: omega, pi, pi_over_bound, coeff, arg
        double complex :: phi_k
        double precision :: sum

        if (.not. present(n_terms)) then
            n = size(phis)
        else
            n = n_terms
        end if

        if (n > size(phis)) then
            print *, "Error: n_terms exceeds the size of phis."
            pdf_val = 0.d0
            return
        end if

        pi = 3.1415926535897932384626433832795d0
        pi_over_bound = pi / (upper_bound - lower_bound)
        
        sum = 0.d0
        do k = 0, n_terms - 1
            omega = k * pi_over_bound
            phi_k = phis(k + 1)
            coeff = 2.d0 / (upper_bound - lower_bound) * real(phi_k * exp(-1.d0 * (0.d0, 1.d0) * omega * lower_bound))
            if (k == 0) coeff = coeff / 2.d0
            arg = omega * (x - lower_bound)
            sum = sum + coeff * cos(arg)
        end do

        pdf_val = sum
    end function cos_approximate


    ! -----------------------------------------------------------------
    ! Approximate a function at multiple points using its characteristic
    ! function.
    !
    ! Arguments:
    !   xs: array of points at which to evaluate the function.
    !   phis: array of characteristic function values. Should be of size
    !       >= `n_terms`.
    !   lower_bound: lower bound of the interval.
    !   upper_bound: upper bound of the interval.
    !   n_terms: number of terms to use in the approximation. If not
    !       provided, defaults to the size of `phis`.
    !
    ! Returns:
    !   pdf_vals: array of approximate values of the function at each
    !       point in `xs`.
    ! -----------------------------------------------------------------
    subroutine cos_approximate_array(xs, phis, lower_bound, upper_bound, n_terms, pdf_vals)
        implicit none
        double precision, intent(in) :: xs(:)
        double complex, intent(in) :: phis(:)
        double precision, intent(in) :: lower_bound, upper_bound
        integer, intent(in), optional :: n_terms
        double precision, intent(out) :: pdf_vals(size(xs))

        integer :: k, i, n
        double precision :: omega, pi, pi_over_bound, coeff, arg
        double complex :: phi_k

        if (.not. present(n_terms)) then
            n = size(phis)
        else
            n = n_terms
        end if

        if (n_terms > size(phis)) then
            print *, "Error: n_terms exceeds the size of phis."
            pdf_vals = 0.d0
            return
        end if

        pi = 3.1415926535897932384626433832795d0
        pi_over_bound = pi / (upper_bound - lower_bound)

        pdf_vals = 0.d0
        do k = 0, n - 1
            omega = k * pi_over_bound
            phi_k = phis(k + 1)
            coeff = 2.d0 / (upper_bound - lower_bound) * real(phi_k * exp(-1.d0 * (0.d0, 1.d0) * omega * lower_bound))
            if (k == 0) coeff = coeff / 2.d0

            do i = 1, size(xs)
                arg = omega * (xs(i) - lower_bound)
                pdf_vals(i) = pdf_vals(i) + coeff * cos(arg)
            end do
        end do
    end subroutine cos_approximate_array

end module cos_approx_mod
