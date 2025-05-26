module cos_approx_mod
    implicit none
    private
    public :: cos_approximate, cos_approximate_array

contains

    ! Single-point COS approximation
    function cos_approximate(x, phis, a, b) result(pdf_val)
        implicit none
        double precision, intent(in) :: x
        double precision, intent(in) :: a, b
        double complex, intent(in) :: phis(32)
        double precision :: pdf_val

        integer :: k, N
        double precision :: omega, pi, coeff, arg
        double complex :: phi_k
        double precision :: sum

        pi = 3.1415926535897932384626433832795d0
        N = 32
        sum = 0.d0

        do k = 0, N - 1
            omega = k * pi / (b - a)
            phi_k = phis(k + 1)
            coeff = 2.d0 / (b - a) * real(phi_k * exp(-1.d0 * (0.d0, 1.d0) * omega * a))
            if (k == 0) coeff = coeff / 2.d0
            arg = omega * (x - a)
            sum = sum + coeff * cos(arg)
        end do

        pdf_val = sum
    end function cos_approximate


    ! Multiple-point COS approximation
    subroutine cos_approximate_array(xs, n_x, phis, a, b, pdf_vals)
        implicit none
        integer, intent(in) :: n_x
        double precision, intent(in) :: xs(n_x)
        double precision, intent(in) :: a, b
        double complex, intent(in) :: phis(32)
        double precision, intent(out) :: pdf_vals(n_x)

        integer :: k, i, N
        double precision :: omega, pi, coeff, arg
        double complex :: phi_k

        pi = 3.1415926535897932384626433832795d0
        N = 32

        pdf_vals = 0.d0

        do k = 0, N - 1
            omega = k * pi / (b - a)
            phi_k = phis(k + 1)
            coeff = 2.d0 / (b - a) * real(phi_k * exp(-1.d0 * (0.d0, 1.d0) * omega * a))
            if (k == 0) coeff = coeff / 2.d0

            do i = 1, n_x
                arg = omega * (xs(i) - a)
                pdf_vals(i) = pdf_vals(i) + coeff * cos(arg)
            end do
        end do
    end subroutine cos_approximate_array

end module cos_approx_mod
