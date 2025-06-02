module funcs
    implicit none
    private
    public :: gaussian_chf_nd

contains

    complex*16 function gaussian_chf_nd(n, omega, mu, sigma) result(val)
        implicit none
        integer, intent(in) :: n
        double precision, intent(in) :: omega(n), mu(n), sigma(n,n)
        double precision :: dot_mu, quad_form
        complex*16 :: im

        im = (0.d0, 1.d0)

        ! Compute ωᵀμ
        dot_mu = sum(omega(:) * mu(:))

        ! Compute ωᵀΣω
        quad_form = 0.d0
        quad_form = sum( (matmul(sigma, omega)) * omega )

        ! Evaluate φ(ω)
        val = exp(im * dot_mu - 0.5d0 * quad_form)
    end function gaussian_chf_nd

end module funcs
