program test_chf_equal
    use constants
    use funcs
    implicit none

    integer, parameter :: n = 4
    double precision :: mu(n), cov(n, n), t(n)
    double precision :: sigma(n), corr, rate, expiry, X_0
    complex*16 :: phi
    double precision :: real_ref, imag_ref
    double precision :: real_err, imag_err
    integer :: i, j, k1, k2, k3, k4
    character(len=512) :: cmd
    character(len=200) :: output_file
    integer :: ios
    double precision :: cpp_real, cpp_imag

    output_file = 'tmp_chf_out.txt'

    ! Parameters
    X_0 = log(100.0d0)
    sigma = 0.4d0
    corr = 0.5d0
    rate = 0.d0
    expiry = 1.d0

    ! Construct mean vector
    do i = 1, n
        mu(i) = X_0 + (rate - 0.5d0 * sigma(i)**2) * expiry
    end do

    ! Construct covariance matrix
    do i = 1, n
        do j = 1, n
            if (i == j) then
                cov(i, j) = sigma(i) * sigma(j) * expiry
            else
                cov(i, j) = sigma(i) * corr * sigma(j) * expiry
            end if
        end do
    end do

    ! Loop over t in [0, 1]^d using 3 points per dimension
    do k1 = 0, 2
        do k2 = 0, 2
            do k3 = 0, 2
                do k4 = 0, 2
                    t = (/ k1/2.d0, k2/2.d0, k3/2.d0, k4/2.d0 /)

                    ! Compute the characteristic function
                    phi = gaussian_chf_nd(n, t, mu, cov)

                    ! Print results from Fortran
                    write(*,'(A)') 'Computed Gaussian characteristic function (Fortran):'
                    write(*,'(F23.20, a, f23.20)') real(phi), ' ', aimag(phi)

                    ! Call the external C++ program
                    write(cmd, '(A,F12.6,F12.6,F12.6,F12.6)') './build/bin/test_chf_equal ', t(1), t(2), t(3), t(4)
                    call system(trim(cmd))
                end do
            end do
        end do
    end do

end program test_chf_equal
