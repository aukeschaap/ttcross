program test_coefficients
    use coefficients_mod
    use s_vector_mod
    use constants
    implicit none

    integer, parameter :: n_dimensions = 4
    integer, parameter :: m0 = 1, m1 = 1, m2 = 1, m3 = 32
    integer :: i, j, k, l, info, me, nproc
    integer :: ind(n_dimensions), mode_sizes(n_dimensions)
    double precision :: mean(n_dimensions), cov(n_dimensions,n_dimensions), sigma(n_dimensions)
    double precision :: lower, upper, X_0, corr, rate, T
    

    mode_sizes = (/ m0, m1, m2, m3 /)


    !-----------------------------------------------
    ! Parameters
    !-----------------------------------------------

    if (me == 0) then
        write(*,*) 'Parameters:'
    end if

    X_0 = log(100.0d0)

    ! allocate(sigma(n_dimensions), stat=info)
    sigma = 0.4d0
    corr = 0.5d0
    rate = 0.d0
    T = 1.d0

    if (me == 0) then
        write(*,'(3x,a,e12.4)') 'X_0     :', X_0
        write(*,'(3x,a,e12.4)') 'sigma   :', sigma(1)
        write(*,'(3x,a,e12.4)') 'corr    :', corr
        write(*,'(3x,a,e12.4)') 'rate    :', rate
        write(*,'(3x,a,e12.4)') 'T       :', T
    end if

    ! Mean
    ! allocate(mean(n_dimensions), stat=info)
    do i = 1, n_dimensions
        mean(i) = X_0 + (rate - 0.5d0 * sigma(i)**2) * T
    end do
    
    ! Covariance matrix
    ! allocate(cov(n_dimensions, n_dimensions), stat=info)
    ! if (info /= 0) stop 'cannot allocate cov'
    do i = 1, n_dimensions
        do j = 1, n_dimensions
            if (i == j) then
                cov(i, j) = sigma(i) * sigma(j) * T
            else
                cov(i, j) = sigma(i) * corr * sigma(j) * T
            end if
        end do
    end do

    ! Print mean vector and covariance matrix if n < 10
    if (n_dimensions < 10) then
        if (me == 0) then
            write(*,*) 'Mean vector (mu):'
            do i = 1, n_dimensions
                write(*,'(F12.6)', advance='no') mean(i)
            end do
            write(*,*)
            write(*,*) 'Covariance matrix:'
            do i = 1, n_dimensions
                write(*,'(*(F12.6))', advance='no') (cov(i, j), j = 1, n_dimensions)
                write(*,*)
            end do
        end if
    end if

    ! Generate s-vectors
    if (me == 0) then
        write(*,'(a)') 'Generating s-vectors...'
    end if

    call generate_s_vectors(n_dimensions)
    if (me == 0) then
        write(*,'(3x,a,i10)') 'Number of s-vectors:', size(s_vectors, 2)
    end if

    

    ! Initialize coefficients for the Gaussian distribution
    call init_coefficients(n_dimensions, mean, cov, lower=0.525170185988090843d0, upper=8.52517018598809173d0)


    ! Loop over all indices
    do i = 1, m0
        ind(1) = i
        do j = 1, m1
            ind(2) = j
            do k = 1, m2
                ind(3) = k
                do l = 1, m3
                    ind(4) = l
                    ! Print the coefficient for the current index
                    call print_coefficient(n_dimensions, ind, mode_sizes)
                end do
            end do
        end do
    end do

contains

    subroutine print_coefficient(n, ind, mode_sizes)
        integer, intent(in) :: n
        integer, intent(in) :: ind(n)
        integer, intent(in) :: mode_sizes(n)
        double precision :: c

        character(len=200) :: line
        character(len=20)  :: idx_str
        integer :: pos, i

        c = calc_coefficient(n, ind, mode_sizes)

        ! Print coefficient
        if (n < 10) then
            line = 'Index:'
            pos = len_trim(line) + 1

            do i = 1, n
                write(idx_str, '(I5)') ind(i)
                line(pos:) = ' '//trim(idx_str)
                pos = len_trim(line) + 1
            end do

            write(line(pos:), '(A, F25.17)') '  Coeff: ', c
            write(*, '(A)') trim(line)
        end if

    end subroutine print_coefficient

end program test_coefficients
