module quad_lib
    implicit none

    double precision, parameter, private :: tpi = 6.2831853071795864769252867665590057683943387987502116419498891846156328125724179972560696506842341359642961730265646132941876892191011644634507188162569622349005682054038770422111192892458979098607639288576219513318668922569512964675735663305424038182912971338469206972209086532964267872145204982825474491740132126311763497630418419256585081834307287357851807200226610610976409330427682939038830232188661145407315191839061843722347638652235862102370961489247599254991347037715054497824558763660238982596673467248813132861720427898927904494743814043597218874055410784343525863535047693496369353388102640011362542905271216555715426855155792183472743574429368818024499068602930991707421015845593785178470840399122242580439217280688363196272595495426199210374144226999999967459560999021194634656321926371900489189106938166052850446165066893700705238623763420200062756775057731750664167628412343553382946071965069808575109374623191257277647075751875039155637155610643424536132260038557532223918184328403d0

contains

    subroutine quad_rinv1(n, q)
        ! Approximation of 1/t using exponential quadrature
        implicit none
        integer, intent(inout) :: n
        double precision :: q(2, n)

        character(len=*), parameter :: subnam = 'quad_rinv1'
        integer :: nq, i, m
        double precision :: sinh_t, cosh_t, exp_sinh, log_huge
        double precision :: weight, alpha, h_step, t

        t = 0.d0
        log_huge = dlog(huge(t))
        nq = (n - 3) / 2
        m = 1
        h_step = dlog(tpi * nq) / nq
        q = 0.d0

        do i = -nq, nq
            t = dble(i) * h_step
            sinh_t = dsinh(t)
            cosh_t = dcosh(t)

            if (sinh_t < -log_huge) then
                ! skip extreme negative sinh_t
            else if (sinh_t > log_huge) then
                ! skip extreme positive sinh_t
            else
                m = m + 1
                exp_sinh = dexp(-sinh_t)
                weight = 2.d0 * cosh_t * h_step / (dsqrt(tpi / 2.d0) * (1.d0 + exp_sinh))
                alpha = (dlog(1.d0 + 1.d0 / exp_sinh)) ** 2
                q(1, m) = weight
                q(2, m) = alpha
            end if
        end do

        n = m
    end subroutine

    subroutine testquad_rinv(nq, q, a, b, n, err)
        ! Test quadrature accuracy for inverse approximation
        implicit none
        integer, intent(in) :: nq
        double precision, intent(in) :: q(2, nq)
        double precision, intent(in) :: a, b
        integer, intent(in) :: n
        double precision, intent(out), optional :: err

        character(len=*), parameter :: subnam = 'testquad_rinv'
        double precision :: t, x, y, log_t, log_a, log_b
        double precision :: approx_val, true_val, rel_error, max_error
        integer :: i, j

        if (a <= 0 .or. b <= 0) then
            write(*,*) subnam, ': illegal interval: ', a, b
            stop
        end if

        if (b < 1) then
            write(*,*) subnam, ': too few step numbers: ', n
            stop
        end if

        x = min(a, b)
        y = max(a, b)
        log_a = dlog(x)
        log_b = dlog(y)
        max_error = 0.d0

        open(10, file = '_testquad.rinv', access = 'sequential')
        do i = 0, n - 1
            log_t = log_a + dble(i) * (log_b - log_a) / (n - 1)
            t = dexp(log_t)
            approx_val = 0.d0

            do j = 1, nq
                approx_val = approx_val + q(1, j) * dexp(-q(2, j) * t * t)
            end do

            rel_error = t * dabs(1.d0 / t - approx_val)
            max_error = dmax1(rel_error, max_error)
            write(10, '(e12.5,3x,2f25.15,3x,e10.2)') t, 1.d0 / t, approx_val, rel_error
        end do
        close(10)

        if (present(err)) err = max_error
    end subroutine

    subroutine lgwt(n, x, w)
        ! Gauss-Legendre quadrature weights and nodes on [-1,1]
        implicit none
        integer, intent(in) :: n
        double precision, intent(out) :: x(n), w(n)

        double precision :: small
        integer :: i, j, m
        double precision :: p1, p2, p3, pp, z, z1

        small = 5 * epsilon(1.d0)
        m = (n + 1) / 2

        do i = 1, m
            z = dcos((tpi * (4 * i - 1)) / (8 * n + 4))

100         p1 = 1.0d0
            p2 = 0.0d0
            do j = 1, n
                p3 = p2
                p2 = p1
                p1 = ((2 * j - 1) * z * p2 - (j - 1) * p3) / j
            end do

            pp = n * (z * p1 - p2) / (z * z - 1)
            z1 = z
            z = z1 - p1 / pp
            if (dabs(z - z1) > small) goto 100

            x(i)       = -z
            x(n + 1 - i) = z
            w(i)       = 2.d0 / ((1 - z * z) * pp * pp)
            w(n + 1 - i) = w(i)
        end do
    end subroutine

end module
