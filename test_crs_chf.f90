program main
    use mvn_pdf_mod
    use tt_lib
    use dmrgg_lib
    use time_lib
    use quad_lib
    use default_lib
    use mat_lib
    use omp_lib
    implicit none
    include 'mpif.h'

    type(dtt) :: tt
    type(ztt) :: qq, tt_z
    integer :: i, j, k, p, n_dimensions, n_quad_points, max_rank, pivoting_strat
    integer :: info, nproc, me, adj
    integer(kind=8) :: neval
    double precision :: t1, t2, tcrs, acc, val, a, b, omega, pi
    double complex :: im, tru, ans(32)
    double precision, allocatable :: par(:)
    double complex, allocatable :: w_complex(:)
    logical :: rescale
    double precision, external :: integrand, dfunc_stdnorm
    double complex, external :: get_reference_val

    ! Initialize constants
    pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286208998628034825342117067982148086513282306647093844609550582231725359408128481117450284102701938521105559644622948954930381964428810975665933446128475648233786783165271201909145648566923460348610454326648213393607260249141273724587006606315588174881520920962829254091715364367892590360011330530548820466521384146951941511609433057270365759591953092186117381932611793105118548074462379962749567351885752724891227938183011949128831426076960045243166083812168476586578265042207222101

    !-----------------------------------------------
    ! Read input parameters
    !-----------------------------------------------
    call readarg(1, n_dimensions, 6)
    call readarg(2, n_quad_points, 65)
    call readarg(3, max_rank, 20)
    call readarg(4, pivoting_strat, 1)

    !-----------------------------------------------
    ! Initialize MPI
    !-----------------------------------------------
    call mpi_init(info)
    if (info /= 0) stop 'mpi: init fail'
    call mpi_comm_size(MPI_COMM_WORLD, nproc, info)
    if (info /= 0) stop 'mpi: comm_size fail'
    call mpi_comm_rank(MPI_COMM_WORLD, me, info)
    if (info /= 0) stop 'mpi: comm_rank fail'

    ! Adjust quadrature if even
    adj = 0
    if (mod(n_quad_points, 2) == 0) then
        n_quad_points = n_quad_points + 1
        adj = 1
    end if

    !-----------------------------------------------
    ! Display configuration summary
    !-----------------------------------------------
    if (me == 0) then
        write(*,*) 'Hi, this is TT cross interpolation for computing integrals...'
        write(*,'(3x,a,i10)') 'dimension:', n_dimensions
        write(*,'(3x,a,i10)', advance='no') 'quadratur:', n_quad_points
        if (adj == 1) then
            write(*,'(a)') ' (adjusted)'
        else
            write(*,*)
        end if
        write(*,'(3x,a,i10)') 'TT ranks :', max_rank
        write(*,'(3x,a,i10)') 'pivoting :', pivoting_strat
        write(*,'(3x,a,i10)') 'MPI procs:', nproc
!$OMP PARALLEL
        if (omp_get_thread_num() == 0) then
            write(*,'(3x,a,i10)') 'OMP thrds:', omp_get_num_threads()
        end if
!$OMP END PARALLEL
        write(*,'(3x,a,i10)') 'sizeof(d):', storage_size(1.d0)
        write(*,'(3x,a,e10.3)') 'epsilon  :', epsilon(1.d0)
    end if

    acc = 500 * epsilon(1.d0)


    !-----------------------------------------------
    ! Allocate and configure quadrature nodes/weights on [a, b]
    !-----------------------------------------------

    if (me == 0) then
        write(*,'(3x,a)') 'Computing quadrature weights...'
    end if
    
    ! Interval of integration (according to cumulants with L = 10)
    a = 0.525170
    b = 8.525170

    allocate(par(2 * n_quad_points), stat=info)
    if (info /= 0) stop 'cannot allocate par'

    ! Compute Gauss-Legendre quadrature on [-1, 1]
    call lgwt(n_quad_points, par(1), par(n_quad_points+1))

    ! Map nodes to [a, b]
    forall(i=1:n_quad_points) par(i) = 0.5d0 * ((b - a) * par(i) + (a + b))

    ! Scale weights for [a, b]
    call dscal(n_quad_points, 0.5d0 * (b - a), par(n_quad_points+1), 1)


    !-----------------------------------------------
    ! Run TT-cross
    !-----------------------------------------------

    if (me == 0) then
        write(*,'(3x,a)') 'Running TT-cross...'
    end if

    t1 = timef()
    tt%l = 1
    tt%m = n_dimensions
    tt%n = n_quad_points
    tt%r = 1
    call alloc(tt)

    call mvn_init(n_dimensions, 0.d0, 1.d0)

    call dtt_dmrgg(tt, integrand, par, maxrank=max_rank, accuracy=acc, &
                   pivoting=pivoting_strat, neval=neval)
    t2 = timef()
    tcrs = t2 - t1

    if (me == 0) write(*,'(a,i12,a,e12.4,a)') '...with', neval, ' evaluations completed in ', tcrs, ' sec.'

    ! Convert tt from dtt to ztt
    tt_z = tt


    !-----------------------------------------------
    ! Build and evaluate the quadrature tensor
    !-----------------------------------------------

    if (me == 0) then
        write(*,'(3x,a)') 'Preparing quadrature tensor...'
    end if

    ! Compute additional complex weights for quadrature
    allocate(w_complex(n_quad_points), stat=info)
    if (info /= 0) stop 'cannot allocate w_complex'

    ! Allocate quadrature tensor
    qq%l = 1
    qq%m = n_dimensions
    qq%n = n_quad_points
    qq%r = 1
    call alloc(qq)

    do k = 0, 31
        omega = k * pi / (300.d0 - 0.d0)
        im = (0.d0, 1.d0)

        do p = 1, n_quad_points
            w_complex(p) = exp(im * omega * exp(par(p)) / dble(n_dimensions))
        end do

        do i = 1, n_dimensions
            forall(p = 1:n_quad_points)
                qq%u(i)%p(1, p, 1) = dcmplx(par(n_quad_points + p), 0.d0) * w_complex(p)
            end forall
        end do

        ans(k+1) = ztt_quad(tt_z, qq)
    end do
    
    
    ! Print quadrature nodes and weights as doubles and complex numbers (not scientific notation)
    ! if (me == 0) then
    !     write(*,'(3x,a)') 'Quadrature nodes and weights:'
    !     write(*,'(3x,a,i10)') 'n_quad_points:', n_quad_points
    !     write(*,'(3x,a)') 'Nodes:'
    !     do i = 1, n_quad_points
    !         write(*,'(3x,a,i10,2f20.10)') 'Node', i, par(i), par(n_quad_points + i)
    !     end do
    !     write(*,'(3x,a)') 'Weights:'
    !     do i = 1, n_quad_points
    !         write(*,'(3x,a,i10,f20.10)') 'Weight', i, par(n_quad_points + i)
    !     end do
    !     write(*,'(3x,a)') 'Complex weight modification:'
    !     do i = 1, n_quad_points
    !         write(*,'(3x,a,i10,2f20.10)') 'Weight', i, real(w_complex(i)), aimag(w_complex(i))
    !     end do
    !     write(*,'(3x,a)') 'Done.'
    ! end if

    !-----------------------------------------------
    ! Evaluate and report
    !-----------------------------------------------

    if (me == 0) then
        do k = 0, 31
            tru = get_reference_val(k)
            write(*,'(a,2e50.40)') 'computed value: ', real(ans(k+1)), aimag(ans(k+1))
            write(*,'(a,2e50.40)') 'analytic value: ', real(tru), aimag(tru)
            write(*,'(a,f7.2)') 'correct digits:', -dlog(zabs((1.d0, 0.d0) - ans(k+1)/tru)) / dlog(10.d0)
        end do
        write(*,'(a)') 'Good bye.'
    end if
    


    call dealloc(tt)
    call dealloc(tt_z)
    call mpi_finalize(info)
    if (info /= 0) stop 'mpi: finalize fail'
end program


double precision function integrand(n_dimensions, ind, mode_sizes, par) result(f)
    use mvn_pdf_mod
    implicit none
    integer, intent(in) :: n_dimensions
    integer, intent(in) :: ind(n_dimensions)
    integer, intent(in) :: mode_sizes(n_dimensions)
    double precision, intent(inout), optional :: par(*)
    integer :: i
    double precision :: x(n_dimensions)

    ! Extract integration node
    do i = 1, n_dimensions
        x(i) = par(ind(i))
    end do

    f = mvn_pdf(x)
end function


function get_reference_val(k) result(val)
    implicit none
    integer, intent(in) :: k
    double complex :: val

    select case (k)
        case(0); val = cmplx(1.0000000000000004, 0.0)
        case(1); val = cmplx(0.477573486783962, 0.8149466884621303)
        case(2); val = cmplx(-0.36837299289083225, 0.7151922711168992)
        case(3); val = cmplx(-0.6257714195536066, 0.09178553732163314)
        case(4); val = cmplx(-0.34684735452966187, -0.3172332243460625)
        case(5); val = cmplx(-0.0017955158585152564, -0.33646866663071234)
        case(6); val = cmplx(0.162377557153463, -0.17003135351832718)
        case(7); val = cmplx(0.16111757391061166, -0.014897840710532378)
        case(8); val = cmplx(0.09249068083644305, 0.06033659639176453)
        case(9); val = cmplx(0.026622324264449382, 0.07016730893534587)
        case(10); val = cmplx(-0.011662817341915883, 0.04952550109973056)
        case(11); val = cmplx(-0.02453611388348422, 0.024273967930164675)
        case(12); val = cmplx(-0.02283078002233625, 0.005817209101882007)
        case(13); val = cmplx(-0.015560865654613152, -0.004063037376331763)
        case(14); val = cmplx(-0.008092447372180414, -0.007401121012863568)
        case(15); val = cmplx(-0.0026887989524065486, -0.007093450692951045)
        case(16); val = cmplx(0.00035180694746948793, -0.005222442957982434)
        case(17); val = cmplx(0.001702458119621754, -0.0031713885853042268)
        case(18); val = cmplx(0.0019787600882656027, -0.0016294446166633578)
        case(19); val = cmplx(0.0017002968172530539, -0.0004927083206423032)
        case(20); val = cmplx(0.0012683915479485622, 0.00011351837686527592)
        case(21); val = cmplx(0.0008009851216763085, 0.0004300613182277521)
        case(22); val = cmplx(0.0004241540880317125, 0.0004423028147576922)
        case(23); val = cmplx(0.00017861244723283787, 0.0004195598508893528)
        case(24); val = cmplx(2.8445454564976957e-05, 0.00031333373400801133)
        case(25); val = cmplx(-4.530757450752984e-05, 0.00022356345058376992)
        case(26); val = cmplx(-8.234702200071308e-05, 0.00016196042183622313)
        case(27); val = cmplx(-9.948295545041033e-05, 8.23964841318507e-05)
        case(28); val = cmplx(-8.455184484541963e-05, 4.940045542157806e-05)
        case(29); val = cmplx(-6.974592247472916e-05, 1.6767171278621708e-05)
        case(30); val = cmplx(-5.3515566669730744e-05, -1.2023150233163912e-05)
        case(31); val = cmplx(-3.397596073667283e-05, -2.6028176412912153e-05)
    end select
end function get_reference_val