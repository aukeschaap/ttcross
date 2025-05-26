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
    integer :: i, j, p, n_dimensions, n_quad_points, max_rank, pivoting_strat
    integer :: info, nproc, me, adj
    integer(kind=8) :: neval
    double precision :: t1, t2, tcrs, acc, val, tru, a, b
    double complex :: cval
    double precision, allocatable :: par(:)
    double complex, allocatable :: w_complex(:)
    logical :: rescale
    double precision, external :: integrand, dfunc_stdnorm

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

    ! True value of the integral
    tru = (1.d0, 0.d0)

    allocate(par(2 * n_quad_points), stat=info)
    if (info /= 0) stop 'cannot allocate par'

    ! Compute Gauss-Legendre quadrature on [-1, 1]
    call lgwt(n_quad_points, par(1), par(n_quad_points+1))

    ! Map nodes to [a, b]
    forall(i=1:n_quad_points) par(i) = 0.5d0 * ((b - a) * par(i) + (a + b))

    ! Scale weights for [a, b]
    call dscal(n_quad_points, 0.5d0 * (b - a), par(n_quad_points+1), 1)

    ! Compute additional complex weights for quadrature
    allocate(w_complex(n_quad_points), stat=info)
    if (info /= 0) stop 'cannot allocate w_complex'

    w_complex = (1.0d0, 0.0d0)

    !-----------------------------------------------
    ! Prepare quadrature tensor
    !-----------------------------------------------

    if (me == 0) then
        write(*,'(3x,a)') 'Preparing quadrature tensor...'
    end if

    qq%l = 1
    qq%m = n_dimensions
    qq%n = n_quad_points
    qq%r = 1
    call alloc(qq)
    
    do i = 1, n_dimensions
        forall(p = 1:n_quad_points)
            qq%u(i)%p(1, p, 1) = dcmplx(par(n_quad_points + p), 0.d0) * w_complex(p)
        end forall
    end do

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
                   pivoting=pivoting_strat, neval=neval, tru=tru)
    t2 = timef()
    tcrs = t2 - t1

    if (me == 0) write(*,'(a,i12,a,e12.4,a)') '...with', neval, ' evaluations completed in ', tcrs, ' sec.'

    ! Convert tt from dtt to ztt
    tt_z = tt

    !-----------------------------------------------
    ! Evaluate and report
    !-----------------------------------------------
    cval = ztt_quad(tt_z, qq)
    if (me == 0) then
        write(*,'(a,2e50.40)') 'computed value:', real(cval)
        write(*,'(a,e50.40)') 'analytic value:', tru
        ! write(*,'(a,f7.2)') 'correct digits:', -dlog(dabs(1.d0 - val/tru)) / dlog(10.d0)
        write(*,'(a,f7.2)') 'correct digits:', -dlog(zabs((1.d0, 0.d0) - cval/tru)) / dlog(10.d0)
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

