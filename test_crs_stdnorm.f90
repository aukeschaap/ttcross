program main
    use tt_lib
    use dmrgg_lib
    use time_lib
    use quad_lib
    use default_lib
    use mat_lib
    use omp_lib
    implicit none
    include 'mpif.h'

    type(dtt) :: tt, qq
    integer :: i, j, n_dimensions, n_quad_points, max_rank, pivoting_strat
    integer :: info, nproc, me, adj
    integer(kind=8) :: neval
    double precision :: t1, t2, tcrs, acc, val, tru, a, b
    double precision, allocatable :: par(:), tmp_nodes(:), tmp_weights(:)
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

    acc = 5 * epsilon(1.d0)


    !-----------------------------------------------
    ! Allocate and configure quadrature nodes/weights on [a, b]
    !-----------------------------------------------
    
    ! Interval of integration
    a = -10.d0
    b = 10.d0

    ! True value of the integral
    ! tru = (b - a) ** n_dimensions
    tru = sqrt(3.141592653589793238d0) ** n_dimensions

    allocate(par(2 * n_quad_points), stat=info)
    if (info /= 0) stop 'cannot allocate par'

    ! Compute Gauss-Legendre quadrature on [-1, 1]
    call lgwt(n_quad_points, par(1), par(n_quad_points+1))

    ! Map nodes to [a, b]
    forall(i=1:n_quad_points) par(i) = 0.5d0 * ((b - a) * par(i) + (a + b))

    ! Scale weights for [a, b]
    call dscal(n_quad_points, 0.5d0 * (b - a), par(n_quad_points+1), 1)

    !-----------------------------------------------
    ! Prepare quadrature tensor
    !-----------------------------------------------
    qq%l = 1
    qq%m = n_dimensions
    qq%n = n_quad_points
    qq%r = 1
    call alloc(qq)
    do i = 1, n_dimensions
        call dcopy(n_quad_points, par(n_quad_points+1), 1, qq%u(i)%p, 1)
    end do
    call ones(qq)

    !-----------------------------------------------
    ! (Optional) rescaling to avoid underflow
    !-----------------------------------------------
    rescale = .false.
    if (rescale) then
        call dscal(n_quad_points, dble(n_quad_points), par(n_quad_points + 1), 1)
        do i = 1, n_dimensions
            qq%u(i)%p = 1.d0 / dble(n_quad_points)
        end do
    end if

    !-----------------------------------------------
    ! Run TT-cross
    !-----------------------------------------------
    t1 = timef()
    tt%l = 1
    tt%m = n_dimensions
    tt%n = n_quad_points
    tt%r = 1
    call alloc(tt)

    call dtt_dmrgg(tt, integrand, par, maxrank=max_rank, accuracy=acc, &
                   pivoting=pivoting_strat, neval=neval, quad=qq, tru=tru)
    t2 = timef()
    tcrs = t2 - t1

    if (me == 0) write(*,'(a,i12,a,e12.4,a)') '...with', neval, ' evaluations completed in ', tcrs, ' sec.'

    !-----------------------------------------------
    ! Evaluate and report
    !-----------------------------------------------
    val = dtt_quad(tt, qq)
    if (me == 0) then
        write(*,'(a,e50.40)') 'computed value:', val
        write(*,'(a,e50.40)') 'analytic value:', tru
        write(*,'(a,f7.2)') 'correct digits:', -dlog(dabs(1.d0 - val/tru)) / dlog(10.d0)
        write(*,'(a)') 'Good bye.'
    end if

    call dealloc(tt)
    call mpi_finalize(info)
    if (info /= 0) stop 'mpi: finalize fail'
end program


double precision function integrand(n_dimensions, ind, n, par) result(f)
    implicit none
    integer, intent(in) :: n_dimensions
    integer, intent(in) :: ind(n_dimensions), n(n_dimensions)
    double precision, intent(inout), optional :: par(*)
    integer :: i
    double precision :: x(n_dimensions)

    ! Extract the actual integration nodes
    do i = 1, n_dimensions
        x(i) = par(ind(i))  ! assumes nodes start at par(1)
    end do

    ! Your integrand goes here (e.g., Gaussian function)
    f = exp(-sum(x**2))  ! Example: multivariate Gaussian

    ! Multiply the corresponding quadrature weights
    do i = 1, n_dimensions
        f = f * par(n(i) + ind(i))
    end do
end function
