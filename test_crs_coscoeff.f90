program main
    use coefficients_mod
    use cos_approx_mod
    use mvn_pdf_mod
    use tt_lib
    use dmrgg_lib
    use time_lib
    use quad_lib
    use default_lib
    use mat_lib
    use omp_lib
    use utils
    use constants
    implicit none
    include 'mpif.h'

    type(dtt) :: tt
    type(ztt) :: qq, tt_z
    integer :: i, j, k, p, n_dimensions, n_quad_points, max_rank, pivoting_strat
    integer :: info, nproc, me, adj
    integer(kind=8) :: neval
    double precision :: t1, t2, tcrs, acc, a, b, val, omega, X_0, corr, rate, T
    double complex :: im, tru, phis(32)
    double precision, allocatable :: par(:), cov(:,:), mean(:), sigma(:)
    double complex, allocatable :: w_complex(:)
    logical :: rescale
    double precision :: x
    integer :: n_pts
    double precision, allocatable :: xs(:), pdf_vals(:)
    character(len=100) :: filename
    integer :: ios, unit_out

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

    !-----------------------------------------------
    ! Display TT-Cross configuration summary
    !-----------------------------------------------
    if (me == 0) then
        write(*,*) 'TT-Cross Configuration:'
        write(*,'(3x,a,i10)') 'dimension:', n_dimensions
        write(*,'(3x,a,i10)') 'quadratur:', n_quad_points
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
    ! Parameters
    !-----------------------------------------------

    if (me == 0) then
        write(*,*) 'Parameters:'
    end if

    X_0 = log(100.0d0)

    allocate(sigma(n_dimensions), stat=info)
    do i = 1, n_dimensions
        sigma(i) = 0.4d0
    end do
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
    allocate(mean(n_dimensions), stat=info)
    do i = 1, n_dimensions
        mean(i) = X_0 + (rate - 0.5d0 * sigma(i)**2) * T
    end do
    
    ! Covariance matrix
    allocate(cov(n_dimensions, n_dimensions), stat=info)
    if (info /= 0) stop 'cannot allocate cov'
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
    if (me ==0 .and. n_dimensions < 10) then
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
        write(*,'(3x,a)') 'Generating s-vectors...'
    end if

    call generate_s_vectors(n_dimensions)
    if (me == 0) then
        write(*,'(3x,a,i10)') 'Number of s-vectors:', size(s_vectors, 2)
    end if

    ! Initialize coefficients for the Gaussian distribution
    call init_coefficients(n_dimensions, mean, cov, lower=0.525170185988090843d0, upper=8.52517018598809173d0)

    !-----------------------------------------------
    ! Allocate and configure quadrature nodes/weights on [a, b]
    !-----------------------------------------------

    if (me == 0) then
        write(*,'(3x,a)') 'Computing quadrature weights...'
    end if
    
    ! Interval of integration (according to cumulants with L = 10)
    a = 0.525170185988090843d0
    b = 8.52517018598809173d0

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

    call dtt_dmrgg(tt, calc_coefficient, maxrank=max_rank, accuracy=acc, &
                   pivoting=pivoting_strat, neval=neval)
    t2 = timef()
    tcrs = t2 - t1

    if (me == 0) write(*,'(a,i12,a,e12.4,a)') '...with', neval, ' evaluations completed in ', tcrs, ' sec.'

    ! Save the TT tensor to HDF5
    call save_dtt_to_hdf5(tt, "out/cos-coeff-tt-4-32-10-0.h5")

    call dealloc(tt)
    call mpi_finalize(info)
    if (info /= 0) stop 'mpi: finalize fail'
end program


