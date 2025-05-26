module mvn_pdf_mod
  implicit none

  type mvn_data_t
    integer :: n
    real(8), allocatable :: mu(:)
    real(8), allocatable :: inv_cov(:,:)
    real(8) :: det_cov
  end type mvn_data_t

  type(mvn_data_t), save :: mvn_data

  contains

  subroutine mvn_init(n, r, T)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: r, T
    real(8), parameter :: sigma = 0.4d0, corr = 0.5d0
    real(8) :: X0, cov(n,n)
    integer :: i, j

    X0 = log(100.0d0)

    allocate(mvn_data%mu(n), mvn_data%inv_cov(n,n))
    mvn_data%n = n

    ! Compute mean
    do i = 1, n
      mvn_data%mu(i) = X0 + (r - 0.5d0 * sigma**2) * T
    end do

    ! Compute covariance matrix
    do i = 1, n
      do j = 1, n
        if (i == j) then
          cov(i,j) = sigma * sigma
        else
          cov(i,j) = sigma * corr * sigma
        end if
        cov(i,j) = cov(i,j) * T
      end do
    end do

    ! Print mean vector and covariance matrix if n < 10
    if (n < 10) then
      print *, 'Mean vector (mu):'
      do i = 1, n
        print '(F12.6)', mvn_data%mu(i)
      end do

      print *, 'Covariance matrix:'
      do i = 1, n
        print '(*(F12.6))', (cov(i,j), j = 1, n)
      end do
    end if

    ! Compute inverse and determinant
    call invert_and_determinant(cov, mvn_data%inv_cov, mvn_data%det_cov, n)
  end subroutine mvn_init


  function mvn_pdf(x) result(pdf)
    implicit none
    real(8), intent(in) :: x(:)
    real(8) :: pdf
    integer :: n, i, j
    real(8), parameter :: pi = 3.141592653589793d0
    real(8) :: diff(size(x)), exponent

    n = mvn_data%n
    diff = x - mvn_data%mu

    ! Mahalanobis distance
    exponent = 0.0d0
    do i = 1, n
      do j = 1, n
        exponent = exponent + diff(i) * mvn_data%inv_cov(i,j) * diff(j)
      end do
    end do

    pdf = exp(-0.5d0 * exponent) / sqrt((2.0d0 * pi)**n * mvn_data%det_cov)
  end function mvn_pdf

  subroutine invert_and_determinant(A, Ainv, detA, n)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: A(n,n)
    real(8), intent(out) :: Ainv(n,n), detA
    integer :: ipiv(n), info, i
    real(8) :: work(n,n)

    work = A
    call dgetrf(n, n, work, n, ipiv, info)

    detA = 1.0d0
    do i = 1, n
      if (ipiv(i) /= i) detA = -detA
      detA = detA * work(i,i)
    end do

    ! Initialize Ainv as identity
    Ainv = 0.0d0
    do i = 1, n
      Ainv(i,i) = 1.0d0
    end do

    ! Invert
    call dgetri(n, work, n, ipiv, Ainv, n, info)
    Ainv = work
  end subroutine invert_and_determinant

end module mvn_pdf_mod
