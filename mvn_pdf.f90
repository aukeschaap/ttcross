module mvn_pdf_mod
  implicit none
  contains

  function mvn_pdf(x, n, r, T) result(pdf)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: x(n), r, T
    real(8) :: pdf
    real(8), parameter :: sigma = 0.4d0, pi = 3.141592653589793d0, corr = 0.5d0
    real(8) :: mu(n), diff(n), chol(n,n), inv_cov(n,n), det_cov, exponent
    integer :: i, j
    real(8) :: X0

    X0 = log(100.0d0)

    ! Compute the mean vector
    do i = 1, n
      mu(i) = X0 + (r - 0.5d0 * sigma**2) * T
    end do

    ! Construct the covariance matrix (correlation matrix scaled)
    do i = 1, n
      do j = 1, n
        if (i == j) then
          chol(i,j) = sigma * sigma
        else
          chol(i,j) = sigma * corr * sigma
        end if
        chol(i,j) = chol(i,j) * T
      end do
    end do

    ! Use Cholesky or LU decomposition to compute inverse and determinant
    call invert_and_determinant(chol, inv_cov, det_cov, n)

    ! Compute x - mu
    diff = x - mu

    ! Compute the Mahalanobis distance
    exponent = 0.0d0
    do i = 1, n
      do j = 1, n
        exponent = exponent + diff(i) * inv_cov(i,j) * diff(j)
      end do
    end do

    pdf = exp(-0.5d0 * exponent) / sqrt((2.0d0 * pi)**n * det_cov)
  end function mvn_pdf

  subroutine invert_and_determinant(A, Ainv, detA, n)
    implicit none
    integer, intent(in) :: n
    real(8), intent(in) :: A(n,n)
    real(8), intent(out) :: Ainv(n,n), detA
    integer :: ipiv(n), info, i
    real(8) :: work(n,n), det_log

    work = A
    call dgetrf(n, n, work, n, ipiv, info)
    det_log = 0.0d0
    detA = 1.0d0
    do i = 1, n
      if (ipiv(i) /= i) detA = -detA
      detA = detA * work(i,i)
    end do
    Ainv = 0.0d0
    do i = 1, n
      Ainv(i,i) = 1.0d0
    end do
    call dgetri(n, work, n, ipiv, Ainv, n, info)
    Ainv = work
  end subroutine invert_and_determinant

end module mvn_pdf_mod
