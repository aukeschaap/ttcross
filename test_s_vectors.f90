program test_s_vectors
    use s_vector_mod
    use default_lib
    implicit none
    integer :: n_dimensions, n_vectors, i, j

    ! Set dimensions
    call readarg(1, n_dimensions, 6)

    ! Generate s_vectors
    call generate_s_vectors(n_dimensions)

    ! Compute number of vectors
    n_vectors = 2**(n_dimensions - 1)

    ! Print s_vectors
    print *, 's_vectors:'
    do i = 1, n_vectors
        write(*,'(3x,*(i3))') (s_vectors(j,i), j=1,n_dimensions)
    end do
end program test_s_vectors
