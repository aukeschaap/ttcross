module s_vector_mod
    implicit none
    integer, allocatable :: s_vectors(:,:)

contains

    subroutine generate_s_vectors(n_dimensions)
        implicit none
        integer, intent(in) :: n_dimensions
        integer :: n_vectors, i, j, bitmask

        n_vectors = 2**(n_dimensions - 1)

        if (.not. allocated(s_vectors)) then
            allocate(s_vectors(n_dimensions, n_vectors))
        end if

        do i = 0, n_vectors - 1
            s_vectors(1, i+1) = 1  ! First component is always 1
            do j = 2, n_dimensions
                bitmask = ibits(i, j-2, 1)  ! Get bit j-2 of i
                if (bitmask == 0) then
                    s_vectors(j, i+1) = 1
                else
                    s_vectors(j, i+1) = -1
                end if
            end do
        end do
    end subroutine generate_s_vectors

end module s_vector_mod
