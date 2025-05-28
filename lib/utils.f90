module utils
   use hdf5
   use tt_lib
   implicit none

contains

    subroutine save_dtt_to_hdf5(tt, filename)
        implicit none
        type(dtt), intent(in) :: tt
        character(len=*), intent(in) :: filename

        integer(hid_t) :: file_id, group_id, dataspace_id, dataset_id, plist_id
        integer :: hdferr
        integer(hsize_t), dimension(1) :: dims1
        integer(hsize_t), dimension(3) :: dims3
        integer :: k, l, m, rk1, nk, rk2
        character(len=100) :: dsetname

        call h5open_f(hdferr)
        call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, hdferr)
        call h5gcreate_f(file_id, "TT", group_id, hdferr)

        ! Save mode sizes and ranks
        dims1(1) = tt%m - tt%l + 1
        call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
        call h5dcreate_f(group_id, "modes", H5T_NATIVE_INTEGER, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, tt%n(tt%l:tt%m), dims1, hdferr)
        call h5dclose_f(dataset_id, hdferr)
        call h5sclose_f(dataspace_id, hdferr)

        dims1(1) = tt%m - tt%l + 2
        call h5screate_simple_f(1, dims1, dataspace_id, hdferr)
        call h5dcreate_f(group_id, "ranks", H5T_NATIVE_INTEGER, dataspace_id, dataset_id, hdferr)
        call h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, tt%r(tt%l-1:tt%m), dims1, hdferr)
        call h5dclose_f(dataset_id, hdferr)
        call h5sclose_f(dataspace_id, hdferr)

        ! Save each core
        do k = tt%l, tt%m
            rk1 = tt%r(k - 1)
            nk  = tt%n(k)
            rk2 = tt%r(k)
            dims3 = (/ rk1, nk, rk2 /)
            call h5screate_simple_f(3, dims3, dataspace_id, hdferr)
            write(dsetname, "(A,I0)") "core_", k - tt%l
            call h5dcreate_f(group_id, trim(dsetname), H5T_NATIVE_DOUBLE, dataspace_id, dataset_id, hdferr)
            call h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, tt%u(k)%p, dims3, hdferr)
            call h5dclose_f(dataset_id, hdferr)
            call h5sclose_f(dataspace_id, hdferr)
        end do

        call h5gclose_f(group_id, hdferr)
        call h5fclose_f(file_id, hdferr)
        call h5close_f(hdferr)

    end subroutine
end module
