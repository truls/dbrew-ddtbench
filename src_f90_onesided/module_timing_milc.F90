!     Copyright (c) 2012 The Trustees of University of Illinois. All 
!     rights reserved. Use of this source code is governed by a 
!     BSD-style license that can be found in the LICENSE file.

      module timing_milc

      private

      public :: timing_milc_su3_zdown_ddt
      public :: timing_milc_su3_zdown_manual
      public :: timing_milc_su3_zdown_mpi_pack_ddt

      contains

      subroutine timing_milc_su3_zdown_ddt( DIM2, DIM3, DIM4, DIM5, &
        outer_loop, inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_5D_complex

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM2, DIM3, DIM4, DIM5
      integer, intent(in) :: outer_loop, inner_loop

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug

      integer, intent(in) :: local_communicator

!> local variables
      complex, dimension(3, DIM2, DIM3, DIM4, DIM5) :: array

      integer :: myrank, ier
      integer :: i, j
      integer :: base, bytes

      integer :: dtype_su3_vector_t, dtype_temp_t, &
        dtype_su3_zdown_t
      integer :: typesize
      integer(kind=MPI_ADDRESS_KIND) :: stride

      character(50) :: method

      integer :: win
      integer(kind=MPI_ADDRESS_KIND) :: win_size
      integer(kind=MPI_ADDRESS_KIND) :: target_disp = 0

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

! assume that one complex element has a size of 8 byte
      win_size = 3 * DIM2 * DIM3 * DIM4 * DIM5 * 8
      call MPI_Win_create( array(1,1,1,1,1), win_size, 8, MPI_INFO_NULL, &
        local_communicator, win, ier )

! initial fence to open epoch
      call MPI_Win_fence( 0, win, ier )

      base = myrank * 3 * DIM2 * DIM3 * DIM4 * DIM5 * 2 + 1
      call fill_unique_array_5D_complex( array, 3, DIM2, DIM3, DIM4, &
        DIM5, base )

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "mpi_ddt"

        call MPI_Type_size( MPI_COMPLEX, typesize, ier )
        bytes = 2 * DIM5 * DIM2*DIM3/2 * 3 * typesize

        call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop
!> modelling the zdown direction
        call MPI_Type_contiguous( 3, MPI_COMPLEX, dtype_su3_vector_t, &
          ier )

        call MPI_Type_vector( DIM5, DIM2*DIM3/2, DIM2*DIM3*DIM4/2, &
          dtype_su3_vector_t, dtype_temp_t, ier )

        call MPI_Type_size( dtype_su3_vector_t, typesize, ier )
        stride = typesize * DIM2 * DIM3 * DIM4 * DIM5/2
        call MPI_Type_create_hvector( 2, 1, stride, dtype_temp_t, &
          dtype_su3_zdown_t, ier )
        call MPI_Type_commit( dtype_su3_zdown_t, ier )

        call MPI_Type_free( dtype_su3_vector_t, ier )
        call MPI_Type_free( dtype_temp_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop

          if ( myrank .EQ. 0 ) then
            call MPI_Put( array(1,1,1,1,1), 1, dtype_su3_zdown_t, 1, &
              target_disp, 1, dtype_su3_zdown_t, win, ier )
            call MPI_Win_fence( 0, win, ier )
            call MPI_Win_fence( 0, win, ier )
            call timing_record(3)
          else
            call MPI_Win_fence( 0, win, ier )
            call MPI_Put( array(1,1,1,1,1), 1, dtype_su3_zdown_t, 0, &
              target_disp, 1, dtype_su3_zdown_t, win, ier )
            call MPI_Win_fence( 0, win, ier )
          endif

        enddo !> inner loop
      
        call MPI_Type_free( dtype_su3_zdown_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo !> outer loop

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      call MPI_Win_free( win, ier )

      end subroutine timing_milc_su3_zdown_ddt

      subroutine timing_milc_su3_zdown_manual( DIM2, DIM3, DIM4, DIM5, &
        outer_loop, inner_loop, correct_flag, ptypesize, testname, &
        filehandle_debug, local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_5D_complex

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM2, DIM3, DIM4, DIM5
      integer, intent(in) :: outer_loop, inner_loop

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug

      integer, intent(in) :: local_communicator

      complex, dimension(3, DIM2, DIM3, DIM4, DIM5) :: array
      complex, dimension(:), allocatable :: buffer

      integer :: myrank, ier
      integer :: i, j, k, l, m, io, ii
      integer :: base, bytes, pos, dsize

      integer :: typesize

      character(50) :: method

      integer :: win
      integer(kind=MPI_ADDRESS_KIND) :: win_size
      integer(kind=MPI_ADDRESS_KIND) :: target_disp = 0

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

      dsize = 2 * DIM5 * DIM2*DIM3/2 * 3
      allocate( buffer(dsize), stat=ier )
! assume that one complex element has a size of 8 byte
      win_size = dsize * 8
      call MPI_Win_create( buffer(1), win_size, 8, MPI_INFO_NULL, &
        local_communicator, win, ier )

! initial fence to open epoch
      call MPI_Win_fence( 0, win, ier )

      base = myrank * 3 * DIM2 * DIM3 * DIM4 * DIM5 * 2 + 1
      call fill_unique_array_5D_complex( array, 3, DIM2, DIM3, DIM4, &
        DIM5, base )

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "manual"

        call MPI_Type_size( MPI_COMPLEX, typesize, ier )
        bytes = 2 * DIM5 * DIM2*DIM3/2 * 3 * typesize

        call timing_init( testname, method, bytes )
      endif

      do io=1,outer_loop

!> modelling the zdown direction
        do ii=1,inner_loop

          if ( myrank .EQ. 0 ) then
            pos = 1
            do i=1,DIM5
              do j =1,DIM4,DIM4/2
                do k=1,DIM3/2
                  do l=1,DIM2
                    do m=1,3
                      buffer(pos) = array(m,l,k,j,i)
                      pos = pos + 1
                    enddo
                  enddo
                enddo
              enddo
            enddo
            call timing_record(2)
            call MPI_Put( buffer(1), dsize, MPI_COMPLEX, 1, &
              target_disp, dsize, MPI_COMPLEX, win, ier )
            call MPI_Win_fence( 0, win, ier )
            call MPI_Win_fence( 0, win, ier )
            call timing_record(3)
            pos = 1
            do i=1,DIM5
              do j =1,DIM4,DIM4/2
                do k=1,DIM3/2
                  do l=1,DIM2
                    do m=1,3
                      array(m,l,k,j,i) = buffer(pos)
                      pos = pos + 1
                    enddo
                  enddo
                enddo
              enddo
            enddo
            call timing_record(4)
          else
            call MPI_Win_fence( 0, win, ier )
            pos = 1
            do i=1,DIM5
              do j =1,DIM4,DIM4/2
                do k=1,DIM3/2
                  do l=1,DIM2
                    do m=1,3
                      array(m,l,k,j,i) = buffer(pos)
                      pos = pos + 1
                    enddo
                  enddo
                enddo
              enddo
            enddo
            pos = 1
            do i=1,DIM5
              do j =1,DIM4,DIM4/2
                do k=1,DIM3/2
                  do l=1,DIM2
                    do m=1,3
                      buffer(pos) = array(m,l,k,j,i)
                      pos = pos + 1
                    enddo
                  enddo
                enddo
              enddo
            enddo
            call MPI_Put( buffer(1), dsize, MPI_COMPLEX, 0, &
              target_disp, dsize, MPI_COMPLEX, win, ier )
            call MPI_Win_fence( 0, win, ier )
          endif

        enddo !> inner loop
      
      enddo !> outer loop

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      call MPI_Win_free( win, ier )

      deallocate( buffer, stat=ier )

      end subroutine timing_milc_su3_zdown_manual

      subroutine timing_milc_su3_zdown_mpi_pack_ddt( DIM2, DIM3, DIM4, &
        DIM5, outer_loop, inner_loop, correct_flag, ptypesize, &
        testname, filehandle_debug, local_communicator )

!      use timing, only: timing_init, timing_record, timing_print
      use utilities, only: fill_unique_array_5D_complex

      implicit none

      include 'mpif.h'

      integer, intent(in) :: DIM2, DIM3, DIM4, DIM5
      integer, intent(in) :: outer_loop, inner_loop

      logical, intent(out) :: correct_flag
      integer, intent(out) :: ptypesize

      character(50), intent(in) :: testname

      integer, intent(in) :: filehandle_debug

      integer, intent(in) :: local_communicator

      complex, dimension(3, DIM2, DIM3, DIM4, DIM5) :: array
      complex, dimension(:), allocatable :: buffer

      integer :: myrank, ier
      integer :: i, j 
      integer :: base, bytes, pos, dsize

      integer :: dtype_su3_vector_t, dtype_temp_t, dtype_su3_zdown_t
      integer :: typesize
      integer(kind=MPI_ADDRESS_KIND) :: stride

      character(50) :: method

      integer :: win
      integer(kind=MPI_ADDRESS_KIND) :: win_size
      integer(kind=MPI_ADDRESS_KIND) :: target_disp = 0

!> just some statements to prevent compiler warnings of unused variables
!> those parameter are included for future features
      correct_flag = .false.
      ptypesize = 0
      typesize = filehandle_debug

      call MPI_Comm_rank( local_communicator, myrank, ier )

      dsize = 2 * DIM5 * DIM2*DIM3/2 * 3
      allocate( buffer(dsize), stat=ier )
      call MPI_Type_size( MPI_COMPLEX, typesize, ier )
      bytes = dsize * typesize

      win_size = bytes
      call MPI_Win_create( buffer(1), win_size, typesize, MPI_INFO_NULL, &
        local_communicator, win, ier )

! initial fence to open epoch
      call MPI_Win_fence( 0, win, ier )

      base = myrank * 3 * DIM2 * DIM3 * DIM4 * DIM5 * 2 + 1
      call fill_unique_array_5D_complex( array, 3, DIM2, DIM3, DIM4, &
        DIM5, base )

      if ( myrank .EQ. 0 ) then
        write (method,'(A)') "mpi_pack_ddt"

        call timing_init( testname, method, bytes )
      endif

      do i=1,outer_loop

!> modelling the zdown direction
        call MPI_Type_contiguous( 3, MPI_COMPLEX, dtype_su3_vector_t, &
          ier )

        call MPI_Type_vector( DIM5, DIM2*DIM3/2, DIM2*DIM3*DIM4/2, &
          dtype_su3_vector_t, dtype_temp_t, ier )

        call MPI_Type_size( dtype_su3_vector_t, typesize, ier )
        stride = typesize*DIM2*DIM3*DIM4*DIM5/2
        call MPI_Type_create_hvector( 2, 1, stride, dtype_temp_t, &
          dtype_su3_zdown_t, ier )
        call MPI_Type_commit( dtype_su3_zdown_t, ier )

        call MPI_Type_free( dtype_su3_vector_t, ier )
        call MPI_Type_free( dtype_temp_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(1)
        endif

        do j=1,inner_loop

          if ( myrank .EQ. 0 ) then
            pos = 0
            call MPI_Pack( array(1,1,1,1,1), 1, dtype_su3_zdown_t, &
              buffer, bytes, pos, local_communicator, ier )
            call timing_record(2)
            call MPI_Put( buffer(1), pos, MPI_PACKED, 1, target_disp, &
              pos, MPI_PACKED, win, ier )
            call MPI_Win_fence( 0, win, ier )
            call MPI_Win_fence( 0, win, ier )
            call timing_record(3)
            pos = 0
            call MPI_Unpack( buffer, bytes, pos, array(1,1,1,1,1), 1, &
              dtype_su3_zdown_t, local_communicator, ier )
            call timing_record(4)
          else
            call MPI_Win_fence( 0, win, ier )
            pos = 0
            call MPI_Unpack( buffer, bytes, pos, array(1,1,1,1,1), 1, &
              dtype_su3_zdown_t, local_communicator, ier )
            pos = 0
            call MPI_Pack( array(1,1,1,1,1), 1, dtype_su3_zdown_t, &
              buffer, bytes, pos, local_communicator, ier )
            call MPI_Put( buffer(1), pos, MPI_PACKED, 0, target_disp, &
              pos, MPI_PACKED, win, ier )
            call MPI_Win_fence( 0, win, ier )
          endif

        enddo !> inner loop
      
        call MPI_Type_free( dtype_su3_zdown_t, ier )

        if ( myrank .EQ. 0 ) then
          call timing_record(5)
        endif

      enddo !> outer loop

      if ( myrank .EQ. 0 ) then
        call timing_print( .true. )
      endif

      call MPI_Win_free( win, ier )
      
      deallocate( buffer, stat=ier )

      end subroutine timing_milc_su3_zdown_mpi_pack_ddt

      end module timing_milc
