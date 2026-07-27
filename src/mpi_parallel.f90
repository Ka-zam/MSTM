module parallel_runtime
   use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
   use, intrinsic :: iso_fortran_env, only: real32, real64
   use mpi, only: mpi_allreduce, mpi_barrier, mpi_bcast, mpi_comm_create, mpi_comm_group, &
                  mpi_comm_null, mpi_comm_rank, mpi_comm_size, mpi_comm_split, mpi_comm_world, &
                  mpi_complex, mpi_double_complex, mpi_double_precision, mpi_finalize, mpi_group_incl, &
                  mpi_init, mpi_integer, mpi_max, mpi_real, mpi_recv, mpi_reduce, mpi_send, mpi_status_size, &
                  mpi_sum, mpi_wtime
   implicit none(type, external)
   private

   integer, public :: mstm_global_rank = 0
   integer, public :: mstm_global_numprocs = 1

   public :: mpi_comm_world, mpi_comm_null
   public :: parallel_wall_time
   public :: parallel_initialize, parallel_finalize
   public :: parallel_rank, parallel_size, parallel_barrier
   public :: parallel_split, parallel_group, parallel_group_include, parallel_communicator_create
   public :: parallel_broadcast, parallel_reduce_sum, parallel_allreduce_sum, parallel_allreduce_max
   public :: parallel_allreduce_sum_complex64_sequence
   public :: parallel_send, parallel_receive

   interface parallel_broadcast
      module procedure broadcast_integer
      module procedure broadcast_real64
      module procedure broadcast_complex64
   end interface parallel_broadcast

   interface parallel_reduce_sum
      module procedure reduce_sum_integer
      module procedure reduce_sum_complex32
      module procedure reduce_sum_real64
      module procedure reduce_sum_complex64
   end interface parallel_reduce_sum

   interface parallel_allreduce_sum
      module procedure allreduce_sum_real64
      module procedure allreduce_sum_complex64
   end interface parallel_allreduce_sum

   interface parallel_send
      module procedure send_complex64
   end interface parallel_send

   interface parallel_receive
      module procedure receive_complex64
   end interface parallel_receive

contains

   subroutine parallel_initialize()
      integer :: error_code

      call mpi_init(error_code)
      call mpi_comm_rank(mpi_comm_world, mstm_global_rank, error_code)
      call mpi_comm_size(mpi_comm_world, mstm_global_numprocs, error_code)
   end subroutine parallel_initialize

   subroutine parallel_finalize()
      integer :: error_code

      call mpi_finalize(error_code)
   end subroutine parallel_finalize

   real(real64) function parallel_wall_time()
      parallel_wall_time = mpi_wtime()
   end function parallel_wall_time

   subroutine parallel_rank(mpi_rank, mpi_comm)
      integer, intent(out) :: mpi_rank
      integer, optional, intent(in) :: mpi_comm
      integer :: communicator, error_code

      communicator = mpi_comm_world
      if (present(mpi_comm)) communicator = mpi_comm
      call mpi_comm_rank(communicator, mpi_rank, error_code)
   end subroutine parallel_rank

   subroutine parallel_size(mpi_size, mpi_comm)
      integer, intent(out) :: mpi_size
      integer, optional, intent(in) :: mpi_comm
      integer :: communicator, error_code

      communicator = mpi_comm_world
      if (present(mpi_comm)) communicator = mpi_comm
      call mpi_comm_size(communicator, mpi_size, error_code)
   end subroutine parallel_size

   subroutine parallel_barrier(mpi_comm)
      integer, optional, intent(in) :: mpi_comm
      integer :: communicator, error_code

      communicator = mpi_comm_world
      if (present(mpi_comm)) communicator = mpi_comm
      call mpi_barrier(communicator, error_code)
   end subroutine parallel_barrier

   subroutine parallel_split(mpi_color, mpi_key, mpi_new_comm, mpi_comm)
      integer, intent(in) :: mpi_color, mpi_key
      integer, intent(out) :: mpi_new_comm
      integer, optional, intent(in) :: mpi_comm
      integer :: communicator, error_code

      communicator = mpi_comm_world
      if (present(mpi_comm)) communicator = mpi_comm
      call mpi_comm_split(communicator, mpi_color, mpi_key, mpi_new_comm, error_code)
   end subroutine parallel_split

   subroutine parallel_group(mpi_group, mpi_comm)
      integer, intent(out) :: mpi_group
      integer, optional, intent(in) :: mpi_comm
      integer :: communicator, error_code

      communicator = mpi_comm_world
      if (present(mpi_comm)) communicator = mpi_comm
      call mpi_comm_group(communicator, mpi_group, error_code)
   end subroutine parallel_group

   subroutine parallel_group_include(mpi_group, mpi_size, mpi_new_group_list, mpi_new_group, mpi_comm)
      integer, intent(in) :: mpi_group, mpi_size, mpi_new_group_list(*)
      integer, intent(out) :: mpi_new_group
      integer, optional, intent(in) :: mpi_comm
      integer :: error_code

      call mpi_group_incl(mpi_group, mpi_size, mpi_new_group_list, mpi_new_group, error_code)
   end subroutine parallel_group_include

   subroutine parallel_communicator_create(mpi_group, mpi_new_comm, mpi_comm)
      integer, intent(in) :: mpi_group
      integer, intent(out) :: mpi_new_comm
      integer, optional, intent(in) :: mpi_comm
      integer :: communicator, error_code

      communicator = mpi_comm_world
      if (present(mpi_comm)) communicator = mpi_comm
      call mpi_comm_create(communicator, mpi_group, mpi_new_comm, error_code)
   end subroutine parallel_communicator_create

   subroutine broadcast_integer(send_buffer, mpi_number, mpi_rank, mpi_comm)
      integer, contiguous, target, intent(inout) :: send_buffer(..)
      integer, intent(in) :: mpi_number, mpi_rank
      integer, optional, intent(in) :: mpi_comm
      integer :: communicator, error_code
      integer, pointer :: buffer(:)

      communicator = mpi_comm_world
      if (present(mpi_comm)) communicator = mpi_comm
      call c_f_pointer(c_loc(send_buffer), buffer, [mpi_number])
      call mpi_bcast(buffer, mpi_number, mpi_integer, mpi_rank, communicator, error_code)
   end subroutine broadcast_integer

   subroutine broadcast_real64(send_buffer, mpi_number, mpi_rank, mpi_comm)
      real(real64), contiguous, target, intent(inout) :: send_buffer(..)
      integer, intent(in) :: mpi_number, mpi_rank
      integer, optional, intent(in) :: mpi_comm
      integer :: communicator, error_code
      real(real64), pointer :: buffer(:)

      communicator = mpi_comm_world
      if (present(mpi_comm)) communicator = mpi_comm
      call c_f_pointer(c_loc(send_buffer), buffer, [mpi_number])
      call mpi_bcast(buffer, mpi_number, mpi_double_precision, mpi_rank, communicator, error_code)
   end subroutine broadcast_real64

   subroutine broadcast_complex64(send_buffer, mpi_number, mpi_rank, mpi_comm)
      complex(real64), contiguous, target, intent(inout) :: send_buffer(..)
      integer, intent(in) :: mpi_number, mpi_rank
      integer, optional, intent(in) :: mpi_comm
      integer :: communicator, error_code
      complex(real64), pointer :: buffer(:)

      communicator = mpi_comm_world
      if (present(mpi_comm)) communicator = mpi_comm
      call c_f_pointer(c_loc(send_buffer), buffer, [mpi_number])
      call mpi_bcast(buffer, mpi_number, mpi_double_complex, mpi_rank, communicator, error_code)
   end subroutine broadcast_complex64

   subroutine reduce_sum_integer(receive_buffer, mpi_number, mpi_rank, mpi_comm, send_buffer)
      integer, contiguous, target, intent(inout) :: receive_buffer(..)
      integer, contiguous, optional, target, intent(in) :: send_buffer(..)
      integer, intent(in) :: mpi_number, mpi_rank
      integer, optional, intent(in) :: mpi_comm
      integer :: communicator, error_code
      integer, allocatable :: local_buffer(:)
      integer, pointer :: receive_values(:), send_values(:)

      communicator = mpi_comm_world
      if (present(mpi_comm)) communicator = mpi_comm
      call c_f_pointer(c_loc(receive_buffer), receive_values, [mpi_number])
      if (present(send_buffer)) then
         call c_f_pointer(c_loc(send_buffer), send_values, [mpi_number])
         call mpi_reduce(send_values, receive_values, mpi_number, mpi_integer, mpi_sum, mpi_rank, communicator, error_code)
      else
         allocate (local_buffer(mpi_number))
         local_buffer = receive_values
         call mpi_reduce(local_buffer, receive_values, mpi_number, mpi_integer, mpi_sum, mpi_rank, communicator, error_code)
      end if
   end subroutine reduce_sum_integer

   subroutine reduce_sum_complex32(receive_buffer, mpi_number, mpi_rank, mpi_comm, send_buffer)
      complex(real32), contiguous, target, intent(inout) :: receive_buffer(..)
      complex(real32), contiguous, optional, target, intent(in) :: send_buffer(..)
      integer, intent(in) :: mpi_number, mpi_rank
      integer, optional, intent(in) :: mpi_comm
      integer :: communicator, error_code
      complex(real32), allocatable :: local_buffer(:)
      complex(real32), pointer :: receive_values(:), send_values(:)

      communicator = mpi_comm_world
      if (present(mpi_comm)) communicator = mpi_comm
      call c_f_pointer(c_loc(receive_buffer), receive_values, [mpi_number])
      if (present(send_buffer)) then
         call c_f_pointer(c_loc(send_buffer), send_values, [mpi_number])
         call mpi_reduce(send_values, receive_values, mpi_number, mpi_complex, mpi_sum, mpi_rank, communicator, error_code)
      else
         allocate (local_buffer(mpi_number))
         local_buffer = receive_values
         call mpi_reduce(local_buffer, receive_values, mpi_number, mpi_complex, mpi_sum, mpi_rank, communicator, error_code)
      end if
   end subroutine reduce_sum_complex32

   subroutine reduce_sum_real64(receive_buffer, mpi_number, mpi_rank, mpi_comm, send_buffer)
      real(real64), contiguous, target, intent(inout) :: receive_buffer(..)
      real(real64), contiguous, optional, target, intent(in) :: send_buffer(..)
      integer, intent(in) :: mpi_number, mpi_rank
      integer, optional, intent(in) :: mpi_comm
      integer :: communicator, error_code
      real(real64), allocatable :: local_buffer(:)
      real(real64), pointer :: receive_values(:), send_values(:)

      communicator = mpi_comm_world
      if (present(mpi_comm)) communicator = mpi_comm
      call c_f_pointer(c_loc(receive_buffer), receive_values, [mpi_number])
      if (present(send_buffer)) then
         call c_f_pointer(c_loc(send_buffer), send_values, [mpi_number])
         call mpi_reduce(send_values, receive_values, mpi_number, mpi_double_precision, mpi_sum, &
                         mpi_rank, communicator, error_code)
      else
         allocate (local_buffer(mpi_number))
         local_buffer = receive_values
         call mpi_reduce(local_buffer, receive_values, mpi_number, mpi_double_precision, mpi_sum, &
                         mpi_rank, communicator, error_code)
      end if
   end subroutine reduce_sum_real64

   subroutine reduce_sum_complex64(receive_buffer, mpi_number, mpi_rank, mpi_comm, send_buffer)
      complex(real64), contiguous, target, intent(inout) :: receive_buffer(..)
      complex(real64), contiguous, optional, target, intent(in) :: send_buffer(..)
      integer, intent(in) :: mpi_number, mpi_rank
      integer, optional, intent(in) :: mpi_comm
      integer :: communicator, error_code
      complex(real64), allocatable :: local_buffer(:)
      complex(real64), pointer :: receive_values(:), send_values(:)

      communicator = mpi_comm_world
      if (present(mpi_comm)) communicator = mpi_comm
      call c_f_pointer(c_loc(receive_buffer), receive_values, [mpi_number])
      if (present(send_buffer)) then
         call c_f_pointer(c_loc(send_buffer), send_values, [mpi_number])
         call mpi_reduce(send_values, receive_values, mpi_number, mpi_double_complex, mpi_sum, &
                         mpi_rank, communicator, error_code)
      else
         allocate (local_buffer(mpi_number))
         local_buffer = receive_values
         call mpi_reduce(local_buffer, receive_values, mpi_number, mpi_double_complex, mpi_sum, &
                         mpi_rank, communicator, error_code)
      end if
   end subroutine reduce_sum_complex64

   subroutine allreduce_sum_real64(receive_buffer, mpi_number, mpi_comm, send_buffer)
      real(real64), contiguous, target, intent(inout) :: receive_buffer(..)
      real(real64), contiguous, optional, target, intent(in) :: send_buffer(..)
      integer, intent(in) :: mpi_number
      integer, optional, intent(in) :: mpi_comm
      integer :: communicator, error_code
      real(real64), allocatable :: local_buffer(:)
      real(real64), pointer :: receive_values(:), send_values(:)

      communicator = mpi_comm_world
      if (present(mpi_comm)) communicator = mpi_comm
      call c_f_pointer(c_loc(receive_buffer), receive_values, [mpi_number])
      allocate (local_buffer(mpi_number))
      if (present(send_buffer)) then
         call c_f_pointer(c_loc(send_buffer), send_values, [mpi_number])
         local_buffer = send_values
      else
         local_buffer = receive_values
      end if
      call mpi_allreduce(local_buffer, receive_values, mpi_number, mpi_double_precision, mpi_sum, communicator, error_code)
   end subroutine allreduce_sum_real64

   subroutine allreduce_sum_complex64(receive_buffer, mpi_number, mpi_comm, send_buffer)
      complex(real64), contiguous, target, intent(inout) :: receive_buffer(..)
      complex(real64), contiguous, optional, target, intent(in) :: send_buffer(..)
      integer, intent(in) :: mpi_number
      integer, optional, intent(in) :: mpi_comm
      integer :: communicator, error_code
      complex(real64), allocatable :: local_buffer(:)
      complex(real64), pointer :: receive_values(:), send_values(:)

      communicator = mpi_comm_world
      if (present(mpi_comm)) communicator = mpi_comm
      call c_f_pointer(c_loc(receive_buffer), receive_values, [mpi_number])
      allocate (local_buffer(mpi_number))
      if (present(send_buffer)) then
         call c_f_pointer(c_loc(send_buffer), send_values, [mpi_number])
         local_buffer = send_values
      else
         local_buffer = receive_values
      end if
      call mpi_allreduce(local_buffer, receive_values, mpi_number, mpi_double_complex, mpi_sum, communicator, error_code)
   end subroutine allreduce_sum_complex64

   subroutine parallel_allreduce_sum_complex64_sequence(receive_buffer, mpi_number, mpi_comm)
      complex(real64), intent(inout) :: receive_buffer(*)
      integer, intent(in) :: mpi_number
      integer, optional, intent(in) :: mpi_comm
      integer :: communicator, error_code
      complex(real64), allocatable :: local_buffer(:)

      communicator = mpi_comm_world
      if (present(mpi_comm)) communicator = mpi_comm
      allocate (local_buffer(mpi_number))
      local_buffer = receive_buffer(1:mpi_number)
      call mpi_allreduce(local_buffer, receive_buffer, mpi_number, mpi_double_complex, mpi_sum, communicator, error_code)
   end subroutine parallel_allreduce_sum_complex64_sequence

   subroutine parallel_allreduce_max(receive_buffer, mpi_number, mpi_comm, send_buffer)
      integer, contiguous, target, intent(inout) :: receive_buffer(..)
      integer, contiguous, optional, target, intent(in) :: send_buffer(..)
      integer, intent(in) :: mpi_number
      integer, optional, intent(in) :: mpi_comm
      integer :: communicator, error_code
      integer, allocatable :: local_buffer(:)
      integer, pointer :: receive_values(:), send_values(:)

      communicator = mpi_comm_world
      if (present(mpi_comm)) communicator = mpi_comm
      call c_f_pointer(c_loc(receive_buffer), receive_values, [mpi_number])
      allocate (local_buffer(mpi_number))
      if (present(send_buffer)) then
         call c_f_pointer(c_loc(send_buffer), send_values, [mpi_number])
         local_buffer = send_values
      else
         local_buffer = receive_values
      end if
      call mpi_allreduce(local_buffer, receive_values, mpi_number, mpi_integer, mpi_max, communicator, error_code)
   end subroutine parallel_allreduce_max

   subroutine send_complex64(send_buffer, mpi_number, mpi_rank, mpi_comm)
      complex(real64), contiguous, target, intent(in) :: send_buffer(..)
      integer, intent(in) :: mpi_number, mpi_rank
      integer, optional, intent(in) :: mpi_comm
      integer :: communicator, error_code
      complex(real64), pointer :: send_values(:)

      communicator = mpi_comm_world
      if (present(mpi_comm)) communicator = mpi_comm
      call c_f_pointer(c_loc(send_buffer), send_values, [mpi_number])
      call mpi_send(send_values, mpi_number, mpi_double_complex, mpi_rank, 1, communicator, error_code)
   end subroutine send_complex64

   subroutine receive_complex64(receive_buffer, mpi_number, mpi_rank, mpi_comm)
      complex(real64), contiguous, target, intent(out) :: receive_buffer(..)
      integer, intent(in) :: mpi_number, mpi_rank
      integer, optional, intent(in) :: mpi_comm
      integer :: communicator, error_code, status(mpi_status_size)
      complex(real64), pointer :: receive_values(:)

      communicator = mpi_comm_world
      if (present(mpi_comm)) communicator = mpi_comm
      call c_f_pointer(c_loc(receive_buffer), receive_values, [mpi_number])
      call mpi_recv(receive_values, mpi_number, mpi_double_complex, mpi_rank, 1, communicator, status, error_code)
   end subroutine receive_complex64

end module parallel_runtime
