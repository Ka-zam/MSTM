module parallel_runtime
   use, intrinsic :: iso_c_binding, only: c_f_pointer, c_loc
   use, intrinsic :: iso_fortran_env, only: real32, real64
   implicit none(type, external)
   private

   integer, public :: mpi_comm_world = 1
   integer, public :: mpi_comm_null = 0
   integer, public :: mstm_global_rank = 0
   integer, public :: mstm_global_numprocs = 1

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
      mstm_global_rank = 0
      mstm_global_numprocs = 1
   end subroutine parallel_initialize

   subroutine parallel_finalize()
   end subroutine parallel_finalize

   real(real64) function parallel_wall_time()
      call cpu_time(parallel_wall_time)
   end function parallel_wall_time

   subroutine parallel_rank(mpi_rank, mpi_comm)
      integer, intent(out) :: mpi_rank
      integer, optional, intent(in) :: mpi_comm

      mpi_rank = 0
   end subroutine parallel_rank

   subroutine parallel_size(mpi_size, mpi_comm)
      integer, intent(out) :: mpi_size
      integer, optional, intent(in) :: mpi_comm

      mpi_size = 1
   end subroutine parallel_size

   subroutine parallel_barrier(mpi_comm)
      integer, optional, intent(in) :: mpi_comm
   end subroutine parallel_barrier

   subroutine parallel_split(mpi_color, mpi_key, mpi_new_comm, mpi_comm)
      integer, intent(in) :: mpi_color, mpi_key
      integer, intent(out) :: mpi_new_comm
      integer, optional, intent(in) :: mpi_comm

      mpi_new_comm = mpi_comm_world
   end subroutine parallel_split

   subroutine parallel_group(mpi_group, mpi_comm)
      integer, intent(out) :: mpi_group
      integer, optional, intent(in) :: mpi_comm

      mpi_group = 0
   end subroutine parallel_group

   subroutine parallel_group_include(mpi_group, mpi_size, mpi_new_group_list, mpi_new_group, mpi_comm)
      integer, intent(in) :: mpi_group, mpi_size, mpi_new_group_list(*)
      integer, intent(out) :: mpi_new_group
      integer, optional, intent(in) :: mpi_comm

      mpi_new_group = 0
   end subroutine parallel_group_include

   subroutine parallel_communicator_create(mpi_group, mpi_new_comm, mpi_comm)
      integer, intent(in) :: mpi_group
      integer, intent(out) :: mpi_new_comm
      integer, optional, intent(in) :: mpi_comm

      mpi_new_comm = mpi_comm_world
   end subroutine parallel_communicator_create

   subroutine broadcast_integer(send_buffer, mpi_number, mpi_rank, mpi_comm)
      integer, contiguous, target, intent(inout) :: send_buffer(..)
      integer, intent(in) :: mpi_number, mpi_rank
      integer, optional, intent(in) :: mpi_comm
   end subroutine broadcast_integer

   subroutine broadcast_real64(send_buffer, mpi_number, mpi_rank, mpi_comm)
      real(real64), contiguous, target, intent(inout) :: send_buffer(..)
      integer, intent(in) :: mpi_number, mpi_rank
      integer, optional, intent(in) :: mpi_comm
   end subroutine broadcast_real64

   subroutine broadcast_complex64(send_buffer, mpi_number, mpi_rank, mpi_comm)
      complex(real64), contiguous, target, intent(inout) :: send_buffer(..)
      integer, intent(in) :: mpi_number, mpi_rank
      integer, optional, intent(in) :: mpi_comm
   end subroutine broadcast_complex64

   subroutine reduce_sum_integer(receive_buffer, mpi_number, mpi_rank, mpi_comm, send_buffer)
      integer, contiguous, target, intent(inout) :: receive_buffer(..)
      integer, contiguous, optional, target, intent(in) :: send_buffer(..)
      integer, intent(in) :: mpi_number, mpi_rank
      integer, optional, intent(in) :: mpi_comm
      integer, pointer :: receive_values(:), send_values(:)

      if (present(send_buffer)) then
         call c_f_pointer(c_loc(receive_buffer), receive_values, [mpi_number])
         call c_f_pointer(c_loc(send_buffer), send_values, [mpi_number])
         receive_values = send_values
      end if
   end subroutine reduce_sum_integer

   subroutine reduce_sum_complex32(receive_buffer, mpi_number, mpi_rank, mpi_comm, send_buffer)
      complex(real32), contiguous, target, intent(inout) :: receive_buffer(..)
      complex(real32), contiguous, optional, target, intent(in) :: send_buffer(..)
      integer, intent(in) :: mpi_number, mpi_rank
      integer, optional, intent(in) :: mpi_comm
      complex(real32), pointer :: receive_values(:), send_values(:)

      if (present(send_buffer)) then
         call c_f_pointer(c_loc(receive_buffer), receive_values, [mpi_number])
         call c_f_pointer(c_loc(send_buffer), send_values, [mpi_number])
         receive_values = send_values
      end if
   end subroutine reduce_sum_complex32

   subroutine reduce_sum_real64(receive_buffer, mpi_number, mpi_rank, mpi_comm, send_buffer)
      real(real64), contiguous, target, intent(inout) :: receive_buffer(..)
      real(real64), contiguous, optional, target, intent(in) :: send_buffer(..)
      integer, intent(in) :: mpi_number, mpi_rank
      integer, optional, intent(in) :: mpi_comm
      real(real64), pointer :: receive_values(:), send_values(:)

      if (present(send_buffer)) then
         call c_f_pointer(c_loc(receive_buffer), receive_values, [mpi_number])
         call c_f_pointer(c_loc(send_buffer), send_values, [mpi_number])
         receive_values = send_values
      end if
   end subroutine reduce_sum_real64

   subroutine reduce_sum_complex64(receive_buffer, mpi_number, mpi_rank, mpi_comm, send_buffer)
      complex(real64), contiguous, target, intent(inout) :: receive_buffer(..)
      complex(real64), contiguous, optional, target, intent(in) :: send_buffer(..)
      integer, intent(in) :: mpi_number, mpi_rank
      integer, optional, intent(in) :: mpi_comm
      complex(real64), pointer :: receive_values(:), send_values(:)

      if (present(send_buffer)) then
         call c_f_pointer(c_loc(receive_buffer), receive_values, [mpi_number])
         call c_f_pointer(c_loc(send_buffer), send_values, [mpi_number])
         receive_values = send_values
      end if
   end subroutine reduce_sum_complex64

   subroutine allreduce_sum_real64(receive_buffer, mpi_number, mpi_comm, send_buffer)
      real(real64), contiguous, target, intent(inout) :: receive_buffer(..)
      real(real64), contiguous, optional, target, intent(in) :: send_buffer(..)
      integer, intent(in) :: mpi_number
      integer, optional, intent(in) :: mpi_comm
      real(real64), pointer :: receive_values(:), send_values(:)

      if (present(send_buffer)) then
         call c_f_pointer(c_loc(receive_buffer), receive_values, [mpi_number])
         call c_f_pointer(c_loc(send_buffer), send_values, [mpi_number])
         receive_values = send_values
      end if
   end subroutine allreduce_sum_real64

   subroutine allreduce_sum_complex64(receive_buffer, mpi_number, mpi_comm, send_buffer)
      complex(real64), contiguous, target, intent(inout) :: receive_buffer(..)
      complex(real64), contiguous, optional, target, intent(in) :: send_buffer(..)
      integer, intent(in) :: mpi_number
      integer, optional, intent(in) :: mpi_comm
      complex(real64), pointer :: receive_values(:), send_values(:)

      if (present(send_buffer)) then
         call c_f_pointer(c_loc(receive_buffer), receive_values, [mpi_number])
         call c_f_pointer(c_loc(send_buffer), send_values, [mpi_number])
         receive_values = send_values
      end if
   end subroutine allreduce_sum_complex64

   subroutine parallel_allreduce_sum_complex64_sequence(receive_buffer, mpi_number, mpi_comm)
      complex(real64), intent(inout) :: receive_buffer(*)
      integer, intent(in) :: mpi_number
      integer, optional, intent(in) :: mpi_comm
   end subroutine parallel_allreduce_sum_complex64_sequence

   subroutine parallel_allreduce_max(receive_buffer, mpi_number, mpi_comm, send_buffer)
      integer, contiguous, target, intent(inout) :: receive_buffer(..)
      integer, contiguous, optional, target, intent(in) :: send_buffer(..)
      integer, intent(in) :: mpi_number
      integer, optional, intent(in) :: mpi_comm
      integer, pointer :: receive_values(:), send_values(:)

      if (present(send_buffer)) then
         call c_f_pointer(c_loc(receive_buffer), receive_values, [mpi_number])
         call c_f_pointer(c_loc(send_buffer), send_values, [mpi_number])
         receive_values = send_values
      end if
   end subroutine parallel_allreduce_max

   subroutine send_complex64(send_buffer, mpi_number, mpi_rank, mpi_comm)
      complex(real64), contiguous, target, intent(in) :: send_buffer(..)
      integer, intent(in) :: mpi_number, mpi_rank
      integer, optional, intent(in) :: mpi_comm
   end subroutine send_complex64

   subroutine receive_complex64(receive_buffer, mpi_number, mpi_rank, mpi_comm)
      complex(real64), contiguous, target, intent(out) :: receive_buffer(..)
      integer, intent(in) :: mpi_number, mpi_rank
      integer, optional, intent(in) :: mpi_comm
   end subroutine receive_complex64

end module parallel_runtime
