!> @brief A subroutine that test the binded version
!!   of mpi in python.
!! 
!!
subroutine simple
  
  use mpi
  !include 'mpif.h'
  integer numtasks, rank, ierr, rc

  call MPI_INIT(ierr)
  if (ierr .ne. MPI_SUCCESS) then
     print *,'Error starting MPI program. Terminating.'
     call MPI_ABORT(MPI_COMM_WORLD, rc, ierr)
  end if

  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  
  if(rank==0) then
     call MPI_COMM_SIZE(MPI_COMM_WORLD, numtasks, ierr)
     write(*,'(A,1x,i2)')'Number of tasks:',numtasks
  endif
  write(*,'(A,1x,i2)') 'Rank: ',rank
  
  call MPI_FINALIZE(ierr)

end subroutine simple
