program main
  
  implicit none
  include 'mpif.h'
  integer rank, num_procs, ierror, tag ,status(MPI_STATUS_SIZE);
  integer source,dest;
  integer iproc;
  character(len=25) :: message;

  call MPI_INIT(ierror);

  !Current rank:
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror);


  if(rank /= 0 ) then
     !Create a stringstream:
     write (message, "(A18,i1)") "hello-word from:",rank;
     message = adjustl(message);
     message = adjustr(message);

     call MPI_SEND(message,25,MPI_CHAR,0, &
          0, MPI_COMM_WORLD,ierror);
  else
     !Here we have the 'base' node
     !Number of processes:
     call MPI_COMM_SIZE(MPI_COMM_WORLD,num_procs,ierror);

     
     do iproc = 1, num_procs-1
        
        call MPI_RECV(message,25, MPI_CHAR, iproc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror);
        print*,message;
     enddo

  endif

  call MPI_FINALIZE(ierror);

end program main
