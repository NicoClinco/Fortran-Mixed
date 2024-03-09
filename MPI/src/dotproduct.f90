!!! Examples of dotproduct in parallel



program main
  implicit none
  real(8) ::y(20),x(20);
  real(8) :: res,loc_res;
  integer :: i;

  integer :: ierror,rank,n,nprocs;
  integer :: loc_start;
  integer :: loc_end;
  include 'mpif.h'
  
  do i =1,20
     y(i) = real(i)/20.0;
  enddo
  x = y*10;
  
  call MPI_INIT(ierror);

  loc_res = 0.0;
  res= 0.0;
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror);

  call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierror);
  n = 20/nprocs;
  loc_start = 1 + rank*n;
  loc_end =   loc_start + (n-1);
  call dotprod(x(loc_start:loc_end),y(loc_start:loc_end),loc_res);
  
  call MPI_REDUCE(loc_res,res,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD,ierror);

  if(rank == 0) print '(A,1X,f8.3)','Result: ',res;

  call MPI_FINALIZE(ierror);
  
contains

  subroutine dotprod(x,y,res)
    implicit none;
    real(8),intent(in) :: x(:),y(:);
    real(8),intent(out) :: res;

    !...Assume that everything is correct...

    !easy:
    res = sum(x*y);
  end subroutine dotprod
end program main
