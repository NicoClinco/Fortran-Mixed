! In this file it is shown
! a different version for the trapezoidal rule:

!We use the reduce operation for summing all the integrals
function fun(x)
  real(kind=8) :: f !The result
  real(kind=8),intent(in) :: x;
  !Definition of the function
  fun = sin(x);  
end function fun


function fvec(x)
  real(kind=8) :: fvec;
  real(kind=8):: x(:);
  fvec = sin(sum(x)); !sin(x^2+y^2+z^2) !For example
end function fvec





!>@brief subroutine for computing the trapezoidal rule
!@param [in] a_loc local starting position
!@param [in] b_loc local end position on the interval
!@param [in] h_loc The step-size considered
!@param [in] n_loc The number of local points
!@param [in] f [function] It is the function to integrate
subroutine trpInt(a_loc,b_loc,n_loc,h_loc,f,res_loc)
  implicit none
  !real(kind=8),intent(in),dimension(:) :: xVec_loc;
  !real(kind=8) h;

  real(kind=8),intent(in) :: a_loc,b_loc,h_loc;
  integer,intent(in) :: n_loc;
  interface 
     function f(x)
       real(kind=8) :: f;
       real(kind=8),intent(in) :: x;
     end function f
  end interface

  real(kind=8),intent(out) :: res_loc;
  real(kind=8) :: xi;
  integer :: i;

  res_loc = (f(a_loc) + f(b_loc))/2.0;
  xi = a_loc;
  do i=1,n_loc-1
     xi = xi + h_loc;
     res_loc = res_loc+f(xi);
  enddo

  res_loc = res_loc*h_loc;

end subroutine trpInt

! Note: In this example we are filling the array xi
! in every process, this is a mess. Don't do it.
program main
  implicit none  
  external trpInt; !External routine for computing the integrak
  real(8), external :: fun;
  include 'mpif.h'

  integer :: n,n_loc,rem; 
  real(8) :: I = 0.0;
  real(8) :: Iloc = 0.0; !Integral value on each rank
  real(8) :: a,b,a_loc,b_loc,h_loc; !Global and local extrema
  
  ! MPI declarations:
  integer :: source,dest,nprocs,ierror,tag,status(MPI_STATUS_SIZE),rank
  integer :: n_loc_rk0;
  integer :: iproc; !counter for processes
  real(8) :: other_rank;

  integer :: bcast_sign;
  
  ! MPI PART START HERE !
  call MPI_INIT(ierror); !!Initialize MPI

  !Get the current rank
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror);


  n = 13;
  a = 0.0;
  b = 1.0;

  !Number of processes:
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierror);

  n_loc = int(n/nprocs);
  rem = mod(n,nprocs); !Compute the remainder

  n_loc_rk0 = n_loc+rem; !Rank-0 compute also the remainder.

  h_loc = (b - a)/n;

  
  if(rank == 0) then
     n_loc = n_loc_rk0;
     a_loc = a + rank*n_loc_rk0*h_loc;
     b_loc = a_loc + n_loc_rk0*h_loc;
     bcast_sign = 31;
  else

     a_loc = a+n_loc_rk0*h_loc + (rank-1)*n_loc*h_loc;
     b_loc = a_loc+ n_loc*h_loc;
     !print '(A,1X,i1,1X,A,1X,i2)','rank:',rank,'bcast_sign-0',bcast_sign;
     ! This print does not work: this line of code
     ! is placed before MPI_BCAST, thus, we cannot determine its behaviour.
  endif

  !Suppose to send to all the processes bcast_sign
  call MPI_BCAST(bcast_sign,1,MPI_INT,0,MPI_COMM_WORLD,ierror);

  print '(A,1X,i1,1X,A,1X,i2)','rank:',rank,'bcast_sign-1',bcast_sign;
  
  
  !print '(A,1X,2(i1,1x),2(f5.2))','rank :',rank,n_loc,a_loc,b_loc
  !Depending where we are, we can decide if perform
  !the summation on proceeding further:
  
  !if(rank == 0)  then
  !
  !   do iproc =1,nprocs-1
  !      other_rank = 0.0;
  !      call MPI_RECV(other_rank,1,MPI_DOUBLE,iproc,0, &
  !          MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierror);
  !      I = I + other_rank;
  !       print '("Value of the integral",1X,f8.4)',I;
  !   enddo
  ! else
  !    Send the value of the integral to the rank 0
  !    call MPI_SEND(I,1,MPI_DOUBLE,0, &
  !         0,MPI_COMM_WORLD,ierror);
  ! endif
  

  !Compute the integral for each rank:
  call trpInt(a_loc,b_loc,n_loc,h_loc,fun,Iloc)
  !The previous block is substituted by this function:
 
  call MPI_REDUCE(Iloc,I,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD,ierror);
  if(rank == 0) print '("Value of the integral",1X,f8.4)',I;
  
  call MPI_FINALIZE(ierror);

end program main
