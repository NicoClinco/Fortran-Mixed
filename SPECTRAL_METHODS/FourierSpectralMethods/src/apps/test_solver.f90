  !! A 1D example of a solver of a an advection equation
  !! with the Fourier collocation method.

! -------------------------------------------------!
!>
!! \begin{equation}
!! 
!! \frac{\partial{u}}{\partial{t}} +
!!    c(x) \frac{\partial{u}}{\partial{x}} = 0 
!!
!! \end{equation}
!!
!! The derivative is computed by the pseudospectral

!! method, and the time step is EXP. EULER 1 order
!--------------------------------------------------!


program main
 
  use fftw3_wrapper
  use fftpack, only : fftfreq
  implicit none
  
  !Compute the rhs 
  real(8),allocatable :: u(:,:)
  real(8),allocatable :: c(:)
  real(8),allocatable :: x(:)

  integer :: i,j,nsteps,nspace
  real(8) :: t0,tfin,dt,curr_time,h

  real(8),parameter :: pi = 3.14159265358979323846


  nsteps = 10000
  nspace = 160
  h = 2.0*pi/real(nspace,8)
  
  allocate(x(nspace),c(nspace),u(nspace,nsteps+1))

  do i=1,nspace
     x(i) = (i-1)*h
  enddo

  c = Conv(x) !Convective
  u = 0.0
  
  u(:,1) = u0(x)   !Initial condition

  t0 = 0.0
  tfin = 8.0
  dt =(tfin-t0)/nsteps
  j = 1
  do while (curr_time < tfin .and. j<=nsteps)

     u(:,j+1) = u(:,j) - dt*rhs(c,u(:,j))
     print '("Current time step: ",f5.3)',curr_time+dt
     curr_time = curr_time + dt
     j=j+1
  end do

  u = transpose(u)
  open(unit=999,file='Solution.csv',status='replace', action="write")
  do i=1,size(u,2)
     do j=1,size(u,1),10
        write(999,'((G0.8,","))',advance='no') u(j,i)
     enddo
     write(999,*)
  enddo
  close(unit=999)
  
  
  
  
  
  
contains
  function rhs(c,u)
    real(8), allocatable :: rhs(:)
    real(8),intent(in) :: c(:); !The convection velocity
    real(8),intent(in) :: u(:); !The solution
    !real(8),allocatable :: D(:,:) !The matrix of the derivatives
    integer :: n
    complex(8),allocatable :: vHat(:);
    complex(8),allocatable :: dudx(:);
    complex,parameter :: ju = (0.0,1.0);

    n = size(c)
    allocate(vHat(n),rhs(n),dudx(n))

    call fft(vHat,u)
    
    vHat = ju*fftfreq(n)*vHat

    call ifft(dudx,vHat)
    
    rhs = c*real(dudx,8)
    
    
  end function rhs

  !Initial solution
  function u0(x)
    implicit none
    real,allocatable :: u0(:);
    real(8) :: x(:)
    allocate(u0(size(x)))
    u0 = exp(-100.0*(x-1.0)**2)
  end function u0

  !Convective velocity:
  function Conv(x)
    implicit none
    real(8) :: x(:)
    real(8),allocatable :: Conv(:);
    allocate(Conv(size(x)))

    !Convective velocity:
    Conv = 0.1+(sin(x-1.0))**2

  end function Conv


end program main
