!! A 1D example of a solver that advance the Fourier modes in time
program test1dFourier
  use fftw3_wrapper
  use SCIFOR, only: print_matrix
  
  implicit none

  integer :: numModes
  real(kind=8) :: h
  real(kind=8),allocatable :: x(:)
  real(kind=8) :: alpha,epsilon !The convection coefficients.
  complex(kind=8),allocatable :: convMatrix(:,:)
  complex(kind=8),allocatable :: diffMatrix(:,:)
  real(8),parameter :: pi = 3.14159265358979323846
  complex(kind=8),allocatable :: uHatStages(:,:) !Solution at every stage [explicit stage]
  complex(kind=8),allocatable :: uHatTimes(:,:) !Solution at every time instance
  real(kind=8),allocatable :: uPhysTimes(:,:) !Solution [physical] at every time instance

  complex(kind=8),allocatable :: tmpuPhys(:);
  complex(kind=8),allocatable :: tmpuHat(:);

  !Time-stepping scheme tableu:
  integer :: numStages = 4
  real(kind=8) :: A(4,4)
  real(kind=8) :: b(4)
  real(kind=8) :: c(4)
  integer :: numSteps,step,stage,i,j
  real(kind=8) :: tstart,tend,dt,curr_time
  !Time advancing stuffs:
  tstart = 0.0
  tend   = 8.0
  numSteps= 200
  dt = (tend-tstart)/real(numSteps)

  numModes = 24
  allocate( uHatTimes(1:numSteps,1:numModes), uPhysTimes(1:numSteps,1:numModes), &
       uHatStages(1:numStages,1:numModes))

  !Allocate temporary vectors:
  allocate(tmpuPhys(1:numModes),tmpuHat(1:numModes))

  alpha = 0.5
  epsilon = 0.0
  
  convMatrix = FourierAdvectionMatrix(numModes,alpha);
  diffMatrix = FourierDiffusionMatrix(numModes,epsilon);

  !write(*,'(A)')'-------------------------------------'
  !call print_matrix(convMatrix)
  !write(*,'(A)')'-------------------------------------'
  !call print_matrix(diffMatrix)

  write(*,'(A)')'=========> STARTING THE TIME-STEPPING SCHEME========>'

  !Construct the butcher tableu:
  call ConstructRK4(A,b,c,4)



  !Mesh:
  allocate(x(numModes))
  h = 2.0*pi/real(numModes,8)
  do i=1,numModes
     x(i) = (i-1)*h
  enddo
  curr_time = tstart
  step=1
  uPhysTimes(1,:) = u0(x)
  !print*,uPhysTimes(step,:)
  !Initial condition
  !call fft(uHatTimes(step,:),uPhysTimes(step,:))
  tmpuPhys = u0(x)


  call fft(tmpuHat,tmpuPhys,num=[numModes])
  uHatTimes(1,:) = tmpuHat
  
  do while(curr_time<tend)
     uHatStages(1,:) = uHatTimes(step,:)
     do stage = 2,numStages
        uHatStages(stage,:) = uHatTimes(step,:)
        do j = 1,stage-1
           uHatStages(stage,:) = uHatStages(stage,:)+dt*A(stage,j)*matmul(convMatrix,uHatStages(j,:))
        enddo
     enddo !end stages
     uHatTimes(step+1,:) = uHatTimes(step,:)
     do stage = 1,numStages
        uHatTimes(step+1,:) = uHatTimes(step+1,:)+dt*b(stage)*matmul(convMatrix,uHatStages(stage,:))
     enddo

     
     tmpuHat = uHatTimes(step+1,:)

     deallocate(tmpuPhys)
     call ifft(tmpuPhys,tmpuHat,num=[numModes])
     uPhysTimes(step+1,:) = real(tmpuPhys,8) !copy
     step=step+1
     curr_time = curr_time + dt
     
  enddo
  
  
  write(*,'(A)')'=========> ENDING THE TIME-STEPPING SCHEME==========>'

  !Here we write the solution:
  open(unit=999,file='sol.csv',status='replace',action='write')
  
  do i=1,size(uPhysTimes,2)
     do j=1,6
        
        if(j.eq.1) then
           write(999,'((G0.8,","))',advance='no') uPhysTimes(1,i)
        else
           write(999,'((G0.8,","))',advance='no') uPhysTimes((j-1)*40,i) 
        endif
 
     enddo
     write(999,*)
  enddo
  close(999)
  
contains

  !> The initial solution here:
  !!
  !!
  !!
  function u0(x)
    implicit none
    real,allocatable :: u0(:);
    real(8) :: x(:)
    allocate(u0(size(x)))
    u0 = exp(-1.0*(x-pi)**2)
  end function u0
  !> Assemble the advection matrix given
  !! the number of fourier modes
  !! 
  !! 
  !!
  function FourierAdvectionMatrix(modes,alpha) result(F)
    use SCIFOR, only : diag
    use fftpack, only : fftfreq
    implicit none
    real(kind=8),intent(in) :: alpha;
    integer,intent(in) :: modes;
    complex(kind=8) :: ks(modes) !The frequencies
    complex(kind=8),allocatable :: F(:,:)
    complex(kind=8) :: ju = (0.0,1.0) !Unitary vector
    ks = cmplx(fftfreq(modes))!Get the frequencies
    F = diag(-ju*alpha*ks);
  end function FourierAdvectionMatrix

  !> @brief Evaluate the diffusion matrix
  !!
  !! @param[in] modes The modes considered
  !! @param[in] The viscosity coefficient
  function FourierDiffusionMatrix(modes,epsilon) result(Fsquare)
    use SCIFOR, only : diag
    use fftpack, only : fftfreq
    implicit none
    real(kind=8),intent(in) :: epsilon
    integer,intent(in) :: modes;
    complex(kind=8) :: ks(modes) !The frequencies
    complex(kind=8),allocatable :: Fsquare(:,:)

    !allocate(Fsquare(modes,modes))
    ks = cmplx(fftfreq(modes))!Get the frequencies
    Fsquare = diag(-epsilon*ks*ks);
  end function FourierDiffusionMatrix

  !> @brief Timestepping scheme
  !! Here we construct the Matrixes correspondant
  !! to the runge-kutta tableau.
  subroutine ConstructRK4(A,b,c,numStages)
    integer,intent(in) :: numStages
    real(kind=8),intent(inout) :: A(numStages,numStages)
    real(kind=8),intent(inout):: b(numStages)
    real(kind=8),intent(inout):: c(numStages)

    A(:,:) = 0.0
    A(2,1) = 0.5
    A(3,2) = 0.5
    A(4,3) = 1.0

    b(1) = 1.0/6.0
    b(2) = 1.0/3.0
    b(3) = 1.0/3.0
    b(4) = 1.0/6.0

    c(1) = 0.0
    c(2) = 0.5
    c(3) = 0.5 
    c(4) = 1.0
  end subroutine ConstructRK4
end program test1dFourier
