!! A 1D example of a solver that advance the Fourier modes in time
!!
!!
!!
!!
program test1dFourier
  use fftw3_wrapper
  use SCIFOR, only: print_matrix
  implicit none

  integer :: numModes
  real(kind=8) :: a,c
  complex(kind=8),allocatable :: convMatrix(:,:)
  complex(kind=8),allocatable :: diffMatrix(:,:)
  real(8),parameter :: pi = 3.14159265358979323846

  convMatrix = FourierAdvectionMatrix(6);
  diffMatrix = FourierDiffusionMatrix(6);

  write(*,'(A)')'-------------------------------------'
  call print_matrix(convMatrix)
  write(*,'(A)')'-------------------------------------'
  call print_matrix(diffMatrix)

  write(*,'(A)')'=========> STARTING THE TIME-STEPPING SCHEME========>'

  

  write(*,'(A)')'=========> ENDING THE TIME-STEPPING SCHEME==========>'
  
contains
  !> Assemble the advection matix given
  !! the number of fourier modes
  !! 
  !!
  !!
  function FourierAdvectionMatrix(modes) result(F)
    use SCIFOR, only : diag
    use fftpack, only : fftfreq
    implicit none
    integer,intent(in) :: modes;
    complex(kind=8) :: ks(modes) !The frequencies
    complex(kind=8),allocatable :: F(:,:)
    complex(kind=8) :: ju = (0.0,1.0) !Unitary vector
    ks = cmplx(fftfreq(modes))!Get the frequencies
    !allocate(F(modes,modes))
    F = diag(ks);

  end function FourierAdvectionMatrix

  !> Evaluate the diffusion matrix
  !!
  !!
  function FourierDiffusionMatrix(modes) result(Fsquare)
    use SCIFOR, only : diag
    use fftpack, only : fftfreq
    implicit none
    integer,intent(in) :: modes;
    complex(kind=8) :: ks(modes) !The frequencies
    complex(kind=8),allocatable :: Fsquare(:,:)

    !allocate(Fsquare(modes,modes))
    ks = cmplx(fftfreq(modes))!Get the frequencies
    Fsquare = diag(-ks*ks);
  end function FourierDiffusionMatrix

  !> Timestepping scheme
  !!
  !!
  !!
  subroutine timeSteppingScheme()
    
  end subroutine timeSteppingScheme
end program test1dFourier
