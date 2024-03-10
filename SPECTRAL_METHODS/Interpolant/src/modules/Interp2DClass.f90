module FourierInterp2dClass
  use fftw3_wrapper, only : fft, ifft !Use the wrapper
  use fftpack_kind, only: rk
  implicit none

  !Basic type for interpolate a 2d function.
  type FourierInterp2d
     !private
     integer  :: m_dim = 2
     real(rk),allocatable :: y(:,:);
     complex(rk),allocatable ::yHat(:,:);
     logical :: m_transformed = .false.;
  end type FourierInterp2d

  real(rk),parameter :: pi =  3.14159265358979323846;!pi
  complex,parameter :: ju = (0.0,1.0);   !img unity

contains

  !>@brief Basic constructor
  !@param Nx: The number of points along x
  !@param Ny: The number of points along y
  !@param xy: The matrix which contains the fitted points for
  !           the function.
  function FourierInterp2d_(Nx,Ny,fxy) result (Iobj)
    implicit none
    integer,intent(in) :: Nx,Ny; !Number of discrete values.
    real(rk),intent(in) :: fxy(Nx,Ny); !Vector of discrete values.
    type (FourierInterp2d) :: Iobj; !Interpolation object

    !Construct the object by assuming
    !that the fourier interpolant is y.

    allocate(Iobj%y(Nx,Ny),Iobj%yHat(Nx,Ny));

    Iobj%y = fxy
    Iobj%m_transformed = .false.
  end function FourierInterp2d_
!----------------------------------------------------
  !>@brief Performs the 2d transform for the object
  !and store it in yHat
  subroutine FwdTransform2d(Iobj)
    implicit none
    class (FourierInterp2d) :: Iobj
    call fft(Iobj%yHat,Iobj%y,shape(Iobj%y));
  end subroutine FwdTransform2d

  !>@brief Interpolate in a specified point
  !>@param [in] xy: An array which contains
  ! the coordinates where we want to interpolate
  !>@param [in] 
  subroutine FourierInterpolate2D_scalar(Iobj,xy,yInt)
    use fftpack, only : fftfreq !Use the frequencies
    implicit none
    class (FourierInterp2d),intent(in) :: Iobj;
    real(rk),dimension(:) :: xy !The points where we want to interpolate
    real(rk),intent(out) :: yInt;
    integer :: i,n(2)
    !integer,allocatable :: kfreqs(:,:),lfreqs(:,:);
    complex(rk),allocatable :: expK(:,:),expL(:,:);
    n = shape(Iobj%y) !Get the dimensions
    

    allocate(expK(n(1),n(2)),expL(n(1),n(2))); !allocate
    !Get the discrete exponential in the specific points:
    do i =1,n(2)
       expK(:,i) = exp(ju*xy(1)*fftfreq(n(1)))

    enddo
    do i=1,n(1)
       expL(i,:) = exp(ju*xy(2)*fftfreq(n(2)))
    enddo
    yInt = real(sum(Iobj%yHat*expK*expL))/(n(1)*n(2));
  end subroutine FourierInterpolate2D_scalar
  
  
end module FourierInterp2dClass
