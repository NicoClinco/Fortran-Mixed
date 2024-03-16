module FourierInterp2dClass
  use fftw3_wrapper, only : fft, ifft !Use the wrapper
  use fftpack_kind, only: rk
  implicit none

  
  interface Interpolate2d
     module procedure Fourierinterpolate2D_scalar,Fourierinterpolate2D_vector
  end interface Interpolate2d
  !Basic type for interpolate a 2d function.
  type FourierInterp2d
     !private
     integer  :: m_dim = 2
     real(rk),allocatable :: y(:,:);
     complex(rk),allocatable ::yHat(:,:);
     logical :: m_transformed = .false.;
   contains
     procedure FwdTransform2d
     procedure :: interpolate2ds => FourierInterpolate2D_scalar
     procedure :: interpolate2dv => FourierInterpolate2D_vector
  end type FourierInterp2d


  type, extends(FourierInterp2d) :: FourierDerivative2d
     !private
     complex(rk),allocatable :: dydxHat(:,:);
     complex(rk),allocatable :: dydyHat(:,:);
     logical :: m_derivtransformed;
   contains
     !Note that the procedure are inherited from the base class
     procedure :: grad2dXs => FourierGradX_scalar
     procedure :: grad2dYs => FourierGradY_scalar
  end type FourierDerivative2d

  real(rk),parameter :: pi =  3.14159265358979323846;!pi
  complex,parameter :: ju = (0.0,1.0);   !img unity

contains

  !>@brief Basic constructor for the derivative object
  !@param Nx: The number of points along x
  !@param Ny: The number of points along y
  !@param xy: The matrix which contains the fitted points for
  !           the function.
  function DerFourierInterp2d_(Nx,Ny,fxy) result (Dobj)
    implicit none
    integer,intent(in) :: Nx,Ny; !Number of discrete values.
    real(rk),intent(in) :: fxy(Nx,Ny); !Vector of discrete values.
    type (FourierDerivative2d) :: Dobj; !Interpolation object

    !Construct the object by assuming
    !that the fourier interpolant is y.

    allocate(Dobj%y(Nx,Ny),Dobj%yHat(Nx,Ny) &
         ,Dobj%dydxHat(Nx,Ny),Dobj%dydyHat(Nx,Ny));

    Dobj%y = fxy
    Dobj%m_transformed = .false.
    Dobj%m_derivtransformed = .false.
  end function DerFourierInterp2d_
  !---------------------------------------------------------

  !>@brief Basic constructor for the interpolation object
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
!-----------------------------------------------------------

  !>@brief Performs the 2d transform for the object
  !and store it in yHat
  subroutine FwdTransform2d(Iobj)
    implicit none
    class (FourierInterp2d) :: Iobj
    call fft(Iobj%yHat,Iobj%y,shape(Iobj%y));
    Iobj%m_transformed = .true.
  end subroutine FwdTransform2d
!------------------------------------------------------------

  !>@brief Initialize the coefficients for computing
  !the derivative along the x direction:
  subroutine derFwdTransform2d_x(Dobj)
    use fftpack, only : fftfreq !Use the frequencies
    implicit none
    integer :: i,nx
    class (FourierDerivative2d),intent(inout) :: Dobj
    if (Dobj%m_transformed .eqv. .false.) then
       call FwdTransform2d(Dobj)
    endif
    
    nx = size(Dobj%y,1)
    do i=1,size(Dobj%y,2)
       Dobj%dydxHat(:,i) = ju*fftfreq(nx)*Dobj%yHat(:,i)
    enddo
      
  end subroutine derFwdTransform2d_x
  !----------------------------------------------------

  subroutine derFwdTransform2d_y(Dobj)
    use fftpack, only : fftfreq !Use the frequencies
    implicit none
    integer :: i,ny
    class (FourierDerivative2d),intent(inout) :: Dobj
    if (Dobj%m_transformed .eqv. .false.) then
       call FwdTransform2d(Dobj)
    endif
    
    ny = size(Dobj%y,2)
    do i=1,size(Dobj%y,1)
       Dobj%dydyHat(i,:) = ju*fftfreq(ny)*Dobj%yHat(i,:)
    enddo
      
  end subroutine derFwdTransform2d_y
  
  !--------------------------------------------------------
  !>@brief Get the derivative in a specified 2d point
  !along the x direction
  !@param[in] Dobj The derivative object
  !@param[in] xy Specific point where we want the derivative
  !@param[out] dyInt The derivative value
  subroutine FourierGradX_scalar(Dobj,xy,dyInt)
    use fftpack, only : fftfreq !Use the frequencies
    implicit none
    class (FourierDerivative2d),intent(inout) :: Dobj
    real(rk),intent(in) :: xy(:)
    real(rk),intent(out) :: dyInt

    integer :: i,n(2)
    complex(rk),allocatable :: expK(:,:),expL(:,:);

    n = shape(Dobj%y)
     
    allocate(expK(n(1),n(2)),expL(n(1),n(2))); !allocate
    
    !Get the discrete exponential in the specific points:
    do i =1,n(2)
       expK(:,i) = exp(ju*xy(1)*fftfreq(n(1)))
    enddo
    do i=1,n(1)
       expL(i,:) = exp(ju*xy(2)*fftfreq(n(2)))
    enddo
    dyInt = real(sum(Dobj%dydxHat*expK*expL))/(n(1)*n(2));
  end subroutine FourierGradX_scalar
!-------------------------------------------------------------------
  !>@brief Get the derivative in a specified 2d point
  !along the y direction
  !@param[in] Dobj the derivative object
  !@param[in] xy the specific point where we want the derivative
  subroutine FourierGradY_scalar(Dobj,xy,dyInt)
    use fftpack, only : fftfreq !Use the frequencies
    implicit none
    class (FourierDerivative2d),intent(inout) :: Dobj
    real(rk),intent(in) :: xy(:)
    real(rk),intent(out) :: dyInt
    
    integer :: i,n(2)
    complex(rk),allocatable :: expK(:,:),expL(:,:);
    
    n = shape(Dobj%y)

    allocate(expK(n(1),n(2)),expL(n(1),n(2))); !allocate the exponentials

    do i =1,n(2) !x
       expK(:,i) = exp(ju*xy(1)*fftfreq(n(1)))
    enddo
    do i=1,n(1) !y
       expL(i,:) = exp(ju*xy(2)*fftfreq(n(2)))
    enddo
    !Derivative with respect y:
    dyInt = real(sum(Dobj%dydyHat*expK*expL))/(n(1)*n(2));
 !---------------------------------------------------------------------   
  end subroutine FourierGradY_scalar
  
  !>@brief Interpolate in a specified point
  !@param [in] xy: An array which contains
  ! the coordinates where we want to interpolate
  !@param [in] 
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
  !----------------------------------------------------------------

  

  
  !>@brief Interpolate to a set of interpolation points
  !@param[in] Iobj
  !@param[in] xy(nPoints,2) a row vector which contains
  !the coordinates
  !@out param yInt [out] The interpolated values.
  subroutine FourierInterpolate2D_vector(Iobj,xy,yInt)
    implicit none;
    class (FourierInterp2d),intent(in) :: Iobj;
    real(rk),intent(in) :: xy(:,:); !A row vector 
    real(rk),intent(out) ::yInt(:);
    integer :: nPoints,ipnts;

    nPoints = size(xy,1) !Get the numOfPoints
    
    do ipnts = 1,nPoints
       call FourierInterpolate2d_scalar(Iobj,xy(ipnts,:),yInt(ipnts))
    enddo
  end subroutine FourierInterpolate2D_vector
  
end module FourierInterp2dClass
