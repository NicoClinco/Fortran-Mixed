module FourierInterp1dClass
  use fftpack, only : fftfreq,fft,ifft;
  use fftpack_kind, only: rk
  implicit none;

  type FourierInterp1d
     !private 
     integer :: m_dim = 1;
     real(rk),allocatable :: y(:);
     complex(rk),allocatable ::yHat(:);
     logical :: m_transformed;
   contains
     
     !procedure :: Interpolate
     procedure :: Fwdtransform
     procedure :: Fourierinterpolate1D_scalar 
  end type FourierInterp1d

  interface Interpolate
     module procedure Fourierinterpolate1D_scalar,Fourierinterpolate1D_vector
  end interface Interpolate
  
  !>@brief Derived class
  type, extends(FourierInterp1d):: FourierDerivative1d
     !private
     complex(rk),allocatable :: dyHat(:);
     logical :: m_derivtransformed;
   contains
     procedure :: derFwdtransform
     procedure :: FourierDerivative1D_scalar
  end type FourierDerivative1d

  real(rk),parameter :: pi =  3.14159265358979323846;!pi
  complex,parameter :: ju = (0.0,1.0);   !img unity

  !interface interpolate1D

  !end interface interpolate1D
contains

  !>@brief Constructor for 1d interpolation
  ! @param [in] The number of discrete values
  ! @param [in] The discrete values
  ! @return The Interpolation-object.
  function FourierInterp1d_(N,y) result(Iobj)
    integer,intent(in) :: N; !Number of discrete values.
    real(rk),intent(in) :: y(N); !Vector of discrete values.
    type (FourierInterp1d) :: Iobj; !Interpolated object

    !Performs the check on N.
    !Allocate the various entries.
    allocate(Iobj%y(N));
    Iobj%y = y;
    allocate(Iobj%yHat(N));
    Iobj%m_transformed = .false.;
  end function FourierInterp1d_
!---------------------------------------------------
  !>@brief Constructor for the 1d derivative.
  function DerFourierInterp1d_(N,y) result(Dobj)
    integer,intent(in) :: N; !Number of discrete values.
    real(rk),intent(in) :: y(N); !Vector of discrete values.
    type (FourierDerivative1d) :: Dobj;

    !Explicit cast:
    allocate(Dobj%y(N));
    Dobj%y = y;
    !Allocate the complex coefficients:
    allocate(Dobj%yHat(N));
    allocate(Dobj%dyHat(N));
    Dobj%m_transformed = .false.;
    Dobj%m_derivtransformed = .false.;

  end function DerFourierInterp1d_
  
  !>@brief Perform the transformation
  !in the fourier space and obtain yHat.
  subroutine Fwdtransform(Iobj)
    class (FourierInterp1d) :: Iobj;
    complex(rk) :: tmp(size(Iobj%y));

    tmp = Iobj%y;

    Iobj%yHat = fft(tmp)
    Iobj%m_transformed = .true.; !Allows us to save the state.
  end subroutine Fwdtransform
!--------------------------------------------------
  subroutine derFwdtransform(Dobj)
    class (FourierDerivative1d) :: Dobj;
    complex(rk) :: tmp(size(Dobj%y));
    
    if(Dobj%m_transformed .eqv. .false.) then
       !Dobj%yHat = fft((Dobj%y));
       tmp = Dobj%y;
       Dobj%yHat = fft(tmp);
       Dobj%m_transformed = .true.;
    endif
    
    Dobj%dyHat = ju * (Dobj%yHat) * fftfreq(size(Dobj%yHat));
    Dobj%m_derivtransformed = .true.;
  end subroutine derFwdtransform
!--------------------------------------------------
  !>@brief Interpolate to a specific point
  ! @param [in]  x The point where we interpolate
  ! @param [out] y the returned value
  subroutine Fourierinterpolate1D_scalar(Iobj,x,yInt)
    class (FourierInterp1d),intent(in) :: Iobj;
    real(rk),intent(in) :: x;

    complex(rk),dimension(size(Iobj%y)) :: discrete_exp;
    integer :: n;
    real(rk),intent(out) :: yInt; !Interpolated value.

    if( Iobj%m_transformed .eqv. .false.) then
       print*,'Error, you must initialize the interpolation object before proceeding';
    endif

    n = size(Iobj%y);

    discrete_exp = exp(ju*fftfreq(n)*x);
    yInt = abs(sum(discrete_exp*Iobj%yHat))/real(n); !Normalize

  end subroutine Fourierinterpolate1D_scalar

  !>@brief Interpolate to a set of points
  !specified in the vector 'x'
  subroutine Fourierinterpolate1D_vector(Iobj,x,yInt)
    implicit none;
    class (FourierInterp1d),intent(in) :: Iobj;
    real(rk),intent(in) :: x(:);     !Vector of points where we want to interpolate.
    real(rk),intent(out) :: yInt(:); !Vector of interpolation points
    complex(rk),dimension(size(x),size(Iobj%y)) :: discrete_exp; !The matrix needed for interpolation
    integer :: i;
    if( Iobj%m_transformed .eqv. .false.) then
       print*,'Error, you must initialize the interpolation object before proceeding';
    endif

    do i=1,size(x)
       discrete_exp(i,:) = exp(ju*fftfreq(size(Iobj%y))*x(i));
    enddo

    
    !Matrix multiplication:
    yInt =  real(matmul(discrete_exp,Iobj%yHat))/(size(Iobj%y));
    
  end subroutine Fourierinterpolate1D_vector

  !>@brief Return the derivative of the function in the
  !point x
  subroutine FourierDerivative1D_scalar(Dobj,x,dyInt)
    class (FourierDerivative1d),intent(in) :: Dobj;
    real(rk),intent(in) :: x;
    real(rk),intent(out) :: dyInt; !The interpolated value
    complex(rk),dimension(size(Dobj%y)) :: discrete_exp;
    
    if((Dobj%m_transformed .eqv. .false.) .or. &
         (Dobj%m_derivtransformed .eqv. .false.)) then
       print*,'Error, you must initialize the interpolation object before proceeding';
    endif

    !n = size(Dobj%y);
    discrete_exp = exp(ju*fftfreq(size(Dobj%y))*x);

    dyInt = abs(sum(discrete_exp*Dobj%dyHat))/real(size(Dobj%y));
    
  end subroutine FourierDerivative1D_scalar
  
end module FourierInterp1dClass
