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
  end type FourierInterp1d

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

  !>@brief Perform the transformation
  !in the fourier space and obtain yHat.
  subroutine Fwdtransform(Iobj)
    type (FourierInterp1d) :: Iobj;
    complex(rk) :: tmp(size(Iobj%y));

    tmp = Iobj%y;

    Iobj%yHat = fft(tmp)
    Iobj%m_transformed = .true.;
  end subroutine Fwdtransform

  !>@brief Interpolate to a specific point
  ! @param [in]  x The point where we interpolate
  ! @param [out] y the returned value
  subroutine Fourierinterpolate1D_scalar(Iobj,x,yInt)
    type (FourierInterp1d),intent(in) :: Iobj;
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

  !subroutine interpolate1D_vector(Iobj,x,yInt)
  !end subroutine interpolate1D_vector

end module FourierInterp1dClass
