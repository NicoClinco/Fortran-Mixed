module Interp1dClass
  use fftpack, only : fft,ifft,fftfreq;
  implicit none;
  
  type Interp1d
     private !Make the 
     integer,parameter :: dimension = 1;
     real,allocatable :: y(:);
     complex,allocatable ::yHat(:);
     logical :: m_transformed;
  end type Interp1d

  real,parameter :: pi =  3.14159265358;!pi
  complex,parameter :: ju = (0.0,1.0);   !img unity
  !interface interpolate1D
  !   
  !end interface interpolate1D
  contains:

  !>@brief Constructor for 1d interpolation
  ! @param [in] The number of discrete values
  ! @param [in] The discrete values
  ! @return The Interpolation-object.
  function Interp1d_(N,y) result(Iobj)
    integer,intent(in) :: N; !Number of discrete values.
    real,intent(in) :: y(N); !Vector of discrete values.
    type (Interp1d) :: Iobj; !Interpolated object

    !Performs the check on N.
    !Allocate the various entries.
    allocate(Iobj%y(N));
    Iobj%y = y;
    allocate(Iobj%yHat(N));
    Iobj%m_transformed = .false.;
  end function Interp1d_
  
  !>@brief Perform the transformation
  !in the fourier space and obtain yHat.
  subroutine Fwdtransform(Iobj)
    type (Interp1d) :: Iobj;
    Iobj%yHat = fft(y);
    Iobj%m_transformed = .true.;
  end subroutine Fwdtransform
  
  !>@brief Interpolate to a specific point
  ! @param [in]  x The point where we interpolate
  ! @param [out] y the returned value
  subroutine interpolate1D_scalar(Iobj,x,yInt)
    type (Interp1d),intent(in) :: Iobj;
    real,intent(in) :: x;
    integer, dimension(size(Iobj%y)) :: freqs;
    complex,dimension(size(Iobj%y)) :: discrete_exp;
    integer :: n;
    real,intent(out) :: yInt; !Interpolated value.
    
    if( Iobj%m_transformed == .false.) then
       print*,'Error, you must initialize the interpolation object before proceeding';
    endif

    n = size(Iobj%y);
    discrete_exp = exp(-ju*fftfreq(n)*2*pi/n*x);
    
    
    yInt = real(sum(discrete_exp*Iobj%yHat)); 
       
  end subroutine interpolate1D_scalar

  !subroutine interpolate1D_vector(Iobj,x,yInt)

  
  
  !end subroutine interpolate1D_vector

end module Interp1dClass
