program main
  use fftpack, only: rk,fft, ifft,fftfreq;
  !use FourierInterp1dClass, only: FourierInterp1d,FourierInterp1d_,FwdTransform,Fourierinterpolate1D_scalar,&
  !     FourierDerivative1d,DerFourierInterp1d_,derFwdTransform,FourierDerivative1D_scalar, InterpolateFourier
  use FourierInterp1dClass
  use fftpack_kind, only: rk
  implicit none;
  integer ::i;
  real(rk),allocatable :: x(:);
  complex(rk), allocatable :: y(:);
  real(rk), allocatable :: y_r(:);
  complex(rk),allocatable :: yHat(:);
  type (FourierInterp1d) :: fourierInterpolant;
  type (FourierDerivative1d) :: fourierDerivative;

  
  !real(rk),parameter :: pi =  3.14159265358979323846;
  real(rk):: y_interp,dy_interp;

  real :: h;
  integer :: n;

  real(rk),allocatable :: y_(:);
  real(rk),allocatable :: y_hat(:);

  n=10;
  allocate(x(n),y(n),yHat(n),y_(n),y_hat(n));

  do i=1,n
     x(i) = (2*pi*(i-1))/real(n);
  enddo

  y = sin(x);

  
  
  fourierInterpolant = FourierInterp1d_(n,real(y));
  fourierDerivative = DerFourierInterp1d_(n,real(y));
  
  call FwdTransform(fourierInterpolant);
  call Fourierinterpolate1D_scalar(fourierInterpolant,0.5_rk,y_interp);
  
  print*,'Interpolated value in x=0.5: ',y_interp;

  call fourierDerivative%derFwdTransform();
  call fourierDerivative%FourierDerivative1D_scalar(0.5_rk,dy_interp);

  call Interpolate(fourierInterpolant,x,y_);
  
  print*,'Interpolation error',abs(y-y_);
  
  print*,'Interpolated value [derivative] in x=0.5: ',dy_interp;
  
  !print '("{",ES10.3,",",1X,ES10.3,"}")',yHat(1:5);


  
 
  


  
  

end program main
