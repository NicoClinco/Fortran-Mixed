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

  !Derivative matrix:
  real(rk),allocatable:: Dn(:,:);
  
  !real(rk),parameter :: pi =  3.14159265358979323846;
  real(rk):: y_interp,dy_interp;

  real :: h;
  integer :: n,npnts;

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

  ! Note that here we have exploited the polymorphism:
  ! The functions/subroutines take a base-class object,
  ! but we can use a derived object:
  call Interpolate(fourierDerivative,x,y_);
  call FwdTransform(fourierDerivative);
  
  print '(A,1X,E11.5)','Interpolation error',abs(sum(y_-y));
  
  print*,'Interpolated value [derivative] in x=0.5: ',dy_interp;
  
  !print '("{",ES10.3,",",1X,ES10.3,"}")',yHat(1:5);

  npnts = 5
  allocate(Dn(0:npnts,0:npnts))
  Dn = GetDerivativeMatrix(npnts)

  open(999, file='DerMatrix.csv', status='replace', action='write')
  do i=1,npnts
     write(999,'(*(E12.6,","))') Dn(i,:)
  enddo
  close(999)
  
  
 
  


  
  

end program main
