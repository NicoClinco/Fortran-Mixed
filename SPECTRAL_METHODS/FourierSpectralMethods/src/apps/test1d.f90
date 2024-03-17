program main
  use fftpack, only: rk,fft, ifft,fftfreq;
  !use FourierInterp1dClass, only: FourierInterp1d,FourierInterp1d_,FwdTransform,Fourierinterpolate1D_scalar,&
  !     FourierDerivative1d,DerFourierInterp1d_,derFwdTransform,FourierDerivative1D_scalar, InterpolateFourier
  use FourierInterp1dClass
  use fftpack_kind, only: rk

  USE SCIFOR, only: parse_cmd_variable,parse_input_variable 
  
  implicit none;
  integer ::i;
  real(rk),allocatable :: x(:);
  complex(rk), allocatable :: y(:);
  real(rk), allocatable :: y_r(:);
  complex(rk),allocatable :: yHat(:);
  type (FourierInterp1d) :: fourierInterpolant;
  type (FourierDerivative1d) :: fourierDerivative;

  !Derivative matrix:
  real(rk),allocatable:: Dn(:,:),In(:,:);
  real(rk),allocatable :: Dny(:),Iny(:);
  !real(rk),parameter :: pi =  3.14159265358979323846;
  real(rk):: y_interp,dy_interp;

  real :: h;
  integer :: n,npnts;

  real(rk),allocatable :: y_(:);
  real(rk),allocatable :: y_hat(:);

  character(len=28) :: cmdlinestring
  
  !---------------> Classical method : explicit evaluation:
  call parse_cmd_variable(cmdlinestring,"cmdlinestring")
  call parse_input_variable(npnts,"num_points",cmdlinestring,default=5,comment='Number of points to test')
  write(*,'("Selected number of points: ",i2)')npnts
  !npnts=9;
  n = npnts
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
  !------------------------------------->
  !print '("{",ES10.3,",",1X,ES10.3,"}")',yHat(1:5);

  !npnts =4
  allocate(Dn(npnts,npnts))
  allocate(In(npnts,npnts))
  Dn = GetDerivativeMatrix(npnts)
  In = GetInterpolationMatrix(npnts)
  
  allocate(Dny(npnts),Iny(npnts))
  
  Dny = matmul(Dn,y)
  Iny = matmul(In,y)

  write(*,'("Interpolation with matrixes",*(E12.4))') Iny
  write(*,'("Exact solution y = sin(x)",*(E12.4))') real(y)
  
  write(*,'("Derivative with matrixes",*(E12.4))') Dny
  write(*,'("Exact solution y = cos(x)",*(E12.4))') cos(x)

  
  


  
  

end program main
