program main
 
  use fftw3_wrapper
  use fftpack, only : fftfreq

  use FourierInterp2dClass
  
  implicit none
  complex(kind=8),dimension(:,:),allocatable:: input_array,org_array,error_array,dfdx;
  
  complex(kind=8),dimension(:,:),allocatable :: output_array,dYhat,yHat,dYtdx;
  real(kind=8),allocatable :: x(:,:),y(:,:);
  integer :: i,n;
  integer :: num(2);
  integer,allocatable :: kfreqs(:,:),lfreqs(:,:);
  !real,parameter :: pi = 3.1416

  type (FourierInterp2d) :: Interpolant2d;
  real(kind=8) :: y_int;
  character(20) :: separator_line;
  separator_line(:) = '********************';


  

  n = 11; !random
  allocate( &
       x(n,n),y(n,n),input_array(n,n),&
       org_array(n,n),error_array(n,n),dfdx(n,n),&
       output_array(n,n),dYhat(n,n),yHat(n,n),dYtdx(n,n),&
       kfreqs(n,n),lfreqs(n,n) &
  )
  do i=1,n
     x(i,:) = (2*pi*(i-1))/real(n);
     y(:,i) = (2*pi*(i-1))/real(n);
  enddo

  input_array = sin(x)*sin(y); !Create the input array
  org_array = input_array;
  num(:) = n;
  call fft(output_array,input_array,num)
  input_array =output_array; !copy
  call ifft(output_array,input_array,num)
  error_array = output_array - org_array

  !do i = 1,n
  !   print '(10(f6.3,2x))',output_array(i,:)
  !enddo
  print '(2A)',separator_line,separator_line
  !do i = 1,n
  !   print '(10(e12.6,2x))',error_array(i,:)
  !enddo
  print '(A,1X,e12.6)','Absolute-error-interpolation:',sum(real(output_array)-(org_array));


  print '(2A)',separator_line,separator_line

  
  
  write(*,'(A40)') '...Filling the frequencies vectors...'

  do i =1,n
     kfreqs(i,:) = fftfreq(n)
     lfreqs(:,i) = fftfreq(n)
  enddo

  

  print '(A,A,A)',separator_line,'Computing the spectral derivative',separator_line

  !Original derivative:
  dfdx = cos(x)*sin(y)
  
  call fft(yHat,org_array)

  !Multiply:
  dyHat = yHat*kfreqs

  call ifft(dYtdx,dyHat)

  print '(A,1X,e12.6)','Absolute-error-derivative:',sum(real(dYtdx)-(dfdx));

  print '(2A)',separator_line,separator_line;

  print '(A)','...Creating the interpolant...'
  
  Interpolant2d = FourierInterp2d_(n,n,real(org_array,8))

  call FwdTransform2d(Interpolant2d)

  call FourierInterpolate2D_scalar(Interpolant2d,[0.2_rk,0.2_rk],y_int)

  
  print '(A,1x,e12.6)','Interpolated value at x,y = (0.2,0.2):',y_int
  
end program main
