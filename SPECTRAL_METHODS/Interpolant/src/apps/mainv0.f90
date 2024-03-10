program main
 
  use fftw3_wrapper
  
  implicit none
  complex(kind=8),dimension(10,10) :: input_array,org_array,error_array;
  complex(kind=8),dimension(:,:),allocatable :: output_array;
  real(kind=8) :: x(10,10),y(10,10);
  integer :: i;
  integer :: num(2);
  integer :: kfreqs(10),lfreqs(10);
  real,parameter :: pi = 3.1416
  !Let's suppose to do the fourier transform of the sin

  character(20) :: separator_line;
  separator_line(:) = '********************';
  
  do i=1,10
     x(i,:) = (2*pi*(i-1))/real(10);
     y(:,i) = (2*pi*(i-1))/real(10);
  enddo

  input_array = sin(x)*sin(y); !Create the input array
  org_array = input_array;
  num(:) = 10;
  allocate(output_array(10,10))
  call fft(output_array,input_array,num)
  input_array =output_array; !copy
  call ifft(output_array,input_array,num)
  error_array = output_array - org_array

  print '(2A)',separator_line
  call fftfreq(kfreqs)
  
  print '(2A)',separator_line
  
  do i = 1,10
     print '(10(f6.3,2x))',output_array(i,:)
  enddo
  print '(2A)',separator_line
  do i = 1,10
     print '(10(e12.6,2x))',error_array(i,:)
  enddo
  print '(A,1X,f6.3)','Absolute-error:',sum(real(output_array)-(org_array));

end program main
