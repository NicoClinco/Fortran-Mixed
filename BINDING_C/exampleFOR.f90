program PrintfromC
  use iso_c_binding
  implicit none
  
  interface
     subroutine helloworld() bind(c)

     end subroutine helloworld
  end interface

  interface
     subroutine printvariable(x) bind(c)
       import c_double
       real(kind=c_double),intent(in),value :: x
     end subroutine printvariable
  end interface

  interface
     subroutine printvector(x,size) bind(c)
       import c_double, c_int
       real(kind=c_double),intent(in) :: x(*)
       integer(kind=c_int),intent(in) :: size
     end subroutine printvector
  end interface

  real(kind=8) :: x_value
  real(kind=c_double):: x_vector(5)

  
  x_vector = 0.0
  
  
  call HelloWorld

  x_value = real(0.5,8)
  call printvariable(x_value)

  call printvector(x_vector,5)
  
end program PrintfromC
