program PrintfromC
  USE, INTRINSIC :: ISO_C_BINDING
  implicit none

  interface
     subroutine helloworld() bind(c)
     end subroutine helloworld
  end interface
!-----------------------------------------------------------
  interface
     subroutine printvariable(x) bind(c)
       import c_double
       real(kind=c_double),intent(in),value :: x
     end subroutine printvariable
  end interface
!-----------------------------------------------------------
  interface
     subroutine printvector(x,size) bind(c)
       USE, INTRINSIC :: ISO_C_BINDING
       integer(kind=c_int),intent(in),value:: size
       real(kind=c_double),intent(in) :: x(*)

     end subroutine printvector
  end interface

  real(kind=8) :: x_value
  real(kind=8):: x_vector(5)

  x_vector = 0.0


  call HelloWorld

  x_vector(1) = 35.5
  call printvector(x_vector,5)

end program PrintfromC
