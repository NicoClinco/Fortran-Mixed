  !>@brief The driver for testing our framework (ADJOINT MODE : SPLIT)
!!
PROGRAM main

  use OAD_active
  use OAD_rev
  !use OAD_tape
  
  USE var
  IMPLICIT NONE

  x%v = 0.5_8
  CALL zero_deriv(x)
  CALL zero_deriv(Y)
  QuadParams(1:3) = 1.0_8
  !Set the derivative of the output to zero?
  Z%d = 1.0_8
  CALL OAD_revTape()
  CALL calc_curve()
  CALL OAD_revAdjoint()
  CALL calc_curve()
  
  WRITE(*,*) 'Derivative of Z wrt the x values:',x%d
  WRITE(*,*) 'Derivative of Y wrt the x values:',Y%d
  
END PROGRAM main
