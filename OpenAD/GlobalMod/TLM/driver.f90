  !>@brief The driver for testing our framework
!!
PROGRAM main
  USE var
  IMPLICIT NONE

  x%v = 0.5_8
  x%d = 1.0_8
  QuadParams(1:3) = 1.0_8

  CALL calc_curve()
  
  WRITE(*,*) 'Derivative of Z wrt the x values:',Z%d
  WRITE(*,*) 'Derivative of Y wrt the x values:',Y%d
  
END PROGRAM main
