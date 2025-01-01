MODULE var

  REAL(KIND=8) :: x
  REAL(KIND=8) :: Y
  REAL(KIND=8) :: Z
  !Quadratic parameters below:
  REAL(KIND=8) :: QuadParams(3)
  
END MODULE var

SUBROUTINE calc_curve
  
  USE var
  IMPLICIT NONE
 !$openad INDEPENDENT(x) 
  Y = QuadParams(1)*x**2+ &
       QuadParams(2)*x  + &
       QuadParams(3)
  !$openad DEPENDENT(Y)
  Z = QuadParams(1)*Y 
  
  !$openad DEPENDENT(Z)
  
  
END SUBROUTINE calc_curve
