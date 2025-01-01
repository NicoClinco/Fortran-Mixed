SUBROUTINE calc_curve
  USE VAR
  !$openad INDEPENDENT(x) 
  Y = QuadParams(1)*x**2+ &
       QuadParams(1)*x  + &
       QuadParams(3)

  Z = QuadParams(1)*Y
  !$openad DEPENDENT(Z)
  
END SUBROUTINE calc_curve
