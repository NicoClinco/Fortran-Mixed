SUBROUTINE calc_curve
  USE VAR
  
  Y = QuadParams(1)*x**2+ &
       QuadParams(1)*x  + &
       QuadParams(3)
  
  Z = QuadParams(1)*Y  

END SUBROUTINE calc_curve
