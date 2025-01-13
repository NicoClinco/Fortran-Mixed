  !>@brief Testing multiple automatic differentiation
  !!modules here
PROGRAM driver
  USE OAD_active
  USE ad_vars1, y1 => y, x1 => x
  USE ad_vars2, y2 => y, x2 => x

  IMPLICIT NONE
  external init_differentiation
  

  x1%v = 0.5_8
  x2%v = 0.4_8

  x1%d = 1.0_8
  x2%d = 1.0_8
  CALL init_ad1()
  CALL init_ad2()
  !CALL init_differentation()
  
  
  write(*,*) 'dy1/dx1:',y1%d
  write(*,*) 'dy2/dx2:',y2%d
  
END PROGRAM driver
