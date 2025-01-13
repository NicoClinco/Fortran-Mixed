!> In this example we englobe together all the differentiable
!! stuffs.
!!
MODULE AD_vars1

  REAL(KIND=8) :: x
  REAL(KIND=8) :: y
CONTAINS
  SUBROUTINE Init_ad1
    IMPLICIT NONE
    !$openAD INDEPENDENT(x)
    y = x**2
    !$openAD DEPENDENT(y)    
  END SUBROUTINE Init_ad1
END MODULE AD_vars1


MODULE AD_vars2

  REAL(KIND=8) :: x
  REAL(KIND=8) :: y
CONTAINS
  SUBROUTINE Init_ad2
    IMPLICIT NONE

    !$openAD INDEPENDENT(x)
    y = x**2
    !$openAD DEPENDENT(y)
  END SUBROUTINE Init_ad2
END MODULE AD_vars2

SUBROUTINE Init_differentiation()
  USE AD_vars1, only : Init_ad1
  USE AD_vars2, only : Init_ad2
  IMPLICIT NONE
  CALL init_ad1
  CALL init_ad2
  
END SUBROUTINE Init_differentiation
