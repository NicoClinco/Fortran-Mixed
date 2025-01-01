!>@brief The subroutine for the forward method
!!
!!@param[in] n
!!@param[in] X
!!@param[in] Y
SUBROUTINE forward(n,X,Y)
  USE VAR, ONLY : W,b
  IMPLICIT NONE
  INTEGER,INTENT(IN)         :: n
  REAL(KIND=8),INTENT(IN)    :: X(n)
  REAL(KIND=8),INTENT(INOUT) :: Y(n)
  Y = W@x + b

END SUBROUTINE forward
