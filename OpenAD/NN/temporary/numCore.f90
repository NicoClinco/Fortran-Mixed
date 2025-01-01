MODULE var

  INTEGER,PARAMETER :: nIN  = 10
  INTEGER,PARAMETER :: nOut = 10 
  
  !Parameters below:
  REAL(KIND=8) :: W(nOut,nIN)
  REAL(KIND=8) :: b(nOut)
END MODULE var

SUBROUTINE forward(nb,nInp,nOutp,X,Y)
  USE VAR, ONLY : nIn,nOut,W,b
  IMPLICIT NONE
  INTEGER,INTENT(IN)         :: nb,nInp,nOutp
  REAL(KIND=8),INTENT(IN)    :: X(nInp,nb)
  REAL(KIND=8),INTENT(INOUT) :: Y(nOutp,nb)
  INTEGER :: i,j,ib
  
  !$openad INDEPENDENT(W)
  DO ib=1,nb
     DO j=1,nInp
        DO i=1,nOutp
           y(i,ib) = y(i,ib)+W(i,j)*x(j,ib)+b(i)
        ENDDO
     ENDDO
  ENDDO
  !$openad DEPENDENT(y)
END SUBROUTINE forward
