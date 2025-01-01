SUBROUTINE forward(nb,nIN,nOut,X,Y)
  USE VAR, ONLY : nIn,nOut,W,b
  IMPLICIT NONE
  INTEGER,INTENT(IN)         :: nb,n
  REAL(KIND=8),INTENT(IN)    :: X(nIn,nb)
  REAL(KIND=8),INTENT(INOUT) :: Y(nOut,nb)
  INTEGER :: i,j,ib
  
  !$openad INDEPENDENT(W)
  DO ib=1,nb
     DO j=1,nIn
        DO i=1,nOut
           y(i,ib) = y(i,ib)+W(i,j)*x(j,ib)+b(i)
        ENDDO
     ENDDO
  ENDDO
  !$openad DEPENDENT(y)
END SUBROUTINE forward
