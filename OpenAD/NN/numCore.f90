MODULE var
  
  !Parameters below
  REAL(KIND=8),ALLOCATABLE :: W(:,:)
  REAL(KIND=8),ALLOCATABLE :: b(:)
END MODULE var

SUBROUTINE forward(nb,nInp,nOutp,X,Y,Wg,bg)
  USE VAR, ONLY : W,b
  IMPLICIT NONE
  INTEGER,INTENT(IN)         :: nb,nInp,nOutp
  REAL(KIND=8),INTENT(IN)    :: X(nInp,nb)
  REAL(KIND=8),INTENT(INOUT) :: Y(nOutp,nb)
  
  
  !Parameters below:
  REAL(KIND=8)               :: Wg(nOutp,nInp,nb)
  REAL(KIND=8)               :: bg(nOutp,nb)
  INTEGER :: i,j,ib
  
  !$openad INDEPENDENT(Wg)
  DO ib=1,nb
     DO i=1,nInp
        DO j=1,nOutp
           Y(j,ib) = Y(j,ib)+Wg(j,i,ib)*x(i,ib)+bg(j,ib)
        ENDDO
     ENDDO
  ENDDO
  !$openad DEPENDENT(y)
END SUBROUTINE forward
