!$openad XXX Template ad_template.f
SUBROUTINE vect_LossFunction(nbatches,x,y_target,theta,yV)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nbatches
  REAL(KIND=8),DIMENSION(nbatches)     :: x
  REAL(KIND=8),DIMENSION(nbatches)     :: y_target
  REAL(KIND=8),DIMENSION(2)            :: theta
  REAL(KIND=8),DIMENSION(2*nbatches)   :: yV
  !REAL(KIND=8),DIMENSION(nbatches)     :: Lfun
  !REAL(KIND=8)                     :: gLfun
  INTEGER      :: i
  REAL(KIND=8) :: rdummy
  
  !$openad INDEPENDENT(theta)
  CALL fvalues(nbatches,x,theta,yV)
  !$openad DEPENDENT(yV)
  
END SUBROUTINE Vect_LossFunction

!$openad XXX Template ad_template.f
SUBROUTINE fvalues(nbatches,x,theta,fval)
  IMPLICIT NONE
  INTEGER,INTENT(IN)               :: nbatches
  REAL(KIND=8),DIMENSION(nbatches) :: x
  REAL(KIND=8),DIMENSION(2)        :: theta
  REAL(KIND=8),DIMENSION(2*nbatches) :: fval

  INTEGER :: i
  !Linear extrapolation:
  do i=1,nbatches
     fval(i) = theta(1)*x(i)
  enddo
  
  do i=nbatches+1,nbatches*2
     fval(i) = 2*theta(2)*x(i-nbatches)
  enddo
  
END SUBROUTINE fvalues
