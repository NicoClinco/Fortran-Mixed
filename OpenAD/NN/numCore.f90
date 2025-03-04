!>@brief The module that contains various parameters to
!!optimize
!!
MODULE var

  INTEGER,PARAMETER :: n_in=2
  INTEGER,PARAMETER :: n_out=2


  REAL(KIND=8) :: W(n_out,n_in)
  REAL(KIND=8) :: b(n_out)
  !The cost-function (that will be transformed)
  REAL(KIND=8) :: cost_fun

CONTAINS
  !>@brief The subroutine that serves for initialize the weights
  !!
  !!@param[in] nin
  !!@param[in] nout
   SUBROUTINE InitWeights(nin,nout)
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: nin,nout

  !   n_in=nin;n_out=nout
  !   ALLOCATE(W(nout,nin),b(nout) )
     W(1:n_out,1:n_in) = 0.0_8;
     b(1:n_out) = 0.0_8
   END SUBROUTINE InitWeights
END MODULE var
!>@brief Perform the forward in our case and evaluate the
!!loss function according the specific flag passed.
!!
SUBROUTINE forward(nb,nIN,nOut,X,Y,&
     eval_loss,ytrgt)
  USE VAR, ONLY : n_In,n_Out,W,b,cost_fun
  IMPLICIT NONE
  INTEGER,INTENT(IN)         :: nb,nIn,nOut
  REAL(KIND=8),INTENT(IN)    :: X(nIn,nb)
  REAL(KIND=8),INTENT(INOUT) :: Y(nOut,nb)
  LOGICAL,INTENT(IN)         :: eval_loss
  REAL(KIND=8),optional      :: ytrgt(nout,nb)
  INTEGER :: i,j,ib
  !$openad INDEPENDENT(W)
  DO ib=1,nb
     DO j=1,nIn
        DO i=1,nOut
           y(i,ib) = y(i,ib)+W(i,j)*x(j,ib)+b(i)
        ENDDO
     ENDDO
  ENDDO
  if(eval_loss .and. present(ytrgt) ) then
     CALL LossFun(nout,nb,y,ytrgt,cost_fun)
  endif
  !$openad DEPENDENT(cost_fun)
END SUBROUTINE forward


SUBROUTINE LossFun(nout,nb,Y,ytrgt,lf)
  USE var, only : W,b
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nout,nb
  REAL(KIND=8),DIMENSION(nout,nb) :: Y
  REAL(KIND=8),DIMENSION(nout,nb) :: ytrgt
  INTEGER :: ib,i,j
  REAL(KIND=8) :: rnb

  !The value of the loss-function:
  REAL(KIND=8) :: lf
  lf = 0.0_8
  rnb = 1.0_8/REAL(nb,8)
  
  DO ib=1,nb
     DO j=1,nout
        lf = lf + (Y(j,ib)-ytrgt(j,ib))*(Y(j,ib)-ytrgt(j,ib))
     ENDDO
  ENDDO
  lf = lf*(rnb)
END SUBROUTINE LossFun
