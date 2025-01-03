!>@brief The module that contains various parameters to
!!optimize (allocatable version)
!!
MODULE var

  !Features:
  INTEGER                  :: n_in,n_out

  REAL(KIND=8),ALLOCATABLE :: W(:,:)
  REAL(KIND=8),ALLOCATABLE :: b(:)
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
    REAL(KIND=8)       :: rnd
    INTEGER :: i,j
       n_in=nin;n_out=nout
       ALLOCATE(W(nout,nin),b(nout) )
    DO i=1,nout
       DO j=1,nin
          CALL random_number(rnd)
          W(j,i) = rnd
       ENDDO
       CALL random_number(rnd)
       b(i) = rnd
    ENDDO
  END SUBROUTINE InitWeights
END MODULE var
!>@brief Perform the forward in our case and evaluate the
!!loss function according the specific flag passed.
!!
SUBROUTINE forward(nb,nIn,nOut,X,Y,&
     eval_loss,ytrgt)
  USE VAR, ONLY : n_In,n_Out,W,b,cost_fun
  IMPLICIT NONE
  INTEGER                    :: nb,nIn,nOut
  REAL(KIND=8),INTENT(IN)    :: X(nIn,nb)
  REAL(KIND=8)               :: Y(nOut,nb)
  LOGICAL,INTENT(IN)         :: eval_loss
  LOGICAL                    :: is_present_trgt
  REAL(KIND=8),optional      :: ytrgt(nOut,nb)
  INTEGER :: i,j,ib
  
  !$openad INDEPENDENT(W)
  
  !$openad xxx simple loop
  DO ib=1,nb
     DO j=1,nIn
        DO i=1,nOut
           y(i,ib) = y(i,ib)+W(i,j)*x(j,ib)+b(i)
        ENDDO
     ENDDO
  ENDDO
  
  CALL LossFun(nout,nb,y,ytrgt,cost_fun)
  !$openad DEPENDENT(cost_fun)
END SUBROUTINE forward


SUBROUTINE LossFun(nout_,nb_,Y,ytrgt,lf)
  IMPLICIT NONE
  INTEGER                           :: nout_,nb_
  REAL(KIND=8),DIMENSION(nout_,nb_) :: Y
  REAL(KIND=8),DIMENSION(nout_,nb_) :: ytrgt
  INTEGER :: ib,i,j
  REAL(KIND=8) :: rnb

  !The value of the loss-function:
  REAL(KIND=8) :: lf
  lf = 0.0_8
  rnb = 1.0_8/REAL(nb_,8)

  !$openad xxx simple loop
  DO ib=1,nb_
     DO j=1,nout_
        lf = lf + (Y(j,ib)-ytrgt(j,ib))*(Y(j,ib)-ytrgt(j,ib))
     ENDDO
  ENDDO
  lf = lf*(rnb)
END SUBROUTINE LossFun

