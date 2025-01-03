  !>@brief Testing the usage of the neural network
  !!below
PROGRAM testNN
  USE w2f__types
  USE OAD_active
  USE OAD_rev
  USE OAD_tape
  USE var, ONLY : W,b,InitWeights,COST_FUN
  IMPLICIT NONE

  INTEGER(w2f__i4),PARAMETER :: Nbatches=10
  INTEGER(w2f__i4),PARAMETER :: nFeatures=3
  REAL(KIND=8) :: xInput(nFeatures,nBatches)
  REAL(KIND=8) :: Y_TRGT(nFeatures,nBatches)
  TYPE(ACTIVE) :: Y(nFeatures,nBatches)
  INTEGER      :: i,j
  REAL(KIND=8) :: tv(nFeatures*nFeatures) !True values of the derivatives
  INTEGER      :: gIndxs(nFeatures,nFeatures)
  LOGICAL      :: GET_LOSS_FUNCTION 

  GET_LOSS_FUNCTION = .true.

  
  !Initialize the array of batches
  CALL InitInputs(nFeatures,nBatches,xInput)
  Y_TRGT(:,:) = 0.0_8
  !CALL InitTRGT(nFeatures,nBatches,Y_TRGT)

  CALL tape_init()
  !Set to zero the derivatives:
  CALL OAD_revTape()
  CALL InitWeights(nFeatures,nFeatures)
  Y(:,:)%v   = 0.0_8
  COST_FUN%d = 1.0_8 
  
  CALL forward(nBatches,nFeatures,nFeatures,xInput,Y,GET_LOSS_FUNCTION,Y_TRGT)
  CALL OAD_revAdjoint()
  CALL forward(nBatches,nFeatures,nFeatures,xInput,Y,GET_LOSS_FUNCTION,Y_TRGT)

  
  WRITE(*,*)'Cost function value=',COST_FUN%v
  DO j=1,nBatches
     write(*,"(20(A))")'*'
     write(*,"(3(F10.5))")y(:,j)%v
     write(*,"(20(A))")'*'
  ENDDO

  !Verification of the true derivative of the loss wrt the parameters:
  DO j=1,nFeatures
     DO i=1,nFeatures
        gIndxs(i,j) = (i-1)*nFeatures+j
        tv(gIndxs(i,j)) = (2.0)/REAL(nBatches,8)*SUM(y(i,:)%v*xInput(j,:))
     ENDDO
  ENDDO
  
  DO j=1,nFeatures
     DO i=1,nFeatures
        write(*,'(A,i1,A1,i1,A1,F10.5)')'dL/dW(',i,',',j,')= ',W(i,j)%d
        write(*,'(A,i1,A1,i1,A1,F10.5)')'dL/dWt(',i,',',j,')=',tv(gIndxs(i,j))
     ENDDO
  ENDDO
  
CONTAINS
!===========================================================
  !>@brief A simple initialization of the weights and biases
  !!@param[in]    nf
  !!@param[in]    nb
  !!@param[inout] x
  SUBROUTINE InitInputs(nf,nb,x)
    IMPLICIT NONE
    INTEGER(kind=4),INTENT(IN)          :: nf,nb
    REAL(KIND=8),INTENT(INOUT)  :: x(nf,nb)
    INTEGER                     :: i,ib
    
    DO ib=1,nb
       DO i=1,nf
          CALL random_number(x(i,ib)) 
       ENDDO
    ENDDO

  END SUBROUTINE InitInputs
!===========================================================
  !>@brief A simple initialization of the target values
  !!
  !!@param[in] nf    : Number of features
  !!@param[in] nb    : Number of batches
  !!@param[in] ytrgt : 
  SUBROUTINE InitTRGT(nf,nb,ytrgt)
    IMPLICIT NONE
    INTEGER(kind=4),INTENT(IN)         :: nf,nb
    REAL(KIND=8),INTENT(INOUT) :: ytrgt(nf,nb)
    INTEGER :: i,ib
    
    DO ib=1,nb
       DO i=1,nf
          CALL random_number(ytrgt(i,ib))
          ytrgt(i,ib) = 5.0+ytrgt(i,ib)
       ENDDO
    ENDDO
  END SUBROUTINE InitTRGT
!======================================================  
END PROGRAM testNN
