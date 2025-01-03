  !>@brief Testing the usage of the neural network
  !!below
PROGRAM testNN
  USE w2f__types
  USE OAD_active
  USE OAD_rev
  USE OAD_tape
  USE var !ALL
  IMPLICIT NONE

  INTEGER(w2f__i4),PARAMETER :: Nbatches  = 10
  INTEGER(w2f__i4),PARAMETER :: nFeatures = 3
  INTEGER(w2f__i4),PARAMETER :: nLayers   = 2
  
  REAL(KIND=8) :: xInput(nFeatures,nBatches)
  REAL(KIND=8) :: Y_TRGT(nFeatures,nBatches)
  TYPE(ACTIVE) :: Y(nFeatures,nBatches)
  INTEGER(w2f__i4)     :: i,j,ilayer,wbindxs(2),ib
  REAL(KIND=8) :: tv(nFeatures*nFeatures)    !True values of the derivatives
  INTEGER      :: gIndxs(nFeatures,nFeatures)
  LOGICAL      :: GET_LOSS_FUNCTION

  !Output of the first layer:
  REAL(KIND=8) :: dummy_out(nFeatures,nbatches)
  REAL(KIND=8) :: dummy_w(nFeatures) !W11^2,W21^2,W31^2
  
  GET_LOSS_FUNCTION = .true.
  

  
  !Initialize the array of batches
  CALL InitInputs(nFeatures,Nbatches,xInput)
  Y_TRGT(:,:) = 0.0_8

  CALL tape_init()
  !Set to zero the derivatives:
  CALL OAD_revTape()
  CALL InitLinearLayers(nFeatures,nFeatures,nLayers)
  !CALL InitTRGT(nFeatures,nBatches,Y_TRGT)
  ! !Set to zero the derivatives:
  ! CALL OAD_revTape()
  ! CALL InitWeights(nFeatures,nFeatures)
  Y(:,:)%v   = 0.0_8
  COST_FUN%d = 1.0_8 
  
  CALL forward(nBatches,xInput,Y,GET_LOSS_FUNCTION,Y_TRGT)
  CALL OAD_revAdjoint()
  CALL forward(nBatches,xInput,Y,GET_LOSS_FUNCTION,Y_TRGT)


  CALL OAD_revPlain()
  
  
  WRITE(*,*)'Cost function value=',COST_FUN%v
  DO j=1,nBatches
     write(*,"(20(A))")'*'
     write(*,"(3(F10.5))")y(:,j)%v
  ENDDO
  
  DO ilayer=1,nlayers
     DO j=1,n_in
        DO i=1,n_out
           CALL getLayerIndexes(n_in,n_out,n_layers,i,j,ilayer,wbindxs)
           write(*,*) '*',i,j,wbindxs(1:2)
           write(*,'(A,i2,A1,i2,A1,F10.5,"|",F10.5)')'dL/dW(',i,',',j,')= ',Wb_global(wbindxs(1))%d,&
                Wb_global(wbindxs(2))%d
           !Output of the first layer:
           
        ENDDO
     ENDDO
  ENDDO

  !True value:
  CALL getLayerIndexes(n_in,n_out,n_layers,1,1,2,wbindxs)
  dummy_w(1) = Wb_global(wbindxs(1))%v
  CALL getLayerIndexes(n_in,n_out,n_layers,2,1,2,wbindxs)
  dummy_w(2) = Wb_global(wbindxs(1))%v
  CALL getLayerIndexes(n_in,n_out,n_layers,3,1,2,wbindxs)
  dummy_w(3) = Wb_global(wbindxs(1))%v
  
  DO ib=1,nbatches
     DO i=1,nFeatures
        dummy_out(i,ib) = dummy_w(i)*xInput(1,ib)
     ENDDO
  ENDDO

  

  

  !Testing the true derivative for the first weigth:
  !tv(1) = (2.0)/REAL(nBatches,8)*DOT_PRODUCT(y(:,:)%v,dummy_out)
  tv(2) = (2.0)/REAL(nBatches,8)*SUM(y(2,:)%v)
  write(*,*)'true, dL/dw11,dL/db',tv(1),tv(2)
  ! DO j=1,nFeatures
  !    DO i=1,nFeatures
  !       write(*,'(A,i1,A1,i1,A1,F10.5)')'dL/dW(',i,',',j,')= ',W(i,j)%d
  !       write(*,'(A,i1,A1,i1,A1,F10.5)')'dL/dWt(',i,',',j,')=',tv(gIndxs(i,j))
  !    ENDDO
  ! ENDDO
  
CONTAINS
  !>@brief A simple initialization of the inputs
  !!@param[in]    nf
  !!@param[in]    nb
  !!@param[inout] x
  SUBROUTINE InitInputs(nf,nb,x)
    IMPLICIT NONE
    INTEGER(kind=4),INTENT(IN)  :: nf,nb
    REAL(KIND=8),INTENT(INOUT)  :: x(nf,nb)
    INTEGER                     :: i,ib
    
    DO ib=1,nb
       DO i=1,nf
          CALL random_number(x(i,ib)) 
       ENDDO
    ENDDO

  END SUBROUTINE InitInputs



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
  
END PROGRAM testNN
