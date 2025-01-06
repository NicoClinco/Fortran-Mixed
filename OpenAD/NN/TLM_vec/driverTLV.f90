!>@brief Testing the usage of the neural network
!!in tangent vector mode.
PROGRAM testNN
  USE w2f__types
  USE OAD_active
  USE var 
  IMPLICIT NONE

  INTEGER(w2f__i4),PARAMETER :: Nbatches  = 10
  INTEGER(w2f__i4),PARAMETER :: nFeatures = 3
  INTEGER(w2f__i4),PARAMETER :: nLayers   = 2
  
  REAL(KIND=8) :: xInput(nFeatures,nBatches)
  REAL(KIND=8) :: Y_TRGT(nFeatures,nBatches)
  INTEGER(w2f__i4)     :: i,j,ilayer,wbindxs(2),ib
  INTEGER      :: gIndxs(nFeatures,nFeatures)
  LOGICAL      :: GET_LOSS_FUNCTION

  !derivatives
  REAL(KIND=8) :: d_W(nFeatures,nFeatures,nlayers),d_b(nFeatures,nlayers)
  
  !True derivatives
  REAL(KIND=8) :: d_Wt(nFeatures,nFeatures,nlayers),d_bt(nFeatures,nlayers)
  
  external forward
  
  GET_LOSS_FUNCTION = .true.
  
  !Initialize the array of batches
  CALL InitInputs(nFeatures,Nbatches,xInput)
  ALLOCATE(Y(nFeatures,Nbatches))
  Y_TRGT(:,:) = 0.0_8

  CALL InitLinearLayers(nFeatures,nFeatures,nLayers)
  
  
  !CALL InitTRGT(nFeatures,nBatches,Y_TRGT)
  Y_TRGT(:,:) = 0.0_8


  DO i=1,size(Wb_global%v)
     CALL zero_deriv(Wb_global(i))
     Wb_global(i)%d(i) = 1.0D0
  ENDDO
  CALL forward(nBatches,xInput,GET_LOSS_FUNCTION,Y_TRGT)

  WRITE(*,*)'Cost function value=',COST_FUN%v
  DO j=1,nBatches
     write(*,"(20(A))")'*'
     write(*,"(3(F10.5))")y(:,j)%v
  ENDDO
  
  DO ilayer=1,nlayers
     DO j=1,n_in
        DO i=1,n_out
           CALL getLayerIndexes(n_in,n_out,n_layers,i,j,ilayer,wbindxs)
           CALL UpdateWeightsAndBiases(nfeatures,nfeatures,nlayers,&
                cost_fun%d,d_W,d_b)
        ENDDO
     ENDDO
  ENDDO
  
  write(*,'(A)')'=== TLVM VALUES:'
  
  CALL printWeights(nFeatures,nFeatures,nlayers,d_W,d_b)
  write(*,'(A)')'========================================'
  write(*,'(A)')'=== TRUE VALUES:'
  CALL UpdateWeightsAndBiases(nfeatures,nfeatures,nlayers,Wb_global%v,W_layers,b_layers)
  CALL dLdW1(nfeatures,nfeatures,nBatches,W_layers(:,:,2),xInput,y%v,d_Wt(:,:,1))
  CALL dLdW2(nfeatures,nfeatures,nBatches,W_layers(:,:,1),xInput,y%v,b_layers(:,1),d_Wt(:,:,2))
  CALL dLdb1(nfeatures,nfeatures,nBatches,W_layers(:,:,2),y%v,d_bt(:,1))
  CALL dLdb2(nfeatures,nfeatures,nBatches,y%v,d_bt(:,2))
  !True:
  ilayer=1
  do i=1,n_out
     d_bt(i,2)   = (2.0)/REAL(nBatches,8)*sum(y(i,:)%v)
  enddo
  CALL printWeights(nFeatures,nFeatures,nlayers,d_Wt,d_bt)
  
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

  !>@brief Get the derivative of the loss function with respect
  !!       the weights of the first layer.
  !!
  !!@note
  !! Here, the equation is the following:
  !!  dL/dWij = W^{2}*[\delta_{ij}]x]^T\cdot y
  !!@endnote
  SUBROUTINE dLdW1(nin,nout,nb,W2,x,y,dW)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nin,nout,nb
    REAL(KIND=8) :: W2(nout,nout)
    REAL(KIND=8) :: x(nin,nb)
    REAL(KIND=8) :: y(nout,nb)
    REAL(KIND=8) :: dW(nout,nout)

    REAL(KIND=8) :: delta_mat(nout,nin)
    REAL(KIND=8) :: dmx(nout)
    REAL(KIND=8) :: W2dmx(nout)
    INTEGER      :: i,j,ib


    delta_mat(:,:) = 0.0_8
    DO i=1,nout
       DO j=1,nout
          delta_mat(i,j) = 1.0_8
          DO ib=1,nb
             dmx = MATMUL(delta_mat,x(:,ib))
             W2dmx = MATMUL(W2,dmx)
             dW(i,j) = dW(i,j)+DOT_PRODUCT(W2dmx,y(:,ib))
          ENDDO
          delta_mat(i,j) = 0.0_8
          dW(i,j) = 2.0_8*dW(i,j)/REAL(nb,8)
       ENDDO
    ENDDO
  END SUBROUTINE dLdW1

  !>@brief Get the derivative of the loss function with respect the weights
  !!       of the second layer
  !!
  !!
  SUBROUTINE dLdW2(nin,nout,nb,W1,x,y,b1,dW)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nin,nout,nb
    REAL(KIND=8) :: W1(nout,nin)
    REAL(KIND=8) :: b1(nout)
    REAL(KIND=8) :: x(nin,nb)
    REAL(KIND=8) :: y(nout,nb)
    REAL(KIND=8) :: dW(nout,nout)

    REAL(KIND=8) :: delta_mat(nout,nout)
    REAL(KIND=8) :: dmW1(nout,nin)
    REAL(KIND=8) :: dmW1x(nout)
    REAL(KIND=8) :: dmb1(nout)
    REAL(KIND=8) :: dw1xPLUSdb1(nout)
    INTEGER      :: i,j,ib


    delta_mat(:,:) = 0.0_8
    DO i=1,nout
       DO j=1,nout
          delta_mat(i,j) = 1.0_8
          DO ib=1,nb
             dmW1 = MATMUL(delta_mat,W1)
             dmW1x = MATMUL(dmW1,xInput(:,ib))
             dmb1 = MATMUL(delta_mat,b1)
             dw1xPLUSdb1 = dmW1x + dmb1
             dW(i,j) = dW(i,j)+DOT_PRODUCT(dw1xPLUSdb1,y(:,ib))
          ENDDO
          delta_mat(i,j) = 0.0_8
          dW(i,j) = 2.0_8*dW(i,j)/REAL(nb,8)
       ENDDO
    ENDDO
  END SUBROUTINE dLdW2

  SUBROUTINE dLdb1(nin,nout,nb,W2,y,db)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nin,nout,nb
    REAL(KIND=8) :: W2(nout,nout)
    REAL(KIND=8) :: b1(nout)
    REAL(KIND=8) :: y(nout,nb)
    REAL(KIND=8) :: db(nout)

    REAL(KIND=8) :: delta_vec(nout)
    REAL(KIND=8) :: W2delta(nout)
    INTEGER      :: i,ib


    delta_vec(:) = 0.0_8
    DO i=1,nout
       delta_vec(i) = 1.0_8
       W2delta = MATMUL(W2,delta_vec)
       DO ib=1,nb
          db(i) = db(i)+DOT_PRODUCT(W2delta,y(:,ib))
       ENDDO
       db(i) = 2.0_8*db(i)/REAL(nb,8)
       delta_vec(i) = 0.0_8
    ENDDO
  END SUBROUTINE dLdb1

  SUBROUTINE dLdb2(nin,nout,nb,y,db)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nin,nout,nb

    REAL(KIND=8) :: y(nout,nb)
    REAL(KIND=8) :: db(nout)
    REAL(KIND=8) :: delta_vec(nout)
    
    INTEGER      :: i,ib

    
    delta_vec(:) = 0.0_8
    DO i=1,nout
       delta_vec(i) = 1.0_8
       DO ib=1,nb
          db(i) = db(i)+DOT_PRODUCT(delta_vec,y(:,ib))
       ENDDO
       db(i) = 2.0_8*db(i)/REAL(nb,8)
       delta_vec(i) = 0.0_8
    ENDDO
    
  END SUBROUTINE dLdb2

END PROGRAM testNN
