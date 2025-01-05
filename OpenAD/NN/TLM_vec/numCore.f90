!>@brief The module that contains various parameters to
!!optimize (allocatable version)
!!
!!The module uses a flattened version for the weights and biases
!!and stores multiple layers
MODULE var

  !Features:
  INTEGER                  :: n_in,n_out,n_layers
  REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: Wb_global

  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: Y
  

  ![nout,nin,nlayers]
  REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: W_layers
  ![nout,nlayers]
  REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: b_layers

  !The cost-function (that will be transformed)
  REAL(KIND=8) :: cost_fun

CONTAINS
  !>@brief The subroutine that serves for initialize the weights
  !!
  !!@param[in] nin
  !!@param[in] nout
  SUBROUTINE InitLinearLayers(nin,nout,nlayers)
    IMPLICIT NONE
    INTEGER,INTENT(IN) :: nin,nout,nlayers
    REAL(KIND=8)       :: rnd
    INTEGER :: i,j
    n_in=nin;n_out=nout;n_layers=nlayers
    ALLOCATE(Wb_global(nout*nin*n_layers+nout*n_layers))
    ALLOCATE(W_layers(n_out,n_in,n_layers),b_layers(n_out,n_layers))
    DO i=1,nout*nin*n_layers+nout*n_layers
       CALL random_number(rnd)
       Wb_global(i) = rnd
    ENDDO
    !ALLOCATE(Y(nout,10))
  END SUBROUTINE InitLinearLayers


  !>@brief Get the global index for our weights and biases
  !!from the global flattened array
  !!
  !!
  SUBROUTINE getLayerIndexes(nin,nout,nlayers,i,j,ilayer,gindxs)
    IMPLICIT NONE
    INTEGER :: nin,nout,nlayers,i,j,ilayer
    INTEGER,INTENT(INOUT) :: gindxs(2)
    gindxs(1) = (ilayer-1)*nin*nout + (i-1)*nin+j
    gindxs(2) = (nlayers)*nin*nout+(ilayer-1)*nout + i
  END SUBROUTINE getLayerIndexes

END MODULE var

SUBROUTINE forward(nb,x,eval_loss,ytrgt)
  USE var, only : wb_global,cost_fun,n_in,n_out,n_layers,&
       W_layers,b_layers,getLayerIndexes,y
  IMPLICIT NONE


  INTEGER                :: nb
  REAL(KIND=8)           :: x(n_in,nb)
  REAL(KIND=8)           :: yi,yj
  LOGICAL,INTENT(IN)     :: eval_loss
  REAL(KIND=8),optional  :: ytrgt(n_out,nb)
  REAL(KIND=8)           :: zBuf(n_out,nb) 
  INTEGER                :: ib,i,j,ilayer,ista,i1,i2
  INTEGER                :: wbindxs(2)

  !$openad INDEPENDENT(Wb_global)
  
  !=========================================================
  ! -Perform the forward method:
  ! -After each layer a buffer zBuff is stored
  !=========================================================

  zBuf(1:n_out,1:nb) = 0.0_8
  
  !To store an intermediate vector.
  DO ib=1,nb
     DO ilayer=1,n_layers
        if(ilayer.eq.1) then
           
           DO i=1,n_out
              
              DO j=1,n_in
                 !CALL getLayerIndexes(n_in,n_out,n_layers,i,j,ilayer,wbindxs)
                 i1 = (ilayer-1)*n_in*n_out + (i-1)*n_in+j
                 y(i,ib) = y(i,ib)+Wb_global(i1)*x(j,ib)
                 
              ENDDO
              i2 = (n_layers)*n_in*n_out+(ilayer-1)*n_out + i
              y(i,ib) = y(i,ib) + Wb_global(i2)
           ENDDO
        else
           y(1:n_out,ib) = 0.0_8
           DO i=1,n_out
              DO j=1,n_in
                 !CALL getLayerIndexes(n_in,n_out,n_layers,i,j,ilayer,wbindxs)
                 i1 = (ilayer-1)*n_in*n_out + (i-1)*n_in+j
                 y(i,ib) = y(i,ib)+Wb_global(i1)*zBuf(j,ib)
              ENDDO
              i2 = (n_layers)*n_in*n_out+(ilayer-1)*n_out + i
              y(i,ib) = y(i,ib) + Wb_global(i2)
           ENDDO
        endif
        zBuf(1:n_out,ib) = y(1:n_out,ib)
        
     ENDDO !END-LAYER-LOOP
  ENDDO
  
  CALL LossFun(nb,y,ytrgt,cost_fun)
  
  !$openad DEPENDENT(cost_fun)
END SUBROUTINE forward
  

!>@brief Computation of the loss function for our neural network
!!
!!@param[in]    nb
!!@param[in]    Y
!!@param[in]    ytrgt
!!@param[inout] cost_fun
SUBROUTINE LossFun(nb,y,ytrgt,lf)
  USE var, only : n_out
  IMPLICIT NONE
  INTEGER                          :: nb
  REAL(KIND=8),DIMENSION(n_out,nb) :: Y
  REAL(KIND=8),DIMENSION(n_out,nb) :: ytrgt
  INTEGER :: ib,i,j
  REAL(KIND=8) :: rnb
  
  !The value of the loss-function:
  REAL(KIND=8) :: lf
  lf = 0.0_8
  rnb = 1.0_8/REAL(nb,8)

  !**openad xxx simple loop
  DO ib=1,nb
     DO j=1,n_out
        lf = lf + (Y(j,ib)-ytrgt(j,ib))*(Y(j,ib)-ytrgt(j,ib))
     ENDDO
  ENDDO
  lf = lf*(rnb)
  
END SUBROUTINE LossFun
