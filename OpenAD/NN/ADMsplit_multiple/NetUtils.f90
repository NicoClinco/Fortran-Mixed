!>@brief Update weights and biases:
!!       Copy the weights and biases from the global-storage
!!       to the W,b.
SUBROUTINE UpdateWeightsAndBiases(nin,nout,nlayers,Wflattened,W,b)
  IMPLICIT NONE
  INTEGER                      :: nin,nout,nlayers
  INTEGER                      :: ilayer,ista,iend
  REAL(KIND=8),INTENT(IN)      :: Wflattened(nin*nout*nlayers+nout*nlayers)
  REAL(KIND=8),INTENT(INOUT)   :: b(nout,nlayers)
  REAL(KIND=8),INTENT(INOUT)   :: W(nout,nin,nlayers)
  INTEGER                      :: i,j,gindx


  !Weigths
  DO ilayer=1,nlayers
     !Weights:
     ista = (ilayer-1)*nin*nout
     iend = ilayer*nin*nout
     DO i=1,nout
        DO j=1,nin
           gindx = (i-1)*nin+j
           W(i,j,ilayer) = Wflattened(ista+gindx)
        ENDDO
     ENDDO

     !biases: weigths-offset+biases index
     ista = (nlayers)*nin*nout+(ilayer-1)*nout
     iend = (nlayers)*nin*nout+ilayer*nout
     DO i=1,nout
        b(i,ilayer) = Wflattened(ista+i)
     ENDDO
  ENDDO
END SUBROUTINE UpdateWeightsAndBiases


!>@brief Print the weights of the neural network
!!       to the std output
!!
SUBROUTINE printWeights(nin,nout,nlayers,W,b)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nin,nout,nlayers
  REAL(KIND=8),INTENT(IN)   :: b(nout,nlayers)
  REAL(KIND=8),INTENT(IN)   :: W(nout,nin,nlayers)

  CHARACTER(40) :: chararr
  INTEGER :: i,j,ilayer
  
  chararr='****************************************'
  
  write(*,'(A40)') chararr
  do ilayer=1,nlayers
     write(*,'(A6,1x,i1," ",A30)') 'LAYER=',ilayer,chararr(1:30)
     DO i=1,nout
        write(*,'(3(F12.8,","),"|",F12.8)')W(i,:,ilayer),b(i,ilayer)
     ENDDO
  enddo
  write(*,'(A40)') chararr
END SUBROUTINE printWeights
