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
     ista = 1+(ilayer-1)*nin*nout
     iend = ilayer*nin*nout
     DO i=1,nout
        DO j=1,nin
           gindx = (i-1)*nin+j
           W(i,j,ilayer) = Wflattened(ista+gindx)
        ENDDO
     ENDDO

     !biases: weigths-offset+biases index
     ista = (nlayers)*nin*nout+1+(ilayer-1)*nout
     iend = (nlayers)*nin*nout+ilayer*nout
     DO i=1,nout
        b(i,ilayer) = Wflattened(ista+i)
     ENDDO
  ENDDO
END SUBROUTINE UpdateWeightsAndBiases
