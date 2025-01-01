  !>@brief Testing the usage of the neural network
  !!below
PROGRAM testNN
  USE OAD_active
  USE OAD_rev
  USE var, ONLY : W,b
  IMPLICIT NONE

  INTEGER,PARAMETER :: Nbatches=10
  INTEGER,PARAMETER :: nFeatures=2
  REAL(KIND=8) :: xInput(nFeatures,nBatches)
  TYPE(ACTIVE) ::Y(nFeatures,nBatches)

  TYPE(ACTIVE) :: Wg(nFeatures,nFeatures,nBatches)
  REAL(KIND=8) :: bg(nFeatures,nBatches)

  INTEGER      :: ib,i

  CALL zero_deriv(Wg)
  xInput(:,:) = 0.5_8
  y(:,:)%d      =1.0_8

  CALL InitInputs(nFeatures,nBatches,xInput)
  
  CALL OAD_revTape()
  CALL forward(Nbatches,nFeatures,nFeatures,xInput,Y,Wg,bg)
  CALL OAD_revAdjoint()
  CALL forward(Nbatches,nFeatures,nFeatures,xInput,Y,Wg,bg)

  DO ib=1,nbatches
     write(*,"(20(A))")'*'
     DO i=1,nFeatures 
        write(*,"(10(F8.5,' '))") Wg(i,:,ib)%d
     ENDDO
     write(*,"(20(A))")'*'
  ENDDO
  
CONTAINS
  SUBROUTINE InitInputs(nf,nb,x)
    IMPLICIT NONE
    INTEGER,INTENT(IN)         :: nf,nb
    REAL(KIND=8),INTENT(INOUT) :: x(:,:)
    INTEGER                    :: i,ib

    DO ib=1,nb
       DO i=1,nf
          CALL random_number(x(i,ib)) 
       ENDDO
    ENDDO

  END SUBROUTINE InitInputs
  
END PROGRAM testNN
