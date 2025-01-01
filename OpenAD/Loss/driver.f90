program driver
  USE OAD_active
  IMPLICIT NONE
  external vect_LossFunction
  !Number of batches:
  INTEGER, PARAMETER :: N = 5
  type(active) :: theta(2),Yvalues(2*N)
  REAL(kind=8) :: x(N),y_target(N) 
  INTEGER      :: i
  
  x(1:N) = 2.0
  theta(1)%v     = 0.1_8
  theta(2)%v     = 0.05_8
  theta%d        = 1.0D0
  y_target = 0.0_8
  CALL vect_LossFunction(N,x,y_target,theta,yvalues)
  write(*,*) 'LOSS-FUNCTION VALUE:',yvalues(1:N)%v
  write(*,*) 'LOSS-FUNCTION VALUE:',yvalues(N+1:2*N)%v

  write(*,*) 'DERIVATIVE-VALUE-1:',yvalues(1:N)%d
  write(*,*) 'DERIVATIVE-VALUE-2:',yvalues(N+1:)%d
  write(*,*) '** Loss function value **'
  !write(*,*) LossValue%v
  
  write(*,*) '** Derivative of the loss function wrt the parameters **'
  !write(*,*) LossValue%d
  write(*,*)'*************************************'
CONTAINS
  SUBROUTINE calc_der_loss(N,y,ytrgt,dLdtheta)
    USE OAD_active
    IMPLICIT NONE
    INTEGER      :: N
    type(active) :: y(2*N)
    REAL(KIND=8) :: ytrgt(N)
    REAL(KIND=8) :: dydtheta(N,2) !Jacobian
    REAL(KIND=8) :: dLdtheta(2)
    INTEGER      :: i
    REAL(KIND=8) :: rBatches

    rBatches = REAL(2.0,8)/N
    dLdtheta(:) = 0.0_8
    DO i=1,N
       dydtheta(i,1) = y(i)%d
       dydtheta(i,2) = y(N+i)%d
    ENDDO
    
    do i=1,N
       dLdtheta(1) = dLdtheta(1)+rBatches*(y(i)%v-ytrgt(i))*dydtheta(i,1)
       dLdtheta(2) = dLdtheta(2)+rBatches*(y(i+N)%v-ytrgt(i))*dydtheta(i,2)
    enddo

    
    
  END SUBROUTINE calc_der_loss

  
end program driver
