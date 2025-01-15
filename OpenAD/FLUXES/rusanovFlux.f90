MODULE SIM_INFO
  
  INTEGER, parameter :: WP = kind(1.d0)
  INTEGER,parameter  :: nEQ = 5
  
  LOGICAL :: run1D = .true.,run2D =.false.
  REAL(KIND=WP) :: gam = 1.4_WP
  REAL(KIND=WP) :: gam_m1 = 0.4_WP
  
  
END MODULE SIM_INFO

!>@brief Computation of the rusanovFlux between left and right
!!state
!!
!!@param[in] Qfl: Left state
!!@param[in] Qfr: Right state
!!@param[in] Fn : Normal flux
!!@param[in] Sf : Surface vectors
!!
SUBROUTINE rusanovFlux(Qfl,Qfr,Fn,Sf)
  USE SIM_INFO
  IMPLICIT NONE


  real(kind=WP)   :: Qfl(nEQ),Qfr(nEQ),Sf(3)
  real(kind=WP)   :: Fn(nEQ)
  real(kind=WP)   :: unL,unR,pL,pR,am,normS,lambda

  REAL(kind=WP) :: Qflr(2*nEQ)

  !$openAD INDEPENDENT(Qflr)
  Qflr(1:nEQ) = Qfl(1:nEQ)
  Qflr(nEQ+1:2*nEQ) = Qfr(1:nEQ)

  CALL pressure(Qflr(1:nEQ),pL)
  CALL pressure(Qflr(nEQ+1:2*nEQ),pR)

  if(pL+pR.lt.0.0_WP) then
     !Rmnn_err = 1
     Fn(:) = 0.0_WP
     !return
  end if

  if(run1D) then !----------------------------------------------
     normS = abs(Sf(1))
     unL = (Qflr(2)*Sf(1))/Qflr(1)
     unR = (Qflr(nEQ+2)*Sf(1))/Qflr(1+nEQ)

     ! Calculate the normal flux component
     Fn(1)   = unL*Qflr(1) + unR*Qflr(1+nEQ)
     Fn(2)   = unL*Qflr(2) + unR*Qflr(2+nEQ) + (pL+pR)*Sf(1)
     Fn(3:4) = 0.0_WP
  else if(run2D) then !-----------------------------------------
     normS = sqrt(Sf(1)**2+Sf(2)**2)
     unL = (Qflr(2)*Sf(1)+Qfl(3)*Sf(2))/Qflr(1)
     unR = (Qflr(2+nEQ)*Sf(1)+Qflr(3+nEQ)*Sf(2))/Qflr(1+nEQ)

     ! Calculate the normal flux component
     Fn(1)   = unL*Qflr(1)   + unR*Qflr(1+nEQ)
     Fn(2:3) = unL*Qfl(2:3) + unR*QfLr(2+nEQ:3+nEQ) + (pL+pR)*Sf(1:2)
     Fn(4)   = 0.0_WP
  else !3D run case --------------------------------------------
     normS = sqrt(Sf(1)**2+Sf(2)**2+Sf(3)**2)
     unL = (Qflr(2)*Sf(1)+Qflr(3)*Sf(2)+Qflr(4)*Sf(3))/Qflr(1)
     unR = (QfLr(2+nEQ)*Sf(1)+Qflr(3+nEQ)*Sf(2)+Qflr(4+nEQ)*Sf(3))/Qflr(1+nEQ)

     ! Calculate the normal flux component
     Fn(1)   = unL*Qflr(1)   + unR*Qflr(1+nEQ)
     Fn(2:4) = unL*Qflr(2:4) + unR*Qflr(2+nEQ:4+nEQ) + (pL+pR)*Sf(1:3)
  end if !------------------------------------------------------
  Fn(5) = unL*(Qflr(5)+pL) + unR*(Qflr(5+nEQ)+pR)

  am  = sqrt(gam*(pL+pR)/(Qflr(1)+Qflr(1+nEQ)))
  lambda = 0.5_WP*abs(unL+unR) + am*normS

  Fn(1:nEQ) = 0.5_WP*(Fn(1:nEQ) - lambda*(Qflr(nEQ+1:2*nEQ)-Qflr(1:nEQ)))

  !$openAD DEPENDENT(Fn)
END SUBROUTINE rusanovFlux

  
SUBROUTINE pressure(Qv,p)
  USE SIM_INFO
  implicit none
  REAL(kind=WP)             :: Qv(nEQ)
  REAL(kind=WP)             :: p
  REAL(kind=WP)             :: rhousq

  if(run1D) then !----------------------------------------------
     rhousq = (Qv(2)**2)/Qv(1)
  else if(run2D) then !-----------------------------------------
     rhousq = (Qv(2)**2+Qv(3)**2)/Qv(1)
  else !3D run case --------------------------------------------
     rhousq = (Qv(2)**2+Qv(3)**2+Qv(4)**2)/Qv(1)
  end if !------------------------------------------------------
  p = gam_m1*(Qv(5)-0.5_WP*rhousq)  
END SUBROUTINE pressure
