!>@brief Testing the main program for AD of the numeric-
!!fluxes subroutines
PROGRAM driver
  USE SIM_INFO
  IMPLICIT NONE
  external rusanovFlux

  REAL(KIND=WP) :: leftState(nEQ)
  REAL(KIND=WP) :: rightState(nEQ)
  REAL(KIND=WP) :: FD_DER(2,nEQ)
  
  TYPE(ACTIVE) :: Fn(nEQ)
  REAL(KIND=WP) :: Sf(3)

  run1D = .true.
  
  leftState(1) = 1.0_WP
  leftState(2) = 1.0_WP
  leftState(5) = 2.0_WP
  
  rightState(1) = 1.0_WP
  rightState(2) = 0.0_WP
  rightState(5) = 1.5_WP
  
  Sf(1) = 1.0_WP

  CALL rusanovFlux(leftState,rightState,Fn,Sf)

  CALL FD_calculation(leftState,rightState,FD_DER(1,1))
  write(*,'("AD: dF/drhoL: ",F13.10)')Fn(1)%d(1)
  write(*,'("FD: dF/drhoL: ",F13.10)')FD_DER(1,1)
  


CONTAINS
  !>@brief Compute the derivative of the flux with the finite
  !!differences. This serves for comparison
  !!
  SUBROUTINE FD_calculation(left,right,dFn_dVar)
    USE OAD_active
    USE SIM_INFO
    IMPLICIT NONE

    REAL(KIND=WP) :: left(nEQ),right(nEQ)
    REAL(KIND=WP) :: ls(nEQ)
    REAL(KIND=WP) :: ls_p1(nEQ)
    REAL(KIND=WP) :: ls_m1(nEQ)
    REAL(KIND=WP) :: rs(nEQ)
    REAL(KIND=WP) :: rs_p1(nEQ)
    TYPE(ACTIVE) :: fn(nEQ),fn_p(nEQ)
    REAL(KIND=WP) :: Sf(3)
    REAL(KIND=WP),INTENT(OUT) :: dFn_dVar
    REAL(KIND=WP) :: dVar

    dVar = 0.00001_WP
    !Computation of dFn/drhoL:
    ls(1:nEQ) = left(1:nEQ)
    rs(1:nEQ) = right(1:nEQ)
    ls_p1(:) = ls(:)
    ls_m1(:) = ls(:)
    ls_p1(1) = ls_p1(1)+dVar
    ls_m1(1) = ls_m1(1)-dVar
    Sf(1) = 1.0_WP
    
    
    CALL rusanovFlux(ls_m1,rs,Fn,Sf)
    CALL rusanovFlux(ls_p1,rs,Fn_p,Sf)
    dFn_dVar = (Fn_p(1)%v - Fn(1)%v)/(2.0_WP*dVar)
    
    
  END SUBROUTINE FD_calculation
  
END PROGRAM driver
