!>@brief Testing the main program for AD of the numeric-
!!fluxes subroutines
PROGRAM driver
  USE SIM_INFO
  IMPLICIT NONE
  external rusanovFlux

  REAL(KIND=WP) :: leftState(nEQ)
  REAL(KIND=WP) :: rightState(nEQ)
  REAL(KIND=WP) :: FD_DER(5,2,nEQ)
  
  TYPE(ACTIVE) :: Fn(nEQ)
  REAL(KIND=WP) :: Sf(3)
  INTEGER       :: runIndxs(3),iF

  run1D = .true.
  
  leftState(1) = 1.0_WP
  leftState(2) = 1.0_WP
  leftState(5) = 2.0_WP
  
  rightState(1) = 1.0_WP
  rightState(2) = 0.0_WP
  rightState(5) = 1.5_WP
  
  Sf(1) = 1.0_WP

  CALL rusanovFlux(leftState,rightState,Fn,Sf)

  CALL FD_calculation(leftState,rightState,FD_DER)

  runIndxs = (/1,2,5/)
  !Printing the left-state:
  DO iF=1,5
     IF((iF .ne.3) .and.(iF .ne.4).and. (iF .ne.3)) then
        write(*,"(A4,3(F11.7,','))") 'FD: ',FD_DER(iF,1,runIndxs)
        write(*,"(A4,3(F11.7,','))") 'AD: ',Fn(iF)%d(runIndxs)
     ENDIF
  ENDDO
  


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
    REAL(KIND=WP) :: rs_m1(nEQ)
    TYPE(ACTIVE) :: fn(nEQ),fn_p(nEQ)
    REAL(KIND=WP) :: Sf(3)
    ![OUT,{LEFT/RIGHT},NEQ]
    REAL(KIND=WP),INTENT(OUT) :: dFn_dVar(5,2,nEQ)
    REAL(KIND=WP) :: dVar
    INTEGER       :: StateIndex

    dVar = 0.001_WP
    !Computation of dFn/drhoL:
    ls(1:nEQ) = left(1:nEQ)
    rs(1:nEQ) = right(1:nEQ)
    ls_p1(:) = ls(:)
    ls_m1(:) = ls(:)
    Sf(1) = 1.0_WP

    ! left state ================================================
    StateIndex=1
    !dFn(1)/drhoL
    ls_p1(1) = ls_p1(1)+dVar
    ls_m1(1) = ls_m1(1)-dVar
    CALL rusanovFlux(ls_m1,rs,Fn,Sf)
    CALL rusanovFlux(ls_p1,rs,Fn_p,Sf)
    dFn_dVar(1,stateIndex,1) = (Fn_p(1)%v - Fn(1)%v)/(2.0_WP*dVar)
    dFn_dVar(2,stateIndex,1) = (Fn_p(2)%v - Fn(2)%v)/(2.0_WP*dVar)
    dFn_dVar(5,stateIndex,1) = (Fn_p(5)%v - Fn(5)%v)/(2.0_WP*dVar)
    !dFn(1)/drhouL
    ls_m1 = ls(:)
    ls_p1 = ls(:)
    ls_p1(2) = ls_p1(2)+dVar
    ls_m1(2) = ls_m1(2)-dVar
    CALL rusanovFlux(ls_m1,rs,Fn,Sf)
    CALL rusanovFlux(ls_p1,rs,Fn_p,Sf)
    dFn_dVar(1,stateIndex,2) = (Fn_p(1)%v - Fn(1)%v)/(2.0_WP*dVar)
    dFn_dVar(2,stateIndex,2) = (Fn_p(2)%v - Fn(2)%v)/(2.0_WP*dVar)
    dFn_dVar(5,stateIndex,2) = (Fn_p(5)%v - Fn(5)%v)/(2.0_WP*dVar)

    !dFn(1)/drhoEl
    ls_m1 = ls(:)
    ls_p1 = ls(:)
    ls_p1(5) = ls_p1(5)+dVar
    ls_m1(5) = ls_m1(5)-dVar    
    CALL rusanovFlux(ls_m1,rs,Fn,Sf)
    CALL rusanovFlux(ls_p1,rs,Fn_p,Sf)
    dFn_dVar(1,stateIndex,5) = (Fn_p(1)%v - Fn(1)%v)/(2.0_WP*dVar)
    dFn_dVar(2,stateIndex,5) = (Fn_p(2)%v - Fn(2)%v)/(2.0_WP*dVar)
    dFn_dVar(5,stateIndex,5) = (Fn_p(5)%v - Fn(5)%v)/(2.0_WP*dVar)

    !right state:
    StateIndex=2
    !dFn(1)/drhoR
    rs_p1    = rs(:)
    rs_m1    = rs(:)
    rs_p1(1) = rs(1)+dVar
    rs_m1(1) = rs(1)-dVar
    CALL rusanovFlux(ls,rs_m1,Fn,Sf)
    CALL rusanovFlux(ls,rs_p1,Fn_p,Sf)
    dFn_dVar(1,stateIndex,1) = (Fn_p(1)%v - Fn(1)%v)/(2.0_WP*dVar)
    dFn_dVar(2,stateIndex,1) = (Fn_p(2)%v - Fn(2)%v)/(2.0_WP*dVar)
    dFn_dVar(5,stateIndex,1) = (Fn_p(5)%v - Fn(5)%v)/(2.0_WP*dVar)
    !dFn(1)/drhouR
    rs_m1 = rs(:)
    rs_p1 = rs(:)
    rs_p1(2) = rs_p1(2)+dVar
    rs_m1(2) = rs_m1(2)-dVar
    CALL rusanovFlux(ls,rs_m1,Fn,Sf)
    CALL rusanovFlux(ls,rs_p1,Fn_p,Sf)
    dFn_dVar(1,stateIndex,2) = (Fn_p(1)%v - Fn(1)%v)/(2.0_WP*dVar)
    dFn_dVar(2,stateIndex,2) = (Fn_p(2)%v - Fn(2)%v)/(2.0_WP*dVar)
    dFn_dVar(5,stateIndex,2) = (Fn_p(5)%v - Fn(5)%v)/(2.0_WP*dVar)

    !dFn(1)/drhoER
    rs_m1 = rs(:)
    rs_p1 = rs(:)
    rs_p1(5) = rs_p1(5)+dVar
    rs_m1(5) = rs_m1(5)-dVar   
    CALL rusanovFlux(ls,rs_m1,Fn,Sf)
    CALL rusanovFlux(ls,rs_p1,Fn_p,Sf)
    dFn_dVar(1,stateIndex,5) = (Fn_p(1)%v - Fn(1)%v)/(2.0_WP*dVar)
    dFn_dVar(2,stateIndex,5) = (Fn_p(2)%v - Fn(2)%v)/(2.0_WP*dVar)
    dFn_dVar(5,stateIndex,5) = (Fn_p(5)%v - Fn(5)%v)/(2.0_WP*dVar)
    
  END SUBROUTINE FD_calculation
  
END PROGRAM driver
