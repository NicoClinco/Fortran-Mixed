!>@brief Testing the main program for AD of the numeric-
!!fluxes subroutines
PROGRAM driver
  USE SIM_INFO
  IMPLICIT NONE
  external rusanovFlux

  REAL(KIND=WP) :: leftState(nEQ)
  REAL(KIND=WP) :: rightState(nEQ)
  
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

  write(*,*) Fn(1)%d(1:10)

  
END PROGRAM driver
