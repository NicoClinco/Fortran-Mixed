  !> Program that computes the derivatives of subroutines in forward vector mode.
  !!
PROGRAM driver

  USE OAD_active
  USE var
  USE AD_vars
  IMPLICIT NONE

  
  INTEGER :: gisp2ic(Nvar*Ns*3),gifp2ic(Nvar*(Ns+1)*3)
  INTEGER :: i,j,ic_sp,ic_fp

  !Solution points in the stencil:
  REAL(KIND=8) :: qs_stncl(Nvar,Ns,3)
  !Flattened
  REAL(KIND=8) :: fqs_stncl(Nvar*Ns*3)


  character(len=30)          :: frmt

  !Initialize:
  Qs(:,:,:) = 1.0_8
  CALL sliceQs(2,qs_stncl)
  fqs_stncl = RESHAPE(qs_stncl,(/Nvar*Ns*3/))
  
  Lmat(:,:) = 0.0_8
  Mmat(:,:) = 0.0_8

  !Create Xs,Xf:
  CALL make_dpoints
  
  !Create the interpolation matrixes (Lmat)
  CALL create_interpolation_matrixes
  frmt = '(1p,1X,   E12.4)'
  write(frmt(8:10),'(I3)') Ns
  do i=1,Ns+1
     write(*,frmt) Lmat(i,1:Ns)
  end do

  
  CALL gsp2_cell(gisp2ic)
  CALL gfp2_cell(gifp2ic)
  
  DO i=1,size(fqs_stncl)
     CALL zero_deriv(flatstencilQS(i))
     flatstencilQS(i)%d(i) = 1.0_8
  ENDDO

  CALL calc_resid(fqs_stncl)

  !Print the derivative of the flux-points with respect to the
  !sol.points
  DO i=1,size(flatstencilQfp)
     DO j=1,size(flatstencilQs)
        ic_sp = gisp2ic(j)
        ic_fp = gifp2ic(i)
        if(ic_sp .eq. ic_fp .and. ic_fp .eq. 1) then
           write(*,*) ic_sp,' ',ic_fp,' ',flatstencilQfp(i)%d(j),flat_fluxes(i)%d(j)
        endif
     ENDDO
  ENDDO

  

CONTAINS

  SUBROUTINE gsp2_cell(g2ic_sp)
    USE var
    USE AD_vars
    IMPLICIT NONE
    INTEGER :: g2ic_sp(Nvar*Ns*3)
    INTEGER :: ic,is,ivar,igsp
    
    DO ic=1,3
       DO is=1,Ns
          DO ivar=1,Nvar
             CALL GLOBALINDX_SP(ivar,is,ic,igsp)
             g2ic_sp(igsp) = ic
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE gsp2_cell

  SUBROUTINE gfp2_cell(g2ic_fp)
    USE var
    USE AD_vars
    IMPLICIT NONE
    
    INTEGER :: g2ic_fp(Nvar*(Ns+1)*3)
    INTEGER :: ic,is,ivar,igfp
    
    DO ic=1,3
       DO is=1,Ns+1
          DO ivar=1,Nvar
             CALL GLOBALINDX_FP(ivar,is,ic,igfp)
             g2ic_fp(igfp) = ic
          ENDDO
       ENDDO
    ENDDO
    
  END SUBROUTINE gfp2_cell

    
  
END PROGRAM driver
