!>@brief The program for testing the 3d interpolation
PROGRAM driver

  USE OAD_active
  USE var
  USE AD_vars, ONLY : flatstencilQs,flatstencilQfp,globalIndx_fp
  IMPLICIT NONE

  REAL(KIND=8) :: Qs_stncl(Nvar,Ns,Ns,Ns)
  INTEGER      :: i,is,js,ks,ivar,gIndx
  INTEGER      :: gindxs2_ijkc_sp(Nvar*Ns*Ns*Ns,4)
  INTEGER      :: gindxs2_ijkc_fp(Nvar*(Ns+1)*Ns*Ns,4)  
  character(len=30)          :: frmt
  REAL(KIND=8)               :: eps_tol = 1.e-8

  !Data-structure which contains the derivatives:
  TYPE :: AD_Qf
     INTEGER :: num_sp
     INTEGER,ALLOCATABLE      :: SP_ijk(:,:) !(num_sp,3)
     REAL(KIND=8),ALLOCATABLE :: DF_SP(:)    !(num_sp)
  END TYPE AD_Qf

  !In each flux point we have the derivative of the fp wrt the
  !sol.points
  TYPE(AD_Qf) :: AD_Qf1(Nvar,Ns+1,Ns,Ns)
  INTEGER     :: num_sp,j,isp,jsp,ksp,ivarsp
  
  Lmat(:,:) = 0.0_8
  Mmat(:,:) = 0.0_8

  CALL make_dpoints
  CALL create_interpolation_matrixes

  frmt = '(1p,1X,   E12.4)'
  write(frmt(8:10),'(I3)') Ns
  do i=1,Ns+1
     write(*,frmt) Lmat(i,1:Ns)
  end do

  
  !Initialize the values to 1.0 for simplicity
  !Initialize the derivatives for the forward mode:
  DO i=1,size(flatstencilQS)
     CALL zero_deriv(flatstencilQS(i))
     flatstencilQs(i)%d(i) = 1.0_8
     flatstencilQs(i)%v    = 1.0_8
  ENDDO

  !Evaluate the discrete points and the derivative also
  CALL calc_qfp()
  CALL build_maps(gindxs2_ijkc_sp,gindxs2_ijkc_fp)

  !For each fp we count the number of non-zero entries
  !in the derivative vector
  ivar=1
  DO ks=1,Ns
     DO js=1,Ns
        DO is=1,Ns+1
           DO ivar=1,Nvar
              CALL globalIndx_fp(ivar,is,js,ks,gIndx)
              AD_Qf1(ivar,is,js,ks)%num_sp = 0
              DO i=1,Ns**3
                 if(abs(flatstencilQfp(gIndx)%d(i)) .gt. eps_tol) then
                    AD_Qf1(ivar,is,js,ks)%num_sp = AD_Qf1(ivar,is,js,ks)%num_sp+1
                 endif
              ENDDO
              num_sp = AD_Qf1(ivar,is,js,ks)%num_sp
              ALLOCATE(AD_Qf1(ivar,is,js,ks)%DF_SP(num_sp))
              ALLOCATE(AD_Qf1(ivar,is,js,ks)%SP_ijk(num_sp,3))
              j=0
              DO i=1,Ns**3
                 if(abs(flatstencilQfp(gIndx)%d(i)) .gt. eps_tol) then
                    j=j+1
                    AD_Qf1(ivar,is,js,ks)%DF_SP(j)  = flatstencilQfp(gIndx)%d(i)
                    AD_Qf1(ivar,is,js,ks)%SP_ijk(j,1:3) = gindxs2_ijkc_sp(i,2:4)
                    ivarsp = gindxs2_ijkc_sp(i,1)
                    isp = gindxs2_ijkc_sp(i,2);jsp = gindxs2_ijkc_sp(i,3);ksp = gindxs2_ijkc_sp(i,4)
                    
                    write(*,'("dQ_",i1,i1,i1,i1,"/","dx_",i1,i1,i1,i1,"=",F11.7)')ivar,is,js,ks,ivarsp,isp,jsp,ksp,&
                         AD_Qf1(ivar,is,js,ks)%DF_SP(j)
                 endif
              END DO
           ENDDO !var
        ENDDO !is
     ENDDO !js
  ENDDO !ks
  

CONTAINS

  !>@brief From global to local indexes below
  !!
  SUBROUTINE build_maps(sp_indxs,fp_indxs)
    USE var
    USE ad_vars
    IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: sp_indxs(Nvar*Ns*Ns*Ns,4)
    INTEGER,INTENT(INOUT) :: fp_indxs(Nvar*(Ns+1)*Ns*Ns,4)
    INTEGER :: is,js,ks,ivar,gIndx
    
    DO ks=1,Ns
       DO js=1,Ns
          DO is=1,Ns+1
             DO ivar=1,Nvar
                if(is .lt. Ns+1) then
                   CALL globalIndx_sp(ivar,is,js,ks,gIndx)
                   sp_indxs(gIndx,1:4) = (/ivar,is,js,ks/)
                endif
                CALL globalIndx_fp(ivar,is,js,ks,gIndx)
                fp_indxs(gIndx,1:4) = (/ivar,is,js,ks/)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE build_maps
  
  
END PROGRAM driver
