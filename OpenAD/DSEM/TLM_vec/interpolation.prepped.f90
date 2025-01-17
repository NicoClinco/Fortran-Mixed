!> @brief The module that stores the matrixes for interpolation
!! and the global variables
MODULE var

  INTEGER,PARAMETER :: Nvar = 3
  INTEGER,PARAMETER :: Nc = 50
  INTEGER,PARAMETER :: Ns = 3

  !SP,FP
  REAL(KIND=8) :: Xs(Ns)
  REAL(KIND=8) :: Xf(Ns+1)

  !MAT (Created in create_interpolation_matrixes)
  REAL(KIND=8) :: Lmat(Ns+1,Ns)
  REAL(KIND=8) :: Mmat(Ns+1,Ns)


  ![3d vector] stored here
  REAL(KIND=8) :: Qs(Nvar,Ns,Ns,Ns,Nc)


CONTAINS

  !>@brief Create the discrete solution points for the SD scheme
  !! 
  SUBROUTINE make_dpoints
    IMPLICIT NONE

    CALL gauss_xs(Ns,Xs,0)
    CALL gauss_xf(Ns+1,Xf,0)

  END SUBROUTINE make_dpoints

  SUBROUTINE create_interpolation_matrixes
    IMPLICIT NONE
    INTEGER, parameter :: RP = kind(1.d0)

    integer :: i,j,jj,ifp,s,k,is
    real(kind=RP) :: num,den,xval,summ,norm

    ! Matrix Lmat is of size N+1 by N
    ! it multiplies the solution at the N solution points
    ! to give solution at the N+1 flux points

    ! ------------------------
    ! Computing Lmat
    ! ------------------------

    do i=1,(Ns+1)/2
       norm = 0.0_RP
       do j=1,Ns-1
          num = 1.0_RP
          den = 1.0_RP
          do jj=1,Ns
             if(jj.ne.j) then
                num = num*(Xf(i)-Xs(jj))
                den = den*(Xs(j)-Xs(jj))
             end if
          end do

          Lmat(i,j) = num/den
          Lmat(Ns+2-i,Ns+1-j) = Lmat(i,j)
          norm = norm + Lmat(i,j)
       end do
       Lmat(i,Ns) = 1.0_RP-norm
       Lmat(Ns+2-i,1) = Lmat(i,Ns)
    end do
    if(mod(Ns+1,2).ne.0) then
       i=(Ns+1)/2+1
       norm = 0.0_RP
       do j=1,Ns/2-1
          num = 1.0_RP
          den = 1.0_RP
          do jj=1,Ns
             if(jj.ne.j) then
                num = num*(Xf(i)-Xs(jj))
                den = den*(Xs(j)-Xs(jj))
             end if
          end do

          Lmat(i,j) = num/den
          Lmat(i,Ns+1-j) = Lmat(i,j)
          norm = norm + Lmat(i,j)
       end do
       Lmat(i,Ns/2) = 0.5_RP-norm
       Lmat(i,Ns/2+1) = Lmat(i,Ns/2)
    end if
    !do i=1,Ns+1
    !do j=1,Ns
    !  call Lagrange(Xf(i),j,Xs,Ns,Lmat(i,j))
    !end do
    !end do

    ! ------------------------
    ! Computing Mmat
    ! ------------------------
    do is=1,Ns/2
       norm = 0.0_RP
       do ifp=1,Ns

          xval = Xs(is)
          summ = 0.0_RP
          do k=1,Ns+1
             if(k.ne.ifp) then

                num = 1.0_RP
                den = 1.0_RP
                do s=1,Ns+1
                   if(s.ne.ifp.and.s.ne.k) then
                      num = num*(xval-Xf(s))
                   end if
                   if(s.ne.ifp) then
                      den = den*(Xf(ifp)-Xf(s))
                   end if
                end do
                summ = summ+num/den
             end if
          end do

          Mmat(ifp,is) = summ
          Mmat(Ns+2-ifp,Ns+1-is) = -Mmat(ifp,is)
          norm = norm + Mmat(ifp,is)

       end do
       Mmat(Ns+1,is) = -norm
       Mmat(1,Ns+1-is) = -Mmat(Ns+1,is)
    end do
    if(mod(Ns,2).ne.0) then
       is = Ns/2+1
       do ifp=1,(Ns+1)/2

          xval = Xs(is)
          summ = 0.0_RP
          do k=1,Ns+1
             if(k.ne.ifp) then

                num = 1.0_RP
                den = 1.0_RP
                do s=1,Ns+1
                   if(s.ne.ifp.and.s.ne.k) then
                      num = num*(xval-Xf(s))
                   end if
                   if(s.ne.ifp) then
                      den = den*(Xf(ifp)-Xf(s))
                   end if
                end do
                summ = summ+num/den
             end if
          end do

          Mmat(ifp,is) = summ
          Mmat(Ns+2-ifp,is) = -Mmat(ifp,is)
       end do
    end if
  END SUBROUTINE create_interpolation_matrixes



END MODULE var

!>@brief The module that contains the variable to be differentiated
!!   
!!
MODULE AD_vars
  USE var, only : Nvar,Ns,Nc
  !Flattened 3D-vectors
  REAL(KIND=8) :: FLATstencilQs(Nvar*Ns*Ns*Ns)
  REAL(KIND=8) :: FLATstencilQfp(Nvar*(Ns+1)*(Ns)*(Ns))
  
CONTAINS
  
  !@brief Computation of the global index for the solution
  !points when we flatten a vector:
  SUBROUTINE globalIndx_sp(ivar,is,js,ks,gisp)
    USE var, only : Nvar,Ns,Nc
    IMPLICIT NONE
    INTEGER :: ivar,is,js,ks
    INTEGER :: gisp
    gisp = (Ns)*Ns*Nvar*(ks-1)+&
           (Ns)*(js-1)*Nvar+&
           (is-1)*Nvar+ivar
    return

    !ns*ns*nvar*(Ns-1)+Ns*(ns-1)*nvar+(NS-1)*NVAR+NVAR
  END SUBROUTINE globalIndx_sp


  !@brief Computation of the global index for the flux points
  !when we flatten a flux-point vector of type 1:
  SUBROUTINE globalIndx_fp(ivar,ifp,js,ks,gifp)
    USE var, only : Nvar,Ns,Nc
    IMPLICIT NONE
    INTEGER :: ivar,ifp,js,ks
    INTEGER :: gifp

    
    gifp = (Ns+1)*Ns*Nvar*(ks-1)+&
           (Ns+1)*(js-1)*Nvar+&
           (ifp-1)*Nvar+ivar
    return 
  END SUBROUTINE globalIndx_fp
END MODULE AD_vars

!>@brief Interpolate the solution to the flux-points
!!@note
!!  This subroutine will be transformed such that the inputs become
!!  differentiable
!!@endnote
SUBROUTINE calc_Qfp
  USE var
  USE AD_vars
  IMPLICIT NONE

  INTEGER :: ivar,is,ifp,jfp,kfp
  INTEGER :: igsp,igfp


  !$openAD INDEPENDENT(FLATstencilQs)
  !Put to zero the flux-points:
  flatstencilQfp(:) = 0.0_8
  
  ivar=1
  !Interpolate from solution points to the flux-points:
  do kfp=1,Ns
     do jfp=1,Ns
        do ifp=1,Ns+1
           do is=1,Ns  
              do ivar=1,Nvar
                 CALL globalIndx_fp(ivar,ifp,jfp,kfp,igfp)
                 CALL globalIndx_sp(ivar,is,jfp,kfp,igsp)
                 flatstencilQfp(igfp) = flatstencilQfp(igfp)+Lmat(ifp,is)*flatstencilQs(igsp)
              enddo
           enddo
        enddo
     enddo
  enddo
  
  !$openAD DEPENDENT(FLATstencilQfp)
END SUBROUTINE calc_Qfp
!===========================================================
subroutine gauss_xs(N,Xs,sp_type)
  implicit none
  integer, parameter :: RP = kind(1.d0);
  real(kind=RP), parameter :: pi  = 3.1415926

  integer,       intent(in)  :: N,sp_type
  real(kind=RP), intent(out) :: Xs(N)
  integer                    :: i

  if(sp_type.eq.1) then   ! Chebyshev-Gauss points

     do i=1,N/2
        Xs(i) = 0.5_RP*(1.0_RP-cos(real(2*i-1,RP)/real(2*N,RP)*pi))
        Xs(N-i+1) = 1.0_RP-Xs(i)
     end do
     if(mod(N,2).ne.0) Xs(N/2+1) = 0.5_RP

  else                   ! Gauss-Legendre points

     if(N.eq.1) then
        Xs(1) = 0.0_RP
     else if(N.eq.2) then
        Xs(1) = -sqrt(1.0_RP/3.0_RP)
     else if(N.eq.3) then
        Xs(1) = -sqrt(3.0_RP/5.0_RP)
        Xs(2) = 0.0_RP
     else if(N.eq.4) then
        Xs(1) = -sqrt((3.0_RP+2.0_RP*sqrt(6.0_RP/5.0_RP))/7.0_RP)
        Xs(2) = -sqrt((3.0_RP-2.0_RP*sqrt(6.0_RP/5.0_RP))/7.0_RP)
     else if(N.eq.5) then
        Xs(1) = -sqrt(5.0_RP+2.0_RP*sqrt(10.0_RP/7.0_RP))/3.0_RP
        Xs(2) = -sqrt(5.0_RP-2.0_RP*sqrt(10.0_RP/7.0_RP))/3.0_RP
        Xs(3) = 0.0_RP
     else if(N.eq.6) then
        Xs(1) = -0.9324695142031520278123016_RP
        Xs(2) = -0.6612093864662645136613996_RP
        Xs(3) = -0.2386191860831969086305017_RP
     else if(N.eq.7) then
        Xs(1) = -0.9491079123427585245261897_RP
        Xs(2) = -0.7415311855993944398638648_RP
        Xs(3) = -0.4058451513773971669066064_RP
        Xs(4) = 0.0_RP
     else if(N.eq.8) then
        Xs(1) = -0.9602898564975362316835609_RP
        Xs(2) = -0.7966664774136267395915539_RP
        Xs(3) = -0.5255324099163289858177390_RP
        Xs(4) = -0.1834346424956498049394761_RP
     else if(N.eq.9) then
        Xs(1) = -0.9681602395076260898355762_RP
        Xs(2) = -0.8360311073266357942994298_RP
        Xs(3) = -0.6133714327005903973087020_RP
        Xs(4) = -0.3242534234038089290385380_RP
        Xs(5) = 0.0_RP
     else if(N.eq.10) then
        Xs(1) = -0.9739065285171717200779640_RP
        Xs(2) = -0.8650633666889845107320967_RP
        Xs(3) = -0.6794095682990244062343274_RP
        Xs(4) = -0.4333953941292471907992659_RP
        Xs(5) = -0.1488743389816312108848260_RP
     end if
     do i=1,N/2
        !Exploit the symmetry
        Xs(N-i+1) = -Xs(i)
     end do

     ! Map to 0 1, obtain the nodes:
     do  i=1,N
        Xs(i) = Xs(i)*0.5_RP+0.5_RP
     end do

  end if
end subroutine gauss_xs
!===========================================================

subroutine gauss_xf(N,Xf,fp_type)
  implicit none
  integer, parameter :: RP = kind(1.d0)
  integer, intent(in) :: N,fp_type
  real(kind=RP), intent(out) :: Xf(N+1)
  real(kind=RP), parameter :: pi  = 3.1415926
  integer :: ifp

  if(fp_type.eq.1) then   ! Chebyshev-Lobatto points

     do ifp=1,(N+1)/2     
        Xf(ifp) = 0.5_RP*(1.0_RP-cos(real(ifp-1,RP)*pi/real(N,RP)));
        Xf(N+2-ifp) = 1.0_RP-Xf(ifp)
     end do
     if(mod(N+1,2).ne.0) Xf((N+1)/2+1) = 0.5_RP
     Xf(1)   = 0.0_RP
     Xf(N+1) = 1.0_RP

  else                   ! Gauss-Legendre points

     Xf(1) = -1.0_RP
     Xf(N+1) = 1.0_RP
     if(N.eq.2) then
        Xf(2) = 0.0_RP
     else if(N.eq.3) then
        Xf(2) = -sqrt(1.0_RP/3.0_RP)
        Xf(3) = -Xf(2)
     else if(N.eq.4) then
        Xf(2) = -sqrt(3.0_RP/5.0_RP)
        Xf(3) = 0.0_RP
        Xf(4) = -Xf(2)
     else if(N.eq.5) then
        Xf(2) = -sqrt((3.0_RP+2.0_RP*sqrt(6.0_RP/5.0_RP))/7.0_RP)
        Xf(3) = -sqrt((3.0_RP-2.0_RP*sqrt(6.0_RP/5.0_RP))/7.0_RP)
        Xf(4) = -Xf(3)
        Xf(5) = -Xf(2)
     else if(N.eq.6) then
        Xf(2) = -sqrt(5.0_RP+2.0_RP*sqrt(10.0_RP/7.0_RP))/3.0_RP
        Xf(3) = -sqrt(5.0_RP-2.0_RP*sqrt(10.0_RP/7.0_RP))/3.0_RP
        Xf(4) = 0.0_RP
        Xf(5) = -Xf(3)
        Xf(6) = -Xf(2)
     else if(N.eq.7) then
        Xf(2) = -0.9324695142031520278123016_RP
        Xf(3) = -0.6612093864662645136613996_RP
        Xf(4) = -0.2386191860831969086305017_RP
        Xf(5) = -Xf(4)
        Xf(6) = -Xf(3)
        Xf(7) = -Xf(2)
     else if(N.eq.8) then
        Xf(2) = -0.9491079123427585245261897_RP
        Xf(3) = -0.7415311855993944398638648_RP
        Xf(4) = -0.4058451513773971669066064_RP
        Xf(5) = 0.0_RP
        Xf(6) = -Xf(4)
        Xf(7) = -Xf(3)
        Xf(8) = -Xf(2)
     else if(N.eq.9) then
        Xf(2) = -0.9602898564975362316835609_RP
        Xf(3) = -0.7966664774136267395915539_RP
        Xf(4) = -0.5255324099163289858177390_RP
        Xf(5) = -0.1834346424956498049394761_RP
        Xf(6) = -Xf(5)
        Xf(7) = -Xf(4)
        Xf(8) = -Xf(3)
        Xf(9) = -Xf(2)
     else if(N.eq.10) then
        Xf(2) = -0.9681602395076260898355762_RP
        Xf(3) = -0.8360311073266357942994298_RP
        Xf(4) = -0.6133714327005903973087020_RP
        Xf(5) = -0.3242534234038089290385380_RP
        Xf(6) = 0.0_RP
        Xf(7) = -Xf(5)
        Xf(8) = -Xf(4)
        Xf(9) = -Xf(3)
        Xf(10)= -Xf(2)
     end if

     ! Map to 0 1
     do  ifp=1,N+1
        Xf(ifp) = Xf(ifp)*0.5_RP+0.5_RP
     end do

  end if

  return
end subroutine gauss_xf



