      subroutine gauss_xf(N,Xf,fp_type)

      
      implicit none
      integer, parameter :: RP = kind(1.d0)
      integer, intent(in) :: N,fp_type
      real(kind=RP), intent(out) :: Xf(N+1)
      real(kind=RP), parameter :: pi  = 4.0_RP*atan(1.0_RP), &
           macheps = epsilon(1.0_RP), &
           dbleeps = 10.0_RP*macheps;
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
      end subroutine
