
subroutine gauss_xs(N,Xs,sp_type)


  implicit none
  integer, parameter :: RP = kind(1.d0);
  real(kind=RP), parameter :: pi  = 4.0_RP*atan(1.0_RP), &
       macheps = epsilon(1.0_RP), &
       dbleeps = 10.0_RP*macheps;

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
