!>@brief Get the interpolation matrix
!!
!!
subroutine get_Interpolation_matrix(nInt,nBasis,Lmat)
  implicit none

  integer,intent(in):: nBasis,nInt
  real(kind=8),intent(inout) :: LMat(nBasis,nInt)
  real(kind=8) :: xBasis(nBasis)
  real(kind=8) :: xInterp(nInt)
  integer :: ib,ii
  real(kind=8) :: low_factor,up_factor

  !Get the interpolation points[Gauss-Legeandre]
  call gauss_xs(nInt,xInterp,0)
  call gauss_xs(nBasis,xBasis,0)
  !Create the interpolation matrix:
  do ii=1,nInt
     do ib=1,nBasis-1

        low_factor = product(( xInterp(ii) - xBasis(1:ib-1))/(xBasis(ib) - xBasis(1:ib-1)))
        up_factor  = product ( ( xInterp(ii) - xBasis(ib+1:nBasis) ) / ( xBasis(ib) - xBasis(ib+1:nBasis) ) )
        LMat(ib,ii) = low_factor*up_factor
     enddo
  enddo

  ib=nBasis
  do ii=1,nInt
     LMat(ib,ii) = product(( xInterp(ii) - xBasis(1:ib-1))/(xBasis(ib) - xBasis(1:ib-1)))
  enddo

contains
  !> @brief The subroutine used for constructing
  !! the solution points, that can be of two types:
  !! Gauss-Legeandre, Gauss Chebychev.
  !!
  !!@param[in] N The number of nodes
  !!@param[out] Xs The solution points
  !!@param[in] sp_typ 0: GLegeandre, 1:Chebychev
  subroutine gauss_xs(N,Xs,sp_type)
    implicit none

    integer,       intent(in)  :: N,sp_type
    real(kind=8), intent(out) ::  Xs(N)
    real(kind=8) :: pi = 3.141592654
    integer                    :: i

    if(sp_type.eq.1) then   ! Chebyshev-Gauss points

       do i=1,N/2
          Xs(i) = 0.5*(1.0-cos(real(2*i-1,8)/real(2*N,8)*pi))
          Xs(N-i+1) = 1.0-Xs(i)
       end do
       if(mod(N,2).ne.0) Xs(N/2+1) = 0.5

    else                   ! Gauss-Legendre points

       if(N.eq.1) then
          Xs(1) = 0.0
       else if(N.eq.2) then
          Xs(1) = -sqrt(1.0/3.0)
       else if(N.eq.3) then
          Xs(1) = -sqrt(3.0/5.0)
          Xs(2) = 0.0
       else if(N.eq.4) then
          Xs(1) = -sqrt((3.0+2.0*sqrt(6.0/5.0))/7.0)
          Xs(2) = -sqrt((3.0-2.0*sqrt(6.0/5.0))/7.0)
       else if(N.eq.5) then
          Xs(1) = -sqrt(5.0+2.0*sqrt(10.0/7.0))/3.0
          Xs(2) = -sqrt(5.0-2.0*sqrt(10.0/7.0))/3.0
          Xs(3) = 0.0
       else if(N.eq.6) then
          Xs(1) = -0.9324695142031520278123016
          Xs(2) = -0.6612093864662645136613996
          Xs(3) = -0.2386191860831969086305017
       else if(N.eq.7) then
          Xs(1) = -0.9491079123427585245261897
          Xs(2) = -0.7415311855993944398638648
          Xs(3) = -0.4058451513773971669066064
          Xs(4) = 0.0
       else if(N.eq.8) then
          Xs(1) = -0.9602898564975362316835609
          Xs(2) = -0.7966664774136267395915539
          Xs(3) = -0.5255324099163289858177390
          Xs(4) = -0.1834346424956498049394761
       else if(N.eq.9) then
          Xs(1) = -0.9681602395076260898355762
          Xs(2) = -0.8360311073266357942994298
          Xs(3) = -0.6133714327005903973087020
          Xs(4) = -0.3242534234038089290385380
          Xs(5) = 0.0
       else if(N.eq.10) then
          Xs(1) = -0.9739065285171717200779640
          Xs(2) = -0.8650633666889845107320967
          Xs(3) = -0.6794095682990244062343274
          Xs(4) = -0.4333953941292471907992659
          Xs(5) = -0.1488743389816312108848260
       end if
       do i=1,N/2
          !Exploit the symmetry
          Xs(N-i+1) = -Xs(i)
       end do

       ! Map to 0 1, obtain the nodes:
       do  i=1,N
          Xs(i) = Xs(i)*0.5+0.5
       end do

    end if

  end subroutine gauss_xs
end subroutine get_Interpolation_matrix
