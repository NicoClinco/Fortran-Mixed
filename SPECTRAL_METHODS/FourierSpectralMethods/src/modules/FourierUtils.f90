!> A module that contains utilities for the spectral solver
!!
!!
MODULE FourierUtils
  PUBLIC


CONTAINS
  !> Transform our solution from the physical to the frequency space
  function fwd_transform_1d(u) result(uh)
    use fftw3_wrapper
    implicit none
    real(kind=8),intent(in) :: u(:)
    complex(kind=8),allocatable :: uh(:)
    integer :: nummodes
    nummodes = size(u)
    allocate(uh(nummodes))
    CALL fft(uh,u,num=[nummodes])
  end function fwd_transform_1d
  !-----------------------------------------------
  !-----------------------------------------------

  !> Perform the backward transform from the complex to the physical space
  function bwd_transform_1d(uhat) result(u)
    use fftw3_wrapper
    implicit none
    complex(kind=8),intent(in) :: uhat(:)
    real(kind=8),allocatable :: u(:)
    integer :: nummodes
    nummodes = size(uhat)
    CALL ifft(u,uhat,num=[nummodes])
  end function bwd_transform_1d
  !-----------------------------------------------
  !-----------------------------------------------

END MODULE FourierUtils
