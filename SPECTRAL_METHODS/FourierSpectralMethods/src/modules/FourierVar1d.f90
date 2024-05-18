!! The module that contains the global
!! solution for our system of Fourier
module FourierVar1d

  public !Default
  
  !! The global fourier solution [1D]
  !!
  type :: fourier_solution
     integer :: m_num_modes
     real(kind=8),allocatable :: m_gu(:)      !global solution
     complex(kind=8),allocatable :: m_uHat(:) !Complex
     real(kind=8),allocatable :: m_ruHat(:)   !Extended [real]
   contains
     procedure :: update_real
  end type fourier_solution

  type(fourier_solution) :: m_sol

  real(kind=8) :: adv_vel   !The advection velocity
  real(kind=8) :: nu        !The diffusion term
   
contains

  !!Initialize the solution with the complex one:
  subroutine InitSol(uhat,nummodes)
    implicit none
    complex(kind=8),intent(in) :: uhat(:)
    integer,intent(in) :: nummodes
    m_sol%m_num_modes = nummodes
    allocate(m_sol%m_uHat(nummodes))
    m_sol%m_uHat(:) = uhat(:)

    !Update the real vector(extended) with the hat one:
    allocate(m_sol%m_ruHat(2*nummodes))
    call m_sol%update_real(uhat)
  end subroutine InitSol
!----------------------------------------------------------------
  !!Update the solution vector of the real variables by
  !!uhat
  subroutine update_real(msol,uhat)
    implicit none
    class(fourier_solution),intent(inout) :: msol
    complex(kind=8),intent(in) :: uhat(:)
    integer :: nummodes
    nummodes = msol%m_num_modes
    if(allocated(msol%m_ruHat)) then
       msol%m_ruHat(1:nummodes) = real(uhat(:))
       msol%m_ruHat(nummodes+1:2*nummodes) = aimag(uhat(:))
    endif
  end subroutine update_real
!---------------------------------------------------------------- 
end module FourierVar1d

