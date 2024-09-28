!! The module that contains the global solution of our solver
MODULE FourierVar1d

  PUBLIC !Default
  
  !! The Fourier solution stored in the MODULE.
  TYPE :: fourier_solution
     INTEGER :: m_num_modes
     REAL(kind=8),allocatable    :: m_gu(:)      !global solution
     COMPLEX(kind=8),allocatable :: m_uHat(:) !Complex
     REAL(kind=8),allocatable    :: m_ruHat(:)   !Extended [real]
   CONTAINS
     PROCEDURE :: update_real
  END TYPE fourier_solution

  ! The global solution stored in the module:
  type(fourier_solution) :: m_sol

  real(kind=8) :: adv_vel   !The advection velocity
  real(kind=8) :: nu        !The diffusion term
   
CONTAINS

  !>@brief Initialize the solution
  !!
  !!@param[in] uhat:     The complex solution
  !!@param[in] nummodes: The number of modes
  SUBROUTINE InitSol(uhat,nummodes)
    implicit none
    complex(kind=8),intent(in) :: uhat(:)
    integer,intent(in) :: nummodes
    m_sol%m_num_modes = nummodes
    allocate(m_sol%m_uHat(nummodes))
    m_sol%m_uHat(:) = uhat(:)

    !Update the real vector(extended) with the hat one:
    allocate(m_sol%m_ruHat(2*nummodes))
    call m_sol%update_real(uhat)
  end SUBROUTINE InitSol
!----------------------------------------------------------------
  !> @brief Copy the complex coefficients to a vector 2 times long
  !!
  !! @param[in] uhat: The complex coefficients.
  SUBROUTINE update_real(msol,uhat)
    implicit none
    class(fourier_solution),intent(inout) :: msol
    complex(kind=8),intent(in) :: uhat(:)
    integer :: nummodes
    nummodes = msol%m_num_modes
    if(allocated(msol%m_ruHat)) then
       msol%m_ruHat(1:nummodes) = real(uhat(:))
       msol%m_ruHat(nummodes+1:2*nummodes) = aimag(uhat(:))
    endif
  end SUBROUTINE update_real
!---------------------------------------------------------------- 
END MODULE FourierVar1d

