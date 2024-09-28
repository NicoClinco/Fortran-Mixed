!> A Program that solvea basic linear advection equation
!!
!!
!! The program explain a basic example in which we solve a linear
!! advection equation with an implicit time-stepping scheme.
!! For the moment, EULER-BDF1 is chosen.
program testImplicitSolver
  USE nlesolver_wrapper
  USE nlesolver_module, wp => nlesolver_rk
  USE FourierVar1d
  USE FourierUtils
  USE SCIFOR, only : print_matrix
  USE SF_FONTS
  USE SF_PARSE_INPUT
  implicit none
  real(kind=8),parameter :: pi = 3.14159265358979323846
  complex(kind=8) :: junit = (0.0,1.0)
  
  real(kind=8),allocatable :: x(:)
  real(kind=8) :: h,kin_vis
  real(kind=8),allocatable :: tmp_u(:)
  complex(kind=8),allocatable :: tmp_uh(:)

  
  real(kind=8),allocatable :: tmp_u_ext(:)
  complex(kind=8),allocatable :: utime(:,:)   !The solution at each time-step (complex)
  real(kind=8),allocatable :: u_phy_time(:,:) !The solution at each time-step (physical)
  real(kind=8) :: t_init,t_fin,dt,time
  integer :: n_time_steps,cur_step
  
  
  integer :: gnum_modes,i
  real(kind=8),allocatable :: solver_params(:)
  character(len=:),allocatable :: buf_title
  character(len=:),allocatable :: arg_name
  
  !-------NL_SYSTEM_SOLVER_SETTINGS----------
  real(kind=8) :: tol = 1.0e-7
  integer :: max_iter = 200
  logical :: verbose = .false.
  logical :: USE_broyden = .true.
  integer :: step_mode = 1
  integer :: n_intervals
  character(len=:),allocatable :: description
  real(kind=8) :: fmin_tol
  !------------------------------------------

  !Parse input/output:
  arg_name = 'num_modes'
  CALL parse_cmd_variable(gnum_modes,trim(arg_name),default=4)
  
  
  buf_title = '--------- IMPLICIT-EULER SOLVER FOR THE ADVECTION EQUATION ---------'
  buf_title = bold(buf_title)
  print*,buf_title
  buf_title = bold('NUM-MODES SELECTED=')
  print*,trim(buf_title),gnum_modes
  !---------------------------------------------------------

  !gnum_modes = 6
  !------------------ MESH ---------------------------------
  allocate(x(gnum_modes))
  h= 2.0*pi/gnum_modes
  do i=1,gnum_modes
     x(i) = (i-1)*h
  enddo
  !------------------ MESH ---------------------------------

  !------------------ INITIAL CONDITIONS -------------------
  allocate(tmp_u(gnum_modes),&
       tmp_uh(gnum_modes))
  allocate(solver_params(2))
  allocate(tmp_u_ext(2*gnum_modes))
  tmp_u = u0(x)
  
  tmp_uh = fwd_transform_1d(tmp_u)
  adv_vel = real(pi,8) !Advection velocity (constant)
  kin_vis = 0.001;
  solver_params(1) = adv_vel;solver_params(2) = kin_vis
  
  !------------------ INITIAL CONDITIONS --------------------
  
  !------------------ global time-stepping -----------------
  t_init = 0.0
  t_fin  = 1.0
  n_time_steps = 100
  dt = (t_fin - t_init)/real(n_time_steps)
  allocate(utime(0:n_time_steps,gnum_modes)) !The solution for each time step
  allocate(u_phy_time(0:n_time_steps,gnum_modes))
  utime(0,:) = tmp_uh(:) 
  u_phy_time(0,:) = bwd_transform_1d(utime(0,:))
  
  fmin_tol = 1.0e-2
  n_intervals = 2
  description = 'Non linear system : assembling'
  description = bold(description)
 
  
  CALL set_parameters(solver_params)
  CALL setParametricJac(jac)
  CALL setParametricFunc(resid)
  CALL initNLsys_global(2*gnum_modes,2*gnum_modes,&
       step_mode =step_mode,&
       USE_broyden=USE_broyden,&
       n_intervals = n_intervals,&
       fmin_tol = fmin_tol,&
       description = description,&
       max_iter = max_iter,&
       tol = tol,&
       verbose=verbose)
  
  CALL INITSOL(tmp_uh,gnum_modes)

  !TIME-STEPPING--------------------------------------------
  cur_step = 0
  do while(time<t_fin)
     cur_step = cur_step + 1
     time = time + dt
     
     write(*,'("Timestep: ",i2)') cur_step
     write(*,'("Time: ",e8.3)') time

     !Extended solution [real,imag]
     tmp_u_ext(:) = m_sol%m_ruHat
     
     CALL NLSYS_SOLVE(tmp_u_ext)
     
     !Update
     m_sol%m_ruHat = tmp_u_ext(:)
     utime(cur_step,:) = tmp_u_ext(1:gnum_modes) + &
          junit*tmp_u_ext(gnum_modes+1:2*gnum_modes)
     u_phy_time(cur_step,:) = bwd_transform_1d(utime(cur_step,:))
     write(*,'("Current-solution= ",*(G9.3,","))') u_phy_time(cur_step,:)
  enddo
  !---------------------------------------------------------
  print*,trim(bold('END TIME STEPPING SCHEME'))
  print*,trim(bold('Saving the global solution in "sol_euler.dat"'))
  CALL print_matrix(u_phy_time,'sol_euler.dat',w=8,d=3)

CONTAINS
  !-----------------------------------------------
  !-----------------------------------------------
  !!The initial solution for our system
  !!@param[in] x: The discrete mesh
  function u0(x)
    implicit none
    real(kind=8),intent(in) :: x(:)
    real(kind=8),allocatable :: u0(:)
    allocate(u0(size(x)))
    u0(:) = cos(x(:))
  end function u0
  !-----------------------------------------------
  !-----------------------------------------------  

  !> Advection operator along the x-direction.
  !! 
  !! @param [in]   uhat  : The solution vector    
  !! @return[out]  op    : The resultant operator
  !! @param [in]   alpha : The
  !!
  !! @note
  !!  -Since the system is complex, the advection operator
  !!   is of length 2*num_modes. In the first part we have the
  !!   real part, while in the second part we have the complex
  !!   part
  !! 
  !! - The Advection operator is the following:
  !!   C_{op} = -i*k*u_{k}
  !!@endnote
  SUBROUTINE xAdvectionOp(uhat,op,alpha)
    USE fftpack, only : fftfreq 
    implicit none
    real(kind=8),dimension(:),intent(in)  :: uhat
    real(kind=8),dimension(:),intent(out) :: op
    real(kind=8),intent(in)               :: alpha
    integer :: num_modes
    complex(kind=8),allocatable :: ks(:)
    complex(kind=8),allocatable :: uhatcmplx(:)
    complex(kind=8) :: ju = (0.0,1.0)
    integer :: dim

    dim = int(size(uhat)/2) !For the complex solution
    allocate(uhatcmplx(dim),ks(dim))
    uhatcmplx(:) = uhat(1:dim)+ju*uhat(dim+1:2*dim)
    ks = cmplx(fftfreq(dim))

    uhatcmplx(:) = -ju*alpha*ks*uhatcmplx
    
    op(1:dim) = real(uhatcmplx(:),8)
    op(dim+1:2*dim) = aimag(uhatcmplx(:))
  end SUBROUTINE xAdvectionOp
  !-----------------------------------------------
  !-----------------------------------------------

  !>@brief Assemble the diffusion operator with the
  !!Fourier-Galerkin space discretization.
  !!
  !!@param[in]   uhat: The solution vector extended
  !!@return[out] op  : The operator
  !!@param[in]   nu  : The kinematic viscosity
  !!
  !!@note
  !!  The diffusion operator is the following:
  !!   D_{op} = (k^2)*\hat{u_k}
  !!@endnote
  SUBROUTINE xDiffusionOp(uhat,op,nu)
    USE fftpack, only : fftfreq 
    IMPLICIT NONE
    real(kind=8),dimension(:),intent(in)  :: uhat
    real(kind=8),dimension(:),intent(out) :: op
    real(kind=8),intent(in)               :: nu
    integer :: dim
    complex(kind=8) :: ju = (0.0,1.0)
    complex(kind=8),allocatable :: uhatcmplx(:)
    complex(kind=8),allocatable :: ks(:)

    dim = int(size(uhat)/2)
    allocate(uhatcmplx(dim),ks(dim))
    ks = cmplx(fftfreq(dim))
    uhatcmplx(:) = uhat(1:dim)+ju*uhat(dim+1:2*dim)
    uhatcmplx(:) = (ks**2)*uhatcmplx(:)*nu

    op(1:dim) = -real(uhatcmplx(:),8)
    op(dim+1:2*dim) = -aimag(uhatcmplx(:))
    
  END SUBROUTINE xDiffusionOp
  
  
  !! The residual of our system (non-linear)
  !! [extended version for the non-linear solver]
  !!
  !!@note
  !!The residual is the extended vector that contains
  !!in the first part the real part and in the second
  !!the imaginary values
  !!@endnote
  !!
  !!@param[in]    u:
  !!@param[inout] res:   The residual
  !!@param[in]    alpha: It is a vector that contains the
  !!                     parameters that one can use in the
  !!                     nonlinear solver
  SUBROUTINE resid(u,res,alpha)
    implicit none
    real(kind=8),dimension(:),intent(in) :: u
    real(kind=8),dimension(:),intent(inout) :: res
    real(kind=8),dimension(:),intent(in) :: alpha
    real(kind=8),dimension(:),allocatable :: Fu
    real(kind=8),dimension(:),allocatable :: Du

    integer :: num_params
    
    num_params = size(alpha)
    allocate(Fu(size(u)),Du(size(u)))
    Fu(:) = 0.0; Du(:) = 0.0;
    CALL xAdvectionOp(u,Fu(1:size(u)),alpha(1))

    if(num_params .ge. 2) then
       CALL xDiffusionOp(u,Du(1:size(u)),alpha(2))
    endif

    res(:) = u(:)+ dt*Fu(:) - dt*Du(:)- m_sol%m_ruHat(:)
  END SUBROUTINE resid
  !-----------------------------------------------
  !-----------------------------------------------
  
  !>@brief The Jacobian of the resulting non linear system.
  !!
  !!@param[in]    u:      The unknown variable
  !!@param[inout] J:      The Jacobian matrix
  !!@param[in]    \alpha: The parameters that contains the advection and diffusivity.
  SUBROUTINE jac(u,J,alpha)
    USE fftpack, only : fftfreq
    
    implicit none
    real(kind=8),dimension(:),intent(in) :: u
    real(kind=8),dimension(:,:),intent(inout) :: J
    real(kind=8),dimension(:),intent(in) :: alpha
    integer :: i,num_modes,sub
    real(kind=8),allocatable :: ks(:)
    
    num_modes = int(size(u)/2)
    J(:,:) = 0.0
    allocate(ks(num_modes))
    ks = real(fftfreq(num_modes),8)
    
    do sub=1,2
       do i=0,num_modes-1
          J(1+(sub-1)*num_modes+i,1+(sub-1)*num_modes+i) =  &
               real(1.0,8)+ks(i+1)**2*dt*alpha(2)
       enddo
    enddo
    !Sub-blocks:
    do i =0,num_modes-1
       J(1+i,1+(1)*num_modes+i) = -alpha(1)*real(ks(i+1))*dt
       J(1+(1)*num_modes+i,1+i) = +alpha(1)*real(ks(i+1))*dt
    enddo
  end SUBROUTINE jac
  


end program testImplicitSolver
