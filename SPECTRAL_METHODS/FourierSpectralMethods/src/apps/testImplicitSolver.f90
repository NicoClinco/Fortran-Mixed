program testImplicitSolve
  use nlesolver_wrapper
  use nlesolver_module, wp => nlesolver_rk
  use FourierVar1d
  use SCIFOR, only : print_matrix
  use SF_FONTS
  use SF_PARSE_INPUT
  implicit none
  real(kind=8),parameter :: pi = 3.14159265358979323846
  complex(kind=8) :: junit = (0.0,1.0)
  
  real(kind=8),allocatable :: x(:)
  real(kind=8) :: h
  real(kind=8),allocatable :: tmp_u(:)
  complex(kind=8),allocatable :: tmp_uh(:)

  
  real(kind=8),allocatable :: tmp_u_ext(:)
  complex(kind=8),allocatable :: utime(:,:) !The solution at each time-step (uHat)
  real(kind=8),allocatable :: u_phy_time(:,:) !The soltion at each time-step (physical)
  real(kind=8) :: t_init,t_fin,dt,time
  integer :: n_time_steps,cur_step
  
  
  integer :: gnum_modes
  integer :: i
  real(kind=8),allocatable :: adv_vels(:)
  character(len=:),allocatable :: buf_title
  character(len=:),allocatable :: arg_name
  
  !-------NL_SYSTEM_SOLVER_SETTINGS----------
  real(kind=8) :: tol = 1.0e-7
  integer :: max_iter = 200
  logical :: verbose = .false.
  logical :: use_broyden = .true.
  integer :: step_mode = 1
  integer :: n_intervals
  character(len=:),allocatable :: description
  real(kind=8) :: fmin_tol
  !------------------------------------------

  !Parse input/output:
  arg_name = 'num_modes'
  call parse_cmd_variable(gnum_modes,trim(arg_name),default=4)
  
  
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
  !Temporary:
  allocate(tmp_u(gnum_modes),&
       tmp_uh(gnum_modes))
  allocate(adv_vels(2*gnum_modes))
  allocate(tmp_u_ext(2*gnum_modes))
  tmp_u = u0(x)
  
  tmp_uh = fwd_transform(tmp_u)
  adv_vel = real(2*pi,8) !Advection velocity
  adv_vels(:) = adv_vel

  !------------------ global time-stepping -----------------
  t_init = 0.0
  t_fin  = 1.0
  n_time_steps = 100
  dt = (t_fin - t_init)/real(n_time_steps)
  allocate(utime(0:n_time_steps,gnum_modes)) !The solution for each time step
  allocate(u_phy_time(0:n_time_steps,gnum_modes))
  utime(0,:) = tmp_uh(:)
  u_phy_time(0,:) = bwd_transform(utime(0,:))
  
  fmin_tol = 1.0e-2
  n_intervals = 2
  description = 'Non linear system : assembling'
  description = bold(description)
 
  
  call set_parameters(adv_vels)
  call setParametricJac(jac)
  call setParametricFunc(resid)
  call initNLsys_global(2*gnum_modes,2*gnum_modes,&
       step_mode =step_mode,&
       use_broyden=use_broyden,&
       n_intervals = n_intervals,&
       fmin_tol = fmin_tol,&
       description = description,&
       max_iter = max_iter,&
       tol = tol,&
       verbose=verbose)
  
  call INITSOL(tmp_uh,gnum_modes)

  !TIME-STEPPING--------------------------------------------
  cur_step = 0
  do while(time<t_fin)
     cur_step = cur_step + 1
     time = time + dt
     
     write(*,'("Timestep: ",i2)') cur_step
     write(*,'("Time: ",e8.3)') time
     
     tmp_u_ext(:) = m_sol%m_ruHat
     
     call NLSYS_SOLVE(tmp_u_ext)
     
     !Update
     m_sol%m_ruHat = tmp_u_ext(:)
     utime(cur_step,:) = tmp_u_ext(1:gnum_modes) + &
          junit*tmp_u_ext(gnum_modes+1:2*gnum_modes)
     u_phy_time(cur_step,:) = bwd_transform(utime(cur_step,:))
     write(*,'("Current-solution= ",*(G9.3,","))') u_phy_time(cur_step,:)
     
     !write(*,'("Current-solution=",*(G9.3,","))') tmp_u_ext
  enddo
  !---------------------------------------------------------
  print*,trim(bold('END TIME STEPPING SCHEME'))
  print*,trim(bold('Saving the global solution in "sol_euler.dat"'))
  CALL print_matrix(u_phy_time,'sol_euler.dat',w=8,d=3)

contains
  !-----------------------------------------------
  !-----------------------------------------------
  !!The initial solution for our system
  function u0(x)
    implicit none
    real(kind=8),intent(in) :: x(:)
    real(kind=8),allocatable :: u0(:)
    allocate(u0(size(x)))
    u0(:) = cos(x(:))
  end function u0
  !-----------------------------------------------
  !-----------------------------------------------

  !! Transform our solution from the physical
  !! to the frequency space
  function fwd_transform(u) result(uh)
    use fftw3_wrapper
    implicit none
    real(kind=8),intent(in) :: u(:)
    complex(kind=8),allocatable :: uh(:)
    integer :: nummodes
    nummodes = size(u)
    allocate(uh(nummodes))
    call fft(uh,u,num=[nummodes])
  end function fwd_transform
  !-----------------------------------------------
  !-----------------------------------------------

  !! Perform the backward transformartion for
  function bwd_transform(uhat) result(u)
    use fftw3_wrapper
    implicit none
    complex(kind=8),intent(in) :: uhat(:)
    real(kind=8),allocatable :: u(:)
    integer :: nummodes
    nummodes = size(uhat)
    call ifft(u,uhat,num=[nummodes])
  end function bwd_transform
  !-----------------------------------------------
  !-----------------------------------------------
  
  !! Advection operator along the x direction
  !! Note that here we pass (real)uhat that is
  !! the extended vector of u
  !!
  !!uhat : The sol. vector (unknown)
  !!op   : The result
  !!alpha: The advection velocity 
  subroutine xAdvectionOp(uhat,op,alpha)
    use fftpack, only : fftfreq 
    implicit none
    real(kind=8),dimension(:),intent(in) :: uhat
    real(kind=8),dimension(:),intent(out) :: op
    real(kind=8),intent(in) :: alpha
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
    
  end subroutine xAdvectionOp
  !-----------------------------------------------
  !-----------------------------------------------
  
  !! The residual of our system (non-linear)
  !! [extended version for the non-linear solver]
  !!
  !!The residual is the extended vector for our
  !!system in the Fourier space
  subroutine resid(u,res,alpha)
    !use FourierVar1d, only: m_sol !uHat,m_ruHat
    implicit none
    real(kind=8),dimension(:),intent(in) :: u
    real(kind=8),dimension(:),intent(inout) :: res
    real(kind=8),dimension(:),intent(in) :: alpha
    real(kind=8),dimension(:),allocatable :: Fu
    
    allocate(Fu(size(u)))
    Fu(:) = 0.0
    call xAdvectionOp(u,Fu(1:size(u)),alpha(1))
    res(:) = u(:) - dt*Fu(:) - m_sol%m_ruHat(:)
    
  end subroutine resid
  !-----------------------------------------------
  !-----------------------------------------------
  
  !! The jacobian matrix for the implicit
  !! Euler first order system
  subroutine jac(u,J,alpha)
    use fftpack, only : fftfreq
    
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
          J(1+(sub-1)*num_modes+i,1+(sub-1)*num_modes+i) = real(1.0,8)
       enddo
    enddo
    do i =0,num_modes-1
       J(1+i,1+(1)*num_modes+i) = -alpha(1)*real(ks(i+1))*dt
       J(1+(1)*num_modes+i,1+i) = +alpha(1)*real(ks(i+1))*dt
    enddo
  end subroutine jac
  


end program testImplicitSolve
