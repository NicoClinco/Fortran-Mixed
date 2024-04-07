  !A basic module that serves as a wrapper
!for nlesolver_module.
module nlesolver_wrapper

  use nlesolver_module
  private 
  real(kind=8),allocatable :: m_params(:)
  
  !Pointers to parametric functions:
  procedure(ParametricFunction),pointer :: pParFun => null()
  procedure(ParametricJacFunc),pointer :: pParJacFun => null()

  !Pointers to standard functions (no-parametric)
  procedure(standardFunction),pointer :: pFun => null()
  procedure(standardJacFunc),pointer :: pJacFun => null()

  logical :: m_is_parametric_func = .false.
  logical :: m_is_parametric_jac = .false.

  logical :: m_is_func = .false.
  logical :: m_is_jac = .false.


  

  !> @brief A basic structure that encapsulate
  !! all the functionalities
  type :: solver_info
     integer :: n
     integer :: m
     integer :: max_iter
     real(kind=8) :: tol
     logical :: verbose

     !TYPE OF SOLVER
     !AND OTHER SETTINGS
     type(nlesolver_type) :: m_solver
     real(kind=8) :: alpha
     logical :: use_broyden
     integer :: step_mode
     integer :: n_intervals
     integer :: istat
     character(len=:), allocatable :: message
     character(len=:),allocatable :: description
     real(kind=8) :: fmin_tol
  end type solver_info
  
  !Here we create the solver object:
  type(solver_info) :: solver_obj

  
  !> @brief interface to specify a parametric function
  !! to be used in conjuction with the nle_solver
  abstract interface
     subroutine ParametricFunction(x,fval,params)
       real(kind=8),intent(in) :: x(:)
       real(kind=8),intent(inout) :: fval(:)
       real(kind=8),intent(in) :: params(:) 
     end subroutine  ParametricFunction

     subroutine ParametricJacFunc(x,gval,params)
       real(kind=8),intent(in) :: x(:)
       real(kind=8),intent(inout) :: gval(:,:)
       real(kind=8),intent(in) :: params(:) 
     end subroutine ParametricJacFunc

     subroutine standardFunction(x,fval)
       real(kind=8),intent(in) :: x(:)
       real(kind=8),intent(inout) :: fval(:)
     end subroutine  StandardFunction

     subroutine standardJacFunc(x,gval)
       real(kind=8),intent(in) :: x(:)
       real(kind=8),intent(inout) :: gval(:,:)
     end subroutine StandardJacFunc
  end interface


  !Set to public the procedures:
  public setParametricFunc, setParametricJac, setStandardFunc, setStandardJac
  public set_parameters
  public wrapped_func, wrapped_jac

  public setNLsysSettings, setNLsysOptions
  public initNLsys_global
  public NLsys_solve
  
 contains
  
  !> @brief Set the stored pointer to a function
  !! passed as argument
  subroutine setParametricFunc(fpointer)
    implicit none
    procedure(ParametricFunction) ::fpointer

    pParFun => fpointer
    m_is_parametric_func = .true.
  end subroutine setParametricFunc

  !> @brief Set the stored pointer to the jacobian
  !! function passed as argument.
  subroutine setParametricJac(jacpointer)
    implicit none
    procedure(ParametricJacFunc) ::jacpointer
    pParJacFun => jacpointer
    m_is_parametric_jac = .true.
  end subroutine setParametricJac
  
  subroutine set_parameters(params)
    implicit none
    real(kind=8) :: params(:)
    if (allocated(m_params)) then
       deallocate(m_params)
    endif
    allocate(m_params(size(params)))
    m_params = params
  end subroutine set_parameters
  
  !> @brief Set the stored pointer to a function
  !! passed as argument
  subroutine setStandardFunc(fpointer)
    implicit none
    procedure(StandardFunction) ::fpointer

    pFun => fpointer
    m_is_func = .true.
  end subroutine setStandardFunc

  !> @brief Set the stored pointer to the jacobian
  !! function passed as argument.
  subroutine setStandardJac(jacpointer)
    implicit none
    procedure(StandardJacFunc) ::jacpointer
    pJacFun => jacpointer
    m_is_jac = .true.
  end subroutine setStandardJac
  
  !> @brief Create a wrapper function needed for the
  !! non linear system solution.
  !!
  !!@note
  !! The user must specify 
  !!@endnote
  subroutine wrapped_func(me,x,f)
    implicit none
    class(nlesolver_type),intent(inout) :: me
    real(8),dimension(:),intent(in)    :: x
    real(8),dimension(:),intent(out)   :: f

    
    if(m_is_parametric_func) then
       call pParFun(x,f,m_params)
    else
       call pFun(x,f)
    end if
  end subroutine wrapped_func

  subroutine wrapped_jac(me,x,g)
    implicit none
    class(nlesolver_type),intent(inout) :: me
    real(8),dimension(:),intent(in)    :: x
    real(8),dimension(:,:),intent(out) :: g

    if(m_is_parametric_jac) then
       call pParJacFun(x,g,m_params)
    else
       call pJacFun(x,g)
    endif
  end subroutine wrapped_jac

  !---------- HERE WE SET THE INFORMATION USEFUL FOR THE NON LINEAR SOLVER ----------- !
  !> @brief Set the settings for our solver.
  !! Basically, the entries in solver_obj are updated
  !!
  !! @param[in] n
  !! @param[out] m
  subroutine setNLsysSettings(n,m,max_iter,tol,verbose)
    implicit none
    integer,intent(in) :: n,m
    integer,intent(inout),optional :: max_iter
    real(kind=8),intent(inout),optional :: tol
    logical,intent(inout),optional :: verbose

    if(.not. present(max_iter)) max_iter = 500
    if(.not. present(tol)) tol = 1.0e-8
    if(.not. present(verbose)) verbose = .false.

    solver_obj%n = n
    solver_obj%m = m
    solver_obj%max_iter = max_iter
    solver_obj%tol = tol
    solver_obj%verbose =verbose
  end subroutine setNLsysSettings

  !> @brief Set the options for the optmizer
  !!
  !!
  subroutine setNLsysOptions( &
    step_mode, &
    use_broyden, &
    n_intervals, &
    fmin_tol, &
    description)
    integer,intent(in)      :: step_mode
    logical,intent(in)      :: use_broyden
    integer,intent(in)      :: n_intervals
    real(kind=8),intent(in) :: fmin_tol
    character(len=:),intent(in),allocatable :: description

    solver_obj%step_mode = step_mode
    solver_obj%use_broyden = use_broyden
    solver_obj%n_intervals = n_intervals
    solver_obj%fmin_tol = fmin_tol
    solver_obj%description = description
  end subroutine setNLsysOptions

  !> @brief Initialize the non linear system
  !! 
  !!@note
  !! It is assumed that we will set the function
  !! and the jacobian of the system previously,
  !! otherwise an error will be raised.
  !!@endnote
  subroutine initNLsys()
    implicit none

    if((m_is_func .eqv. .false.) .and. &
       (m_is_parametric_func).eqv. .false.) then
       error stop 'initNLsys: PLEASE, initialize the nl system!!!'
    endif

    !Print the message info:
    write(*,'(A)') solver_obj%description
    write(*,'(A)')'[initNLsys:]-->Initializing the internal solver object<--'
    call solver_obj%m_solver%initialize( &
    n=solver_obj%n, &
    m=solver_obj%m, &
    max_iter=solver_obj%max_iter, &
    tol=solver_obj%tol, &
    func = wrapped_func, &
    grad = wrapped_jac, &
    step_mode = solver_obj%step_mode, &
    use_broyden = solver_obj%use_broyden, &
    n_intervals = solver_obj%n_intervals, &
    fmin_tol = solver_obj%fmin_tol, &
    verbose = solver_obj%verbose &
    )
    !export_iteration = export <- TODO
    
  end subroutine initNLsys

  !!> @brief Set the global non linear system in our case.
  !!
  !!
  !!
  subroutine initNLsys_global( &
       n, &
       m, &
       step_mode,&
       use_broyden,&
       n_intervals, &
       fmin_tol, &
       description,&
       max_iter, &
       tol, &
       verbose)
    implicit none
    integer, intent(in) :: n,m
    integer,intent(inout),optional :: max_iter
    real(kind=8),intent(inout),optional :: tol
    logical,intent(inout),optional :: verbose
    integer,intent(in) :: step_mode,n_intervals
    logical,intent(in) :: use_broyden
    real(kind=8),intent(in) :: fmin_tol
    character(len=:),allocatable :: description
    !NOTE:
    ! here we should check if the entries are all present!!
    call setNLsysSettings(n,m,max_iter,tol,verbose)
    call setNLsysOptions(step_mode,use_broyden,n_intervals,fmin_tol,description)
    call initNLsys
  end subroutine initNLsys_global

  !> @brief Solve the non linear system
  !! @param[in] x0 The initial condition.
  subroutine NLsys_solve(x0)
    implicit none
    real(kind=8) :: x0(:)
    write(*,'(A,1X,*(E10.3,","))')'Initial condition: ',x0(:)
    call solver_obj%m_solver%solve(x0)
  end subroutine NLsys_solve
  
end module nlesolver_wrapper
