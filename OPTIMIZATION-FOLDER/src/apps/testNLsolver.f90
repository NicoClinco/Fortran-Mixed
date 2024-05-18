
!>
!  Test of a small, square (`n=m`) problem.
program nlesolver_test_1

  use nlesolver_wrapper!, only : pParFun, pParJacFun
  use nlesolver_module, wp => nlesolver_rk

  implicit none

  !Settings of the system:
  integer,parameter :: n = 2
  integer,parameter :: m = 2
  integer,parameter :: max_iter = 100
  real(wp),parameter :: tol = 1.0e-8_wp
  logical,parameter :: verbose = .false.

  integer :: max_iter_ = 200
  real(wp) :: tol_ = 1.0e-8
  logical :: verbose_= .false.
  !Solver parameters:
  type(nlesolver_type) :: solver
  real(wp) :: alpha
  logical :: use_broyden
  integer :: step_mode
  integer :: n_intervals
  integer :: istat !! Integer status code.
  character(len=:),allocatable :: message  !! Text status message
  real(wp),dimension(n) :: x
  integer :: f_evals
  integer :: i
  character(len=:),allocatable :: description
  real(wp) :: fmin_tol
  real(wp) :: parameters_vector(3) !Set the parameters vector.
  !We can use the alpha_v as a vector of parameters:
  !real(wp) :: alpha_v(2)
  !alpha_v(:) = 1.0

    
  fmin_tol = 1.0e-2_wp ! don't need a tight tol for this
  n_intervals = 2
  alpha = 1.0_wp
  print '(A)','------->TEST THE MODULE CAPABILITIES<------'

  !Here we test the module that allows us to create parametric
  !functions for solving the non-linear system:
  
  !call setStandardFunc(func_org)
  !call setStandardJac(jac_func_org)

  !--- Parametric version here ---!
  parameters_vector(:)= 1.0
  call set_parameters(parameters_vector)
  call setParametricJac(jac_func_p)
  call setParametricFunc(func_p)
  !Initialize the internal non linear solver:
  step_mode = 1
  use_broyden = .false.
  description= 'Optimizer: CONSTANT-ALPHA'
  call initNLsys_global( n, &
       m, &
       step_mode = step_mode,&
       use_broyden = use_broyden,&
       n_intervals = n_intervals, &
       fmin_tol = fmin_tol, &
       description=description, &
       tol = tol_, &
       verbose = verbose_)!, &
       !max_iter = max_iter_)
       !tol = tol_)
  x = [1.0_wp, 2.0_wp]
  !Solve the nonlinear system
  call NLsys_solve(x)
  write(*,*)'Final root solution',x
  
  print '(A)','------->TEST THE MODULE<------'

  write(*,*) ''
  write(*,*) '***********************'
  write(*,*) '* nlesolver_test_1    *'
  write(*,*) '***********************'
  write(*,*) ''
  do i = 1, 8
     select case (i)
     case(1)
        step_mode = 1
        use_broyden = .false.
        f_evals = 0
        description = 'Constant alpha'
     case(2)
        step_mode = 1
        use_broyden = .true.
        f_evals = 0
        description = 'Constant alpha + broyden'
     case(3)
        step_mode = 2
        use_broyden = .false.
        f_evals = 0
        description = 'Backtracking line search'
     case(4)
        step_mode = 2
        use_broyden = .true.
        f_evals = 0
        description = 'Backtracking line search + broyden'
     case(5)
        step_mode = 3
        use_broyden = .false.
        f_evals = 0
        description = 'Exact line search'
     case(6)
        step_mode = 3
        use_broyden = .true.
        f_evals = 0
        description = 'Exact line search + broyden'
     case(7)
        step_mode = 4
        use_broyden = .false.
        f_evals = 0
        description = 'Fixed point search'
     case(8)
        step_mode = 4
        use_broyden = .true.
        f_evals = 0
        description = 'Fixed point search + broyden'
     case default
        error stop 'invalid case'
     end select

     write(*,*) '-------------------------------------------------------'
     write(*,'(A,I3,A,A)') 'Case ', i, ' : ', description
     write(*,*) ''

     call initNLsys_global( n, &
          m, &
          !max_iter = max_iter_, &
          tol = tol_, &
          step_mode = step_mode,&
          use_broyden = use_broyden,&
          n_intervals = n_intervals, &
          fmin_tol = fmin_tol, &
          verbose = verbose_, &
          description=description)
     
     
     call NLsys_status(istat,message)
     write(*,'(I3,1X,A)') istat, message
     if(istat /= 0) error stop
     
     !Solve the nonlinear system
     x = [1.0_wp, 2.0_wp]
     call NLsys_solve(x)
     write(*,'(*(E10.3,1X))')x 
     !The first step is the initialization:
     !call solver%initialize( n = n, &
     !     m = m, &
     !     max_iter = max_iter, &
     !     tol = tol, &
     !     func = wrapped_func, &
     !     grad = wrapped_jac, &
     !     step_mode = step_mode,&
     !     use_broyden = use_broyden,&
     !     export_iteration = export,&
     !     n_intervals = n_intervals, &
     !     fmin_tol = fmin_tol, &
     !     verbose = verbose)
     !call solver%status(istat, message)
     !write(*,'(I3,1X,A)') istat, message
     !if (istat /= 0) error stop
     
     
     !x = [1.0_wp, 2.0_wp]
     !call solver%solve(x)
     !write(*,*)'---',x

     !call solver%status(istat, message)
     !write(*,'(I3,1X,A)') istat, message
     !write(*,*) ''

  end do

  
contains

  !> @brief The original function [non-parametric]
  subroutine func_org(x,f)
    implicit none
    real(wp),dimension(:),intent(in)    :: x
    real(wp),dimension(:),intent(inout)   :: f
    f(1) = x(1)**2 + x(2) - 0.1_wp
    f(2) = x(2) + 0.2_wp
  end subroutine func_org

  !> @brief The jacobian function [non parametric]
  subroutine jac_func_org(x,g)
    implicit none
    real(wp),dimension(:),intent(in)    :: x
    real(wp),dimension(:,:),intent(inout) :: g

    f_evals = f_evals + 2   ! to approximate forward diff derivatives
    g(1,1) = 2.0_wp * x(1)  !df(1)/dx
    g(2,1) = 0.0_wp         !df(2)/dx

    g(1,2) = 1.0_wp         !df(1)/dy
    g(2,2) = 1.0_wp         !df(2)/dy
    
  end subroutine jac_func_org

  !> @brief The original function [parametric]
  subroutine func_p(x,f,p)
    implicit none
    real(wp),dimension(:),intent(in)    :: x
    real(wp),dimension(:),intent(inout)   :: f
    real(wp),dimension(:),intent(in) :: p
    f(1) = x(1)**2 + p(1)*x(2) - 0.1_wp
    f(2) = p(2)*x(2) + 0.2_wp
    print*,'prova'
  end subroutine func_p

  !> @brief The jacobian function [parametric]
  !!
  !!@note
  !! The user can set the parameter array on the right
  !!@endnote 
  subroutine jac_func_p(x,g,p)
    implicit none
    real(wp),dimension(:),intent(in)    :: x
    real(wp),dimension(:,:),intent(inout) :: g
    real(wp),dimension(:),intent(in) :: p
    
    f_evals = f_evals + 2   ! to approximate forward diff derivatives
    g(1,1) = 2.0_wp * x(1)  !df(1)/dx
    g(2,1) = 0.0_wp         !df(2)/dx

    g(1,2) = p(1)        !df(1)/dy
    g(2,2) = p(2)        !df(2)/dy

  end subroutine jac_func_p
  
  
  subroutine func_wrapper(me,x,f)
    implicit none
    class(nlesolver_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)    :: x
    real(wp),dimension(:),intent(out)   :: f

    !alpha_v is a global variable
    !(that can be contained in a module)
    call func_org(x,f)
    f_evals = f_evals + 1
   
  end subroutine func_wrapper
  
  subroutine func(me,x,f)
    !! compute the function
    implicit none
    class(nlesolver_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)    :: x
    real(wp),dimension(:),intent(out)   :: f

    f_evals = f_evals + 1
    
    
    f(1) = x(1)**2 + x(2) - 0.1_wp
    f(2) = x(2) + 0.2_wp

    ! root is 5.477226E-01  -2.000000E-01

  end subroutine func

  !> @brief Compute the gradient
  !!
  !!@param[inout] me: The non-linear solver object
  !!@param [in] x: The specific point
  !!@param [inout] g: The gradient we want to specify
  subroutine grad(me,x,g)
    !! compute the gradient of the function (Jacobian):
    implicit none
    class(nlesolver_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)    :: x
    real(wp),dimension(:,:),intent(out) :: g

    f_evals = f_evals + 2   ! to approximate forward diff derivatives
    
    g(1,1) = 2.0_wp * x(1)  !df(1)/dx
    g(2,1) = 0.0_wp         !df(2)/dx

    g(1,2) = 1.0_wp         !df(1)/dy
    g(2,2) = 1.0_wp         !df(2)/dy

  end subroutine grad

  !!> @brief Print the value of an iteration
  !!
  subroutine export(me,x,f,iter)
    
    implicit none
    class(nlesolver_type),intent(inout) :: me
    real(wp),dimension(:),intent(in)    :: x
    real(wp),dimension(:),intent(in)    :: f
    integer,intent(in)                  :: iter !! iteration number

    write(*,'(1P,I3,1X,A,I3,A,*(E15.6))') iter, '(',f_evals,')', x, norm2(f)

  end subroutine export


  subroutine print_status
    implicit none
    write(*,*) f_evals
  end subroutine print_status
  
end program nlesolver_test_1

