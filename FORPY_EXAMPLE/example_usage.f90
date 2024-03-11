
program mymodule_example
  use forpy_mod
  implicit none

  integer :: ierror
  
  
  type(module_py) :: readDataModule !Declaration of read-data module
  type(object) :: return_value
  type(tuple) :: args
  type(list) :: paths
  
  real,allocatable :: fortran_array(:) !fortran array that must be passed to python.
  type(ndarray) :: np_array !The numpy array-object
  type(ndarray) :: array_ret

  allocate(fortran_array(10))

  fortran_array = 44.0 ! Set the fortran array to 44.0
  !Copy to numpy array (or move)
  
  
  ierror = forpy_initialize()

  
  ! Instead of setting the environment variable PYTHONPATH,
  ! we can add the current directory "." to sys.path
  ierror = get_sys_path(paths)
  ierror = paths%append("./python_modules")
  
  ierror = import_py(readDataModule,"readModule")
  ierror = ndarray_create(np_array, fortran_array)
  !ierror = print_py(np_array)
  ! Python: 
  !return_value = readDataModule.multiplyArrayBy2(np_array)
  ierror = tuple_create(args, 0)
  ierror = args%setitem(0,np_array)
  ierror = call_py(return_value,readDataModule, "multiplyArrayBy2",args)
  
  !ierror = ndarray_create_empty(array_ret,10, dtype="float32")
  ierror = cast(array_ret,return_value)
  
  !ierror = print_py(return_value)

  ! For call_py, args and kwargs are optional
  ! use call_py_noret to ignore the return value
  ! E. g.:
  ! ierror = call_py_noret(mymodule, "print_args")

  
  !call np_array_ret%destroy
  call np_array%destroy
  call readDataModule%destroy
  call return_value%destroy
  call paths%destroy
  
  call forpy_finalize

end program


