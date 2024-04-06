! -------------------------------------------------!
!>
!! A program that serves as test
!! for the FTorch library.
!--------------------------------------------------!


program main
 
  use ftorch
  implicit none
  integer :: j
  integer,parameter :: num_inputs = 1
  type(torch_tensor), dimension(num_inputs) :: in_tensor     !The input
  type(torch_tensor) :: out_tensor                           !The output
  type(torch_module) :: mlp_model
  
  real(kind=4),dimension(1000,216), target  :: input_array
  real(kind=4), dimension(1000,5), target :: out_array
  integer, parameter :: in_dims = 2
  integer :: in_layout(in_dims) = [1,2]

  do j=1,size(input_array,1)
     input_array(j,:) = real(j)*0.00001
  enddo


  in_tensor(1) = torch_tensor_from_array(input_array, in_layout, torch_kCPU) !With my CPU
  out_tensor = torch_tensor_from_array(out_array, in_layout,torch_kCPU)
  
  mlp_model = torch_module_load('mlp.pt')
  call torch_module_forward(mlp_model, in_tensor, num_inputs, out_tensor)

  open(unit=998,file='test.csv',status='new',action='write')

  do j=1,size(out_array,1)
     write(998,'(5(G0.8,","))') out_array(j,:)
     write(998,*)''
  enddo

  close(998)
  

  
  call torch_module_delete(mlp_model)
  call torch_tensor_delete(in_tensor(1))
  call torch_tensor_delete(out_tensor)

  
end program main
