! -------------------------------------------------!
!>
!! A program that serves as test
!! for the FTorch library.
!!
!! In this example, we have a model which accept only one tensor
!! in input.
!! 
!!@note
!!  The user must be consistent with the model defined in torch
!!  during the training. In our case, we have a model that uses
!!  float32 numbers, thus our inputs-outputs must be kind=4 in
!!  fortran.
!!@endnote
!!--------------------------------------------------!
program main
 
  use ftorch
  implicit none
  integer :: j
  integer,parameter :: num_inputs = 1
  type(torch_tensor), dimension(num_inputs) :: in_tensor     
  type(torch_tensor) :: out_tensor                           
  type(torch_module) :: mlp_model

  !Define the inputs with [batches,num_features]:
  !real(kind=4),dimension(79507,216), target  :: input_array
  !real(kind=4), dimension(79507,6), target :: out_array

  integer, parameter :: in_dims = 2
  integer :: in_layout(in_dims) = [1,2]

  real(kind=4),allocatable,target :: input_array(:,:)
  real(kind=4),allocatable,target :: out_array(:,:)

  integer :: num_batches = 80000
  integer :: num_features_in = 216
  integer :: num_features_out = 6

  allocate( input_array(num_batches,num_features_in))
  allocate( out_array(num_batches,num_features_out))

  !Create dummy inputs (to test the neural network here)
  do j=1,size(input_array,1)
     input_array(j,:) = real(1.05,4)
  enddo

  
  in_tensor(1) = torch_tensor_from_array(input_array, in_layout, torch_kCPU) 
  out_tensor = torch_tensor_from_array(out_array, in_layout,torch_kCPU)
  
  mlp_model = torch_module_load('21p5_scripted.pt')
  call torch_module_forward(mlp_model, in_tensor, num_inputs, out_tensor)

  open(unit=998,file='test.csv',status='new',action='write')

  do j=1,size(out_array,1)
     write(998,'(6(G0.8,","))') out_array(j,:)
     !write(998,*)''
  enddo

  close(998)
  

  
  call torch_module_delete(mlp_model)
  call torch_tensor_delete(in_tensor(1))
  call torch_tensor_delete(out_tensor)

  
end program main
