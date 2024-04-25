!This module is just a container of the global variables
!that are stored internally useful for the whole program
module global_variables

  real(kind=8),allocatable :: Lmat(:,:) !The interpolation matrix
  real(kind=8),allocatable :: Qs(:,:,:,:,:) !The solution variables

  !cells and number of cells global [old mesh]
  integer :: ncell_glob_old
  integer :: nvert_glob_old

  !cells and number of cells global [new mesh]
  integer :: ncell_glob_new
  integer :: nvert_glob_new

  !cells and number of vertexes [new mesh] stored
  !in each processor.
  integer :: ncell_new
  !integer :: nvert_new
  
  !The mesh cells and vertexes:
  integer,allocatable :: glob_vert_old(:,:) ![Nvertex,3]
  integer,allocatable :: glob_cells_old(:,:) ![Ncells,8]

  !Global vertexes and cells:
  integer,allocatable :: glob_vert_new(:,:) ![nvert_glob_new,3]
  integer,allocatable :: glob_cells_new(:,:) ![ncell_glob_new,8]

  !Local vertexes and cells:
  !real(kind=8),allocatable :: vert_new(:,:) ![nvert_new,3]
  integer,allocatable :: cells_new(:,:) ![ncell_new,8]

  !Conversion:
  integer,allocatable :: ic2icg(:)
  
contains

  subroutine Initialize_old_mesh(gvert,gcells,nv,nc)
    use mpi
    implicit none
    integer,intent(in) :: nv,nc
    integer,intent(in) :: gvert(nv,3) 
    integer,intent(in) :: gcells(nc,8)
    integer :: rank,ierr

    !Debuggining
    !call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    !write(*,'(A,1x,i2)')'Hello from:',rank

    !Allocate and copy:
    allocate(glob_vert_old(nv,3),&
         glob_cells_old(nc,8))
    glob_vert_old(:,:) = gvert(:,:)
    glob_cells_old(:,:) = gcells(:,:)

    ncell_glob_old = nc
    nvert_glob_old = nv
  end subroutine Initialize_old_mesh

  !> @brief Initialize the new mesh and put it
  !! into the rank=0.
  !!
  !!
  !!@note
  !! Note that the global vector of vertexes is
  !! allocated in the all processors in order
  !! to avoid the creation of a local vector of vertexes.
  !!@endnote
  subroutine Initialize_new_mesh(gvert,gcells,nv,nc)
    use mpi
    implicit none
    
    integer,intent(in) :: nv,nc
    integer,intent(in) :: gvert(nv,3) 
    integer,intent(in) :: gcells(nc,8)
    integer :: rank,ierr
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    if(rank==0) then
       write(*,'(A,1x,i1)')'====> Initializing the new mesh on rank ',rank
       allocate(glob_cells_new(nc,8))
       glob_cells_new(:,:) = gcells(:,:)
       ncell_glob_new = nc
       nvert_glob_new = nv
    endif
    !The global vector of the vertexes is allocated
    !in the all processors:
    allocate(glob_vert_new(nv,3))
    glob_vert_new(:,:) = gvert(:,:)
    
  end subroutine Initialize_new_mesh

  !> @brief Partition the mesh and initialize
  !! the arrays contained in each processor.
  !!
  !! We 'partition' in the rank=0: We divide the 
  !! 
  !! @note
  !!   ic2icg is required for converting the local index
  !!   of the cell to the global index.
  !! @endnote
  !!
  subroutine PartitionMesh()
    use mpi
    implicit none
    integer :: iproc,num_procs,rank
    integer :: ierr
    integer,allocatable :: buf_cell_loc(:)
    integer,allocatable :: buf_ic2icg(:)
    integer :: icg_start,icg_end,ic_proc
    integer :: ncell_new_dummy

    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_procs, ierr)
    
    if(rank == 0) then
       ncell_new = ncell_glob_new/num_procs
       
       allocate(ic2icg(ncell_new))
       allocate(cells_new(ncell_new,8))
       allocate(buf_cell_loc(ncell_new*8))
       allocate(buf_ic2icg(ncell_new))
       buf_cell_loc(:) = 0
       ic2icg(:) = 0
       !Loop for each rank:
       
       do iproc=0,num_procs-1
          icg_start = iproc*ncell_new
          icg_end   = (iproc+1)*ncell_new
          do ic_proc=1,ncell_new
             buf_ic2icg(ic_proc) = ic_proc + icg_start
             if(iproc == 0) then
                ic2icg(ic_proc) = ic_proc 
                cells_new(ic_proc,:) = glob_cells_new(ic2icg(ic_proc),:)
             else
                buf_cell_loc(1+8*(ic_proc-1):(8*(ic_proc))) = glob_cells_new(ic2icg(ic_proc),:)
             endif
          enddo
          write(*,*) ic2icg(:)
          if(iproc/=0) then
             call MPI_SEND(ncell_new,1,MPI_INT,iproc,1,&
                  MPI_COMM_WORLD,ierr)
             
             call MPI_SEND(buf_ic2icg(1:ncell_new),ncell_new,MPI_INT,iproc,2,&
                  MPI_COMM_WORLD,ierr)
             
             call MPI_SEND(buf_cell_loc(1:ncell_new*8),ncell_new*8,MPI_INT,iproc,3,&
                  MPI_COMM_WORLD,ierr)
          endif
       enddo!end proc.
    else
       !Receive the number of cells and allocate the buffer:
       call MPI_RECV(ncell_new,1,MPI_INT,0,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
       allocate(ic2icg(ncell_new))
       allocate(buf_cell_loc(ncell_new*8))
       call MPI_RECV(ic2icg(1:ncell_new),ncell_new,MPI_INT,0,2,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
       call MPI_RECV(buf_cell_loc(1:ncell_new*8),ncell_new,MPI_INT,0,3,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)

       !Store the local cells that contains the vertexes indexes:
       allocate(cells_new(ncell_new,8))
       do ic_proc=1,ncell_new
          cells_new(ic_proc,:) = buf_cell_loc(1+8*(ic_proc-1):(8*(ic_proc)))
       enddo
    endif
    
  end subroutine PartitionMesh

  !> [  DEBUGGING  ] 
  !! [  DEBUGGING  ]
  !!@brief Print the new mesh stored in each processor
  subroutine printNewLocalMesh()
    use mpi
    implicit none
    integer :: ic_proc
    integer :: rank,ierr
    !Get the rank:

    call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
    write(*,'("Rank: ",i1)') rank
    do ic_proc = 1,ncell_new
       write(*,'(*(i2,1x))') cells_new(ic_proc,:)
    enddo
    write(*,*) ic2icg(:)
  end subroutine printNewLocalMesh
  
end module global_variables
