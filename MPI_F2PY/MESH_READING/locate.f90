!> Given the coordinate x, we return the
!! index ic_old that contains the x coordinate
!! of the cell.
!!
!!
!!@note
!! This subroutines is run by every processor and basically it loops
!! all over the cells of the previous mesh and get the coordinates of the
!! cells.
!!@endnote
subroutine locate_from_old_mesh(x,ic_old)
  use global_variables, only: ncell_glob_old, nvert_glob_old, &
       glob_cells_old, glob_vert_old
  
  integer,intent(out) :: ic_old
  real(kind=8),intent(in) :: x(3) !The coordinate to localize in the
  ! old mesh.
  real(kind=8) :: xvertexes_cell(3,8)
  real(kind=8) :: xmin,ymin,zmin,ymin,ymax,zmax


  ic_old = -1
  !Loop all over the cells of the old mesh
  do i=1,ncell_glob_old
     !Get the coordinates of the vertexes: [x,y,z]
     do j=1,3
        xvertexes_cell(j,:) = glob_vert_old(glob_cells_old(i,:),j)
     enddo
     xmin = min(xvertexes_cell(1,:))
     ymin = min(xvertexes_cell(2,:))
     zmin = min(xvertexes_cell(3,:))
     xmax = max(xvertexes_cell(1,:))
     ymax = max(xvertexes_cell(2,:))
     zmax = max(xvertexes_cell(3,:))
     if((x(1) .lt. xmax .and. x(1) .gt. xmin) .and. &
          (x(2) .lt. ymax .and. x(2) .gt. ymin) .and. &
          (x(3) .lt. xmax .and. x(3) .gt. zmin) ) then
        ic_old = i
        return
     endif
  enddo

  !Check:
  if(ic_old == -1) then
     write(*,'(A,1x,*(F8.3))') 'Error, unable to locate:'x(:)
     stop
  endif

  

end subroutine locate_from_old_mesh
