!NOTE THAT IN FORTRAN THE 'MODULE' IS
!THE CLASS ITSELF

module Class_Rectangle
  use Class_Shape
  implicit none;
  public :: Rectangle
  
  !The rectangle type inherits
  !from the shape class. In fortran
  !90 the 'is-a' is achieved by
  !storing the base class itself:
  type Rectangle
     private
     real :: x,y;
     type (ShapeM) :: shape;
     
  end type Rectangle


contains
  !>@brief Basic Constructor
  function Rectangle_(x,y) result(rect)
    type (Rectangle) :: rect;
    type (ShapeM) :: base_shape;
    real,intent(in) :: x,y;
    base_shape = Shape_('Rectangle');
    rect = Rectangle(x,y,base_shape);
    return;
  end function Rectangle_

  !>@brief Compute the area:
  function Area(rect)
    type (Rectangle) :: rect;
    real :: Area;
    Area = (rect%x)*(rect%y);
  end function Area
end module Class_Rectangle
