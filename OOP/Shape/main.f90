  ! Example of the usage of the classes above:

program main
  use Class_Shape;
  use Class_Rectangle;
  implicit none
 
  type(ShapeM) :: rect_shape;
  type(Rectangle) :: rect_obj;
  !rect_shape = ShapeM('bho'); PRIVATE!
  rect_shape = Shape_('Rectangular')
  call printShapeType(rect_shape)

  Rect_obj = Rectangle_(2.0,3.0);
  print *,'Area rectangle: ',Area(Rect_obj);
end program main
