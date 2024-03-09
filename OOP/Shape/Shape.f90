module Class_Shape
  implicit none;
  public :: ShapeM;
  
  !> @brief The abstract base type
  !! @ param[in] m_name : The name of the shape
  type ShapeM
     private
     character (10) :: m_name;
  end type ShapeM

contains
  !> @brief Public constructor
  function Shape_(name) result (s)
    type (ShapeM) :: s;
    character (10) :: name;
    s = ShapeM(name);
  end function Shape_

  subroutine printShapeType(s)
    type (ShapeM) :: s;
    print*, 'Shape type: ',s%m_name;
  end subroutine printShapeType
  
end module Class_Shape
