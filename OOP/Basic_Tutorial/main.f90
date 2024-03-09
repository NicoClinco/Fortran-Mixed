
program main
  use RationalClass
  implicit none;

  type (Rational) :: r1,r2,r3
  real :: x,y;
  !Construct:

  x = 0.4;y=0.2;
  r1 = Rational_(x,y);
  r2 = Rational_(x,y);

  r3 = r1+r2;

  call printRational(r3);
  
end program main
