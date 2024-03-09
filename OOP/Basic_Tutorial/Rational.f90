!! This is the module which contains the class
!! 'rational'.

!! This file explains the basic usage of classes 
!! in Fortran90 by using the OOP paradigm
!!
!! The derived type Rational
!! has overloaded operators.
!!

module RationalClass

  type Rational
     private !Make the properties private and
             !and accessible only in the module
     integer :: m_num,m_den;
  end type Rational

  !The overloaded operators
  interface assignment (+)
    module procedure sumRational;
  end interface assignment (+)

  interface operator (-)
     module procedure diffRational;
  end interface operator (-)

  interface operator (**)
     module procedure powRational;
  end interface operator (**)

contains

  !> @brief Basic constructor:
  function Rational_(num,den) result(frac)
    real, intent(in) :: num,den;
    type (Rational) :: frac;
    if(den .le. 0.0) then
       print*,'Error, attempting to have a denominator<0';
       frac = Rational(0.0,0.0); ! Return a basic fraction
    endif
    frac = Rational(num,den); !Explicitly construct
  end function Rational_

  !! Here we place the function/methods
  ! that will be  overloaded
  function sumRational(r1,r2) result(x)
    implicit none;
    type (Rational),intent(in) :: r1,r2;
    type (Rational) :: x;

    x%m_num = r1%m_num*r2%m_den + r2%m_num*r1%m_den;
    x%m_den = r1%m_den*r2%m_den;

  end function sumRational

  function diffRational(r1,r2) result(x)
    implicit none;
    type (Rational),intent(in) :: r1,r2;
    type (Rational) :: x;

    x%m_num = r1%m_num*r2%m_den - r2%m_num*r1%m_den;
    x%m_den = r1%m_den*r2%m_den;

  end function diffRational

  function powRational(r1,exp) result(x)
    implicit none;
    type (Rational),intent(in):: r1;
    type (Rational) :: x;
    real,intent(in) :: exp;

    x%m_num = r1%m_num**exp;
    x%m_den = r1%m_den**exp;

  end function powRational
  
  subroutine printRational(x)
    implicit none;
    type (Rational) :: x;
    print *, x%m_num;
    print *, '---'
    print*, x%m_den;
  end subroutine printRational
  
end module RationalClass
