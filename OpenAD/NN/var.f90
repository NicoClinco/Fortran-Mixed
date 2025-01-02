!>@brief The module that contains various parameters to
!!optimize
!!
MODULE var

  INTEGER,PARAMETER :: n_in=2
  INTEGER,PARAMETER :: n_out=2


  REAL(KIND=8) :: W(n_out,n_in)
  REAL(KIND=8) :: b(n_out)
  !The cost-function (that will be transformed)
  REAL(KIND=8) :: cost_fun

CONTAINS
  !>@brief The subroutine that serves for initialize the weights
  !!
  !!@param[in] nin
  !!@param[in] nout
   SUBROUTINE InitWeights(nin,nout)
     IMPLICIT NONE
     INTEGER,INTENT(IN) :: nin,nout

  !   n_in=nin;n_out=nout
  !   ALLOCATE(W(nout,nin),b(nout) )
     W(1:n_out,1:n_in) = 0.0_8;
     b(1:n_out) = 0.0_8
   END SUBROUTINE InitWeights
END MODULE var
