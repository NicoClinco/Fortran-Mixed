
MODULE oad_intrinsics
use OAD_active
use w2f__types
IMPLICIT NONE
SAVE
!
!     **** Statements ****
!
END MODULE

SUBROUTINE vect_lossfunction(NBATCHES, X, Y_TARGET, THETA, YV)
use OAD_active
use w2f__types
use oad_intrinsics
IMPLICIT NONE
!
!     **** Parameters and Result ****
!
INTEGER(w2f__i4) NBATCHES
INTENT(IN) NBATCHES
REAL(w2f__8) X(1 : NBATCHES)
REAL(w2f__8) Y_TARGET(1 : NBATCHES)
type(active) :: THETA(1:2)
type(active) :: YV(1:INT((NBATCHES*2)))
!
!     **** Local Variables and Functions ****
!
EXTERNAL fvalues
!
!     **** Top Level Pragmas ****
!
!$OPENAD INDEPENDENT(THETA)
!$OPENAD DEPENDENT(YV)
!
!     **** Statements ****
!
!$OPENAD XXX Template ad_template.f
CALL fvalues(NBATCHES,X,THETA,YV)
END SUBROUTINE

SUBROUTINE fvalues(NBATCHES, X, THETA, FVAL)
use OAD_active
use w2f__types
use oad_intrinsics
IMPLICIT NONE
!
!     **** Parameters and Result ****
!
INTEGER(w2f__i4) NBATCHES
INTENT(IN) NBATCHES
REAL(w2f__8) X(1 : NBATCHES)
type(active) :: THETA(1:2)
type(active) :: FVAL(1:INT((NBATCHES*2)))
!
!     **** Local Variables and Functions ****
!
INTEGER(w2f__i4) I
REAL(w2f__8) OpenAD_acc_0
REAL(w2f__8) OpenAD_aux_0
REAL(w2f__8) OpenAD_lin_0
REAL(w2f__8) OpenAD_lin_1
!
!     **** Statements ****
!
!$OPENAD XXX Template ad_template.f
DO I = 1, NBATCHES, 1
  OpenAD_lin_0 = X(I)
  FVAL(INT(I))%v = (X(I)*THETA(1)%v)
  CALL sax(OpenAD_lin_0,THETA(1),FVAL(I))
END DO
DO I = NBATCHES+1,(NBATCHES*2),1
  OpenAD_aux_0 = (THETA(2)%v*2.0D00)
  OpenAD_lin_1 = X(I-NBATCHES)
  FVAL(INT(I))%v = (X(I-NBATCHES)*OpenAD_aux_0)
  OpenAD_acc_0 = (2.0D00*OpenAD_lin_1)
  CALL sax(OpenAD_acc_0,THETA(2),FVAL(I))
END DO
END SUBROUTINE
