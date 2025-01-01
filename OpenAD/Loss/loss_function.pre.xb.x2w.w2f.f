
      MODULE oad_intrinsics
      use w2f__types
      IMPLICIT NONE
      SAVE
C
C     **** Statements ****
C
      END MODULE

      SUBROUTINE vect_lossfunction(NBATCHES, X, Y_TARGET, THETA, YV)
      use w2f__types
      use oad_intrinsics
      IMPLICIT NONE
C
C     **** Parameters and Result ****
C
      INTEGER(w2f__i4) NBATCHES
      INTENT(IN)  NBATCHES
      REAL(w2f__8) X(1 : NBATCHES)
      REAL(w2f__8) Y_TARGET(1 : NBATCHES)
      TYPE (oadactive) THETA(1 : 2)
      TYPE (oadactive) YV(1 : INT((NBATCHES * 2)))
C
C     **** Local Variables and Functions ****
C
      EXTERNAL fvalues
C
C     **** Top Level Pragmas ****
C
C$OPENAD INDEPENDENT(THETA)
C$OPENAD DEPENDENT(YV)
C
C     **** Statements ****
C
C$OPENAD XXX Template ad_template.f
      CALL fvalues(NBATCHES, X, __deriv__(THETA), __deriv__(YV))
      END SUBROUTINE

      SUBROUTINE fvalues(NBATCHES, X, THETA, FVAL)
      use w2f__types
      use oad_intrinsics
      IMPLICIT NONE
C
C     **** Parameters and Result ****
C
      INTEGER(w2f__i4) NBATCHES
      INTENT(IN)  NBATCHES
      REAL(w2f__8) X(1 : NBATCHES)
      TYPE (oadactive) THETA(1 : 2)
      TYPE (oadactive) FVAL(1 : INT((NBATCHES * 2)))
C
C     **** Local Variables and Functions ****
C
      INTEGER(w2f__i4) I
      REAL(w2f__8) OpenAD_acc_0
      REAL(w2f__8) OpenAD_aux_0
      REAL(w2f__8) OpenAD_lin_0
      REAL(w2f__8) OpenAD_lin_1
C
C     **** Statements ****
C
C$OPENAD XXX Template ad_template.f
      DO I = 1, NBATCHES, 1
        OpenAD_lin_0 = X(I)
        __value__(FVAL(INT(I))) = (X(I) * __value__(THETA(1)))
        CALL sax(OpenAD_lin_0, __deriv__(THETA(1)), __deriv__(FVAL(I)))
      END DO
      DO I = NBATCHES + 1, (NBATCHES * 2), 1
        OpenAD_aux_0 = (__value__(THETA(2)) * 2.0D00)
        OpenAD_lin_1 = X(I - NBATCHES)
        __value__(FVAL(INT(I))) = (X(I - NBATCHES) * OpenAD_aux_0)
        OpenAD_acc_0 = (2.0D00 * OpenAD_lin_1)
        CALL sax(OpenAD_acc_0, __deriv__(THETA(2)), __deriv__(FVAL(I)))
      END DO
      END SUBROUTINE
