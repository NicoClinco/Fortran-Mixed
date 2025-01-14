
MODULE oad_intrinsics
use OAD_active
use w2f__types
IMPLICIT NONE
SAVE
!
!     **** Statements ****
!
END MODULE

MODULE sim_info
use OAD_active
use w2f__types
IMPLICIT NONE
SAVE
!
!     **** Global Variables & Derived Type Definitions ****
!
REAL(w2f__8) GAM
REAL(w2f__8) GAM_M1
LOGICAL(w2f__i4) RUN1D
LOGICAL(w2f__i4) RUN2D
!
!     **** Local Variables and Functions ****
!
INTEGER(w2f__i4) NEQ
PARAMETER ( NEQ = 5)
INTEGER(w2f__i4) WP
PARAMETER ( WP = 8)
!
!     **** Initializers ****
!
DATA GAM / 1.39999999999999991118D00 /
DATA GAM_M1 / 4.00000000000000022204D-01 /
DATA RUN1D / .TRUE. /
DATA RUN2D / .FALSE. /
!
!     **** Statements ****
!
END MODULE

SUBROUTINE rusanovflux(QFL, QFR, FN, SF)
use OAD_active
use w2f__types
use oad_intrinsics
use sim_info
IMPLICIT NONE
!
!     **** Parameters and Result ****
!
REAL(w2f__8) QFL(1 : 5)
REAL(w2f__8) QFR(1 : 5)
type(active) :: FN(1:5)
REAL(w2f__8) SF(1 : 3)
!
!     **** Local Variables and Functions ****
!
type(active) :: AM
type(active) :: LAMBDA
REAL(w2f__8) NORMS
type(active) :: PL
type(active) :: PR
EXTERNAL pressure
type(active) :: QFLR(1:10)
type(active) :: UNL
type(active) :: UNR
REAL(w2f__8) OpenAD_acc_0
REAL(w2f__8) OpenAD_acc_1
REAL(w2f__8) OpenAD_acc_10
REAL(w2f__8) OpenAD_acc_11
REAL(w2f__8) OpenAD_acc_12
REAL(w2f__8) OpenAD_acc_13
REAL(w2f__8) OpenAD_acc_14
REAL(w2f__8) OpenAD_acc_15
REAL(w2f__8) OpenAD_acc_16
REAL(w2f__8) OpenAD_acc_2
REAL(w2f__8) OpenAD_acc_3
REAL(w2f__8) OpenAD_acc_4
REAL(w2f__8) OpenAD_acc_5
REAL(w2f__8) OpenAD_acc_6(1 : 5)
REAL(w2f__8) OpenAD_acc_7
REAL(w2f__8) OpenAD_acc_8
REAL(w2f__8) OpenAD_acc_9
REAL(w2f__8) OpenAD_aux_0
REAL(w2f__8) OpenAD_aux_1
REAL(w2f__8) OpenAD_aux_10
REAL(w2f__8) OpenAD_aux_11(1 : 5)
REAL(w2f__8) OpenAD_aux_12(1 : 5)
REAL(w2f__8) OpenAD_aux_13
REAL(w2f__8) OpenAD_aux_14
REAL(w2f__8) OpenAD_aux_15
REAL(w2f__8) OpenAD_aux_16
REAL(w2f__8) OpenAD_aux_17
REAL(w2f__8) OpenAD_aux_18
REAL(w2f__8) OpenAD_aux_2
REAL(w2f__8) OpenAD_aux_3
REAL(w2f__8) OpenAD_aux_4
REAL(w2f__8) OpenAD_aux_5
REAL(w2f__8) OpenAD_aux_6
REAL(w2f__8) OpenAD_aux_7
REAL(w2f__8) OpenAD_aux_8
REAL(w2f__8) OpenAD_aux_9
REAL(w2f__8) OpenAD_lin_0
REAL(w2f__8) OpenAD_lin_1
REAL(w2f__8) OpenAD_lin_10
REAL(w2f__8) OpenAD_lin_11
REAL(w2f__8) OpenAD_lin_12
REAL(w2f__8) OpenAD_lin_13
REAL(w2f__8) OpenAD_lin_14
REAL(w2f__8) OpenAD_lin_15
REAL(w2f__8) OpenAD_lin_16
REAL(w2f__8) OpenAD_lin_17
REAL(w2f__8) OpenAD_lin_18
REAL(w2f__8) OpenAD_lin_19
REAL(w2f__8) OpenAD_lin_2
REAL(w2f__8) OpenAD_lin_20
REAL(w2f__8) OpenAD_lin_21
REAL(w2f__8) OpenAD_lin_22
REAL(w2f__8) OpenAD_lin_23
REAL(w2f__8) OpenAD_lin_24
REAL(w2f__8) OpenAD_lin_25(1 : 5)
REAL(w2f__8) OpenAD_lin_26
REAL(w2f__8) OpenAD_lin_27
REAL(w2f__8) OpenAD_lin_28
REAL(w2f__8) OpenAD_lin_29
REAL(w2f__8) OpenAD_lin_3
REAL(w2f__8) OpenAD_lin_30
REAL(w2f__8) OpenAD_lin_31
REAL(w2f__8) OpenAD_lin_32
REAL(w2f__8) OpenAD_lin_33
REAL(w2f__8) OpenAD_lin_34
REAL(w2f__8) OpenAD_lin_35
REAL(w2f__8) OpenAD_lin_36
REAL(w2f__8) OpenAD_lin_37
REAL(w2f__8) OpenAD_lin_38(1 : 2)
REAL(w2f__8) OpenAD_lin_39
REAL(w2f__8) OpenAD_lin_4
REAL(w2f__8) OpenAD_lin_40(1 : 2)
REAL(w2f__8) OpenAD_lin_41(1 : 2)
REAL(w2f__8) OpenAD_lin_42
REAL(w2f__8) OpenAD_lin_43
REAL(w2f__8) OpenAD_lin_44
REAL(w2f__8) OpenAD_lin_45
REAL(w2f__8) OpenAD_lin_46
REAL(w2f__8) OpenAD_lin_47
REAL(w2f__8) OpenAD_lin_48
REAL(w2f__8) OpenAD_lin_49
REAL(w2f__8) OpenAD_lin_5
REAL(w2f__8) OpenAD_lin_50
REAL(w2f__8) OpenAD_lin_51
REAL(w2f__8) OpenAD_lin_52
REAL(w2f__8) OpenAD_lin_53
REAL(w2f__8) OpenAD_lin_54
REAL(w2f__8) OpenAD_lin_55
REAL(w2f__8) OpenAD_lin_56
REAL(w2f__8) OpenAD_lin_57(1 : 3)
REAL(w2f__8) OpenAD_lin_58
REAL(w2f__8) OpenAD_lin_59(1 : 3)
REAL(w2f__8) OpenAD_lin_6
REAL(w2f__8) OpenAD_lin_60(1 : 3)
REAL(w2f__8) OpenAD_lin_7
REAL(w2f__8) OpenAD_lin_8
REAL(w2f__8) OpenAD_lin_9
type(active) :: OpenAD_prp_0
type(active) :: OpenAD_prp_1
type(active) :: OpenAD_prp_2
type(active) :: OpenAD_prp_3
type(active) :: OpenAD_prp_4(1:5)
type(active) :: OpenAD_prp_5(1:5)
type(active) :: OpenAD_prp_6
type(active) :: OpenAD_prp_7

REAL(w2f__8) EKL,EKR,IEL,IER
INTEGER :: i
!
!     **** Top Level Pragmas ****
!
!$OPENAD INDEPENDENT(QFLR)
!$OPENAD DEPENDENT(FN)
!
!     **** Statements ****
!
QFLR(1:5)%v = QFL(1:5)
QFLR(6:10)%v = QFR(1:5)
CALL zero_deriv(QFLR(1:5))
CALL zero_deriv(QFLR(6:10))

do i=1,10
   QFLR(i)%d(i) = 1.0_WP
enddo

!Left:
QFLR(2)%d(1) = QFL(2)/QFL(1)
QFLR(3)%d(1) = QFL(3)/QFL(1)
QFLR(4)%d(1) = QFL(4)/QFL(1)

EKL=0.0_WP
EKR=0.0_WP
EKL = 0.5_WP*SUM(QFL(2:4)**2)/QFL(1)**2
IEL = QFL(5)/QFL(1)-EKL
QFLR(5)%d(2) = 2.0_WP*(QFL(2)/QFL(1))
QFLR(5)%d(3) = 2.0_WP*(QFL(3)/QFL(1))
QFLR(5)%d(4) = 2.0_WP*(QFL(4)/QFL(1))
QFLR(5)%d(1) = IEL+EKL


!Right

QFLR(7)%d(6) = QFR(2)/QFR(1)
QFLR(8)%d(6) = QFR(3)/QFR(1)
QFLR(9)%d(6) = QFR(4)/QFR(1)
EKR = 0.5_WP*SUM(QFR(2:4)**2)/QFR(1)**2
IER = QFR(5)/QFR(1)-EKR
QFLR(10)%d(7) = 2.0_WP*(QFR(2)/QFR(1))
QFLR(10)%d(8) = 2.0_WP*(QFR(3)/QFR(1))
QFLR(10)%d(9) = 2.0_WP*(QFR(4)/QFR(1))
QFLR(10)%d(6) = IER+EKR




CALL pressure(QFLR(1:5),PL)
CALL pressure(QFLR(6:10),PR)
IF ((PL%v+PR%v).LT.0.0D00) THEN
  FN(1:5)%v = 0.0D00
  CALL zero_deriv(FN(1:5))
ENDIF
IF(RUN1D) THEN
  NORMS = ABS(SF(1))
  OpenAD_aux_0 = (SF(1)*QFLR(2)%v)
  OpenAD_lin_2 = SF(1)
  OpenAD_lin_0 = (INT(1_w2f__i8)/QFLR(1)%v)
  OpenAD_lin_1 = (-(OpenAD_aux_0/(QFLR(1)%v*QFLR(1)%v)))
  UNL%v = (OpenAD_aux_0/QFLR(1)%v)
  OpenAD_aux_1 = (SF(1)*QFLR(7)%v)
  OpenAD_lin_5 = SF(1)
  OpenAD_lin_3 = (INT(1_w2f__i8)/QFLR(6)%v)
  OpenAD_lin_4 = (-(OpenAD_aux_1/(QFLR(6)%v*QFLR(6)%v)))
  UNR%v = (OpenAD_aux_1/QFLR(6)%v)
  OpenAD_lin_6 = UNL%v
  OpenAD_lin_7 = QFLR(1)%v
  OpenAD_lin_8 = UNR%v
  OpenAD_lin_9 = QFLR(6)%v
  FN(1)%v = (QFLR(1)%v*UNL%v+QFLR(6)%v*UNR%v)
  OpenAD_acc_0 = (OpenAD_lin_5*OpenAD_lin_3)
  OpenAD_acc_1 = (OpenAD_lin_2*OpenAD_lin_0)
  CALL sax(OpenAD_lin_1,QFLR(1),UNL)
  CALL saxpy(OpenAD_acc_1,QFLR(2),UNL)
  CALL sax(OpenAD_lin_4,QFLR(6),UNR)
  CALL saxpy(OpenAD_acc_0,QFLR(7),UNR)
  CALL sax(OpenAD_lin_6,QFLR(1),FN(1))
  CALL saxpy(OpenAD_lin_7,UNL,FN(1))
  CALL saxpy(OpenAD_lin_8,QFLR(6),FN(1))
  CALL saxpy(OpenAD_lin_9,UNR,FN(1))
  OpenAD_aux_2 = (PL%v+PR%v)
  OpenAD_lin_10 = UNL%v
  OpenAD_lin_11 = QFLR(2)%v
  OpenAD_lin_12 = UNR%v
  OpenAD_lin_13 = QFLR(7)%v
  OpenAD_lin_14 = SF(1)
  FN(2)%v = (QFLR(2)%v*UNL%v+QFLR(7)%v*UNR%v+SF(1)*OpenAD_aux_2)
  FN(3:4)%v = 0.0D00
  CALL zero_deriv(FN(3:4))
  CALL setderiv(OpenAD_prp_0,PL)
  CALL inc_deriv(OpenAD_prp_0,PR)
  CALL sax(OpenAD_lin_14,OpenAD_prp_0,FN(2))
  CALL saxpy(OpenAD_lin_10,QFLR(2),FN(2))
  CALL saxpy(OpenAD_lin_11,UNL,FN(2))
  CALL saxpy(OpenAD_lin_12,QFLR(7),FN(2))
  CALL saxpy(OpenAD_lin_13,UNR,FN(2))
ELSE
  IF(RUN2D) THEN
    NORMS = SQRT((SF(1)**2)+(SF(2)**2))
    OpenAD_aux_13 = (QFL(3)*SF(2)+SF(1)*QFLR(2)%v)
    OpenAD_lin_29 = SF(1)
    OpenAD_lin_27 = (INT(1_w2f__i8)/QFLR(1)%v)
    OpenAD_lin_28 = (-(OpenAD_aux_13/(QFLR(1)%v*QFLR(1)%v)))
    UNL%v = (OpenAD_aux_13/QFLR(1)%v)
    OpenAD_aux_14 = (SF(1)*QFLR(7)%v+SF(2)*QFLR(8)%v)
    OpenAD_lin_32 = SF(1)
    OpenAD_lin_33 = SF(2)
    OpenAD_lin_30 = (INT(1_w2f__i8)/QFLR(6)%v)
    OpenAD_lin_31 = (-(OpenAD_aux_14/(QFLR(6)%v*QFLR(6)%v)))
    UNR%v = (OpenAD_aux_14/QFLR(6)%v)
    OpenAD_lin_34 = UNL%v
    OpenAD_lin_35 = QFLR(1)%v
    OpenAD_lin_36 = UNR%v
    OpenAD_lin_37 = QFLR(6)%v
    FN(1)%v = (QFLR(1)%v*UNL%v+QFLR(6)%v*UNR%v)
    OpenAD_acc_8 = (OpenAD_lin_33*OpenAD_lin_30)
    OpenAD_acc_9 = (OpenAD_lin_32*OpenAD_lin_30)
    OpenAD_acc_10 = (OpenAD_lin_29*OpenAD_lin_27)
    CALL sax(OpenAD_lin_28,QFLR(1),UNL)
    CALL saxpy(OpenAD_acc_10,QFLR(2),UNL)
    CALL sax(OpenAD_lin_31,QFLR(6),UNR)
    CALL saxpy(OpenAD_acc_8,QFLR(8),UNR)
    CALL saxpy(OpenAD_acc_9,QFLR(7),UNR)
    CALL sax(OpenAD_lin_37,UNR,FN(1))
    CALL saxpy(OpenAD_lin_36,QFLR(6),FN(1))
    CALL saxpy(OpenAD_lin_35,UNL,FN(1))
    CALL saxpy(OpenAD_lin_34,QFLR(1),FN(1))
    OpenAD_aux_15 = (PL%v+PR%v)
    OpenAD_lin_38 = QFL(2:3)
    OpenAD_lin_39 = UNR%v
    OpenAD_lin_40 = QFLR(7:8)%v
    OpenAD_lin_41 = SF(1:2)
    FN(2:3)%v = (QFL(2:3)*UNL%v+QFLR(7:8)%v*UNR%v+SF(1:2)*OpenAD_aux_15)
    FN(4)%v = 0.0D00
    CALL zero_deriv(FN(4))
    CALL setderiv(OpenAD_prp_6,PL)
    CALL inc_deriv(OpenAD_prp_6,PR)
    CALL sax(OpenAD_lin_41,OpenAD_prp_6,FN(2:3))
    CALL saxpy(OpenAD_lin_38,UNL,FN(2:3))
    CALL saxpy(OpenAD_lin_39,QFLR(7:8),FN(2:3))
    CALL saxpy(OpenAD_lin_40,UNR,FN(2:3))
  ELSE
    NORMS = SQRT((SF(1)**2)+(SF(2)**2)+(SF(3)**2))
    OpenAD_aux_16 = (SF(1)*QFLR(2)%v+SF(2)*QFLR(3)%v+SF(3)*QFLR(4)%v)
    OpenAD_lin_44 = SF(1)
    OpenAD_lin_45 = SF(2)
    OpenAD_lin_46 = SF(3)
    OpenAD_lin_42 = (INT(1_w2f__i8)/QFLR(1)%v)
    OpenAD_lin_43 = (-(OpenAD_aux_16/(QFLR(1)%v*QFLR(1)%v)))
    UNL%v = (OpenAD_aux_16/QFLR(1)%v)
    OpenAD_aux_17 = (SF(1)*QFLR(7)%v+SF(2)*QFLR(8)%v+SF(3)*QFLR(9)%v)
    OpenAD_lin_49 = SF(1)
    OpenAD_lin_50 = SF(2)
    OpenAD_lin_51 = SF(3)
    OpenAD_lin_47 = (INT(1_w2f__i8)/QFLR(6)%v)
    OpenAD_lin_48 = (-(OpenAD_aux_17/(QFLR(6)%v*QFLR(6)%v)))
    UNR%v = (OpenAD_aux_17/QFLR(6)%v)
    OpenAD_lin_52 = UNL%v
    OpenAD_lin_53 = QFLR(1)%v
    OpenAD_lin_54 = UNR%v
    OpenAD_lin_55 = QFLR(6)%v
    FN(1)%v = (QFLR(1)%v*UNL%v+QFLR(6)%v*UNR%v)
    OpenAD_acc_11 = (OpenAD_lin_51*OpenAD_lin_47)
    OpenAD_acc_12 = (OpenAD_lin_46*OpenAD_lin_42)
    OpenAD_acc_13 = (OpenAD_lin_50*OpenAD_lin_47)
    OpenAD_acc_14 = (OpenAD_lin_49*OpenAD_lin_47)
    OpenAD_acc_15 = (OpenAD_lin_45*OpenAD_lin_42)
    OpenAD_acc_16 = (OpenAD_lin_44*OpenAD_lin_42)
    CALL sax(OpenAD_lin_43,QFLR(1),UNL)
    CALL saxpy(OpenAD_acc_12,QFLR(4),UNL)
    CALL saxpy(OpenAD_acc_15,QFLR(3),UNL)
    CALL saxpy(OpenAD_acc_16,QFLR(2),UNL)
    CALL sax(OpenAD_lin_48,QFLR(6),UNR)
    CALL saxpy(OpenAD_acc_11,QFLR(9),UNR)
    CALL saxpy(OpenAD_acc_13,QFLR(8),UNR)
    CALL saxpy(OpenAD_acc_14,QFLR(7),UNR)
    CALL sax(OpenAD_lin_55,UNR,FN(1))
    CALL saxpy(OpenAD_lin_54,QFLR(6),FN(1))
    CALL saxpy(OpenAD_lin_53,UNL,FN(1))
    CALL saxpy(OpenAD_lin_52,QFLR(1),FN(1))
    OpenAD_aux_18 = (PL%v+PR%v)
    OpenAD_lin_56 = UNL%v
    OpenAD_lin_57 = QFLR(2:4)%v
    OpenAD_lin_58 = UNR%v
    OpenAD_lin_59 = QFLR(7:9)%v
    OpenAD_lin_60 = SF(1:3)
    FN(2:4)%v = (QFLR(2:4)%v*UNL%v+QFLR(7:9)%v*UNR%v+SF(1:3)*OpenAD_aux_18)
    CALL setderiv(OpenAD_prp_7,PL)
    CALL inc_deriv(OpenAD_prp_7,PR)
    CALL sax(OpenAD_lin_60,OpenAD_prp_7,FN(2:4))
    CALL saxpy(OpenAD_lin_56,QFLR(2:4),FN(2:4))
    CALL saxpy(OpenAD_lin_57,UNL,FN(2:4))
    CALL saxpy(OpenAD_lin_58,QFLR(7:9),FN(2:4))
    CALL saxpy(OpenAD_lin_59,UNR,FN(2:4))
  ENDIF
ENDIF
OpenAD_aux_3 = (QFL(5)+PL%v)
OpenAD_aux_4 = (QFR(5)+PR%v)
OpenAD_lin_15 = OpenAD_aux_3
OpenAD_lin_16 = UNL%v
OpenAD_lin_17 = OpenAD_aux_4
OpenAD_lin_18 = UNR%v
FN(5)%v = (UNL%v*OpenAD_aux_3+UNR%v*OpenAD_aux_4)
OpenAD_aux_8 = (PL%v+PR%v)
OpenAD_aux_6 = (GAM*OpenAD_aux_8)
OpenAD_aux_7 = (QFLR(1)%v+QFLR(6)%v)
OpenAD_aux_5 = SQRT(OpenAD_aux_6/OpenAD_aux_7)
OpenAD_lin_22 = GAM
OpenAD_lin_20 = (INT(1_w2f__i8)/OpenAD_aux_7)
OpenAD_lin_21 = (-(OpenAD_aux_6/(OpenAD_aux_7*OpenAD_aux_7)))
OpenAD_lin_19 = (5.0D-01/OpenAD_aux_5)
AM%v = OpenAD_aux_5
OpenAD_aux_10 = (UNL%v+UNR%v)
OpenAD_aux_9 = ABS(OpenAD_aux_10)
OpenAD_lin_23 = NORMS
OpenAD_lin_24 = SIGN(1.0D00,OpenAD_aux_10)
LAMBDA%v = (AM%v*NORMS+OpenAD_aux_9*5.0D-01)
OpenAD_acc_2 = (OpenAD_lin_19*OpenAD_lin_23)
OpenAD_acc_3 = (OpenAD_lin_21*OpenAD_acc_2)
OpenAD_acc_4 = (OpenAD_lin_22*OpenAD_lin_20*OpenAD_acc_2)
OpenAD_acc_5 = (OpenAD_lin_24*5.0D-01)
CALL setderiv(OpenAD_prp_1,PL)
CALL inc_deriv(OpenAD_prp_1,PR)
CALL setderiv(OpenAD_prp_2,QFLR(1))
CALL inc_deriv(OpenAD_prp_2,QFLR(6))
CALL setderiv(OpenAD_prp_3,UNL)
CALL inc_deriv(OpenAD_prp_3,UNR)
CALL sax(OpenAD_lin_15,UNL,FN(5))
CALL saxpy(OpenAD_lin_16,PL,FN(5))
CALL saxpy(OpenAD_lin_17,UNR,FN(5))
CALL saxpy(OpenAD_lin_18,PR,FN(5))
CALL sax(OpenAD_acc_3,OpenAD_prp_2,LAMBDA)
CALL saxpy(OpenAD_acc_4,OpenAD_prp_1,LAMBDA)
CALL saxpy(OpenAD_acc_5,OpenAD_prp_3,LAMBDA)
OpenAD_aux_12 = (QFLR(6:10)%v-QFL(1:5))
OpenAD_aux_11 = (FN(1:5)%v-LAMBDA%v*OpenAD_aux_12)
OpenAD_lin_25 = OpenAD_aux_12
OpenAD_lin_26 = LAMBDA%v
FN(1:5)%v = (OpenAD_aux_11*5.0D-01)
OpenAD_acc_6 = (OpenAD_lin_25*INT((-1_w2f__i8)))
OpenAD_acc_7 = (OpenAD_lin_26*INT((-1_w2f__i8)))
CALL setderiv(OpenAD_prp_4,FN(1:5))
CALL setderiv(OpenAD_prp_5,OpenAD_prp_4)
CALL saxpy(OpenAD_acc_6,LAMBDA,OpenAD_prp_5)
CALL saxpy(OpenAD_acc_7,QFLR(6:10),OpenAD_prp_5)
CALL sax(5.0D-01,OpenAD_prp_5,FN(1:5))
END SUBROUTINE

SUBROUTINE pressure(QV, P)
use OAD_active
use w2f__types
use oad_intrinsics
use sim_info
IMPLICIT NONE
!
!     **** Parameters and Result ****
!
type(active) :: QV(1:5)
type(active) :: P
!
!     **** Local Variables and Functions ****
!
type(active) :: RHOUSQ
REAL(w2f__8) OpenAD_acc_17
REAL(w2f__8) OpenAD_acc_18
REAL(w2f__8) OpenAD_acc_19
REAL(w2f__8) OpenAD_acc_20
REAL(w2f__8) OpenAD_acc_21
REAL(w2f__8) OpenAD_acc_22
REAL(w2f__8) OpenAD_aux_19
REAL(w2f__8) OpenAD_aux_20
REAL(w2f__8) OpenAD_aux_21
REAL(w2f__8) OpenAD_aux_22
REAL(w2f__8) OpenAD_lin_61
REAL(w2f__8) OpenAD_lin_62
REAL(w2f__8) OpenAD_lin_63
REAL(w2f__8) OpenAD_lin_64
REAL(w2f__8) OpenAD_lin_65
REAL(w2f__8) OpenAD_lin_66
REAL(w2f__8) OpenAD_lin_67
REAL(w2f__8) OpenAD_lin_68
REAL(w2f__8) OpenAD_lin_69
REAL(w2f__8) OpenAD_lin_70
REAL(w2f__8) OpenAD_lin_71
REAL(w2f__8) OpenAD_lin_72
REAL(w2f__8) OpenAD_lin_73
type(active) :: OpenAD_prp_8
!
!     **** Statements ****
!
IF(RUN1D) THEN
  OpenAD_aux_19 = (QV(2)%v**2)
  OpenAD_lin_63 = (2*(QV(2)%v**(2-INT(1_w2f__i8))))
  OpenAD_lin_61 = (INT(1_w2f__i8)/QV(1)%v)
  OpenAD_lin_62 = (-(OpenAD_aux_19/(QV(1)%v*QV(1)%v)))
  RHOUSQ%v = (OpenAD_aux_19/QV(1)%v)
  OpenAD_acc_17 = (OpenAD_lin_63*OpenAD_lin_61)
  CALL sax(OpenAD_lin_62,QV(1),RHOUSQ)
  CALL saxpy(OpenAD_acc_17,QV(2),RHOUSQ)
ELSE
  IF(RUN2D) THEN
    OpenAD_aux_21 = ((QV(2)%v**2)+(QV(3)%v**2))
    OpenAD_lin_67 = (2*(QV(2)%v**(2-INT(1_w2f__i8))))
    OpenAD_lin_68 = (2*(QV(3)%v**(2-INT(1_w2f__i8))))
    OpenAD_lin_65 = (INT(1_w2f__i8)/QV(1)%v)
    OpenAD_lin_66 = (-(OpenAD_aux_21/(QV(1)%v*QV(1)%v)))
    RHOUSQ%v = (OpenAD_aux_21/QV(1)%v)
    OpenAD_acc_18 = (OpenAD_lin_68*OpenAD_lin_65)
    OpenAD_acc_19 = (OpenAD_lin_67*OpenAD_lin_65)
    CALL sax(OpenAD_lin_66,QV(1),RHOUSQ)
    CALL saxpy(OpenAD_acc_18,QV(3),RHOUSQ)
    CALL saxpy(OpenAD_acc_19,QV(2),RHOUSQ)
  ELSE
    OpenAD_aux_22 = ((QV(2)%v**2)+(QV(3)%v**2)+(QV(4)%v**2))
    OpenAD_lin_71 = (2*(QV(2)%v**(2-INT(1_w2f__i8))))
    OpenAD_lin_72 = (2*(QV(3)%v**(2-INT(1_w2f__i8))))
    OpenAD_lin_73 = (2*(QV(4)%v**(2-INT(1_w2f__i8))))
    OpenAD_lin_69 = (INT(1_w2f__i8)/QV(1)%v)
    OpenAD_lin_70 = (-(OpenAD_aux_22/(QV(1)%v*QV(1)%v)))
    RHOUSQ%v = (OpenAD_aux_22/QV(1)%v)
    OpenAD_acc_20 = (OpenAD_lin_73*OpenAD_lin_69)
    OpenAD_acc_21 = (OpenAD_lin_72*OpenAD_lin_69)
    OpenAD_acc_22 = (OpenAD_lin_71*OpenAD_lin_69)
    CALL sax(OpenAD_lin_70,QV(1),RHOUSQ)
    CALL saxpy(OpenAD_acc_20,QV(4),RHOUSQ)
    CALL saxpy(OpenAD_acc_21,QV(3),RHOUSQ)
    CALL saxpy(OpenAD_acc_22,QV(2),RHOUSQ)
  ENDIF
ENDIF
OpenAD_aux_20 = (QV(5)%v-RHOUSQ%v*5.0D-01)
OpenAD_lin_64 = GAM_M1
P%v = (GAM_M1*OpenAD_aux_20)
CALL setderiv(OpenAD_prp_8,QV(5))
CALL saxpy(-5.0D-01,RHOUSQ,OpenAD_prp_8)
CALL sax(OpenAD_lin_64,OpenAD_prp_8,P)
END SUBROUTINE
