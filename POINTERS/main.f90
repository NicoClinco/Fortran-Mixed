!---------------------------------------!
!  This file contains an example        !
!  of usage for pointers                !
!                                       !
!---------------------------------------!

MODULE LIN_ALGEBRA
  IMPLICIT NONE;

  !This data structure encapsulate 
  !pointers:
  TYPE Row
     REAL,DIMENSION(:), POINTER :: row
  END TYPE Row

  TYPE Triangular_matrix
     !Dimension
     INTEGER :: m_nrows; 
     !Array of pointers:
     TYPE(Row),DIMENSION(:),ALLOCATABLE :: m_rows;
     !TO DO: INSERT THE SUBROUTINE TO GET THE ELEMENTS:
     
  END TYPE Triangular_matrix

  !> @brief Dense_Matrix_class
  TYPE Dense_Matrix
     REAL,DIMENSION(:,:),ALLOCATABLE :: m_data;
  END TYPE Dense_Matrix

CONTAINS
  SUBROUTINE INITIALIZE_MATRIX_ONES(TM,nrows)
    IMPLICIT NONE;
    TYPE(Triangular_matrix),INTENT(INOUT) :: TM;
    INTEGER,INTENT(IN) :: nrows;
    INTEGER :: i;
    !CHARACTER (*),INTENT(IN) :: matrix_type; 
    IF(nrows .le. 0 ) THEN
       PRINT *,'ERROR, TRIED TO CREATE AN EMPTY MATRIX';
       RETURN;
    ELSE
       TM%m_nrows = nrows;
       ALLOCATE(TM%m_rows(nrows));

       DO i = 1,TM%m_nrows
          !Allocate every single row:
          ALLOCATE(TM%m_rows(i)%row(1:i));
          TM%m_rows(i)%row(1:i) = 1.0; !Initialize to 1 by default;
       END DO
    END IF
  END SUBROUTINE INITIALIZE_MATRIX_ONES

  SUBROUTINE DENSE_MUL(DM,DM1,DMR)
    implicit none;
    integer :: i,j;
    TYPE(Dense_Matrix),INTENT(IN) :: DM,DM1;
    TYPE(Dense_Matrix),INTENT(OUT):: DMR;
    
    !REAL,DIMENSION(size(DM%m_data,1),size(DM%m_data,2)),TARGET:: l_data;
    REAL,DIMENSION(:,:),ALLOCATABLE,TARGET:: l_data;
    REAL,DIMENSION(:,:),ALLOCATABLE,TARGET:: r_data;
    !REAL,DIMENSION(size(DM1%m_data,1),size(DM1%m_data,2)),TARGET:: r_data;
    
    REAL,DIMENSION(:,:),POINTER :: left;
    REAL,DIMENSION(:,:),POINTER :: right;

    l_data = DM%m_data;
    r_data = DM1%m_data;

    !NOTE THAT HERE WE ARE POINTING TO l_data
    !AND NOT TO DM%m_data!
    left  => l_data;
    right => r_data;

    left(1,1) = 350.0;

    print*,'VALUE :',DM%m_data(1,1);
    print*,"";

    !Check the rows and cols:
    IF (SIZE(left,2) .NE. SIZE(right,1) ) THEN
       !PRINT *,'Error, tried to multiply matrix with incompatible size';
       ERROR STOP 'Error, tried to multiply matrix with incompatible size';
    ENDIF

    !Deallocate and resize:
    if (allocated(DMR%m_data)) then
       deallocate(DMR%m_data);
    endif
    allocate(DMR%m_data(SIZE(left,1),SIZE(right,2)));
    do i=1,SIZE(left,1)
       do j=1,SIZE(right,2)
          DMR%m_data(i,j) = SUM(left(i,:)*right(:,j))
       enddo
    enddo
    
  END SUBROUTINE DENSE_MUL
END MODULE LIN_ALGEBRA


PROGRAM MAIN
  USE LIN_ALGEBRA;
  IMPLICIT NONE;
  INTEGER :: I,J = 0;
  REAL,TARGET :: R = 13;
  REAL, POINTER :: P;
  REAL, POINTER :: P1;
  TYPE(Triangular_matrix):: my_tm;
  INTEGER :: AlloState;

  TYPE(Dense_matrix) :: dm1,dm2,dm3;
  
  
  CALL MAKE_POINTER(P,R);
  
  P = 32.0; !In fortran P is not the address,
  !but it is the reference itself!
  PRINT*,'Is P associated with the target R?:',ASSOCIATED(P,R);
  PRINT*, P, R

  !We can also allocate memory here:
  ALLOCATE(P1,STAT =AlloState);
  IF(AlloState == 0) P1 = 0.5;
  DEALLOCATE(P1);
  
  !Examples of the usage for allocating a triangular matrix:
  CALL INITIALIZE_MATRIX_ONES(my_tm,4);

  PRINT*,'Matrix number of rows: ',my_tm%m_nrows;
  !PRINT*,'Matrix first row: ',my_tm%m_rows(1)%row(:);
  DO I = 1,my_tm%m_nrows
     !999 FORMAT(' ');
     PRINT *,my_tm%m_rows(i)%row(1:i);
     !DO J = 1,SIZE(my_tm%m_rows(i)%row)
     !   PRINT '(f4.2)',my_tm%m_rows(i)%row(J);
     !END DO
  END DO


  !--Testing the allocation for the dense matrix as below--!
  ALLOCATE( &
       dm1%m_data(4,4),&
       dm2%m_data(4,4))
  dm1%m_data = 1.5;
  dm2%m_data = 2.0;

  CALL DENSE_MUL(dm1,dm2,dm3);

  print*,'Results of the multiplication is:'
  print *,dm3%m_data;
  
  
CONTAINS
  SUBROUTINE MAKE_POINTER(P,R)
    REAL,TARGET,INTENT(INOUT) :: R;
    REAL,POINTER,INTENT(INOUT) :: P;
    NULLIFY(P); !MAKE_UNDEFINED.
    P=>R; !Make pointer
  END SUBROUTINE MAKE_POINTER
  
END PROGRAM MAIN
