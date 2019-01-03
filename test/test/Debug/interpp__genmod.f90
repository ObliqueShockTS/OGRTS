        !COMPILER-GENERATED INTERFACE MODULE: Thu Jan 03 16:49:49 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE INTERPP__genmod
          INTERFACE 
            SUBROUTINE INTERPP(X,Y,Z,XVAL,YVAL,ZVAL,FCN_3D,VAL)
              REAL(KIND=8) :: X(4)
              REAL(KIND=8) :: Y(4)
              REAL(KIND=8) :: Z(4)
              REAL(KIND=8) :: XVAL
              REAL(KIND=8) :: YVAL
              REAL(KIND=8) :: ZVAL
              REAL(KIND=8) :: FCN_3D(4,4,4)
              REAL(KIND=8) :: VAL
            END SUBROUTINE INTERPP
          END INTERFACE 
        END MODULE INTERPP__genmod
