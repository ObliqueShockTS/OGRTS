        !COMPILER-GENERATED INTERFACE MODULE: Thu Jan 03 16:49:48 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE SPLINE__genmod
          INTERFACE 
            SUBROUTINE SPLINE(X,Y,N,YP1,YPN,Y2)
              INTEGER(KIND=4) :: N
              REAL(KIND=4) :: X(N)
              REAL(KIND=4) :: Y(N)
              REAL(KIND=4) :: YP1
              REAL(KIND=4) :: YPN
              REAL(KIND=4) :: Y2(N)
            END SUBROUTINE SPLINE
          END INTERFACE 
        END MODULE SPLINE__genmod
