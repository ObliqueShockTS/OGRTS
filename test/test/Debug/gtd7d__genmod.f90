        !COMPILER-GENERATED INTERFACE MODULE: Thu Jan 03 16:49:48 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE GTD7D__genmod
          INTERFACE 
            SUBROUTINE GTD7D(IYD,SEC,ALT,GLAT,GLONG,STL,F107A,F107,AP,  &
     &MASS,D,T)
              INTEGER(KIND=4) :: IYD
              REAL(KIND=4) :: SEC
              REAL(KIND=4) :: ALT
              REAL(KIND=4) :: GLAT
              REAL(KIND=4) :: GLONG
              REAL(KIND=4) :: STL
              REAL(KIND=4) :: F107A
              REAL(KIND=4) :: F107
              REAL(KIND=4) :: AP(7)
              INTEGER(KIND=4) :: MASS
              REAL(KIND=4) :: D(9)
              REAL(KIND=4) :: T(2)
            END SUBROUTINE GTD7D
          END INTERFACE 
        END MODULE GTD7D__genmod
