        !COMPILER-GENERATED INTERFACE MODULE: Sun Dec 09 17:00:33 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE CELLCENTER__genmod
          INTERFACE 
            SUBROUTINE CELLCENTER(X,Y,Z,XC,YC,ZC,NX,NY,NZ)
              INTEGER(KIND=4), INTENT(IN) :: NZ
              INTEGER(KIND=4), INTENT(IN) :: NY
              INTEGER(KIND=4), INTENT(IN) :: NX
              REAL(KIND=8), INTENT(IN) :: X(NX,NY,NZ)
              REAL(KIND=8), INTENT(IN) :: Y(NX,NY,NZ)
              REAL(KIND=8), INTENT(IN) :: Z(NX,NY,NZ)
              REAL(KIND=8), INTENT(OUT) :: XC(NX-1,NY-1,NZ-1)
              REAL(KIND=8), INTENT(OUT) :: YC(NX-1,NY-1,NZ-1)
              REAL(KIND=8), INTENT(OUT) :: ZC(NX-1,NY-1,NZ-1)
            END SUBROUTINE CELLCENTER
          END INTERFACE 
        END MODULE CELLCENTER__genmod
