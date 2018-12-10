        !COMPILER-GENERATED INTERFACE MODULE: Sun Dec 09 17:01:35 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE WPROF__genmod
          INTERFACE 
            FUNCTION WPROF(Z,ZL,S,UINF,ULB,ULBD,MN1,ZN1,UN1,UGN1,MN2,ZN2&
     &,UN2,UGN2)
              INTEGER(KIND=4) :: MN2
              INTEGER(KIND=4) :: MN1
              REAL(KIND=4) :: Z
              REAL(KIND=4) :: ZL
              REAL(KIND=4) :: S
              REAL(KIND=4) :: UINF
              REAL(KIND=4) :: ULB
              REAL(KIND=4) :: ULBD
              REAL(KIND=4) :: ZN1(MN1)
              REAL(KIND=4) :: UN1(MN1)
              REAL(KIND=4) :: UGN1(2)
              REAL(KIND=4) :: ZN2(MN2)
              REAL(KIND=4) :: UN2(MN2)
              REAL(KIND=4) :: UGN2(2)
              REAL(KIND=4) :: WPROF
            END FUNCTION WPROF
          END INTERFACE 
        END MODULE WPROF__genmod
