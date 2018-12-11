        !COMPILER-GENERATED INTERFACE MODULE: Wed Nov 28 15:12:00 2018
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DENSU__genmod
          INTERFACE 
            FUNCTION DENSU(ALT,DLB,TINF,TLB,XM,ALPHA,TZ,ZLB,S2,MN1,ZN1, &
     &TN1,TGN1)
              INTEGER(KIND=4) :: MN1
              REAL(KIND=4) :: ALT
              REAL(KIND=4) :: DLB
              REAL(KIND=4) :: TINF
              REAL(KIND=4) :: TLB
              REAL(KIND=4) :: XM
              REAL(KIND=4) :: ALPHA
              REAL(KIND=4) :: TZ
              REAL(KIND=4) :: ZLB
              REAL(KIND=4) :: S2
              REAL(KIND=4) :: ZN1(MN1)
              REAL(KIND=4) :: TN1(MN1)
              REAL(KIND=4) :: TGN1(2)
              REAL(KIND=4) :: DENSU
            END FUNCTION DENSU
          END INTERFACE 
        END MODULE DENSU__genmod
