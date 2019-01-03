        !COMPILER-GENERATED INTERFACE MODULE: Thu Jan 03 16:49:48 2019
        ! This source file is for reference only and may not completely
        ! represent the generated interface used by the compiler.
        MODULE DENSM__genmod
          INTERFACE 
            FUNCTION DENSM(ALT,D0,XM,TZ,MN3,ZN3,TN3,TGN3,MN2,ZN2,TN2,   &
     &TGN2)
              INTEGER(KIND=4) :: MN2
              INTEGER(KIND=4) :: MN3
              REAL(KIND=4) :: ALT
              REAL(KIND=4) :: D0
              REAL(KIND=4) :: XM
              REAL(KIND=4) :: TZ
              REAL(KIND=4) :: ZN3(MN3)
              REAL(KIND=4) :: TN3(MN3)
              REAL(KIND=4) :: TGN3(2)
              REAL(KIND=4) :: ZN2(MN2)
              REAL(KIND=4) :: TN2(MN2)
              REAL(KIND=4) :: TGN2(2)
              REAL(KIND=4) :: DENSM
            END FUNCTION DENSM
          END INTERFACE 
        END MODULE DENSM__genmod
