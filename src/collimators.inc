C// Collimator parameters for NESS

      INTEGER*4 NCO
      PARAMETER (NCO=7)  ! number of collimators in the setup
      INTEGER*4 NLAM(NCO),VLAM(NCO),CTYP(NCO)
      LOGICAL*4 OSC(NCO)

      REAL*8 ROH(NCO),ROV(NCO),GAMH(NCO),GAMV(NCO),REFH(NCO),REFV(NCO),
     &       DLAMH(NCO),DLAMV(NCO)

      COMMON /COLLIMATORS/ ROH,ROV,GAMH,GAMV,REFH,REFV,DLAMH,DLAMV,
     &                     NLAM,VLAM,CTYP,OSC




