C $Id: sim_con.f,v 1.1 2009/02/17 23:35:42 saroun Exp $
C//////////////////////////////////////////////////////////////////////
C////
C////  R E S T R A X   4.4
C////
C////  Conversion subroutines between RESCAL and TRAX parameters
C////
C////  * SUBROUTINE RT_CONVRT
C////  * SUBROUTINE VIVF(VI,VF)
C////
C//////////////////////////////////////////////////////////////////////

C----------------------------------------------------
      SUBROUTINE READ_SOU(OBJ,LINE,IER)
C Read parameters of a SOURCE from the CFG file
C-----------------------------------------------------------
      use VCTABLE
      IMPLICIT NONE
      TYPE(TCLASS) :: OBJ
      INTEGER*4 NS,IER,LP
      REAL*8 DIA,W,H,VEC(1:4)
      CHARACTER*128 LINE
      character(64) :: SARG
      SARG=' '


      IER=0
      READ(LINE,*,ERR=99) NS,DIA,W,H

c       write(*,*) 'READ_SOU ',NS,DIA,W,H

C set values to OBJ
90    SELECT CASE (NS)
      CASE (0)
        VEC=(/2.D0,DIA*10.,DIA*10.,1.D-1/) ! disc
      CASE (1)
        VEC=(/3.D0,W*10.,H*10.,1.D-1/)     ! box
      CASE (2)
        VEC=(/2.D0,W*10.,H*10.,1.D-1/)     ! disc
      CASE (3)
        VEC=(/1.D0,W*10.,H*10.,1.D-1/)     ! cylinder
      END SELECT
!       CALL CMDOBJ_INP(OBJ,'SHAPE:SIZE',(/3.D0,W*10.,H*10.,1.D-1/),' ',LR,SARG,LP)
      CALL PCLASS_INP(OBJ,'SHAPE:SIZE',VEC,SARG,LP)
      RETURN

99    IER=1
      END

C-----------------------------------------------------------
      SUBROUTINE READ_DET(OBJ,LINE,IER)
C Read parameters of a DETECTOR from the CFG file
C-----------------------------------------------------------
      use VCTABLE
      IMPLICIT NONE
      TYPE(TCLASS) :: OBJ

      INTEGER*4 NS,NSEG,IER,LP,TYP
      REAL*8 DIA,W,H,EFF,SPACE,VEC(4),RESX,RESY
      CHARACTER*128 LINE
      character(64) :: SARG
      SARG=' '

      IER=0
C full format
      READ(LINE,*,ERR=1) TYP,NS,DIA,W,H,RESX,RESY,EFF,NSEG,SPACE
      GOTO 90
C older format
1     READ(LINE,*,ERR=2) NS,DIA,W,H,EFF,NSEG,SPACE
      IF (EFF.GT.0.AND.NSEG.GT.0) THEN
         TYP=1
      ELSE
         TYP=0
      ENDIF
      RESX=0.D0
      RESY=0.D0
      GOTO 90
C primitive format
2     READ(LINE,*,ERR=99) NS,DIA,W,H
      TYP=0
      RESX=0.D0
      RESY=0.D0
      EFF=0.D0
      NSEG=1
      SPACE=0.D0

C set values to OBJ
90    IF (NSEG.GT.64) NSEG=64
      SELECT CASE (NS)  ! shape, size
      CASE (0)
        VEC=(/2.D0,DIA*10.,DIA*10.,W*10./)  ! disc, circular, thickness = W
      CASE (1)
        VEC=(/3.D0,W*10.,H*10.,DIA*10./)    ! box (W x H), thickness = DIA
      CASE (2)
        VEC=(/2.D0,W*10.,H*10.,DIA*10./)    ! disc, elliptical (W x H) (axis along beam), thickness = DIA
      CASE (3)
        VEC=(/1.D0,W*10.,H*10.,DIA*10./)    ! cylinder elliptical (axis vertical), thickness = DIA
      END SELECT
      CALL PCLASS_INP(OBJ,'SHAPE:SIZE',VEC,SARG,LP)
      CALL PCLASS_INP(OBJ,'RES',(/RESX,RESY/),SARG,LP)
      CALL PCLASS_INP(OBJ,'TYPE:ALPHA:ND:SPACE',(/TYP*1.D0,EFF,NSEG*1.D0,SPACE/),SARG,LP)

      RETURN
99    IER=1
      END

C----------------------------------------------------
      SUBROUTINE READ_MONO(OBJ,LINE,IER)
C Read parameters of a CRYSTAL from the CFG file
C-----------------------------------------------------------
      use VCTABLE
      IMPLICIT NONE
      INCLUDE 'const.inc'
      TYPE(TCLASS) :: OBJ

      REAL*8 CHI,ANIZ,POIS,THICK,HEIGHT,WIDTH,V(3),DV(3)
      INTEGER*4 NH,NV,NB,IER,LP
      CHARACTER*128 LINE
      character(64) :: SARG
      SARG=' '

c      write(*,*) 'READ_MONO'

      IER=0
      V(1:3)=0.D0
      DV(1:3)=0.1D0 ! default spacing between segments in [mm]
      READ(LINE,*,ERR=1) CHI,ANIZ,POIS,THICK,HEIGHT,WIDTH,NH,NV,NB,V,DV ! full version (velocity + segemnt spacings in mm)
      GOTO 90
1     READ(LINE,*,ERR=2) CHI,ANIZ,POIS,THICK,HEIGHT,WIDTH,NH,NV,NB,V ! version with velocity [m/s]
      GOTO 90
2     READ(LINE,*,ERR=3) CHI,ANIZ,POIS,THICK,HEIGHT,WIDTH,NH,NV,NB   ! segments in 3 dimensions
      GOTO 90
3     READ(LINE,*,ERR=99) CHI,ANIZ,POIS,THICK,HEIGHT,WIDTH,NH,NV     ! segments in 2 dimensions
      NB=1

C set values to OBJ
90    CALL PCLASS_INP(OBJ,'SIZE',(/WIDTH*10.,HEIGHT*10.,THICK*10./),SARG,LP)
      CALL PCLASS_INP(OBJ,'VEL',V,SARG,LP)
      CALL PCLASS_INP(OBJ,'N',(/NH*1.D0,NV*1.D0,NB*1.D0/),SARG,LP)
      CALL PCLASS_INP(OBJ,'D',DV,SARG,LP)
      CALL PCLASS_INP(OBJ,'CHI:ANIZ:POISS',(/CHI,ANIZ,POIS/),SARG,LP)
      RETURN
99    IER=1
      END

C-----------------------------------------------------------
      SUBROUTINE READ_COL(OBJ,LINE,DIST,IC,IER)
C Read parameters of a GUIDE/COLLIMATOR from the CFG file
C IC .. collimator index
C-----------------------------------------------------------
      use VCTABLE
      IMPLICIT NONE
      INCLUDE 'collimators.inc'
      TYPE(TCLASS) :: OBJ

      REAL*8 DIST,LENG,H1,H2,V1,V2
      INTEGER*4 IC,IER,LP
      CHARACTER*128 LINE
      character(64) :: SARG
      SARG=' '

10    FORMAT(a15,8(1x,G10.4))
      IER=0

c      write(*,*) 'READ_COL'

C full format
      READ(LINE,*,ERR=1) CTYP(IC),DIST,LENG,H1,H2,V1,V2,ROH(IC),ROV(IC),
     &   GAMH(IC),GAMV(IC),REFH(IC),REFV(IC),NLAM(IC),VLAM(IC),
     &   DLAMH(IC),DLAMV(IC)
c      write(*,10) 'full format',IC,CTYP(IC),ROH(IC),ROV(IC)

      IF (CTYP(IC).GT.4.OR.CTYP(IC).LT.-1) GOTO 1 ! -1 >= CTYP <= 4, else error
      GOTO 90

C old RESTRAX format, extended by NLAM, VLAM, DLAM
1     READ(LINE,*,ERR=2) DIST,LENG,H1,H2,V1,V2,ROH(IC),GAMH(IC),
     &   GAMV(IC),REFH(IC),REFV(IC),NLAM(IC),VLAM(IC),DLAMH(IC)
c      write(*,*) 'format 1!'
      CTYP(IC)=0
      DLAMV(IC)=DLAMH(IC)
      ROV(IC)=0
      IF (GAMH(IC).GT.0.OR.GAMV(IC).GT.0) CTYP(IC)=1
      GOTO 90

C old RESTRAX format, extended by NLAM,DLAM
2     READ(LINE,*,ERR=3) DIST,LENG,H1,H2,V1,V2,ROH(IC),GAMH(IC),
     &                   GAMV(IC),REFH(IC),REFV(IC),NLAM(IC),DLAMH(IC)
c      write(*,*) 'format 2!'
      VLAM(IC)=0
      DLAMV(IC)=DLAMH(IC)
      ROV(IC)=0
      CTYP(IC)=0
      IF (GAMH(IC).GT.0.OR.GAMV(IC).GT.0) CTYP(IC)=1
      GOTO 90

C old RESTRAX format
3     READ(LINE,*,ERR=99)  DIST,LENG,H1,H2,V1,V2,ROH(IC),GAMH(IC),
     &     GAMV(IC),REFH(IC),REFV(IC)

c      write(*,*) 'format 3!'
      NLAM(IC)=0
      VLAM(IC)=0
      DLAMH(IC)=0.0
      DLAMV(IC)=DLAMH(IC)
      ROV(IC)=0
      CTYP(IC)=0
      IF (GAMH(IC).GT.0.OR.GAMV(IC).GT.0) CTYP(IC)=1


C set values to OBJ
90    IF (LENG.LE.0) CTYP(IC)=-1
      CALL PCLASS_INP(OBJ,'SIZE',(/H1*10.,V1*10.,LENG*10./),SARG,LP)
      CALL PCLASS_INP(OBJ,'DIST',(/DIST*10./),SARG,LP)
      CALL PCLASS_INP(OBJ,'EXIT',(/H2*10.,V2*10./),SARG,LP)
      CALL PCLASS_INP(OBJ,'RO',(/ROH(IC),ROV(IC)/),SARG,LP)
      CALL PCLASS_INP(OBJ,'N',(/NLAM(IC)*1.D0,VLAM(IC)*1.D0/),SARG,LP)
      CALL PCLASS_INP(OBJ,'DL',(/DLAMH(IC),DLAMV(IC)/),SARG,LP)
      CALL PCLASS_INP(OBJ,'M',(/GAMH(IC),GAMV(IC)/),SARG,LP)
      CALL PCLASS_INP(OBJ,'REF',(/REFH(IC),REFV(IC)/),SARG,LP)
      CALL PCLASS_INP(OBJ,'TYPE',(/1.D0*CTYP(IC)/),SARG,LP)

      RETURN
99    IER=1
      END


C--------------------------------------------------------------
      SUBROUTINE READ_COL1(OBJ,LINE,NF,DIST,IC,IER)
C Read parameters of a GUIDE/COLLIMATOR from the CFG file, with presence indicator
C--------------------------------------------------------------
      use VCTABLE
      IMPLICIT NONE
      INCLUDE 'collimators.inc'
      TYPE(TCLASS) :: OBJ

      REAL*8 DIST
      INTEGER*4 IC,NF,IS,IL,IER
      CHARACTER*128 LINE

      IER=0
      IS=1
      CALL FINDPAR(LINE,1,IS,IL)
      IF (IL.LE.0) GOTO 99
      READ(LINE(IS:IS+IL-1),*,ERR=99) NF
      CALL READ_COL(OBJ,LINE(IS+IL:128),DIST,IC,IER)

      RETURN
99    IER=1
      END

C----------------------------------------------------------
      SUBROUTINE READCFG(FNAME)
C  Imports TAS configuration from old format (*.cfg) file
C  RESCAL parameters (RES_DAT) are applied, so thay should be imported before
C----------------------------------------------------------
      USE COMPONENTS
      use VCTABLE
      use COMPONENTS_IO
      use INSTCONTROL
      use xmlparse
      use XMLSIMRES
      use FILETOOLS
      use IO
      IMPLICIT NONE

      INCLUDE 'const.inc'
      INCLUDE 'collimators.inc'
      INCLUDE 'rescal.inc'
      INCLUDE 'trax.inc'
      INCLUDE 'ness_common.inc'
      INCLUDE 'ness_cmd.inc'

      INTEGER*4 IU
      PARAMETER (IU=25)
      CHARACTER*(*) FNAME

      CHARACTER*128 LINE
      INTEGER*4 IER,ILINE,IS,IL,I,NCOL(4)
      REAL*8 VL0,VL1,VL2,VL3
      REAL*8 VLCANM,VLCANS,VLCANA,VLCAND,DGUIDE,DGA,DIST2A
      real(KIND(1.D0)) :: LAMI,LAMF,Z,D
  !    type(XML_PARSE) :: info

      type(TCLASS) :: myOBJ
  ! samples
      TYPE(TSCRYST) :: mySCRYST
  ! interfaces
      TYPE(TSPEC)   :: myINST


CxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCxxxxxxxxx
100   FORMAT(A60)
101   FORMAT(a)

c      write(*,*) 'READCFG'

      CALL BOUNDS(FNAME,IS,IL)
      OPEN(UNIT=IU,FILE=FNAME(IS:IS+IL-1),STATUS='OLD',ERR=999,IOSTAT=IER)

! initialize structures
      call SOURCE_DEFAULT(SOU)
      call GUIDE_DEFAULT(GDEA)
      call GUIDE_DEFAULT(GUIDE)
      call GUIDE_DEFAULT(SOL1)
      call GUIDE_DEFAULT(SOL2A)
      call GUIDE_DEFAULT(SOL2)
      call GUIDE_DEFAULT(SOL3)
      call GUIDE_DEFAULT(SOL4)
      call CRYSTAL_DEFAULT(MON)
      call CRYSTAL_DEFAULT(ANA)
      call DETECTOR_DEFAULT(DET)

C ***  READ CFG FILE  ***
      ILINE=0
C* Source
  3   READ(IU,*,END=997)
      READ(IU,100,END=997) CFGTITLE
      READ(IU,*,END=997)
      READ(IU,101,END=997) LINE
      myOBJ=PCLASS(PCOMOBJ(SOU))
      call SetClassAttrib(myOBJ,'SOU','source')
      CALL READ_SOU(myOBJ,LINE,IER)
      IF (IER.NE.0) GOTO 997
      ILINE=ILINE+2

C* Guide A
      READ(IU,*,END=997)
      READ(IU,101,END=997) LINE
      myOBJ=PCLASS(PCOMOBJ(GDEA))
      call SetClassAttrib(myOBJ,'GUIDEA','guide A')
      CALL  READ_COL(myOBJ,LINE,DGA,1,IER)
      IF (IER.NE.0) GOTO 997
      ILINE=ILINE+2

C* Guide B
      READ(IU,*,END=997)
      READ(IU,101,END=997) LINE
      myOBJ=PCLASS(PCOMOBJ(GUIDE))
      call SetClassAttrib(myOBJ,'GUIDEB','guide B')
      CALL  READ_COL1(myOBJ,LINE,NFG,DGUIDE,2,IER)
      IF (IER.NE.0) GOTO 997
      ILINE=ILINE+2
      IF (NFG.EQ.0) THEN ! ignore guides
         CTYP(2)=-1
         CTYP(1)=-1
      ENDIF

C* Monochromator
      READ(IU,*,END=997)
      READ(IU,101,END=997) LINE
      myOBJ=PCLASS(PCOMOBJ(MON))
      call SetClassAttrib(myOBJ,'MON','monochromator')
      CALL  READ_MONO(myOBJ,LINE,IER)
!      call PCLASS_INP(myOBJ,'SGN',RES_DAT(i_SM),SARG,LR)
!      call PCLASS_INP(myOBJ,'MOS',RES_DAT(i_ETAM),SARG,LR)
!      call PCLASS_INP(myOBJ,'DHKL',RES_DAT(i_DM),SARG,LR)
!      call PCLASS_INP(myOBJ,'RO',(/RES_DAT(i_ROMH),RES_DAT(i_ROMV),0.D0/),SARG,LR)
      IF (IER.NE.0) GOTO 997
      ILINE=ILINE+2
! no name for crystals - set by READCRYST

C* Analyzer
      READ(IU,*,END=997)
      READ(IU,101,END=997) LINE
      myOBJ=PCLASS(PCOMOBJ(ANA))
      call SetClassAttrib(myOBJ,'ANA','analyzer')
!      call PCLASS_INP(myOBJ,'SGN',RES_DAT(i_SA),SARG,LR)
!      call PCLASS_INP(myOBJ,'MOS',RES_DAT(i_ETAA),SARG,LR)
!      call PCLASS_INP(myOBJ,'DHKL',RES_DAT(i_DA),SARG,LR)
!      call PCLASS_INP(myOBJ,'RO',(/RES_DAT(i_ROAH),RES_DAT(i_ROAV),0.D0/),SARG,LR)
      CALL  READ_MONO(myOBJ,LINE,IER)
      IF (IER.NE.0) GOTO 997
      ILINE=ILINE+2
! no name for crystals - set by READCRYST

C* Detector
      READ(IU,*,END=997)
      READ(IU,101,END=997) LINE
      myOBJ=PCLASS(PCOMOBJ(DET))
      call SetClassAttrib(myOBJ,'DET','detector')
      call READ_DET(myOBJ,LINE,IER)
      IF (IER.NE.0) GOTO 997
      ILINE=ILINE+2

C* Distances
      READ(IU,*,END=997)
      READ(IU,*,END=997,ERR=997) VL0,VL1,VL2,VL3
      ILINE=ILINE+2
C* COL1
      READ(IU,*,END=997)
      READ(IU,101,END=997) LINE
      myOBJ=PCLASS(PCOMOBJ(SOL1))
      call SetClassAttrib(myOBJ,'COL1','collimator 1')
      CALL  READ_COL(myOBJ,LINE,VLCANM,3,IER)
      IF (IER.NE.0) GOTO 997
      ILINE=ILINE+2

C* COL2 A
      READ(IU,*,END=997)
      READ(IU,101,END=997) LINE
      myOBJ=PCLASS(PCOMOBJ(SOL2A))
      call SetClassAttrib(myOBJ,'COL2A','collimator 2A')
      CALL  READ_COL(myOBJ,LINE,DIST2A,4,IER)
      IF (IER.NE.0) GOTO 997
      ILINE=ILINE+2

C* COL2 B
      READ(IU,*,END=997)
      READ(IU,101,END=997) LINE
      myOBJ=PCLASS(PCOMOBJ(SOL2))
      call SetClassAttrib(myOBJ,'COL2B','collimator 2B')
      CALL  READ_COL(myOBJ,LINE,VLCANS,5,IER)
      IF (IER.NE.0) GOTO 997
      ILINE=ILINE+2

C* COL3
      READ(IU,*,END=997)
      READ(IU,101,END=997) LINE
      myOBJ=PCLASS(PCOMOBJ(SOL3))
      call SetClassAttrib(myOBJ,'COL3','collimator 3')
      CALL  READ_COL(myOBJ,LINE,VLCANA,6,IER)
      IF (IER.NE.0) GOTO 997
      ILINE=ILINE+2

C* COL4
      READ(IU,*,END=997)
      READ(IU,101,END=200) LINE ! last item can be followed by EOF => continue reading
      myOBJ=PCLASS(PCOMOBJ(SOL4))
      call SetClassAttrib(myOBJ,'COL4','collimator 4')
200   CALL  READ_COL(myOBJ,LINE,VLCAND,7,IER)
      IF (IER.NE.0) GOTO 997
      ILINE=ILINE+2

C ***  END OF READ CFG FILE, go to interpretation  ***
      CLOSE(IU)
C* SAMPLE
      call PSAMOBJ_DEFAULT(PSAMOBJ(mySCRYST))
      mySCRYST%SAM%FRAME%SHAPE=FRAME_SHAPE_CYLINDER
      myOBJ=PCLASS(PSAMOBJ(mySCRYST))
      call SetClassAttrib(myOBJ,'SAM','sample')
!      call PCLASS_INP(myOBJ,'MOS',(/RES_DAT(i_ETAS)/),SARG,LR)
!      call PCLASS_INP(myOBJ,'SIZE',(/RES_DAT(i_SDI)*10.,RES_DAT(i_SHI)*10.,RES_DAT(i_SDI)*10./),SARG,LR)
!      call PCLASS_INP(myOBJ,'CELLS',RES_DAT(i_AS),SARG,LR)
!      call PCLASS_INP(myOBJ,'CELLA',RES_DAT(i_AA),SARG,LR)
!      call PCLASS_INP(myOBJ,'VECA',RES_DAT(i_AX),SARG,LR)
!      call PCLASS_INP(myOBJ,'VECB',RES_DAT(i_BX),SARG,LR)
!      TAU=(/NINT(RES_DAT(i_QH)),NINT(RES_DAT(i_QK)),NINT(RES_DAT(i_QL))/)
!      Q=(/RES_DAT(i_QH)-TAU(1),RES_DAT(i_QK)-TAU(2),RES_DAT(i_QL)-TAU(3)/)
!      call PCLASS_INP(myOBJ,'TAU',TAU,SARG,LR)
!      call PCLASS_INP(myOBJ,'QHKL',Q,SARG,LR)
!      call PCLASS_INP(myOBJ,'EN',RES_DAT(i_EN),SARG,LR)
!      call PCLASS_INP(myOBJ,'GHKL',RES_DAT(i_GH),SARG,LR)
!      call PCLASS_INP(myOBJ,'GRAD',RES_DAT(i_GMOD),SARG,LR)

C Distances
C assume source->sample direction of ray-tracing
      SOU%FRAME%DIST=0.D0
      if (NFG.GT.0) then
        SOL1.FRAME.DIST=VLCANM*10.+GUIDE.FRAME.SIZE(3)
      ELSE
        SOL1.FRAME.DIST=(VLCANM-DGUIDE)*10.
      ENDIF
      MON.FRAME.DIST=(VL0-VLCANM)*10.
      GDEA.FRAME.DIST=(DGUIDE-DGA)*10.-GDEA.FRAME.SIZE(3)
      GUIDE.FRAME.DIST=DGUIDE*10.-GDEA.FRAME.DIST
      SOL2A.FRAME.DIST=(VLCANS-DIST2A)*10.-SOL2A.FRAME.SIZE(3)
      SOL2.FRAME.DIST=DIST2A*10.+SOL2A.FRAME.SIZE(3)
      mySCRYST%SAM%FRAME%DIST=VL1*10.-SOL2.FRAME.DIST-SOL2A.FRAME.DIST
      SOL3.FRAME.DIST=VLCANA*10.
      SOL4.FRAME.DIST=VLCAND*10.
      ANA.FRAME.DIST=VL2*10.-SOL3.FRAME.DIST
      DET.FRAME.DIST=VL3*10.-SOL4.FRAME.DIST

  ! add interface
      call PINSOBJ_DEFAULT(PINSOBJ(myINST))
      myOBJ=PCLASS(PINSOBJ(myINST))
      call SetClassAttrib(myOBJ,'TAS','three-axis spectrometer')
!      call PCLASS_INP(myOBJ,'KI:KF:Q0',(/STP%KI,STP%KF,STP%Q/),SARG,LR)
!      call PCLASS_INP(myOBJ,'FIX',(/1.D0*(STP%NFX-1)/),SARG,LR)
!      call PCLASS_INP(myOBJ,'SS',(/1.D0*(STP%SS)/),SARG,LR)

C additional settings dependent on RESCAL parameters
C (copied from NESS_CONV)
      LAMI=2*PI/STP.KI
      LAMF=2*PI/STP.KF
      NCOL=(/NFM,NFS,NFA,NFD/)
      DO I=1,4
        ALPHA(I)=RES_DAT(i_ALF1+I-1)
        BETA(I)=RES_DAT(i_BET1+I-1)
        NCOL(I)=-1
        select case(I)
        case(1)
          D=SOL1%FRAME%DIST
        case(2)
          D=SOL2%FRAME%DIST
        case(3)
          D=SOL3%FRAME%DIST
        case(4)
          D=SOL4%FRAME%DIST
        end select
        IF ((ALPHA(I).GE.500).AND.(D.NE.0)) THEN
          ALPHA(I)=0.
          NCOL(I)=1
        ENDIF
      ENDDO
      NFM=NCOL(1)
      NFS=NCOL(2)
      NFA=NCOL(3)
      NFD=NCOL(4)
      Z=0.D0
      IF(ALPHA(2).GT.0) Z=500.D0 ! no lamellae in SOL2A
      CALL CREATE_COL1(GDEA,NFG,0.D0,0.D0,1,LAMI)
      CALL CREATE_COL1(GUIDE,NFG,0.D0,0.D0,2,LAMI)
      CALL CREATE_COL1(SOL1,NFM,ALPHA(1),BETA(1),3,LAMI)
      CALL CREATE_COL1(SOL2A,NFS,Z,Z,4,LAMI)
      CALL CREATE_COL1(SOL2,NFS,ALPHA(2),BETA(2),5,LAMI)
      CALL CREATE_COL1(SOL3,NFA,ALPHA(3),BETA(3),6,LAMF)
      CALL CREATE_COL1(SOL4,NFD,ALPHA(4),BETA(4),7,LAMF)

      call DeleteComponents
      MONOCHROMATORSID='MON'
      ANALYZERSID='ANA'
  ! interface
      i=AddInterface(PINSOBJ(myINST))
  ! primary beam
      i=AddComponent(PCOMOBJ(SOU),1)
      i=AddComponent(PCOMOBJ(GDEA),1)
      i=AddComponent(PCOMOBJ(GUIDE),1)
      i=AddComponent(PCOMOBJ(SOL1),1)
      i=AddComponent(PCOMOBJ(MON),1)
      i=AddComponent(PCOMOBJ(SOL2A),1)
      i=AddComponent(PCOMOBJ(SOL2),1)
  ! sample
      i=AddSpecimen(PSAMOBJ(mySCRYST))
  ! secondary beam
      i=AddComponent(PCOMOBJ(SOL3),2)
      i=AddComponent(PCOMOBJ(ANA),2)
      i=AddComponent(PCOMOBJ(SOL4),2)
      i=AddComponent(PCOMOBJ(DET),2)

      call RESCAL_IMPORT
      call PCLASS_CLEAR(myOBJ)

      write(SOUT,*) 'Configuration updated from ',FNAME(IS:IS+IL-1)
      RETURN

C ***  ERROR while reading, show message and go to interpretation anyway ***
997   WRITE(SOUT,*) 'ERROR at line ',ILINE+2,' in file '//FNAME
      WRITE(SOUT,*) LINE
      pause
      CLOSE(IU)
      CALL BEAM_INP  ! copy back data from TAS array to objects
      RETURN

899   WRITE(*,898) IER
898   FORMAT('Fatal error: ',I5,' cannot open output file res_setup.xml. ',/,
     * 'Check privileges or disk space !')
      STOP

999   WRITE(*,998) IER
998   FORMAT('Fatal error: ',I5,' cannot open configuration file. ',/,
     * 'Check privileges or disk space !')
      STOP
      END


C-------------------------------------------
      SUBROUTINE RESCAL_IMPORT
C Import rescal parameters
! Works only for TAS configurations
C-------------------------------------------
      USE CONSTANTS
      USE FILETOOLS
      use VCTABLE
      use COMPONENTS
      use COMPONENTS_IO
      use INSTCONTROL
      IMPLICIT NONE
      include 'rescal.inc'
      real(KIND(1.D0)) :: KFIX,KI,KF,Z,TAU(3),Q(3)
      INTEGER :: NFX,LR,ierr
      type(TCLASS) :: myOBJ
      character(64) :: SARG
      SARG=' '
      if (SAMOBJ%SCLS.eq.SCLS_SINGLE) then
        myOBJ=PCLASS(SAMOBJ)
        call PCLASS_INP(myOBJ,'MOS',(/RES_DAT(i_ETAS)/),SARG,LR)
        call PCLASS_INP(myOBJ,'SIZE',(/RES_DAT(i_SDI)*10.,RES_DAT(i_SHI)*10.,RES_DAT(i_SDI)*10./),SARG,LR)
        call PCLASS_INP(myOBJ,'CELLS',RES_DAT(i_AS),SARG,LR)
        call PCLASS_INP(myOBJ,'CELLA',RES_DAT(i_AA),SARG,LR)
        call PCLASS_INP(myOBJ,'VECA',RES_DAT(i_AX),SARG,LR)
        call PCLASS_INP(myOBJ,'VECB',RES_DAT(i_BX),SARG,LR)
        TAU=(/NINT(RES_DAT(i_QH)),NINT(RES_DAT(i_QK)),NINT(RES_DAT(i_QL))/)
        Q=(/RES_DAT(i_QH)-TAU(1),RES_DAT(i_QK)-TAU(2),RES_DAT(i_QL)-TAU(3)/)
        call PCLASS_INP(myOBJ,'TAU',TAU,SARG,LR)
        call PCLASS_INP(myOBJ,'QHKL',Q,SARG,LR)
        call PCLASS_INP(myOBJ,'EN',RES_DAT(i_EN),SARG,LR)
        call PCLASS_INP(myOBJ,'GHKL',RES_DAT(i_GH),SARG,LR)
        call PCLASS_INP(myOBJ,'GRAD',RES_DAT(i_GMOD),SARG,LR)
        RES_DAT(i_A3)=SAMOBJ%P_FRAME%GON(1)/deg
        RES_DAT(i_A4)=SAMOBJ%P_FRAME%AX(1)/deg
      endif

      if (MONOCHROMATORS(1)%CCLS.eq.CCLS_CRYSTAL) then
        myOBJ=PCLASS(MONOCHROMATORS(1))
        call PCLASS_INP(myOBJ,'SGN',RES_DAT(i_SM),SARG,LR)
        call PCLASS_INP(myOBJ,'MOS',RES_DAT(i_ETAM),SARG,LR)
        call PCLASS_INP(myOBJ,'DHKL',RES_DAT(i_DM),SARG,LR)
        call PCLASS_INP(myOBJ,'RO',(/RES_DAT(i_ROMH),RES_DAT(i_ROMV),0.D0/),SARG,LR)
      endif

      if (ANALYZERS(1)%CCLS.eq.CCLS_CRYSTAL) then
        myOBJ=PCLASS(ANALYZERS(1))
        call PCLASS_INP(myOBJ,'SGN',RES_DAT(i_SA),SARG,LR)
        call PCLASS_INP(myOBJ,'MOS',RES_DAT(i_ETAA),SARG,LR)
        call PCLASS_INP(myOBJ,'DHKL',RES_DAT(i_DA),SARG,LR)
        call PCLASS_INP(myOBJ,'RO',(/RES_DAT(i_ROAH),RES_DAT(i_ROAV),0.D0/),SARG,LR)
      endif

      if (INSOBJ%ICLS.gt.0) then
        NFX=NINT(RES_DAT(i_FX))
        KFIX=RES_DAT(i_KFIX)
        KI=KFIX
        KF=KFIX
        IF (NFX.EQ.1) THEN
          Z=KFIX**2-RES_DAT(i_EN)/HSQOV2M
          if (Z.gt.0.D0) KF=SQRT(Z)
        ELSE
          Z=KFIX**2+RES_DAT(i_EN)/HSQOV2M
          if (Z.gt.0.D0) KI=SQRT(Z)
        ENDIF
        myOBJ=PCLASS(INSOBJ)
        call PCLASS_INP(myOBJ,'KI:KF',(/KI,KF/),SARG,LR)
        call PCLASS_INP(myOBJ,'FIX',(/1.D0*(NFX-1)/),SARG,LR)
        call PCLASS_INP(myOBJ,'SS',(/1.D0*(RES_DAT(i_SS))/),SARG,LR)
      endif
      call INST_ADJUST(INSOBJ,ierr)
      CALL BEAM_OUT  ! make a copy from objects to BPAR_ARRAY
      end SUBROUTINE RESCAL_IMPORT

C-------------------------------------------
      SUBROUTINE RESCAL_EXPORT
C Export rescal parameters
! Works only for TAS configurations
C-------------------------------------------
      USE CONSTANTS
      USE FILETOOLS
      use VCTABLE
      use COMPONENTS
      IMPLICIT NONE
      include 'rescal.inc'
      real(KIND(1.D0)) :: TAU(3),Q(3),AUX(64)
      character(64) :: SARG
      INTEGER :: NFX,LR,i
      type(TCLASS) :: myOBJ
      if (SAMOBJ%SCLS.eq.SCLS_SINGLE) then
        myOBJ=PCLASS(SAMOBJ)
        call PCLASS_OUT(myOBJ,'MOS',RES_DAT(i_ETAS),SARG,LR)
        call PCLASS_OUT(myOBJ,'SIZE',AUX,SARG,LR)
        RES_DAT(i_SDI:i_SHI)=AUX(1:2)/10.D0
        call PCLASS_OUT(myOBJ,'CELLS',RES_DAT(i_AS),SARG,LR)
        call PCLASS_OUT(myOBJ,'CELLA',RES_DAT(i_AA),SARG,LR)
        call PCLASS_OUT(myOBJ,'VECA',RES_DAT(i_AX),SARG,LR)
        call PCLASS_OUT(myOBJ,'VECB',RES_DAT(i_BX),SARG,LR)
        call PCLASS_OUT(myOBJ,'TAU',TAU,SARG,LR)
        call PCLASS_OUT(myOBJ,'QHKL',Q,SARG,LR)
        do i=1,3
          RES_DAT(i_QH+i-1)=TAU(i)+Q(i)
        enddo
        call PCLASS_OUT(myOBJ,'EN',RES_DAT(i_EN),SARG,LR)
        call PCLASS_OUT(myOBJ,'GHKL',RES_DAT(i_GH),SARG,LR)
        call PCLASS_OUT(myOBJ,'GRAD',RES_DAT(i_GMOD),SARG,LR)
      endif

      if (MONOCHROMATORS(1)%CCLS.eq.CCLS_CRYSTAL) then
        myOBJ=PCLASS(MONOCHROMATORS(1))
        call PCLASS_OUT(myOBJ,'SGN',RES_DAT(i_SM),SARG,LR)
        call PCLASS_OUT(myOBJ,'MOS',RES_DAT(i_ETAM),SARG,LR)
        call PCLASS_OUT(myOBJ,'DHKL',RES_DAT(i_DM),SARG,LR)
        call PCLASS_OUT(myOBJ,'RO',AUX,SARG,LR); RES_DAT(i_ROMH:i_ROMV)=AUX(1:2)
      endif

      if (ANALYZERS(1)%CCLS.eq.CCLS_CRYSTAL) then
        myOBJ=PCLASS(ANALYZERS(1))
        call PCLASS_OUT(myOBJ,'SGN',RES_DAT(i_SA),SARG,LR)
        call PCLASS_OUT(myOBJ,'MOS',RES_DAT(i_ETAA),SARG,LR)
        call PCLASS_OUT(myOBJ,'DHKL',RES_DAT(i_DA),SARG,LR)
        call PCLASS_OUT(myOBJ,'RO',AUX,SARG,LR); RES_DAT(i_ROAH:i_ROAV)=AUX(1:2)
      endif

      if (INSOBJ%ICLS.gt.0) then
        myOBJ=PCLASS(INSOBJ)
        call PCLASS_OUT(myOBJ,'FIX',AUX,SARG,LR)
        call PCLASS_OUT(myOBJ,'SS',RES_DAT(i_SS),SARG,LR)
        RES_DAT(i_FX)=AUX(1)+1.D0
        NFX=NINT(RES_DAT(i_FX))
        IF (NFX.EQ.1) THEN
          call PCLASS_OUT(myOBJ,'KI',RES_DAT(i_KFIX),SARG,LR)
        ELSE
          call PCLASS_OUT(myOBJ,'KF',RES_DAT(i_KFIX),SARG,LR)
        ENDIF
      endif
      if (associated(SAMOBJ%P_FRAME)) then
        RES_DAT(i_A3)=SAMOBJ%P_FRAME%GON(1)/deg
        RES_DAT(i_A4)=SAMOBJ%P_FRAME%AX(1)/deg
      endif

      end SUBROUTINE RESCAL_EXPORT

C----------------------------------------------------------
      SUBROUTINE TURN_FRAME(OBJ,DIST)
C turn FRAME by 180 deg around vertical axis
C----------------------------------------------------------
      USE COMPONENTS
      IMPLICIT NONE
      TYPE(TFRAME) :: OBJ
      REAL*8 DIST

      OBJ.GON(2)=-OBJ.GON(2)
      OBJ.STA(1)=-OBJ.STA(1)
      OBJ.STA(3)=-OBJ.STA(3)
      OBJ.DIST=DIST
      OBJ.AX(3)=-OBJ.AX(3)
      END


C----------------------------------------------------------
      SUBROUTINE TURN_BENDER(OBJ,DIST)
C turn BENDER by 180 deg around vertical axis
C apply axis parameters from arguments (DIST,AXH,AXV)
C----------------------------------------------------------
      USE COMPONENTS
      IMPLICIT NONE
      TYPE(BENDER) :: OBJ
      REAL*8 DIST,W1,H1

      W1=OBJ.FRAME.SIZE(1)
      H1=OBJ.FRAME.SIZE(2)
      OBJ.FRAME.SIZE(1)=OBJ.W2
      OBJ.FRAME.SIZE(2)=OBJ.H2
      OBJ.W2=W1
      OBJ.H2=H1
      IF (OBJ.TYP.LE.1) OBJ.ROH=-OBJ.ROH
      call TURN_FRAME(OBJ.FRAME,DIST)

      END

C----------------------------------------------------------
      SUBROUTINE TURN_CRYSTAL(OBJ,DIST)
C turn CRYSTAL by 180 deg around vertical axis
C apply axis parameters from arguments (DIST,AXH,AXV)
C----------------------------------------------------------
      USE COMPONENTS
      IMPLICIT NONE
      TYPE (CRYSTAL) :: OBJ
      REAL*8 DIST

      OBJ.FRAME.GON(1)=-OBJ.FRAME.GON(1)-2*OBJ.CHI
      OBJ.FRAME.AX(1)=-OBJ.FRAME.AX(1)
      call TURN_FRAME(OBJ.FRAME,DIST)
      END


C----------------------------------------------------------
      SUBROUTINE TURN_UPSTREAM
C converts the primary spectrometer to up-stream geometry
C exchange entry/exit windows, distances etc.
C NOTE: TAS array always contains parameters in downstream geometry
C----------------------------------------------------------
      USE COMPONENTS
      IMPLICIT NONE
      INCLUDE 'ness_common.inc'

      REAL*8 DIST

c      DIST=GUIDE.FRAME.DIST+GUIDE.FRAME.SIZE(3)
      DIST=GDEA.FRAME.DIST+GDEA.FRAME.SIZE(3)

      call TURN_FRAME(SOU%FRAME,DIST)

      DIST=GUIDE.FRAME.DIST+GUIDE.FRAME.SIZE(3)-GDEA.FRAME.SIZE(3)
      call TURN_BENDER(GDEA,DIST)

      DIST=SOL1.FRAME.DIST+SOL1.FRAME.SIZE(3)-GUIDE.FRAME.SIZE(3)
      call TURN_BENDER(GUIDE,DIST)

      DIST=MON.FRAME.DIST-SOL1.FRAME.SIZE(3)
      call TURN_BENDER(SOL1,DIST)

      DIST=SOL2A.FRAME.DIST+SOL2A.FRAME.SIZE(3)
      call TURN_CRYSTAL(MON,DIST)

      DIST=SOL2.FRAME.DIST+SOL2.FRAME.SIZE(3)-SOL2A.FRAME.SIZE(3)
      call TURN_BENDER(SOL2A,DIST)

      DIST=SAM.FRAME.DIST-SOL2.FRAME.SIZE(3)
      call TURN_BENDER(SOL2,DIST)

      DIST=0.D0
      call TURN_FRAME(SAM%FRAME,DIST)
      END

C----------------------------------------------------------
      SUBROUTINE WRITEDEFCFG
C  Write default config. file to the current directory
C  Use the full format (new in since version 4.9.92)
C TYP=-1 no collimator
C TYP=0  standard collimator (course or soller)
C TYP=1  guide (or bender), can be curved (RO means curvature in [1/m] )
C TYP=2  parabolic guide, equal lengths of the lamellae (RO means focal distance in [cm] !)
C TYP=3  parabolic guide, optimized lengths of the lamellae
C TYP=4  elliptic guide, wider window is the smaller ellipse axis
C
C extended detector format:
C =========================
C if eta>0, then assumes vertical tube(s)
C nseg = number of tubes
C space = space between tubes [mm]
C phi = inclination angle [deg], phi=0 means tube is vertical
C detector efficiency: 1-EXP(-eta*lambda*pathlength)
C----------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*16 C(13)
      CHARACTER*82 VALS(13),HEAD(13)
      INTEGER*4 I
C element names
      DATA C/
     &  'title',
     &  'n-guide A',
     &  'n-guide B',
     &  'monochromator',
     &  'analyzer',
     &  'segments',
     &  'detector',
     &  'distances',
     &  '1st collimator',
     &  '2nd collimator A',
     &  '2nd collimator B',
     &  '3nd collimator',
     &  '4th collimator'/
c some headers

      VALS(1)='default setup (IN14 with PG and Soller)'
      VALS(2)='0   21.   6.   12.'
      VALS(3)='1  7.5  584.5  6.  6. 12. 12.  0.  0.  1.2  1.2  1  1  1  1  0. 0.'
      VALS(4)='1  1  807.5  1050.  6.  6.  12.  12.  0.  0. 1.2  1.2  1  1  1  1  0. 0.'
      VALS(5)='0.  1.    0.3  0.2  12.  15.  1  7  1'
      VALS(6)='0.  1.    0.3  1.0  12.  16.  1  5  1'
      VALS(7)='0   1     4.0  4.0  6.0  0.0  0.0  0.0   1   0.'
      VALS(8)='22.5  270.  155.  64.'
      VALS(9)= '0     1.     5.  10.  10.   15.   15.  0.  0.  0.  0.  0.  0. 1  1  0. 0.'
      VALS(10)='-1    2.     0.   6.   6    12.   12.  0.  0.  0.  0.  0.  0. 1  1  0. 0.'
      VALS(11)='0   100.    20.   4.   4.   12.   12   0.  0.  0.  0.  0.  0. 1  1  0. 0.'
      VALS(12)='0   113.    20.   4.   4.   12.   12.  0.  0.  0.  0.  0.  0. 1  1  0. 0.'
      VALS(13)='0    35.    20.   4.   4.   12.   12.  0.  0.  0.  0.  0.  0. 1  1  0. 0.'

      HEAD(1)='title (max.60 characters)'
      HEAD(2)='shape, dia, width, height'
      HEAD(3)='type, gap, len, H1, H2, V1, V2, roh, rov, mh, mv, refh, refv, nh, nv, dh, dv'
      HEAD(4)='use, type, dist, len, H1, H2, V1, V2, roh, rov, mh, mv, refh, refv, nh, nv, dh, dv'
      HEAD(5)='chi, aniz., poiss., thick., height, length, nh, nv, nt'
      HEAD(6)=HEAD(5)
      HEAD(7)='type, shape, dia, width, height, resx[mm], resy[mm], alpha[1/cm/A], nseg,space[mm]'
      HEAD(8)='sou-mono, mono-sample, sample-anal, anal-det'
      HEAD(9)='typ, dist, len, H1, H2, V1, V2, roh, rov, mh, mv, refh, refv, nh, nv, dh, dv'
      HEAD(10)=HEAD(3)
      HEAD(11)=HEAD(9)
      HEAD(12)=HEAD(9)
      HEAD(13)=HEAD(9)

1     FORMAT(a13,' (',a,')')
2     FORMAT(a)


      OPEN(UNIT=24,FILE='simres00.cfg',STATUS='UNKNOWN',ERR=10)
      DO I=1,13
         write(24,1) C(I),HEAD(I)
         write(24,2) VALS(I)
      ENDDO
      CLOSE(24)
10    continue
      END




