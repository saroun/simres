C//////////////////////////////////////////////////////////////////////
C////
C////  R E S T R A X   4.1
C////
C////  Subroutines from the program TRAX 
C////
C////  * SUBROUTINE THRAX
C////  * SUBROUTINE SSTV
C////  * SUBROUTINE TRCAN(VL0,VLI,VLC,S10,S20,NF,SHI)
C////  * SUBROUTINE TCR
C////  * SUBROUTINE HIRANG(OB,HII,HI)
C////  
C//////////////////////////////////////////////////////////////////////

C************************
      SUBROUTINE THRAX(READCONF)
C************************      
C     COMPUTES THE RESOLUTION MATRIX, THE NORMALIZATION FACTORS AND
C     DIFFERENT SCAN WIDTHS AND ABSOLUTE INTENSITIES AT DETECTOR
C     FOR A THREE-AXIS NEUTRON SPECTROMETER
      
      
      INCLUDE 'const.inc'

      implicit double precision (a-h,o-z)
      DIMENSION  AAA(4,4),EM(4,4)
      INTEGER*4 READCONF
      double precision CM1V(5,5),X33(1,1),AUX3(3,3),AUX4(4,4),
     +AU2(2,2),
     1AHP(3,4),AUX6(6,6),CNH(6,6),CNV(3,3),AU23(2,3),AUX2(2,2),
     2AH(4,6),AV(3,5),AV1(1,3),AU3(3,3),AH8(6,8),CKI0(3,3)
      CHARACTER*8 CRYTI(2)

      COMMON/MATRIX/ AAA
      COMMON/EXTRA/XXC,YYC,F,W,WI,WF      
      COMMON /ERROR/IERR            

      COMMON/CRYTI/CRYTI/CRYD/CRYD(2)/D0KOVD/D0KOVD(2)/SV/SVMO(31)
     1/COMIU/COMIU(10)/POISS/POISS(2)
      COMMON /DIMENS/NSOU,WSOU,HSOU,DIASOU,WMON,THMON,HMON,
     1WSAM,THSAM,HSAM,DIASAM,WANA,THANA,HANA
      COMMON /DETECT/NDET,NSEGDET,WDET,HDET,DIADET,ADET,
     &       SPACEDET,PHIDET
c      COMMON /EM/EMMIN1(4,4),EM(4,4),EMPMIN(4,4),VOLRES,VOLX,CMX(3,3)
      COMMON /EM/EMMIN1(4,4),VOLRES
      COMMON /CTR/DV(5,5),BM,BA,AHS(3,4)
      COMMON /SETP/VL0,VL1,VL2,VL3,VLMIN,VLMAX
      COMMON /ANGLES/OM,HIM,ETAMEF,OS,HIS,OA,HIA,ETAAEF,PHI,PHIM
      COMMON /VOLS/VOLSAM,VOLEF,FRACV,RADET,VOLCKI,VOLCKF,CKI(3,3),
     1CKF(3,3)
      COMMON /MEAN/VKI,SLAMDI,EI0,VKF,SLAMDF,EF0,Q0,HOMEGA,GRAD(3),DIRQ
     *,GRADI(3),TAU,DIRTAU,DIRTD,DIRQD,EPSLN,QSMALL,Q,DIRQT
      COMMON /COVAR/CMH(8,8),CMV(5,5),CMHP(4,4),CMVP(3,3)
      COMMON/STV/SH(10,10),SH1(10,10),SV(5,5),SV1(5,5),SHP(5,5),SVP(3,3)
      COMMON/MONOCH/TETAM,HIMON,ROHM,ROVM,ETAM,NGRADM,PM,PM0,RKINM,RLAMM
      COMMON/ANALYS/TETAA,HIANA,ROHA,ROVA,ETAA,NGRADA,PA,PA0,RKINA,RLAMA
      COMMON /VERTIC/ETVM,ETVMEF,ANRM,ETVA,ETVAEF,ANRA
      COMMON /SAMPLE/TETAS,HISAM,ETAS
      COMMON /INDEX/NDA,NCY,NEFIX,IM,IA,ISC,NCR1,NCR2,NSAM,NOPT,NWR,NCF
      COMMON /ABNORM/R0,R0CN,FLU,R0FLU
      COMMON /COARSE/NFM,VLCANM,VLSM,HDM1,HDM2,VDM1,VDM2,
     1               NFS,VLCANS,VLMS,HDS1,HDS2,VDS1,VDS2,
     2               NFA,VLCANA,VLSA,HDA1,HDA2,VDA1,VDA2,
     3               NFD,VLCAND,VLAD,HDD1,HDD2,VDD1,VDD2
      COMMON /SOLLER/NGUIDE,GAMACR,ALPHA(4),BETA(4)
      COMMON /NORMS/ RELTR,RELMC 
      
      RIGHT=PI/2.
  
      IF (READCONF.GT.0) CALL RT_CONVRT
              
C      write(*,*) 'THRAX - RT_CONVRT O.K.'      


      IF(NEFIX.GT.1)GO TO 14
      VKI=SQRT(2./HSQOVM*EI0)
      SLAMDI=2.*PI/VKI
      OM=ASIN(SLAMDI/2./CRYD(1))*IM
      TETAM=OM/deg
      SIM=SIN(OM)
      VKFM=PI/CRYD(2)
      VKF=VKI*VKI-HOMEGA*2./HSQOVM
      IF(VKF.LE.(VKFM*VKFM))GO TO 7
      VKF=SQRT(VKF)
      SLAMDF=2.*PI/VKF
      OA=ASIN(VKFM/VKF*IA)
      TETAA=OA/deg
      EF0=0.5*HSQOVM*VKF*VKF
      GO TO 15

   14 VKF=SQRT(2./HSQOVM*EF0)
      SLAMDF=2.*PI/VKF
      OA=ASIN(SLAMDF/2./CRYD(2))*IA
      TETAA=OA/deg
      SIA=SIN(OA)
      VKIM=PI/CRYD(1)
      VKI=VKF*VKF+HOMEGA*2./HSQOVM
      IF(VKI.LE.(VKIM*VKIM))GO TO 7
      VKI=SQRT(VKI)
      SLAMDI=2.*PI/VKI
      OM=ASIN(VKIM/VKI)*IM
      TETAM=OM/deg
      EI0=.5*HSQOVM*VKI*VKI

   15 Q0=Q

C      write(*,*) 'THRAX - 15 O.K.'      

     
   36 TWO=.5*(VKI/VKF+VKF/VKI-Q0*Q0/VKI/VKF)
      IF(ABS(TWO).GE.1.) GO TO 9
      TWO=ACOS(TWO)
      OS=ISC*TWO/2.
      TETAS=OS/deg
      PHI=ATAN2(-VKF*SIN(2.*OS),VKI-VKF*COS(2.*OS))
      PHIM=PHI-2.*OS
      CALL HIRANG(OM,HIMON,HIM)
      CALL HIRANG(OA,HIANA,HIA)
C      DIRTAU=DIRTD*deg
      IF(NSAM.EQ.0) GO TO 24

      IF(HIS.LE.(RIGHT-OS)) GO TO 29
      HIS=HIS-PI
C      DIRTAU=DIRTAU-PI
     
   29 IF(HIS.GE.(-RIGHT-OS)) GO TO 24
      HIS=HIS+PI
C      DIRTAU=DIRTAU+PI
      
   24 COSE=1.
      SINE=0.
      EPSLN=0.
C      write(*,*) 'THRAX - 24 O.K.'      

    8 ETAMEF=ETAM*minute
      ETVMEF=ETAM*ANRM*minute
      AR=ABS(OM+HIM)/(OM+HIM)
      CO=COS(HIM)
      SI=SIN(HIM)
      IF(ROHM.EQ.0.) GO TO 18
      TTEMP=POISS(1)/(1-POISS(1))
      AM = 1 - (1+TTEMP)*CO*CO
      BM = (1+TTEMP)*SI*CO
         
 
   18 ETAAEF=ETAA*minute
      ETVAEF=ETAA*ANRA*minute
      BET=ABS(OA+HIA)/(OA+HIA)
      CO=COS(HIA)
      SI=SIN(HIA)
      IF(ROHA.EQ.0.) GO TO 21
      TTEMP=POISS(2)/(1-POISS(2))
      AA = 1 - (1+TTEMP)*CO*CO
      BA = (1+TTEMP)*SI*CO
      
   21 PM=1.                   ! peak reflectivity=1 is supposed
      PA=1. 
   
C      write(*,*) 'THRAX - 21 O.K.'      
      CALL SSTV
C      write(*,*) 'THRAX - SSTV O.K.'                  
      CALL TCR
C      write(*,*) 'THRAX - TCR O.K.'      

      CALL BABT(CMV,DV,5,5,CM1V)
C     CM1V IS THE COVARIANCE MATRIX OF THE VERTICAL VARIABLES 
C     (DELTA0,DELTA1,DELTA2,DELTA3,ZS). 
      DO 1 I=1,6
      DO 1 J=1,8
    1 AH8(I,J)=0. 
      COT=1./TAN(OM)
      DO 30 I=1,3 
      DO 30 J=1,4 
   30 AHP(I,J)=0. 
      CO=COS(OM-HIM)/VL1*VKI
      SI=VKI*SIN(OM-HIM)/VL1
      AH8(1,1)=COT*(VKI*ROHM*AR-SI) 
      AH8(1,2)=-COT*CO+ROHM*AR*VKI*(-AM+COT*BM)
      AH8(2,2)=CO 
      AH8(2,1)=SI 
      CO=COS(OS+HIS)
      SI=SIN(OS+HIS)
      AH8(3,3)=CO 
      AH8(3,4)=SI 
      AH8(4,3)=-SI
      AH8(4,4)=CO 
      CO=VKI*CO/VL1 
      SI=VKI*SI/VL1 
      AH8(1,3)=-COT*CO
      AH8(1,4)=-COT*SI
      AH8(1,7)=VKI*COT
      AHP(1,3)=-AH8(1,7)/VL1
      IF(ETAM.EQ.0.) AH8(1,7)=AH8(1,7)*SIN(OM+HIM)/ABS(SIN(OM-HIM)) 
      AHP(1,4)=AH8(1,7) 
      AHP(1,1)=AH8(1,1) 
      AHP(1,2)=AH8(1,2) 
      AHP(2,1)=AH8(2,1) 
      AHP(2,2)=AH8(2,2) 
      AHP(3,3)=1. 
      AHP(2,3)=VKI/VL1
      AH8(2,3)=CO 
      AH8(2,4)=SI 
      COT=1./TAN(OA)
      CO=VKF*COS(OS-HIS)/VL2
      SI=VKF*SIN(OS-HIS)/VL2
      AH8(5,3)=-CO*COT
      AH8(5,4)=SI*COT 
      AH8(6,3)=-CO
      AH8(6,4)=SI 
      CO=VKF*COS(OA+HIA)/VL2
      SI=VKF*SIN(OA+HIA)/VL2
      AH8(5,5)=COT*(SI-ROHA*VKF*BET)
      AH8(5,6)=-COT*CO-VKF*ROHA*BET*(AA+COT*BA)
      AH8(5,8)=-VKF*COT 
      IF(ETAA.EQ.0.) AH8(5,8)=AH8(5,8)*SIN(OA-HIA)/ABS(SIN(OA-HIA)) 
      AH8(6,5)=SI 
      AH8(6,6)=-CO
      CALL BABT(CMH,AH8,8,6,CNH)
C     THE COVARIANCE MATRIX OF THE VECTOR (DKI,KI*GAMAI,YS,XS,DKF,
C     KF*GAMAF) IS CNH(6,6) AND HAS BEEN COMPUTED 
      CALL BABT(CMHP,AHP,4,3,AUX3)
      DO 31 I=1,2 
      DO 31 J=1,2 
   31 CKI0(I,J)=AUX3(I,J) 
      COM=DETERM(AUX3,3,AU3)
      
C      write(*,*) 'COM1: ', COM
      
      IF(COM.LT.0.) GO TO 10
      V1H=SQRT(COM) 
      DO 32 J=1,4 
   32 AHS(1,J)=AHP(1,J) 
      CALL BABT(CMHP,AHS,4,3,AUX3)
      COM=DETERM(AUX3,3,AU3)
      
C      write(*,*) 'COM2: ', COM
      
      
      IF(COM.LT.0.) GO TO 10
      V0H=SQRT(COM) 
      RMH=V0H/V1H 
      DO 33 I=1,2 
      DO 33 J=1,3 
   33 AU23(I,J)=0.
      AU23(1,1)=-VKI/VL0
      AU23(1,2)=-AU23(1,1)
      AU23(2,1)=1.
      CALL BABT(CMVP,AU23,3,2,AU2)
      COM=DETERM(AU2,2,AUX2)
      
C      write(*,*) 'COM3: ', COM
      
      
      IF(COM.LT.0.) GO TO 10
      V0V=SQRT(COM) 
      AU23(1,1)=0.
      AU23(1,2)=-VKI/VL1
      AU23(1,3)=-AU23(1,2)
      AU23(2,1)=0.
      AU23(2,3)=1.
      CALL BABT(CMVP,AU23,3,2,AU2)
      COM=DETERM(AU2,2,AUX2)
      
C      write(*,*) 'COM4: ', COM
      
      
      IF(COM.LT.0.) GO TO 10
      V1V=SQRT(COM) 
      CKI0(3,3)=AU2(1,1)
      RMV=V0V/V1V 
      VOLPH1=V1H*V1V*2.*PI*SQRT(2.*PI)
C     VOLPH1 IS THE PHASE VOLUME AFTER MONOCHROMATOR FOR AN INFINITELY
C     EXTENDED SAMPLE     
      ERM=RMH*RMV 
      DO 2 I=1,3
      DO 2 J=1,5
    2 AV(I,J)=0.
      AV(1,2)=VKI 
      AV(2,3)=VKF 
      AV(3,5)=1.
      CALL BABT(CM1V,AV,5,3,CNV)
C     THE COVARIANCE MATRIX OF THE VECTOR (KI*DELTAI,KF*DELTAF,ZS) IS 
C     CNV(3,3) AND HAS BEEN COMPUTED
      DO 11 I=1,2 
      DO 11 J=1,2 
   11 CKI(I,J)=CNH(I,J) 
      DO 12 I=1,2 
      CKI(I,3)=0. 
      CKF(I,3)=0. 
      CKF(3,I)=0. 
      CKI0(I,3)=0.
      CKI0(3,I)=0.
   12 CKI(3,I)=0. 
      CKI(3,3)=CNV(1,1) 
C     CKI(3,3) IS THE COVARIANCE MATRIX OF THE (DKI,KI*GAMAI,KI*DELTAI) 
C     VECTOR AND DESCRIBES THE EFFECTIVE MONOCHROMATOR ELLIPSOID
C     CKI0(3,3) IS THE SAME THING FOR AN INFINITELY EXTENDED SAMPLE 
      COM=DETERM(CKI,3,AUX3)
      
C      write(*,*) 'COM5: ', COM
      
      
      IF(COM.LT.0.) GO TO 10
      VOLCKI=2.*PI*SQRT(2.*PI*COM)
      DO 34 I=1,2 
      DO 34 J=1,2 
      CKF(I,J)=CNH(I+4,J+4)
   34 CONTINUE
      CKF(3,3)=CNV(2,2)
C//added/////////////////////
      COM=DETERM(CKF,3,AUX3)
      
C      write(*,*) 'COM1: ', COM
      
      
      IF(COM.LT.0.) GO TO 10
      VOLCKF=2.*PI*SQRT(2.*PI*COM)
C/////////////////////////////////
      WDTHEF=SQRT(12.*CNH(3,3)) 
      RADET=SQRT((SH(4,4)*SH(5,5)-SH(4,5)*SH(5,4))*SV(3,3)) 
      FRACV=(CMH(3,3)*CMH(4,4)-CMH(3,4)*CMH(4,3))*CMV(3,3)
      IF(FRACV.LT.0.) GO TO 10
      FRACV=RADET*SQRT(FRACV)
      IF(NSAM.GT.0) VOLSAM=WSAM*THSAM*HSAM
      IF(NSAM.EQ.0) VOLSAM=.25*PI*DIASAM*DIASAM*HSAM
      VOLEF=FRACV*VOLSAM
      ZEF=SQRT(12.*CNV(3,3))
      AREAEF=WDTHEF*ZEF 
      PATHEF=VOLEF/AREAEF 
      DO 3 I=1,4
      DO 3 J=1,6
    3 AH(I,J)=0.
      COM=COS(PHIM) 
      SIM=SIN(PHIM) 
      CO=COS(PHI) 
      SI=SIN(PHI) 
      AH(1,1)=CO
      AH(1,2)=SI
      AH(1,5)=-COM
      AH(1,6)=-SIM
      AH(2,1)=-SI 
      AH(2,2)=CO
      AH(2,5)=SIM 
      AH(2,6)=-COM
      AH(4,1)=HSQOVM*VKI
      AH(4,5)=-HSQOVM*VKF 
      CALL BABT(CNH,AH,6,4,EMMIN1)
      AV1(1,1)=1. 
      AV1(1,2)=-1.
      AV1(1,3)=0. 
      CALL BABT(CNV,AV1,3,1,X33)
      DETAS=ETAS*minute*Q0/R8LN2
      DETAS=DETAS*DETAS
      EMMIN1(2,2)=EMMIN1(2,2)+DETAS
      EMMIN1(3,3)=X33(1,1)+DETAS
      CALL INV1(EMMIN1,4,AUX4,EM) 
C     AT THIS MOMENT THE RESOLUTION MATRIX EM(4,4) AND ITS INVERSE, THE 
C     COVARIANCE MATRIX EMMIN1(4,4) HAVE BEEN COMPUTED. 

      VOLRES=DETERM(EMMIN1,4,AUX4)
      IF(VOLRES.LT.0.) GO TO 10
      VOLRES=4.*PI*PI*SQRT(VOLRES)
      DO 17 I=1,3
      DO 17 J=1,3
   17 AU3(I,J)=CNH(I,J)
      DET1=DETERM(AU3,3,AUX3) 
      IF(DET1.LT.0.) GO TO 10
      DET2=CNV(1,1)*CNV(3,3)-CNV(1,3)*CNV(3,1)
      FACT=ERM
      CO=4.*PI*PI*PM*FACT 
      FLU=1.D+13     ! FLX(VKI)
      SI=SQRT(DET1)*SQRT(DET2*2.*PI)
      YSAM=FLU*CO*SI
      
C      write(*,*) 'TRAX: ',SI*4.*PI**2*FACT,'  ',FACT

      DET3=DETERM(CKI,3,AUX3)      
      
      YSAM0=FLU*CO*VOLPH1/2./PI 
      DET1=DETERM(CNH,6,AUX6) 
      IF(DET1.LE.0.) GO TO 10
      DET2=DETERM(CNV,3,AUX3) 
      IF(DET2.LE.0.) GO TO 10
      R0=SQRT((2.*PI)**9)*SQRT(DET1)*SQRT(DET2)
      R0FLU=R0*FLU
      
      RELTR=(2*PI)**3*SQRT(DET1*DET2/DET3)
      
C      write(*,*) "RATIO=",RADET*VOLSAM/2/PI/SQRT(2*PI) 
C      write(*,*) "RELTR=",RELTR
C      R0=2.*PI*CO*PA*SQRT(DET1)*SQRT(DET2)*RADET

      DO 100 I=1,4
          DO 100  J=1,4
  100 AAA(I,J)=EM(I,J)
      RETURN 

    7 WRITE(6,98)HOMEGA
   98 FORMAT(5X,'THE REQUIRED ENERGY TRANSFER,',F8.3,
     1' MILLIEV, IS INACCESIBLE')
      IERR=1
      return
    9 WRITE(6,99)Q0
   99 FORMAT(5X,'THE REQUIRED WAVE-VECTOR TRANSFER,',F8.3,
     1' 1/A, IS INACCESIBLE')
      IERR=1
      RETURN
   10 WRITE(6,95)
   95 FORMAT(5X,'ERROR IN INVERTING A MATRIX DUE TO A TOO SMALL'/5X,
     1'MOSAIC SPREAD OR COLLIMATOR DIVERGENCE')
      IERR=1
C      WRITE(*,*) COM,FRACV,VOLRES,DET1,DET2
      
      RETURN
      END
C 

c**************************************************************
c
      SUBROUTINE SSTV
      INCLUDE 'const.inc'

C     THIS SUBROUTINE COMPUTES THE INITIAL PROBABILITY MATRICES OF THE
C     HORIZONTAL AND VERTICAL SPATIAL VARIABLES (Y0,LM,GM,LS,GS,LA,GA,YD, 
C     CSIM,CSIA),AND (Z0,ZM,ZS,ZA,ZD)
      implicit double precision (a-h,o-z) 
      COMMON/STV/SH(10,10),SH1(10,10),SV(5,5),SV1(5,5),SHP(5,5),SVP(3,3)
      COMMON /ANGLES/OM,HIM,ETAMEF,OS,HIS,OA,HIA,ETAAEF,PHI,PHIM
      COMMON /DIMENS/NSOU,WSOU,HSOU,DIASOU,WMON,THMON,HMON, 
     1WSAM,THSAM,HSAM,DIASAM,WANA,THANA,HANA
      COMMON /DETECT/NDET,NSEGDET,WDET,HDET,DIADET,ADET,
     &       SPACEDET,PHIDET
      COMMON /INDEX/NDA,NCY,NEFIX,IM,IA,ISC,NCR1,NCR2,NSAM,NOPT,NWR,NCF 
      DATA TWEL,SIXT/12.,16./
     

      DO 1 I=1,10 
      DO 1 J=1,10 
    1 SH(I,J)=0.
      DO 2 I=1,5
      DO 2 J=1,5
      SHP(I,J)=0. 
    2 SV(I,J)=0.
      DO 3 I=1,3
      DO 3 J=1,3
    3 SVP(I,J)=0. 
      IF(NSOU.GT.0)GO TO 4
      SH(1,1)=SIXT/DIASOU/DIASOU
      SV(1,1)=SH(1,1) 
      GO TO 5 
    4 SH(1,1)=TWEL/WSOU/WSOU
      SV(1,1)=TWEL/HSOU/HSOU
    5 SH(2,2)=TWEL/WMON/WMON
      SH(3,3)=TWEL/THMON/THMON
      SV(2,2)=TWEL/HMON/HMON

      DO 6 I=1,3
      DO 6 J=1,3
    6 SHP(I,J)=SH(I,J)
      SVP(1,1)=SV(1,1)
      SVP(2,2)=SV(2,2)
      IF(NSAM.GT.0)GO TO 7
      W=SIXT/DIASAM/DIASAM
      SH(4,4)=W 
      SH(5,5)=W 
      GO TO 8 
    7 SH(4,4)=TWEL/WSAM/WSAM
      SH(5,5)=TWEL/THSAM/THSAM
    8 SV(3,3)=TWEL/HSAM/HSAM
      SH(6,6)=TWEL/WANA/WANA
      SH(7,7)=TWEL/THANA/THANA
      SV(4,4)=TWEL/HANA/HANA

      W=R8LN2*R8LN2 
      SH(9,9)=W/ETAMEF/ETAMEF 
      SHP(5,5)=SH(9,9)
      SH(10,10)=W/ETAAEF/ETAAEF 
      IF(NDET) 9,9,10 
    9 SH(8,8)=SIXT/DIADET/DIADET
      SV(5,5)=SH(8,8) 

      RETURN
   10 SH(8,8)=TWEL/WDET/WDET
      SV(5,5)=TWEL/HDET/HDET
      RETURN
      END 
C
c**************************************************************
c
      SUBROUTINE TRCAN(VL0,VLI,VLC,S10,S20,NF,SHI)
      implicit double precision (a-h,o-z)
C     THIS SUBROUTINE COMPUTES THE TRANSMISSION MATRIX FOR A COARSE COLL
C     IMATOR OF CIRCULAR (NF=0) OR RECTANGULAR (NF=1) CROSS SECTION 
      double precision SHI(2,2)
      DATA UN,CON1,CON2,CON3/1.,.866,1.3333,.6667/ 
      IF(VLC.EQ.0.) VLC=0.1
      DO 4 I=1,2
      DO 4 J=1,2
    4 SHI(I,J)=0. 
      IF(VL0.LT.(VLI+VLC))GO TO 15
      ALFA=VLI/VL0
      BETA=VLC/VL0
      EPS=ALFA+BETA
      S1=S10
      S2=S20
      IF(NF)10,10,5 
    5 S1=CON1*S10 
      S2=CON1*S20 
   10 IF(VLC.EQ.0.)GO TO 30 
      Y0K=ABS(S1*EPS-S2*ALFA)/BETA
      Y0J=(EPS*S1+S2*ALFA)/BETA 
      AI0=2.*(S1/ALFA+S2/EPS)*(Y0J-Y0K)-(UN/ALFA-UN/EPS)*(Y0J**2-Y0K**2)
     1+4.*Y0K*S2/EPS
      A11=CON2*S2*Y0K**3/EPS+CON3*(S2/EPS+S1/ALFA)*(Y0J**3-Y0K**3)- 
     1.5*(UN/ALFA-UN/EPS)*(Y0J**4-Y0K**4) 
      A22=CON3*((S1/ALFA)**3+(S2/EPS)**3)*(Y0J-Y0K)-(S1**2*(UN-ALFA)/
     2ALFA**3-S2**2*(UN-EPS)/EPS**3)*(Y0J**2-Y0K**2)+CON3*(S1*(UN-ALFA
     3)**2/ALFA**3+S2*(UN-EPS)**2/EPS**3)*(Y0J**3-Y0K**3)-.1667*((UN/
     4ALFA-UN)**3-(UN/EPS-UN)**3)*(Y0J**4-Y0K**4)+CON2*Y0K*(S2/EPS)**
     53+CON2*S2*(UN-EPS)**2*(Y0K/EPS)**3
      A12=.5*((S1/ALFA)**2-(S2/EPS)**2)*(Y0J**2-Y0K**2)-CON3*(S1*(UN-
     1ALFA)/ALFA**2+S2*(UN-EPS)/EPS**2)*(Y0J**3-Y0K**3)+.25*((UN/ALFA-
     2UN)**2-(UN/EPS-UN)**2)*(Y0J**4-Y0K**4)-CON2*S2*(UN-EPS)/EPS**2*
     3Y0K**3
      X=A11/AI0
      Y=A22/AI0 
      Z=A12/AI0 
      D=X*Y-Z**2
      SHI(1,1)=Y/D
      SHI(2,2)=X/D
      SHI(1,2)=-Z/D 
      SHI(2,1)=SHI(1,2) 
      GO TO 25
   30 XMO=S1*S1/12.
      SHI(1,1)=(UN-ALFA)**2/XMO 
      SHI(2,2)=ALFA**2/XMO
      SHI(1,2)=ALFA*(UN-ALFA)/XMO 
      SHI(2,1)=SHI(1,2) 
      GO TO 25
   15 WRITE(6,20) VL0,VLI,VLC
   20 FORMAT(2X,'WRONG INPUT DATA: ',3(G10.4,1x),
     &' THE PRESENCE OF THIS COLLIMATOR IS IGNORED') 
   25 RETURN
      END 
C
                                                      
c**************************************************************

      SUBROUTINE TCR
      INCLUDE 'const.inc'
C     COMPUTES THE PROBABILITY MATRICES OF THE SPATIAL VARIABLES AS 
C     MODIFIED BY THE BRAGG CONSTRAINTS, NEUTRON GUIDE, SOLLER
C     COLLIMATORS AND COARSE COLLIMATORS OR SLITS
      implicit double precision (a-h,o-z) 
      double precision DH(4,10),GH(10,8),TV(2,5),ACO(4,4),
     +AU4(4,4),F(2,2), 
     1AU8(8,8),DV1(4,5),CAN(10,10),VCAN(5,5),SHI(2,2),AH8(8,8),TV1(1,3) 
     2,CANP(5,5),VCANP(3,3),DHP(2,5),ACOP(2,2),DVP(2,3),GHP(5,4),F1(1)
      COMMON/STV/SH(10,10),SH1(10,10),SV(5,5),SV1(5,5),SHP(5,5),SVP(3,3)
      COMMON /COARSE/NFM,VLCANM,VLSM,HDM1,HDM2,VDM1,VDM2, 
     1               NFS,VLCANS,VLMS,HDS1,HDS2,VDS1,VDS2, 
     2               NFA,VLCANA,VLSA,HDA1,HDA2,VDA1,VDA2, 
     3               NFD,VLCAND,VLAD,HDD1,HDD2,VDD1,VDD2
      COMMON /VOLS/VOLSAM,VOLEF,FRACV,RADET,VOLCKI,VOLCKF,CKI(3,3),
     1CKF(3,3)
      COMMON /QEF/PATHEF,WDTHEF,AREAEF,ZEF
      COMMON /SETP/VL0,VL1,VL2,VL3,VLMIN,VLMAX 
      COMMON/MONOCH/TETAM,HIMON,ROHM,ROVM,ETAM,NGRADM,PM,PM0,RKINM,RLAMM
      COMMON/ANALYS/TETAA,HIANA,ROHA,ROVA,ETAA,NGRADA,PA,PA0,RKINA,RLAMA
      COMMON /VERTIC/ETVM,ETVMEF,ANRM,ETVA,ETVAEF,ANRA
      COMMON /ANGLES/OM,HIM,ETAMEF,OS,HIS,OA,HIA,ETAAEF,PHI,PHIM
      COMMON /SOLLER/NGUIDE,GAMACR,ALPHA(4),BETA(4)
      COMMON /CTR/DV(5,5),BM,BA,AHS(3,4)
      COMMON /COVAR/CMH(8,8),CMV(5,5),CMHP(4,4),CMVP(3,3) 
      COMMON /MEAN/VKI,SLAMDI,EI0,VKF,SLAMDF,EF0,Q0,HOMEGA,GRAD(3),DIRQ 
     *,GRADI(3),TAU,DIRTAU,DIRTD,DIRQD,EPSLN,QSMALL,Q,DIRQT 
C     DATA GAMACR/.00167/


  101 FORMAT(' TCR')
 
      DO 111 I=1,10 
      DO 111 J=1,10 
      SH1(I,J)=SH(I,J)
  111 CAN(I,J)=0.
      DO 112 I=1,5
      DO 112 J=1,5
      CANP(I,J)=0.
      SV1(I,J)=SV(I,J)
  112 VCAN(I,J)=0.
      DO 110 I=1,2
      DO 110 J=1,5
  110 DHP(I,J)=0. 
      DO 113 I=1,3
      DO 113 J=1,3
  113 VCANP(I,J)=0. 
      DO 116 I=1,3
      DO 116 J=1,4
  116 AHS(I,J)=0. 

1234  FORMAT (a,10(1x,G10.4))
      IF(NFM.LT.0)GO TO 15
C      write(*,1234) 'TCR ',VL0,VLCANM,VLSM,HDM1,HDM2,NFM,SHI
      CALL TRCAN(VL0,VLCANM,VLSM,HDM1,HDM2,NFM,SHI)
C      write(*,*) 'TCR - TRCAN 1 O.K.'      
       
      CO=COS(OM+HIM)
      SI=SIN(OM+HIM)
      CAN(1,1)=SHI(1,1)
      CAN(1,2)=SHI(1,2)*SI
      CAN(1,3)=-SHI(1,2)*CO 
      CAN(2,2)=SHI(2,2)*SI*SI 
      CAN(2,3)=-SHI(2,2)*CO*SI
      CAN(3,3)=SHI(2,2)*CO*CO 
      CAN(2,1)=CAN(1,2)
      CAN(3,1)=CAN(1,3)
      CAN(3,2)=CAN(2,3)

      CALL TRCAN(VL0,VLCANM,VLSM,VDM1,VDM2,NFM,SHI) 
      VCAN(1,1)=SHI(1,1)
      VCAN(2,2)=SHI(2,2)
      VCAN(1,2)=SHI(1,2)
      VCAN(2,1)=VCAN(1,2)
      DO 114 I=1,3
      DO 114 J=1,3
      VCANP(I,J)=VCAN(I,J)
  114 CANP(I,J)=CAN(I,J)
   15 IF(NFS.LT.0)GO TO 20
      CO=COS(OM-HIM)
      SI=SIN(OM-HIM)
      A=COS(OS+HIS)
      B=SIN(OS+HIS)
C      write(*,1234) 'TCR ',VL1,VLCANS,VLMS,HDS1,HDS2,NFS,SHI
      CALL TRCAN(VL1,VLCANS,VLMS,HDS1,HDS2,NFS,SHI) 
C      write(*,*) 'TCR - TRCAN 2 O.K.'      
      
      
      CAN(2,2)=CAN(2,2)+SHI(1,1)*SI*SI
      CAN(3,3)=CAN(3,3)+SHI(1,1)*CO*CO
      CAN(2,3)=CAN(2,3)+SHI(1,1)*SI*CO
      CAN(3,2)=CAN(2,3)
      CANP(4,4)=SHI(2,2)
      CAN(4,4)=SHI(2,2)*A*A 
      CAN(5,5)=SHI(2,2)*B*B 
      CAN(4,5)=SHI(2,2)*A*B 
      CAN(5,4)=CAN(4,5)
      CAN(2,4)=CAN(2,4)-SHI(1,2)*SI*A 
      CAN(4,2)=CAN(2,4)
      CAN(2,5)=CAN(2,5)-SHI(1,2)*SI*B 
      CAN(5,2)=CAN(2,5)
      CAN(3,4)=CAN(3,4)-SHI(1,2)*CO*A 
      CAN(4,3)=CAN(3,4)
      CAN(3,5)=CAN(3,5)-SHI(1,2)*CO*B 
      CAN(5,3)=CAN(3,5)
      CALL TRCAN(VL1,VLCANS,VLMS,VDS1,VDS2,NFS,SHI) 
      VCAN(2,2)=VCAN(2,2)+SHI(1,1)
      VCAN(3,3)=SHI(2,2)
      VCAN(2,3)=SHI(1,2)
      VCAN(3,2)=VCAN(2,3)
      DO 115 I=1,3
      DO 115 J=1,3
      CANP(I,J)=CAN(I,J)
  115 VCANP(I,J)=VCAN(I,J)
   20 IF(NFA.LT.0)GO TO 25
      CO=COS(OA+HIA)
      SI=SIN(OA+HIA)
      A=COS(OS-HIS)
      B=SIN(OS-HIS)
C      write(*,1234) 'TCR ',VL2,VLCANA,VLSA,HDA1,HDA2,NFA,SHI
      CALL TRCAN(VL2,VLCANA,VLSA,HDA1,HDA2,NFA,SHI) 
C      write(*,*) 'TCR - TRCAN 3 O.K.'      
      
      
      CAN(4,4)=CAN(4,4)+SHI(1,1)*A*A
      CAN(5,5)=CAN(5,5)+SHI(1,1)*B*B
      CAN(4,5)=CAN(4,5)-SHI(1,1)*A*B
      CAN(5,4)=CAN(4,5)
      CAN(6,6)=SHI(2,2)*SI*SI 
      CAN(7,7)=SHI(2,2)*CO*CO 
      CAN(6,7)=-SHI(2,2)*SI*CO
      CAN(7,6)=CAN(6,7)
      CAN(4,6)=CAN(4,6)+SHI(1,2)*SI*A 
      CAN(6,4)=CAN(4,6)
      CAN(4,7)=CAN(4,7)-SHI(1,2)*A*CO 
      CAN(7,4)=CAN(4,7)
      CAN(5,6)=CAN(5,6)-SHI(1,2)*SI*B 
      CAN(6,5)=CAN(5,6) 
      CAN(5,7)=CAN(5,7)+SHI(1,2)*CO*B 
      CAN(7,5)=CAN(5,7)
      CALL TRCAN(VL2,VLCANA,VLSA,VDA1,VDA2,NFA,SHI) 
      VCAN(3,3)=VCAN(3,3)+SHI(1,1)
      VCAN(4,4)=VCAN(4,4)+SHI(2,2)
      VCAN(3,4)=SHI(1,2)
      VCAN(4,3)=VCAN(3,4)
   25 IF(NFD.LT.0)GO TO 30
      CALL TRCAN(VL3,VLCAND,VLAD,HDD1,HDD2,NFD,SHI)
C      write(*,*) 'TCR - TRCAN 4 O.K.'      
      
      
      CO=COS(OA-HIA)
      SI=SIN(OA-HIA)
      CAN(6,6)=CAN(6,6)+SHI(1,1)*SI*SI
      CAN(7,7)=CAN(7,7)+SHI(1,1)*CO*CO
      CAN(6,7)=CAN(6,7)+SHI(1,1)*SI*CO
      CAN(7,6)=CAN(6,7)
      CAN(8,8)=SHI(2,2)
      CAN(6,8)=-SHI(1,2)*SI 
      CAN(8,6)=CAN(6,8)
      CAN(7,8)=-SHI(1,2)*CO 
      CAN(8,7)=CAN(7,8)
      CALL TRCAN(VL3,VLCAND,VLAD,VDD1,VDD2,NFD,SHI)
      VCAN(4,4)=VCAN(4,4)+SHI(1,1)
      VCAN(5,5)=SHI(2,2)
      VCAN(4,5)=SHI(1,2)
      VCAN(5,4)=VCAN(4,5)
   30 CALL SUM(SH,CAN,10,SH1) 
      CALL SUM(SV,VCAN,5,SV1) 
      CALL SUM(SHP,CANP,5,SHP)
      CALL SUM(SVP,VCANP,3,SVP) 
C     THE PRESENCE OF COARSE COLLIMATORS OR SLITS WAS ACCOUNTED FOR

C      write(*,*) 'TCR - COLIM O.K.'      
      DO 1 I=1,4
      DO 1 J=1,10 
    1 DH(I,J)=0.
      DO 2 I=1,5
      DO 2 J=1,5
    2 DV(I,J)=0.
      DO 3 I=1,10 
      DO 3 J=1,8
    3 GH(I,J)=0.
      GH(2,1)=1.
      GH(3,2)=1.
      GH(4,3)=1.
      GH(5,4)=1.
      GH(6,5)=1.
      GH(7,6)=1.
      GH(9,7)=1.
      GH(10,8)=1. 
      CO=COS(OM+HIM)
      SI=SIN(OM+HIM)

      SSS=ABS(SI)/SI
      A=1./VL0
      DH(1,1)=-A
      DH(1,3)=-A*CO 
      DH(1,2)=A*SI
      DV(1,1)=-A
      DV(1,2)=A 
      A=1./VL1
      GH(1,1)=SI-2.*VL0*ROHM*SSS 
      GH(1,2)=-CO-2.*VL0*BM*ROHM*SSS 
      CO=COS(OM-HIM)
      SI=SIN(OM-HIM)
      DH(2,3)=A*CO
      DH(2,2)=A*SI
      DV(2,2)=-A
      DV(2,3)=A 
      GH(1,1)=GH(1,1)+VL0*A*SI
      GH(1,2)=GH(1,2)+VL0*A*CO
      CO=COS(OS+HIS)
      SI=SIN(OS+HIS)
      DH(2,5)=A*SI
      DH(2,4)=A*CO
      DHP(2,4)=A
      A=VL0*A
      GH(1,3)=A*CO
      GH(1,4)=A*SI
      CO=COS(OS-HIS)
      SI=SIN(OS-HIS)
      A=1./VL2
      DH(3,5)=A*SI
      DH(3,4)=-A*CO 
      DV(3,3)=-A
      DV(3,4)=A 
      GH(8,3)=A*VL3*CO
      GH(8,4)=-A*VL3*SI
      CO=COS(OA+HIA)
      SI=SIN(OA+HIA)

      SSS=ABS(SI)/SI
      DH(3,7)=-A*CO 
      DH(3,6)=A*SI
      GH(8,5)=-VL3*A*SI+2.*VL3*ROHA*SSS
      GH(8,6)=VL3*A*CO+2.*VL3*BA*ROHA*SSS
      A=1./VL3
      CO=COS(OA-HIA)
      SI=SIN(OA-HIA)
      DH(4,7)=A*CO
      DH(4,6)=A*SI
      DH(4,8)=A 
      DV(4,4)=-A
      DV(4,5)=A 
      DV(5,3)=1.
      GH(8,5)=GH(8,5)-SI
      GH(8,6)=GH(8,6)-CO
      GH(1,7)=-2.*VL0 
      IF(ETAM.EQ.0.) GH(1,7)=GH(1,7)*SIN(HIM)*COS(OM)/ABS(SIN(OM-HIM))
      GH(8,8)=2.*VL3
      IF(ETAA.EQ.0.) GH(8,8)=GH(8,8)*SIN(HIA)*COS(OA)/ABS(SIN(OA-HIA))
      DO 16 I=1,5 
      DO 16 J=1,4 
   16 GHP(I,J)=GH(I,J)
      GHP(1,3)=VL0/VL1
      GHP(1,4)=GH(1,7)
      A=R8LN2*R8LN2 

C      write(*,*) 'TCR - ETAM',ETVMEF,ETAMEF,SIN(OM+HIM)      
      IF(ETAM.NE.0.)F(1,1)=A/ETVMEF/ETVMEF
      IF(ETAM.EQ.0.)F(1,1)=A/ETAMEF/ETAMEF*ABS(SIN(OM-HIM)/SIN(OM+HIM)) 
C      write(*,*) 'TCR - ETAA',ETVAEF,ETAAEF,SIN(OA+HIA)      
      IF(ETAA.NE.0.)F(2,2)=A/ETVAEF/ETVAEF
      IF(ETAA.EQ.0.)F(2,2)=A/ETAAEF/ETAAEF*ABS(SIN(OA-HIA)/SIN(OA+HIA)) 
C      write(*,*) 'TCR - test OK'      
      F(1,2)=0. 
      F(2,1)=0. 
      DO 4 J=1,5
      TV(1,J)=.5*(DV(1,J)-DV(2,J))/SIN(OM)
    4 TV(2,J)=.5*(DV(3,J)-DV(4,J))/SIN(OA)
C      write(*,*) 'TCR - test OK'      

      TV(1,2)=TV(1,2)-ROVM*ABS(OM)/OM*ABS(COS(HIM)) 
      TV(2,4)=TV(2,4)-ROVA*ABS(OA)/OA*ABS(COS(HIA)) 
      CALL BTAB(F,TV,2,5,VCAN)
      CALL SUM(SV1,VCAN,5,SV1)
C      write(*,*) 'TCR - test OK'      
      DO 18 J=1,3 
   18 TV1(1,J)=TV(1,J)
      F1(1)=F(1,1)
      CALL BTAB(F1,TV1,1,3,VCANP) 
      CALL SUM(SVP,VCANP,3,SVP) 
C      write(*,*) 'TCR - test OK'      
      DO 5 I=1,4
      DO 5 J=1,4
    5 ACO(I,J)=0. 
      CO=0. 
      SI=3./(GAMACR*SLAMDI)**2
      DO 6 I=1,4
      CO=CO+ALPHA(I)
    6 IF(ALPHA(I).NE.0.) 
     &  ACO(I,I)=A/ALPHA(I)/ALPHA(I)/minute/minute 
      IF(NGUIDE.NE.0) ACO(1,1)=ACO(1,1)+SI
      CO=CO+NGUIDE
C      write(*,*) 'TCR - 6 OK'      
      IF(CO.EQ.0.) GO TO 7
      CALL BTAB(ACO,DH,4,10,CAN)
      CALL SUM(SH1,CAN,10,SH1)
      CO=0. 
      ACOP(1,2)=0.
      ACOP(2,1)=0.
      DO 11 I=1,2 
      CO=CO+ALPHA(I)
   11 ACOP(I,I)=ACO(I,I)
      CO=CO+NGUIDE
C      write(*,*) 'TCR - 6 OK'      
      IF(CO.EQ.0.) GO TO 7
      DO 13 I=1,2 
      DO 13 J=1,3 
   13 DHP(I,J)=DH(I,J)
      CALL BTAB(ACOP,DHP,2,5,VCAN)
      CALL SUM(SHP,VCAN,5,SHP)
    7 DO 9 I=1,4
      DO 9 J=1,4
    9 ACO(I,J)=0. 
C      write(*,*) 'TCR - 7 OK'      
      CO=0. 
      DO 10 I=1,4 
      CO=CO+BETA(I) 
   10 IF(BETA(I).NE.0.) 
     &  ACO(I,I)=A/BETA(I)/BETA(I)/minute/minute
      IF(NGUIDE.NE.0) ACO(1,1)=ACO(1,1)+SI
      CO=CO+NGUIDE
      IF(CO.EQ.0.)GO TO 8 
      DO 12 I=1,4
      DO 12 J=1,5
   12 DV1(I,J)=DV(I,J)
      CALL BTAB(ACO,DV1,4,5,VCAN) 
      CALL SUM(SV1,VCAN,5,SV1)
      CO=0. 
      DO 14 I=1,2 
      CO=CO+BETA(I) 
   14 ACOP(I,I)=ACO(I,I)
      CO=CO+NGUIDE
      IF(CO.EQ.0.) GO TO 8
      DO 17 I=1,2 
      DO 17 J=1,3 
   17 DVP(I,J)=DV(I,J)
      CALL BTAB(ACOP,DVP,2,3,VCANP) 
      CALL SUM(SVP,VCANP,3,SVP) 
    8 CONTINUE
C      write(*,*) 'TCR - 8 OK'      
      DO 117 J=1,4
      AHS(3,J)=GHP(1,J) 
  117 AHS(2,J)=-VKI/VL0*GHP(1,J)
      AHS(2,1)=AHS(2,1)+VKI/VL0*SIN(OM+HIM) 
      AHS(2,2)=AHS(2,2)-VKI/VL0*COS(OM+HIM) 
      CALL BTAB(SHP,GHP,5,4,ACO)
C      write(*,*) 'TCR - BTAB'      
      CALL INV1(ACO,4,AU4,CMHP) 
C      write(*,*) 'TCR - INV ACO'      
C     SOLLER COLLIMATORS AND NEUTRON GUIDE HAVE BEEN ACCOUNTED FOR
      CALL BTAB(SH1,GH,10,8,AH8)
      CALL INV1(AH8,8,AU8,CMH)
C      write(*,*) 'TCR - INV AH8'      
      CALL INV1(SV1,5,VCAN,CMV) 
C      write(*,*) 'TCR - INV SV1'      
      CALL INV1(SVP,3,VCANP,CMVP)
C      write(*,*) 'TCR - INV SVP'      
C     BRAGG CONSTRAINTS HAVE BEEN ACCOUNTED FOR
C     CMH AND CMV ARE THE MODIFIED COVARIANCE MATRICES OF THE 
C     INDEPENDENT VARIABLES (HORIZONTAL AND VERTICAL, RESPECTIVELY) 
      RETURN
      END 
         

c**************************************************************
c
      SUBROUTINE HIRANG(OB,HII,HI)
C     MAKES SURE THAT THE CRYSTAL INCLINATION ANGLE IS IN THE CORRECT
C     RANGE (MINUS TETA BRAGG TO PI MINUS TETA BRAGG FOR POSITIVE TETA,
C     MINUS PI MINUS TETA BRAGG TO MINUS TETA BRAGG FOR NEGATIVE TETA)
      implicit double precision (a-h,o-z)
      INCLUDE 'const.inc'
      HI=deg*HII
      IF(OB.LT.0.) GO TO 1
      IF(HI.LT.(-OB)) HI=HI+PI
      IF(HI.GT.(PI-OB)) HI=HI-PI
      GO TO 2
    1 IF(HI.LT.(-PI-OB)) HI=HI+PI
      IF(HI.GT.(-OB)) HI=HI-PI
    2 HII=HI/deg
      RETURN
      END    

C
        SUBROUTINE VANAD(DVN,YVN)
C***********************************************************************
C   returns Vanad scan width and intensity (in rel. units)
C   added by J.S. 3/6/1997
C***********************************************************************
        IMPLICIT REAL*8 (A-H,O-Z)
        DIMENSION A(4,4)        
        COMMON /MATRIX/A
        COMMON /EM/EMMIN1(4,4),VOLRES
        COMMON /ABNORM/R0,R0CN,FLU,R0FLU        

        HUITLOG2=8.D0*LOG(2.D0)        
        DVN=(A(2,4)*A(1,1)-A(1,4)*A(1,2))**2/(A(1,1)*(A(2,2)*A(1,1)
     1        -A(1,2)**2))
        DVN=SQRT(HUITLOG2/(A(4,4)-A(1,4)**2/A(1,1)-DVN))
        
        YVN=R0/SQRT(EMMIN1(4,4))        ! intensity in rel. units
        
        
c        write(*,*) 'R0, EM, YVN : ',R0, EMMIN1(4,4), YVN
        
        END
        
C
        REAL*4 FUNCTION OPTV(RO)
C***********************************************************************
C   returns value to be minimized when optimizing bending radii
C   added by J.S. 3/6/1997
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      
      PARAMETER(EPS=1.D-15)
      REAL*4 RO(4),B(4)        
      COMMON DUMMY(42),ROMH,ROMV,ROAH,ROAV,SDI,SHI      
      COMMON /ABNORM/R0,R0CN,FLU,R0FLU        
      COMMON/MONOCH/TETAM,HIMON,ROHM,ROVM,ETAM,NGRADM,PM,PM0,RKINM,RLAMM
      COMMON/ANALYS/TETAA,HIANA,ROHA,ROVA,ETAA,NGRADA,PA,PA0,RKINA,RLAMA
      INTEGER*4 OPTPAR,OPTMERIT
      REAL*8 OPTEV
      COMMON /MCOPTIM/ OPTPAR,OPTMERIT,OPTEV
        B(1)=ROHM
        B(2)=ROVM
        B(3)=ROHA
        B(4)=ROVA
        ROHM=RO(1)/100
        ROVM=RO(2)/100
        ROHA=RO(3)/100
        ROVA=RO(4)/100

        CALL THRAX(0)
        CALL VANAD(DVN,YVN)
c10      format('OPTV ',6(2x,G12.5)) 
c        write(*,*) (RO(I),I=1,4),DVN,YVN
c        pause
      OPTV=0.
      IF (OPTMERIT.EQ.1) THEN
        IF (ABS(YVN).GE.EPS) OPTV=1./YVN/1E10
      ELSE IF (OPTMERIT.EQ.2) THEN
        OPTV=DVN
      ELSE IF (OPTMERIT.EQ.3) THEN
        IF (ABS(YVN).GE.EPS) OPTV=DVN/YVN/1E10
      ELSE IF (OPTMERIT.EQ.4) THEN
        IF (ABS(YVN).GE.EPS) OPTV=DVN**2/YVN/1E10
      ENDIF  

        ROHM=B(1)
        ROVM=B(2)
        ROHA=B(3)
        ROVA=B(4)        
c        N=4
      END
        
