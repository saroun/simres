!//////////////////////////////////////////////////////////////////////
!////  $Id: ness_crystal_mosaic.f,v 1.14 2019/08/16 17:16:26 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.14 $
!////     $Date: 2019/08/16 17:16:26 $
!////////////////////////////////////////////////////////////////////////
C////
C////  Subroutines describing objects - CRYSTAL
C////  Implements:
C////  - reflectivity by analytical mosaic model (Bacon & Lowde, Sears)
C////  - Debye-Waller factor and absorption (Freund)
C////  - Reflectivity of elastically bent perfect crystals (Kulda)
C////  - mosaicity distribution function
C////    + corresponding cumulative probability and its inverse
C////////////////////////////////////////////////////////////////////////

C------------------------------------------------------------------
      SUBROUTINE MOSAIC_TEST
C test ,osaic model: return parameters of reflectivity and absorption
C------------------------------------------------------------------
      USE CRYSTALS
      USE CONSTANTS
      IMPLICIT NONE

      TYPE (CRYSTAL) :: CR
      REAL*8 GETMI

      REAL*8 LAMBDA,ETA,THICK,CHI,R0,fwhm,Integ,weff,wmos,GETQKIN,reff
      REAL*8 KI,QK,MI,TEMP
      INTEGER*4 ierr
      CHARACTER*8 CRNAME,SDUM
      DATA CRNAME/'        '/
      SAVE CRNAME

1     FORMAT('Crystal name [',a8,']: ',$)
2     FORMAT(a8)
3     FORMAT('ETA['']  D[mm]  chi[deg]  T[K]: ',$)
4     FORMAT('ki = ',$)
5     FORMAT( 2x,'ki          ',
     *        'fwhm        ',
     *        'R           ',
     *        'Integral    ',
     *        'Qkin.       ',
     *        'mi [cm^-1]  ',
     *        'sigma [b]   ')


15    FORMAT(10(F10.6,2x))
6     FORMAT('ETA['']  D[mm]  chi[deg]  T[K]: ',4(2x,F10.6))

50    WRITE(*,*)
      WRITE(*,1) CRNAME
      READ(*,2) SDUM
      CALL MakeUpCase(SDUM)
      if ((SDUM(1:3).EQ.'END').OR.(SDUM(1:3).EQ.'QUI')) THEN
         goto 200
      else if (SDUM.NE.'        ') THEN
         CRNAME=SDUM
         call SET_CRYST_PARAM(CR,CRNAME,.false.)
      endif
      WRITE(*,3)
      READ(*,*) ETA,THICK,CHI,TEMP
      CR%HMOS=ETA*PI/180/60/R8LN2
      CR%FRAME%SIZE(3)=THICK
      CR%CHI=CHI*PI/180
      IERR=0
      WRITE(*,*)
      WRITE(*,*) CRNAME
      WRITE(*,6) ETA,THICK,CHI,TEMP
      WRITE(*,5)
      DO WHILE (IERR.EQ.0)
        write(*,4)
        read(*,*,err=100,iostat=ierr) KI
        if (KI.eq.0) goto 200
        lambda=2*pi/KI
        QK= GETQKIN(CR,lambda)
        MI= GETMI(CR,lambda,TEMP)
        call getrefpar(CR,lambda,QK,MI,R0,Integ,wmos,fwhm)
        if (R0.GT.0) THEN
          weff=Integ/SQRT2PI/R0*R8LN2
          reff=Integ/SQRT2PI/wmos*R8LN2
         else
          weff=0
          reff=0
        endif
        weff=weff*180*60/PI
        wmos=wmos*180*60/PI
        fwhm=fwhm*180*60/PI
        write(*,15) KI,fwhm,R0,reff,QK,MI,MI*CR%VOL
      ENDDO
100   continue
      GOTO 50
200   CONTINUE
      END SUBROUTINE MOSAIC_TEST

CxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCX
C--------------------------------------------------------
      REAL*8 FUNCTION GETQKIN(CR,lambda)
C return kinematical reflectivity
C--------------------------------------------------------
      USE COMPONENTS
      USE CONSTANTS
      IMPLICIT NONE
      TYPE (CRYSTAL) :: CR
      REAL*8 STHB,lambda,CTHB

      STHB=lambda/2/CR%dhkl
      IF (STHB.GE.9.99999D-1) THEN
         GETQKIN=0
        RETURN
      ENDIF
      CTHB=SQRT(1.D0-STHB**2)
c      write(*,*) 'GETQKIN: ',CR%QML,CR%DHKL,STHB
      GETQKIN=CR%QML*CR%DHKL*STHB**2/PI/CTHB

      END FUNCTION GETQKIN

C--------------------------------------------------------
      REAL*8 FUNCTION GETREFDYN(CR,lambda)
C return peak reflectivity of an elastically bent single crystal
C--------------------------------------------------------
      USE CRYSTALS
      USE CONSTANTS
      IMPLICIT NONE
      TYPE (CRYSTAL) :: CR
      REAL*8 Z,B0,BT,grad
      INTEGER*4 I,J
      REAL*8 STHB,lambda,CTHB,THB,QK,sigma,ki(3),kg(3)
      REAL*8 GETDW

      STHB=lambda/2/CR%dhkl
      IF (STHB.GE.9.99999D-1) THEN
        GETREFDYN=0
        RETURN
      ENDIF

      CTHB=SQRT(1.D0-STHB**2)
      QK=CR%QML*CR%DHKL*STHB**2/PI  !  = Qkin*cos(thb)
      call DWFACT(CR,300.D0,lambda,B0,BT)
      QK=QK*GETDW(B0,BT,CR%DHKL)/10.D0
      THB=ASIN(STHB)
      ki(1)=-cos(THB-CR%CHI)
      ki(2)=0.
      ki(3)=-sin(THB-CR%CHI)
      kg(1)=-cos(THB+CR%CHI)
      kg(2)=0.
      kg(3)=+sin(THB+CR%CHI)
      grad=0.
      DO I=1,3
        Z=0.
        DO J=1,3
          Z=Z+CR%DG_DR(I,J)*KI(J)
        ENDDO
        grad=grad+Z*KG(I)
      ENDDO
      grad=abs(grad)/CR%GTOT
      if (grad/QK.lt.1D-2) then
         Z=1.D0
      else
         sigma=QK/grad
         Z=1.D0-exp(-sigma)
      endif
      GETREFDYN=Z
      END FUNCTION GETREFDYN

C-------------------------------------------------------------
      REAL*8 FUNCTION GETDTH(CR)
C Get the spread of Bragg angles due to the elastic deformation
C-------------------------------------------------------------
      USE CRYSTALS
      use TRACELIB
      IMPLICIT NONE
      TYPE (CRYSTAL) :: CR
      REAL*8 Z,grad
      INTEGER*4 I,J
      REAL*8 ki(3),kg(3),R(3),TIN,TOUT,STHB,CTHB,THB
      THB=ASIN(CR%LAMBDA/2.D0/CR%DHKL)
      STHB=sin(abs(THB))
      IF (STHB.GE.9.99999D-1) THEN
         GETDTH=0.D0
        RETURN
      ENDIF

      ki(1)=-cos(THB-CR%CHI)
      ki(2)=0.
      ki(3)=-sin(THB-CR%CHI)
      kg(1)=-cos(THB+CR%CHI)
      kg(2)=0.
      kg(3)=+sin(THB+CR%CHI)
      grad=0.
      DO I=1,3
        Z=0.
        DO J=1,3
          Z=Z+CR%DG_DR(I,J)*KI(J)
        ENDDO
        grad=grad+Z*KG(I)
      ENDDO
c      write(*,*) CR%FRAME%NAME,grad,cos(THB),CR%GTOT
      CTHB=SQRT(1.D0-STHB**2)
      if (CR%GTOT.GE.1.D-10) THEN
         grad=abs(grad)/CTHB/CR%GTOT
      else
         grad=0.D0
      endif
      R(1)=0.
      R(2)=0.
      R(3)=0.
      call BORDER_BOX(R,KI,CR%FRAME%SIZE,TIN,TOUT)
c      write(*,*) TIN,TOUT
      GETDTH=grad*(TOUT-TIN)
      END FUNCTION GETDTH

C--------------------------------------------------------------------
      real(kind(1.D0)) FUNCTION GETEFFMOS(CR)
C Get effective mosaicity
C = change of refl. plane orientation along incident beam path
C--------------------------------------------------------------------
      USE CRYSTALS
      IMPLICIT NONE
      TYPE (CRYSTAL) :: CR
      real(kind(1.D0)) :: Z,G01,G02
      INTEGER :: I,J
      real(kind(1.D0)) :: ri(3),rf(3),G1(3),G2(3),CTHB,SCHI,CCHI,THB
      
      THB=ASIN(CR%LAMBDA/2.D0/CR%DHKL)
      CTHB=cos(THB-CR%CHI)
      ri(1)=0.5*CTHB*CR%FRAME%SIZE(3)
      ri(2)=0.
      ri(3)=0.5*CR%FRAME%SIZE(3)
      rf(1)=-0.5*CTHB*CR%FRAME%SIZE(3)
      rf(2)=0.
      rf(3)=-0.5*CR%FRAME%SIZE(3)
      DO I=1,3
        G1(I)=0.
        G2(I)=0.
        DO J=1,3
          G1(I)=G1(I)+CR%DG_DR(I,J)*ri(J)
          G2(I)=G2(I)+CR%DG_DR(I,J)*rf(J)
        ENDDO
      ENDDO
      SCHI=sin(CR%CHI)
      CCHI=SQRT(1.D0-SCHI**2)
      G1(1)=G1(1)+CR%GTOT*SCHI
      G1(3)=G1(3)-CR%GTOT*CCHI
      G2(1)=G2(1)+CR%GTOT*SCHI
      G2(3)=G2(3)-CR%GTOT*CCHI
      Z=0.
      DO I=1,3
        Z=Z+G1(I)*G2(I)
      ENDDO
      IF (ABS(Z).LT.1.D-10) THEN
            Z=0.D0
      ELSE
            G01=SQRT(G1(1)**2+G1(3)**2)
            G02=SQRT(G2(1)**2+G2(3)**2)
            if (G01*G02>1.D-10) then
                  Z=ACOS(Z/G01/G02)
            endif
      ENDIF
      GETEFFMOS=Z
      END FUNCTION GETEFFMOS

!--------------------------------------------------------
      subroutine DWFACT(CR,T,lambda,B0,BT)
! return  Debye-Waller coefficients B0,BT
! for cubic crystals, Debye model
!--------------------------------------------------------
      USE CRYSTALS
      USE CONSTANTS
      use SIMATH
      IMPLICIT NONE
      TYPE (CRYSTAL) :: CR
      real(kind(1.D0)),intent(in) :: T,lambda
      real(kind(1.D0)),intent(out) :: B0,BT
      real(kind(1.D0)) :: x,phi
      B0=2872.556/CR%A/CR%thetaD
      x=CR%thetaD/max(T,1D-1)
      phi=DW_PHI(x)
      BT=4*B0*phi/x**2
      END subroutine DWFACT

!--------------------------------------------------------
      real(kind(1.D0)) FUNCTION GETDW(B0,BT,dhkl)
! return  Debye-Waller factor, exp(-2W)
! for given B0,BT
!--------------------------------------------------------
      IMPLICIT NONE
      real(kind(1.D0)) :: dhkl,B0,BT
      GETDW=exp(-0.5D0*(B0+BT)/(dhkl**2))
      END FUNCTION GETDW

C--------------------------------------------------------
      REAL*8 FUNCTION GETMI(CR,lambda,T)
C//  calculates absorption according to A.Freund
C//  (Nucl.Inst.Meth. 213 (1983) 495-501)
C//  sig0 .... absorption
C//  sig1 .... multiple-phonon
C//  sig2 .... single-phonon

C//  sigmaa ... absorption for 1A neutrons [barn*A^-1]
C//  sigmab ... bound-atom scattering cross-section [barn]
C//  V0 .... volume [A^3]/atom
C//  A  .... atomic number
C//  thetaD .... Debye temperature (K)
C//  C2 .... constant from the Freund's paper  [A^-2 meV^-1]
C--------------------------------------------------------
      USE CRYSTALS
      USE CONSTANTS
      use SIMATH
      IMPLICIT NONE
      TYPE (CRYSTAL) :: CR
      real(kind(1.D0)) :: lambda,E,T,R,x,Xn,Ifact,B(0:30)
      real(kind(1.D0)) :: sig0,sig1,sig2,B0,BT,sigtot,DWMF
      real(kind(1.D0)) :: GETDW
      INTEGER :: I
1     format(a,' :',6(1x,G12.5))

!      DATA B /1,-0.5,0.166667,0,-0.033333,0,0.0238095,0,-0.033333,
!     *  0,0.0757576,0,-0.253114,0,1.16667,0,-7.09216,0,54.9712,
!     *  0,-529.124,0,6192.12,0,-86580.3,0,1.42551717e6,0,
!     * -2.7298231e7,0,6.01580874e8/

      GETMI=0
      E=HSQOV2M*(2*PI/lambda)**2

      call BERNOULLI(30,B)
      x=CR%thetaD/max(1.D-1,T)
      if (x.le.6) then
        R=0.D0
        Ifact=1.D0
        Xn=1.D0/x
        DO I=0,30
          R=R+B(I)*Xn/(Ifact*(I+2.5D0))
          Xn=Xn*x
          Ifact=Ifact*(I+1)
        ENDDO
      else
        R=3.3/SQRT(x**7)
      endif


      call DWFACT(CR,T,lambda,B0,BT)
      CR%DW=GETDW(B0,BT,CR%DHKL)
      DWMF=1.D0-exp(-(B0+BT)*CR%C2*E)
      sig0=CR%sigmaa*lambda
      sig1=CR%sigmab*(CR%A/(CR%A+1))**2*DWMF
      sig2=3*CR%sigmab/CR%A*SQRT(k_B*CR%thetaD/E)*R
      sigtot=(sig0+sig1+sig2)/CR%VOL
      GETMI=sigtot
      !write(*,1) 'GETMI lambda,x,R',lambda,x,R
      !write(*,1) 'GETMI B0,BT',B0,BT
      !write(*,1) 'GETMI sig0, sig1, sig2',sig0, sig1, sig2
      END FUNCTION GETMI

C-----------------------------------------------------------------
      SUBROUTINE GETREFPAR(CR,lambda,QK,MI0,R0,Integ,wg,fwhm)
C get integral parameters of the reflectivity curve
C-----------------------------------------------------------------
      USE CRYSTALS
      USE CONSTANTS
      IMPLICIT NONE
      TYPE (CRYSTAL) :: CR
      REAL*8 lambda,R0,fwhm,Integ,X(101),Y(101),SUMSQ,MOSREF
      REAL*8 MI0,QK,maxy,X2,X1,DX,WG,XMIN,STHB,THB,CTHB,STMCHI,STPCHI
      INTEGER*4 I,I1,I2

      IF(QK.LE.0) THEN
        R0=0
        Integ=0
        wg=0
        fwhm=0
        return
      ENDIF

      maxy=0
  ! mosaic crystal
      IF(CR%HMOS.GT.SEC) THEN
        DX=CR%HMOS*10./101.
        XMIN=-CR%HMOS*5
  ! bent perfect crystal
      ELSE
        STHB=lambda/2/CR%dhkl
        IF (ABS(STHB).GE.1) THEN
          return
        ENDIF
        THB=ASIN(STHB)
        CTHB=COS(THB)
        STMCHI=ABS(sin(THB-CR%CHI))
        STPCHI=ABS(sin(THB+CR%CHI))
        DX=CR%RH*cos(CR%CHI)/CTHB*(1-(1+CR%POI)*STPCHI*STMCHI)
        DX=DX+CR%DGR*1.D-4*cos(CR%DGA-CR%CHI)/CTHB
        IF(STMCHI.GT.0.05) THEN
          DX=DX*CR%FRAME%SIZE(3)/STMCHI
        ELSE
          DX=DX*CR%FRAME%SIZE(1)/STPCHI
        ENDIF
        DX=ABS(DX)*1.5
        XMIN=-DX/2
        DX=DX/101
      ENDIF

      DO I=1,101
        X(I)=XMIN+DX*I
        Y(I)=MOSREF(CR,X(I),lambda,MI0,QK)
        IF (maxy.lt.Y(i)) maxy=Y(I)
      ENDDO
      I=1
      DO WHILE ((I.LE.101).AND.(Y(I).LE.maxy/2))
         I=I+1
      ENDDO
      I1=I
      I=101
      DO WHILE ((I.GE.1).AND.(Y(I).LE.maxy/2))
         I=I-1
      ENDDO
      I2=I
      X1=X(I1-1)+(maxy/2-Y(I1-1))/(Y(I1)-Y(I1-1))*DX
      X2=X(I2)+(maxy/2-Y(I2))/(Y(I2+1)-Y(I2))*DX
      fwhm=X2-X1

      Integ=Y(1)+Y(101)
      DO I=2,100
         IF (MOD(I,2).EQ.0) THEN
            Integ=Integ+4*Y(I)
         ELSE
            Integ=Integ+2*Y(I)
         ENDIF
      ENDDO
      Integ=Integ/3*(X(2)-X(1))
      SUMSQ=Y(1)*X(1)**2+Y(101)*X(101)**2
      DO I=2,100
         IF (MOD(I,2).EQ.0) THEN
            SUMSQ=SUMSQ+4*Y(I)*X(I)**2
         ELSE
            SUMSQ=SUMSQ+2*Y(I)*X(I)**2
         ENDIF
      ENDDO
      SUMSQ=SUMSQ/3*(X(2)-X(1))
      wg=SQRT(SUMSQ/Integ)*R8LN2
      R0=MOSREF(CR,0.D0,lambda,MI0,QK)
      END SUBROUTINE GETREFPAR



C------------------------------------------------------
      REAL*8 FUNCTION MOSREF(CR,X,LAMBDA,MI0,QK)
C reflecivity function for mosaic crystal according to Bacon, Lowde, Sears
C treats also bent perfect crystals (box-shape)
C-------------------------------------------------------
      USE CRYSTALS
      USE CONSTANTS
      IMPLICIT NONE
      REAL*8 X,LAMBDA
      TYPE (CRYSTAL) :: CR
      REAL*8 ETA1,STPCHI,STMCHI,THB,STHB,CTHB,MI0,CGEOM
      REAL*8 GAMMA1,GAMMA2,QK,SIGMA,RB,RL,DUM


      ETA1=CR%HMOS
      STHB=lambda/2/CR%dhkl
      IF (ABS(STHB).GE.1) THEN
         MOSREF=0
         return
      ENDIF
      THB=ASIN(STHB)
      CTHB=COS(THB)
      STMCHI=ABS(sin(THB-CR%CHI))
      STPCHI=ABS(sin(THB+CR%CHI))

      IF (ETA1.GT.SEC.AND.(STMCHI.LT.1.D-5.OR.
     *    STPCHI.LT.1.D-5)) THEN
        write(*,*) 'Cannot calculate reflectivity '
        MOSREF=0
        RETURN
      ENDIF
  ! bent perfect crystal
      IF (ETA1.LE.SEC) THEN
        DUM=CR%RH*cos(CR%CHI)/CTHB*(1-(1+CR%POI)*STPCHI*STMCHI)
        DUM=DUM+CR%DGR*1.D-4*cos(CR%DGA-CR%CHI)/CTHB
        IF(STMCHI.GT.0.05) THEN
          DUM=DUM*CR%FRAME%SIZE(3)/STMCHI
        ELSE
          DUM=DUM*CR%FRAME%SIZE(1)/STPCHI
        ENDIF
        IF (ABS(DUM)/2.GT.ABS(X)) THEN
          MOSREF=CR%REF
        ELSE
          MOSREF=0
        ENDIF
  ! mosaic crystal
      ELSE
C!!! sign of CHI is oposite to the definition in TRAX
        CGEOM=(STPCHI-STMCHI)/(STMCHI+STPCHI)
        GAMMA1=(1/STMCHI+1/STPCHI)/2.D0
        GAMMA2=(1/STMCHI-1/STPCHI)/2.D0
        IF(MI0.EQ.0) THEN

         SIGMA=exp(-0.5*(X/ETA1)**2)/ETA1/SQRT2PI
         IF(THB.LT.CR%CHI) THEN
           MOSREF=(1+CGEOM)/2*(1.D0-exp(-QK*CR%DW*gamma1*CR%FRAME%SIZE(3)/10.*2*SIGMA))     ! Laue
         ELSE
           IF (GAMMA2.NE.0) THEN                                     ! Bragg
             MOSREF=(1.D0+CGEOM)/
     1      (1+CGEOM/TANH(gamma2*QK*CR%DW*CR%FRAME%SIZE(3)/10.*SIGMA))
           ELSE
             MOSREF=1.D0/(1.D0+1.D0/(gamma1*QK*CR%DW*CR%FRAME%SIZE(3)/10.*SIGMA))
           ENDIF
         ENDIF

        ELSE  ! absorbing case
          SIGMA=QK*CR%DW/MI0*exp(-0.5*(X/ETA1)**2)/ETA1/SQRT2PI
          RB=SQRT(1+2*SIGMA+(SIGMA*CGEOM)**2)
          RL=SQRT((1+2*SIGMA)*CGEOM**2+SIGMA**2)
          IF(THB.LT.CR%CHI) THEN
           DUM=SIGMA*(1+CGEOM)*EXP(-(1+SIGMA)*CR%FRAME%SIZE(3)/
     1        10.*MI0*gamma1)
           MOSREF=DUM*SINH(CR%FRAME%SIZE(3)/10*MI0*gamma1*RL)/RL       ! Laue
          ELSE
           MOSREF= SIGMA*(1+CGEOM)/(1+SIGMA+RB/TANH(CR%FRAME%SIZE(3)/
     1      10.*MI0*gamma1*RB))  ! Bragg
          ENDIF
        ENDIF
      ENDIF
      END FUNCTION MOSREF

C     ---------------------------------------------------
      SUBROUTINE MAKEUPCASE(LINE)
C     converts LINE to uppercase
C     ---------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) LINE
      INTEGER*4 L,I
      L=LEN(LINE)
      DO i=1,L
         IF((LINE(i:i).GE.'a').AND.(LINE(i:i).LE.'z')) THEN
            LINE(i:i)=CHAR(ICHAR(LINE(i:i))-32)
         ENDIF
      END DO
      END SUBROUTINE MAKEUPCASE



C---------------------------------------------------------------
      REAL*8 FUNCTION HMOS_DIST(X)
C Distribution of horizontal angular deviation of mosaic blocks
C --------------------------------------------------------------
      USE CRYSTALS
      USE CONSTANTS
      IMPLICIT NONE
      REAL(KIND(1.D0)),intent(in) :: X
      REAL(KIND(1.D0)), PARAMETER :: SQRT2LN2=1.177410022515474691D0
      REAL(KIND(1.D0)), PARAMETER :: DD0=0.270347525583091532681D0   ! 1/(pi*sqrt(2*ln2))

      IF (MDIST.EQ.1) THEN
C pseudo-Voigt
        HMOS_DIST=0.5*DD0/(1. + (X/SQRT2LN2)**2)+0.5*EXP(-0.5*X**2)/SQRT2PI
      ELSE IF (MDIST.EQ.2) THEN
C Lorenz
        HMOS_DIST=DD0/(1. + (X/SQRT2LN2)**2)
      ELSE IF (MDIST.EQ.3) THEN
C Rectangle
        IF (ABS(X).LE.0.5) THEN
          HMOS_DIST=1.
        ELSE
          HMOS_DIST=0.
        ENDIF
      ELSE
C Gauss
        HMOS_DIST=EXP(-0.5*X**2)/SQRT2PI
      ENDIF
      END FUNCTION HMOS_DIST

C-----------------------------------------------------
      SUBROUTINE ERF_INIT(F,AMIN,AMAX)
C Calculate lookup tables for ERF_CR function
C-----------------------------------------------------
      IMPLICIT NONE
      INTEGER*4 I,J
      REAL*8 SUM,Z1,Z2,Z3,A,B,DET,X1I,XJ,AMIN,AMAX
      INTEGER*4 DIM
      PARAMETER(DIM=1025)
      REAL*8 XMIN,DX,Y(DIM),XMIN1,DX1,Y1(DIM)
      COMMON /ERFCOM/ XMIN,DX,Y,XMIN1,DX1,Y1
      REAL*8 F
      EXTERNAL F

! Generate cumulative function from F(X)
      SUM=0.
      DX=(AMAX-AMIN)/(DIM-1)
      XMIN=AMIN
      Y(1)=0.
      DO I=1,DIM-1
        Z1=XMIN+(I-1)*DX
        Z2=Z1+DX/2
        Z3=Z1+DX
        SUM=SUM+F(Z1)+4*F(Z2)+F(Z3)
        Y(I+1)=SUM
      END DO
      DO I=1,DIM
        Y(I)=Y(I)/Y(DIM)
      ENDDO
! Generate inverse cumulative function
      DX1=1.D+0/(DIM-1)
      XMIN1=Y(1)
      Y1(1)=XMIN
      Y1(DIM)=XMIN+(DIM-1)*DX
      J=1
      DO I=2,DIM-1
        X1I=XMIN1+(I-1)*DX1
        DO WHILE ((J.LT.DIM-1).and.(Y(J).LT.X1I))
          J=J+1
        END DO
10      XJ=XMIN+(J-1)*DX
        A=(Y(J+1)+Y(J-1)-2*Y(J))/2
        B=(Y(J+1)-Y(J-1))/2
!				goto 20
        IF (ABS(A).LT.1D-30) THEN
          if (J.le.2) then
				    write(*,*) 'error in ERF_INIT, J=',J,'  I=',I
          endif
          J=J-1
          GOTO 10
        ELSE
          DET=B**2-4*A*(Y(J)-X1I)
          IF (DET.LE.0.D0) then
            write(*,*) 'error in ERF_INIT: ',DET,A,B
          ENDIF
          Z1=XJ+DX*(-B+SQRT(DET))/2/A
          Z2=XJ+DX*(-B-SQRT(DET))/2/A
          IF (ABS(Z2-XJ).LT.ABS(Z1-XJ)) Z1=Z2
          Y1(I)=Z1
        ENDIF
!20    continue
      END DO
      END SUBROUTINE ERF_INIT

C----------------------------------------------------------------
      REAL*8 FUNCTION ERF_CR(ARG,INV)
c Return cumulative function (or inverse, if INV=0)
c Uses lookup table generated by CUM_INIT
C----------------------------------------------------------------
      IMPLICIT NONE
      INTEGER*4 DIM,INV,J1,J2,J,I
      REAL*8 ARG,A,B,Z,XJ,DET,Z1,Z2,ARG1
      PARAMETER(DIM=1025)
      REAL*8 XMIN,DX,Y(DIM),XMIN1,DX1,Y1(DIM)
      REAL*8 ERF_INTERP
      COMMON /ERFCOM/ XMIN,DX,Y,XMIN1,DX1,Y1

C Cumul. function
      IF (INV.NE.1) THEN
        ERF_CR=ERF_INTERP(XMIN,DX,Y,DIM,ARG)
C else Iverse cumul. function
      ELSE
        IF(ARG.LE.XMIN1) THEN
           ERF_CR=Y1(1) ! left limit
        ELSE IF(ARG.GE.XMIN1+(DIM-1)*DX1) THEN
           ERF_CR=Y1(DIM)  ! right limit
        ELSE
C Find J1,J2 so that Y(J1) > A >= Y(J2) and J2=J1+1
          ARG1=ARG
          IF(ARG.GT.9.D-1) ARG1=1.D0-ARG
          Z=(ARG1-XMIN1)/DX1+1
          I=INT(Z)
          Z1=(Y1(I)-XMIN)/DX+1
          Z2=(Y1(I+1)-XMIN)/DX+1
          J1=INT(Z1)
          J2=INT(Z2)+1
          DO WHILE(J2.GT.J1+1)
            J=(J2+J1)/2
            IF(Y(J).LT.ARG1) THEN
               J1=J
            ELSE
               J2=J
            ENDIF
          ENDDO
C Set J so that Y(J) is close to ARG
          J=J1
          IF(ARG1-Y(J1).GT.Y(J2)-ARG1) J=J2
C but avoid J=1 or J=DIM
          IF(J.LT.2) J=J+1
          IF(J.GT.DIM-1) J=J-1
C interpolate quadratically between Y(J-1) and Y(J+1)
C return inverse value for Y=ARG
10        XJ=XMIN+(J-1)*DX
          A=(Y(J+1)+Y(J-1)-2*Y(J))/2
          B=(Y(J+1)-Y(J-1))/2
          IF(ABS(B).LT.1.D-30) THEN
             ERF_CR=XJ
          ELSE IF (ABS(A).LT.1D-30) THEN
            ERF_CR=XJ+Y(J+1)
            J=J-1
            GOTO 10
          ELSE
            DET=B**2-4*A*(Y(J)-ARG1)
            DET=SQRT(DET)
            Z1=XJ+DX*(-B+DET)/2/A
            Z2=XJ+DX*(-B-DET)/2/A
            IF (ABS(Z2-XJ).LT.ABS(Z1-XJ)) Z1=Z2
            IF(ARG.GT.9.D-1) Z1=-Z1
            ERF_CR=Z1
          ENDIF
        ENDIF
      ENDIF
      END FUNCTION ERF_CR

C-------------------------------------------------------------------
      REAL*8 FUNCTION ERF_INTERP(XMIN,DX,Y,DIM,A)
C Quadratic interpolation in Y array with equidistant X
C given by XMIN,DX. MIN,MAX are values to be returned outside limits
C-------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER*4 DIM,I
      REAL*8 XMIN,DX,A,Z,XMAX
      REAL*8 Y(DIM)

      XMAX=XMIN+(DIM-1)*DX
      if (A.LE.XMIN) THEN
         ERF_INTERP=Y(1)
      elseif (A.GE.XMAX) THEN
         ERF_INTERP=Y(DIM)
      elseif (A.LE.XMIN+DX) THEN
         ERF_INTERP=Y(1)+(Y(2)-Y(2))*(A-XMIN)/DX
      elseif (A.GE.XMAX-DX) THEN
         ERF_INTERP=Y(DIM-1)+(Y(DIM)-Y(DIM-1))*(A-XMAX+DX)/DX
      else
        Z=(A-XMIN)/DX+1
        I=NINT(Z)
        if (I.EQ.DIM-1) I=I-1
        IF (I.EQ.2) I=3
        ERF_INTERP=Y(I)+(Y(I+1)-Y(I-1))/2*(Z-I)+(Y(I+1)+Y(I-1)-2*Y(I))/2*(Z-I)**2
      endif
      end FUNCTION ERF_INTERP

