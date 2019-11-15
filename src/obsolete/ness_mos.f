c      PROGRAM MOSAIC
c      CALL TEST
c      END

      SUBROUTINE TEST
      USE COMPONENTS
      IMPLICIT NONE
      INCLUDE 'const.inc'

      TYPE (CRYSTAL) :: CR
      REAL*8 GETMI,SQRT2PI
      PARAMETER (SQRT2PI=2.506628275)

      REAL*8 LAMBDA,ETA,THICK,CHI,R0,fwhm,Integ,weff,wmos,GETQKIN,reff
      REAL*8 KI,QKIN,MI,TEMP
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
         CALL READCRYST(CR,CRNAME)
      endif
      WRITE(*,3)
      READ(*,*) ETA,THICK,CHI,TEMP
      CR.HMOS=ETA*PI/180/60/R8LN2
      CR.FRAME.SIZE(3)=THICK
      CR.CHI=CHI*PI/180
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
        QKIN= GETQKIN(CR,lambda)
        MI= GETMI(CR,lambda,TEMP)
        call getrefpar(CR,lambda,QKIN,MI,R0,Integ,wmos,fwhm)
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
        write(*,15) KI,fwhm,R0,reff,QKIN,MI,MI*CR.VOL
      ENDDO
100   continue
      GOTO 50
200   CONTINUE
      END

CxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCX
C--------------------------------------------------------
      REAL*8 FUNCTION GETQKIN(CR,lambda)
C--------------------------------------------------------
      USE COMPONENTS
      IMPLICIT NONE
      INCLUDE 'const.inc'
      TYPE (CRYSTAL) :: CR
      REAL*8 STHB,lambda,CTHB

      STHB=lambda/2/CR.dhkl
      IF (STHB.GE.9.99999D-1) THEN
         GETQKIN=0
        RETURN
      ENDIF
      CTHB=SQRT(1.D0-STHB**2)
      GETQKIN=CR.QML*CR.DHKL*STHB**2/PI/CTHB

      END

C--------------------------------------------------------
      REAL*8 FUNCTION GETREFDYN(CR,lambda)
C--------------------------------------------------------
      USE COMPONENTS
      IMPLICIT NONE
      INCLUDE 'const.inc'
      TYPE (CRYSTAL) :: CR
      REAL*8 Z,Z1,Z2,grad
      INTEGER*4 I,J
      REAL*8 STHB,lambda,CTHB,THB,QK,sigma,ki(3),kg(3)
      REAL*8 GETDW

      STHB=lambda/2/CR.dhkl
      IF (STHB.GE.9.99999D-1) THEN
        GETREFDYN=0
        RETURN
      ENDIF

      CTHB=SQRT(1.D0-STHB**2)
      QK=CR.QML*CR.DHKL*STHB**2/PI  !  = Qkin*cos(thb)
      QK=QK*GETDW(CR,300.D0,Z1,Z2,lambda)/10.
      THB=ASIN(STHB)
      ki(1)=-cos(THB-CR.CHI)
      ki(2)=0.
      ki(3)=-sin(THB-CR.CHI)
      kg(1)=-cos(THB+CR.CHI)
      kg(2)=0.
      kg(3)=+sin(THB+CR.CHI)
      grad=0.
      DO I=1,3
        Z=0.
        DO J=1,3
          Z=Z+CR.DG_DR(I,J)*KI(J)
        ENDDO
        grad=grad+Z*KG(I)
      ENDDO
      grad=abs(grad)/CR.GTOT
      if (grad/QK.lt.1D-2) then
         Z=1.D0
      else
         sigma=QK/grad
         Z=1.D0-exp(-sigma)
      endif
      GETREFDYN=Z
      END

C-------------------------------------------------------------
      REAL*8 FUNCTION GETDTH(CR)
C Get the spread of Bragg angles due to the elastic deformation
C-------------------------------------------------------------
      USE COMPONENTS
      IMPLICIT NONE
      TYPE (CRYSTAL) :: CR
      REAL*8 Z,grad
      INTEGER*4 I,J
      REAL*8 ki(3),kg(3),R(3),TIN,TOUT,STHB,CTHB

      STHB=sin(abs(CR.THB))
      IF (STHB.GE.9.99999D-1) THEN
         GETDTH=0.D0
        RETURN
      ENDIF

      ki(1)=-cos(CR.THB-CR.CHI)
      ki(2)=0.
      ki(3)=-sin(CR.THB-CR.CHI)
      kg(1)=-cos(CR.THB+CR.CHI)
      kg(2)=0.
      kg(3)=+sin(CR.THB+CR.CHI)
      grad=0.
      DO I=1,3
        Z=0.
        DO J=1,3
          Z=Z+CR.DG_DR(I,J)*KI(J)
        ENDDO
        grad=grad+Z*KG(I)
      ENDDO
c      write(*,*) CR.FRAME.NAME,grad,cos(CR.THB),CR.GTOT
      CTHB=SQRT(1.D0-STHB**2)
      if (CR.GTOT.GE.1.D-10) THEN
         grad=abs(grad)/CTHB/CR.GTOT
      else
         grad=0.D0
      endif
      R(1)=0.
      R(2)=0.
      R(3)=0.
      call BORDER_BOX(R,KI,CR.FRAME.SIZE,TIN,TOUT)
c      write(*,*) TIN,TOUT
      GETDTH=grad*(TOUT-TIN)
      END

C--------------------------------------------------------------------
      REAL*8 FUNCTION GETEFFMOS(CR)
C Get effective mosaicity
C = change of refl. plane orientation along incident beam path
C--------------------------------------------------------------------
      USE COMPONENTS
      IMPLICIT NONE
      TYPE (CRYSTAL) :: CR
      REAL*8 Z,G01,G02
      INTEGER*4 I,J
      REAL*8 ri(3),rf(3),G1(3),G2(3),CTHB,SCHI,CCHI

      CTHB=cos(CR.THB-CR.CHI)
      ri(1)=0.5*CTHB*CR.FRAME.SIZE(3)
      ri(2)=0.
      ri(3)=0.5*CR.FRAME.SIZE(3)
      rf(1)=-0.5*CTHB*CR.FRAME.SIZE(3)
      rf(2)=0.
      rf(3)=-0.5*CR.FRAME.SIZE(3)
      DO I=1,3
        G1(I)=0.
        G2(I)=0.
        DO J=1,3
          G1(I)=G1(I)+CR.DG_DR(I,J)*ri(J)
          G2(I)=G2(I)+CR.DG_DR(I,J)*rf(J)
        ENDDO
      ENDDO
      SCHI=sin(CR.CHI)
      CCHI=SQRT(1.D0-SCHI**2)
      G1(1)=G1(1)+CR.GTOT*SCHI
      G1(3)=G1(3)-CR.GTOT*CCHI
      G2(1)=G2(1)+CR.GTOT*SCHI
      G2(3)=G2(3)-CR.GTOT*CCHI
      Z=0.
      DO I=1,3
        Z=Z+G1(I)*G2(I)
      ENDDO
      IF (ABS(G01*G02).LT.1.D-10) THEN
         Z=0.D0
      ELSE
        G01=SQRT(G1(1)**2+G1(3)**2)
        G02=SQRT(G2(1)**2+G2(3)**2)
        Z=ACOS(Z/G01/G02)
      ENDIF
      GETEFFMOS=Z
      END

C--------------------------------------------------------
      REAL*8 FUNCTION GETDW(CR,T,B0,BT,lambda)
C///  Debye-Waller factor
C/// cubic crystals, Debye model
C--------------------------------------------------------
      USE COMPONENTS
      IMPLICIT NONE
      INCLUDE 'const.inc'
      TYPE (CRYSTAL) :: CR
      INTEGER*4 P
      REAL*8 B0,BT,X,Z,DX,ksi,T,E,lambda
      INTEGER*4 I
      B0=2873/CR.A/CR.thetaD
      x=CR.thetaD/T
      Z=0
      dx=x/100
      DO i=1,101
        if (i.eq.1.or.i.eq.101) THEN
          P=1
        else if (MOD(i,2).eq.1) THEN
          P=4
        else
          P=2
        endif
        if (i.eq.1) then
          Z=-1
        else
          ksi=(i-1)*dx
          Z=Z+P*ksi/(exp(ksi)-1)
        endif
      ENDDO
      Z=Z*dx/3
      BT=4*B0*Z/x**2
      E=HSQOV2M*(2*PI/lambda)**2
c      GETDW=exp(-(B0+BT)*CR.C2*E)
      GETDW=exp(-2.0*(B0+BT)/(CR.dhkl*2)**2)

      END

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
      USE COMPONENTS
      IMPLICIT NONE
      INCLUDE 'const.inc'
      TYPE (CRYSTAL) :: CR
      REAL*8 KB
      PARAMETER (KB  = 0.086174)   ! Boltzman constant [meV*K^-1]
      REAL*8 lambda,E,T,R,x,Ifact,B(0:30)
      REAL*8 sig0,sig1,sig2,GETDW,B0,BT,sigtot
      INTEGER*4 I

      DATA B /1,-0.5,0.166667,0,-0.033333,0,0.0238095,0,-0.033333,
     *  0,0.0757576,0,-0.253114,0,1.16667,0,-7.09216,0,54.9712,
     *  0,-529.124,0,6192.12,0,-86580.3,0,1.42551717e6,0,
     * -2.7298231e7,0,6.01580874e8/

      GETMI=0
      E=HSQOV2M*(2*PI/lambda)**2


      x=CR.thetaD/T
      R=0
      Ifact=1
      DO I=0,20
        Ifact=Ifact*I
        IF (I.EQ.0) IFact=1
        R=R+B(I)*X**(I-1)/(Ifact*(I+5/2))
      ENDDO
      R=3*(R+B(22)*X**21/(Ifact*21*22*(22+5/2))/2)


      CR.DW=GETDW(CR,T,B0,BT,lambda)
      sig0=CR.sigmaa*lambda
      sig1=CR.sigmab*(CR.A/(CR.A+1))**2*(1-exp(-(B0+BT)*CR.C2*E))
      if (x.le.6) then
        sig2=CR.sigmab/CR.A*SQRT(KB*CR.thetaD/E)*R
      else
        sig2=CR.sigmab/CR.A*SQRT(KB*CR.thetaD/E)*3.3/SQRT(x**7)
      endif
      sigtot=(sig0+sig1+sig2)/CR.VOL
      GETMI=sigtot
c      IF(INDEX(CR.FRAME.NAME,'PG').GT.0) THEN  ! Pyrolytic graphite
c         IF(sigtot*CR.VOL.lt.4) GETMI=4/CR.VOL
c      ENDIF

c10    FORMAT(a13,3(G13.5,2x))
c      write(*,10) 'sigma     : ',sig0,sig1,sig2
c      write(*,10) 'DW,B0,BT  : ',DW,B0,BT
c      write(*,10) 'R,t/E,s/A : ',R,KB*CR.thetaD/E,CR.sigmab/CR.A
c      write(*,10) 'V0, MI    : ',CR.VOL,(sig0+sig1+sig2)/CR.VOL

      END

C-----------------------------------------------------------------
      SUBROUTINE GETREFPAR(CR,lambda,QKIN,MI0,R0,Integ,wg,fwhm)
C-----------------------------------------------------------------
      USE COMPONENTS
      IMPLICIT NONE
      INCLUDE 'const.inc'
      TYPE (CRYSTAL) :: CR
      REAL*8 lambda,R0,fwhm,Integ,X(101),Y(101),SUMSQ,MOSREF
      REAL*8 MI0,QKIN,maxy,X2,X1,DX,WG,XMIN,STHB,THB,CTHB,STMCHI,STPCHI
      INTEGER*4 I,I1,I2

      IF(QKIN.LE.0) THEN
        R0=0
        Integ=0
        wg=0
        fwhm=0
        return
      ENDIF

      maxy=0
      IF(CR.HMOS.GT.SEC) THEN
        DX=CR.HMOS*10./101.
        XMIN=-CR.HMOS*5
      ELSE
        STHB=lambda/2/CR.dhkl
        IF (ABS(STHB).GE.1) THEN
          return
        ENDIF
        THB=ASIN(STHB)
        CTHB=COS(THB)
        STMCHI=ABS(sin(THB-CR.CHI))
        STPCHI=ABS(sin(THB+CR.CHI))
        DX=CR.RH*cos(CR.CHI)/CTHB*(1-(1+CR.POI)*STPCHI*STMCHI)
        DX=DX+CR.DGR*1.D-4*cos(CR.DGA-CR.CHI)/CTHB
        IF(STMCHI.GT.0.05) THEN
          DX=DX*CR.FRAME.SIZE(3)/STMCHI
        ELSE
          DX=DX*CR.FRAME.SIZE(1)/STPCHI
        ENDIF
        DX=ABS(DX)*1.5
        XMIN=-DX/2
        DX=DX/101
      ENDIF

      DO I=1,101
        X(I)=XMIN+DX*I
        Y(I)=MOSREF(CR,X(I),lambda,MI0,QKIN)
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
      R0=MOSREF(CR,0.D0,lambda,MI0,QKIN)
      END



C-------------------------------------------------
      SUBROUTINE READCRYST(CR,CRNAME)
C// QML ... Maier-Leibnitz reflectivity in A^-1 cm^-1
C// QML = 4*PI*(F*dhkl/V0)**2
C-------------------------------------------------
      USE COMPONENTS
      IMPLICIT NONE
      INCLUDE 'inout.inc'
      TYPE (CRYSTAL) :: CR
      CHARACTER*8 CRNAME1
      CHARACTER*(*) CRNAME
      CHARACTER*128 LINE
      INTEGER*4 INDX,IRES,IS,IL,ILINE
1     FORMAT(a)
2     FORMAT('Format error in the crystal library, line ',I5 )
3     FORMAT('Crystal ',a,' not found in the library' )
4     FORMAT('Cannot open crystal library' )
      INDX=0

      CALL BOUNDS(CRNAME,IS,IL)
      CRNAME1=CRNAME(IS:IS+IL-1)
      CALL MakeUpCase(CRNAME1)

      CALL OPENRESFILE('crystal.lib',1,IRES,1)
      IF(IRES.LE.0) GOTO 99
c      OPEN(unit=1,file='crystal.lib',status='OLD',err=98)
      ILINE=0
      DO WHILE (INDX.EQ.0)
         READ(1,1,END=100,ERR=97) LINE
         ILINE=ILINE+1
         CALL MakeUpCase(LINE)
         INDX=INDEX(LINE,CRNAME1)
      ENDDO
100   close(1)
      IF(INDX.EQ.0) GOTO 98
      READ(LINE(INDX+8:LEN(LINE)),*,err=97)
     * CR.DHKL,CR.QML,CR.sigmab,CR.sigmaa,CR.VOL,CR.A,CR.thetaD,CR.C2,
     * CR.poi
      CR.C2=CR.C2/1000  ! from eV^-1 to meV^-1
      RETURN

C format error
97    write(smes,2) ILINE
      CR.VOL=0
      close(1)
      return

C crystal not found
98    write(smes,3) CRNAME1(1:IL)
      CR.VOL=0
      close(1)
      return

C cannot open
99    write(smes,4)
      CR.VOL=0
      return


      END


C-------------------------------------------------
      REAL*8 FUNCTION MOSREF(CR,X,LAMBDA,MI0,QKIN)
C-------------------------------------------------
      USE COMPONENTS
      IMPLICIT NONE
      INCLUDE 'const.inc'
      REAL*8 SQRT2PI
      PARAMETER(SQRT2PI=2.506628275)
      REAL*8 X,LAMBDA
      TYPE (CRYSTAL) :: CR
      REAL*8 ETA1,STPCHI,STMCHI,THB,STHB,CTHB,MI0,CGEOM
      REAL*8 GAMMA1,GAMMA2,QKIN,SIGMA,RB,RL,DUM


      ETA1=CR.HMOS
      STHB=lambda/2/CR.dhkl
      IF (ABS(STHB).GE.1) THEN
         MOSREF=0
         return
      ENDIF
      THB=ASIN(STHB)
      CTHB=COS(THB)
      STMCHI=ABS(sin(THB-CR.CHI))
      STPCHI=ABS(sin(THB+CR.CHI))

      IF (ETA1.GT.SEC.AND.(STMCHI.LT.1.D-5.OR.
     *    STPCHI.LT.1.D-5)) THEN
        write(*,*) 'Cannot calculate reflectivity '
        MOSREF=0
        RETURN
      ENDIF
      IF (ETA1.LE.SEC) THEN     ! bent perfect crystal
        DUM=CR.RH*cos(CR.CHI)/CTHB*(1-(1+CR.POI)*STPCHI*STMCHI)
        DUM=DUM+CR.DGR*1.D-4*cos(CR.DGA-CR.CHI)/CTHB
        IF(STMCHI.GT.0.05) THEN
          DUM=DUM*CR.FRAME.SIZE(3)/STMCHI
        ELSE
          DUM=DUM*CR.FRAME.SIZE(1)/STPCHI
        ENDIF
        IF (ABS(DUM)/2.GT.ABS(X)) THEN
          MOSREF=CR.REF
        ELSE
          MOSREF=0
        ENDIF
      ELSE          ! mosaic crystal
C!!! sign of CHI is oposite to the definition in TRAX
        CGEOM=(STPCHI-STMCHI)/(STMCHI+STPCHI)
        GAMMA1=(1/STMCHI+1/STPCHI)/2
        GAMMA2=(1/STMCHI-1/STPCHI)/2
        IF(MI0.EQ.0) THEN

         SIGMA=exp(-0.5*(X/ETA1)**2)/ETA1/SQRT2PI
         IF(THB.LT.CR.CHI) THEN
           MOSREF=(1+CGEOM)/2*
     1       (1-exp(-QKIN*CR.DW*gamma1*CR.FRAME.SIZE(3)/10.*2*SIGMA))     ! Laue
         ELSE
           IF (GAMMA2.NE.0) THEN                                     ! Bragg
             MOSREF=(1+CGEOM)/
     1      (1+CGEOM/TANH(gamma2*QKIN*CR.DW*CR.FRAME.SIZE(3)/10.*SIGMA))
           ELSE
             MOSREF=
     1       1/(1+1/(gamma1*QKIN*CR.DW*CR.FRAME.SIZE(3)/10.*SIGMA))
           ENDIF
         ENDIF

        ELSE  ! absorbing case
          SIGMA=QKIN*CR.DW/MI0*exp(-0.5*(X/ETA1)**2)/ETA1/SQRT2PI
          RB=SQRT(1+2*SIGMA+(SIGMA*CGEOM)**2)
          RL=SQRT((1+2*SIGMA)*CGEOM**2+SIGMA**2)
          IF(THB.LT.CR.CHI) THEN
           DUM=SIGMA*(1+CGEOM)*EXP(-(1+SIGMA)*CR.FRAME.SIZE(3)/
     1        10.*MI0*gamma1)
           MOSREF=DUM*SINH(CR.FRAME.SIZE(3)/10*MI0*gamma1*RL)/RL       ! Laue
          ELSE
           MOSREF= SIGMA*(1+CGEOM)/(1+SIGMA+RB/TANH(CR.FRAME.SIZE(3)/
     1      10.*MI0*gamma1*RB))  ! Bragg
          ENDIF
        ENDIF
      ENDIF
      END

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
      END



