C $Id: ness_mcryst.f,v 1.1 2007/10/08 08:35:37 saroun Exp $
C//////////////////////////////////////////////////////////////////////
C////                                                              ////
C////  NEutron Scattering Simulation - v.2.0, (c) J.Saroun, 2000   ////
C////                                                              ////
C//////////////////////////////////////////////////////////////////////
C////
C////  Subroutines describing objects - CRYSTAL ARRAY
C////
C//////////////////////////////////////////////////////////////////////
C ---------------------------------------
      SUBROUTINE CRYST_INIT2(CR)
C ---------------------------------------
      USE COMPONENTS
      implicit none
      INCLUDE 'const.inc'
      TYPE (CRYSTAL) :: CR
      REAL*8 Z,GETQKIN,GETREFDYN,GETMI
      INTEGER*4 I,J

C///  THB,DHKL,etc.. must be specified before !

      CALL SLIT_INIT(CR.FRAME)

C// determine crystal type
  ! ignored
      if (CR%NH.LE.0) then
        CR%TYP=ctyp_none
  ! elastically bent
      else if (CR%HMOS.LT.sec) then
        CR%TYP=ctyp_bent
      else
  ! mosaic with gradient
        if ((CR.DGR.NE.0.D0).or.((CR%RH.NE.0.D0).and.(CR%NH.EQ.1))) then
          CR%TYP=ctyp_gradient
        else
  ! mosaic
          CR%TYP=ctyp_mosaic
        endif
      endif

      CR.GTOT=2*PI/ABS(CR.DHKL)
      CR.G(1)=CR.GTOT*SIN(CR.CHI)
      CR.G(2)=0.
      CR.G(3)=CR.GTOT*COS(CR.CHI)
      CR.STMCH=SIN(CR.THB-CR.CHI)
      CR.CTMCH=COS(CR.THB-CR.CHI)
      if (ABS(CR.THB).LT.1E-5) then
        CR.LAMBDA=1.8
      else
        CR.LAMBDA=2.*CR.DHKL*SIN(CR.THB)
      ENDiF
      DO I=1,3
        CR.MAPG(I)=.FALSE.
        DO J=1,3
          CR.DG_DR(I,J)=0.D0
        ENDDO
      ENDDO

c G-gradient
      select case(CR%TYP)
      case(ctyp_bent)
          CR.DG_DR(1,1)=-COS(CR.CHI)*CR.GTOT*CR.RH
          CR.DG_DR(1,3)=SIN(CR.CHI)*CR.GTOT*CR.RH
          CR.DG_DR(3,1)=SIN(CR.CHI)*CR.GTOT*CR.RH
          CR.DG_DR(2,2)=0.D0    ! no vertical bending
          CR.DG_DR(3,3)=-CR.POI*COS(CR.CHI)*CR.GTOT*CR.RH
          CR.MAPG(1)=.TRUE.
          CR.MAPG(3)=.TRUE.
      case(ctyp_gradient)
          Z=1.D-4*CR.GTOT*CR.DGR
          CR.DG_DR(1,1)=CR.DG_DR(1,1)+Z*cos(CR.DGA+CR.CHI)
          CR.DG_DR(1,3)=CR.DG_DR(1,3)-Z*sin(CR.DGA+CR.CHI)
          CR.DG_DR(3,1)=CR.DG_DR(3,1)-Z*sin(CR.DGA+CR.CHI)
          CR.DG_DR(3,3)=CR.DG_DR(3,3)-Z*cos(CR.DGA+CR.CHI)
          CR.MAPG(1)=.TRUE.
          CR.MAPG(3)=.TRUE.
      END select
c unit vector |- to G
      CR.gama(1)=COS(CR.CHI)
      CR.gama(3)=-SIN(CR.CHI)
      CR.gama(2)=0
  ! kin. reflectivity
      CR.QHKL=GETQKIN(CR,CR.lambda)
  ! extinction length
      CR.dext=CR.DHKL/CR.LAMBDA*SQRT(4*PI/CR.QML)
  ! absorption coefficient
      CR.MI=GETMI(CR,CR.lambda,300.D0)
  ! peak reflectivity (no mosaic)
      CR.REF=GETREFDYN(CR,CR.LAMBDA)
  ! Darwin box width
      CR.DELTA=CR.QHKL*CR.Dext*1D-4/PI

  ! special
      select case(CR%TYP)
      case (ctyp_mosaic,ctyp_gradient)
    ! mosaicity is always larger than Darwin box
        CR%HMOS=MAX(CR%HMOS,CR%DELTA)
    ! primary extinction
        Z=CR%dlam/CR%Dext
        IF (Z.GT.1D-5) CR%Ext1=TANH(Z)/Z
      case default
        CR.Ext1=1.D0
      END select

      END


C-------------------------------------------
      LOGICAL FUNCTION CRYST_GO2(CRYST,NEU)
C-------------------------------------------
      USE COMPONENTS
      USE TRACING
      implicit none
      INCLUDE 'const.inc'
      INCLUDE 'crystal.inc'
      INCLUDE 'ness_common.inc'
      INCLUDE 'randvars.inc'

      TYPE (CRYSTAL) :: CRYST
      TYPE(NEUTRON) :: NEU
      REAL*8 V(3),K(3),R(3),KaG(3)
      REAL*8 PP,DT,H,TIN,TOUT,DT1
      INTEGER*4 i

      TOF=0.D0
      CALL SLIT_PRE(CRYST%FRAME,NEU)

C/// if nh=0 then accept neutrons without transformations
      IF (CRYST%NH.EQ.0) then
        DT=-NEU%R(3)/NEU.K(3)
        DO I=1,3
          NEU%R(I)=NEU%R(I)+DT*NEU%K(I)
        END DO
        NEU%T=NEU%T+DT/HOVM
        CRYST%FRAME%COUNT=CRYST%FRAME%COUNT+1
        GOTO 10
      endif

! polarization condition
! quite stupid, should be improved
      if (CRYST%MAG*NEU%S(1).LT.0) GOTO 99

C/// move to the entry point
      call SLIT_ENTRY(CRYST%FRAME,NEU)
! log event at the entry
      call AddEventLog(CRYST%FRAME,NEU)

C work with a local copy of neutron coordinates
      V=NEU%R
      K=NEU%K

C trace through the  array of crystals
      call CRYST_ARRAY(CRYST,V,K,R,KaG,PP)

! reject all events with transmission lower than 0.0001
      if (PP.LE.1.D-4) goto 99
C accept event: update NEU structure with exit values
      do i=1,3
        NEU%R(i)=R(i)
        NEU%K(i)=KaG(i)
      enddo
      NEU%T=NEU%T+TOF/HOVM
      NEU%P=NEU%P*PP
    ! increment event counter
      CRYST%FRAME%COUNT=CRYST%FRAME%COUNT+1

! log event at the exit
      call AddEventLog(CRYST%FRAME,NEU)
! convert to exit coordinates
      CALL SLIT_POST(CRYST%FRAME,NEU)

10    CRYST_GO2=.TRUE.
      RETURN

99    NEU.P=0
      CRYST_GO2=.FALSE.
      END


C---------------------------------------------------------------
      REAL*8 FUNCTION HMOS_DIST(X)
C Distribution of horizontal angular deviation of mosaic blocks
C --------------------------------------------------------------
      USE COMPONENTS
      IMPLICIT NONE
      INCLUDE 'const.inc'
      INCLUDE 'ness_common.inc'
      REAL*8 X,SQRT2PI,SQRT2LN2,DD0
      PARAMETER (SQRT2PI=2.506628274631, SQRT2LN2=1.177410023, DD0=0.270347526)

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
      END

C-----------------------------------------------------------
      subroutine CRYST_ARRAY(CR,R1,K1,R2,K2,PP)
C Transmission function for mosaic crystal
C (eta>")
C simplified version (Darwin width=0)
C-----------------------------------------------------------
      USE COMPONENTS
      USE TRACING, ONLY: AddEventLog,NEUTRON
      implicit none
      INCLUDE 'const.inc'
      INCLUDE 'crystal.inc'
      REAL*8 EPS
      PARAMETER (EPS=1D-5)
      TYPE (CRYSTAL) :: CR
      REAL*8 R1(3),K1(3),R2(3),K2(3),DT
      REAL*8 G(3)
      INTEGER*4 I
      REAL*8 PP,TIN,TOUT
      REAL*8 KK,K0,Qhkl,Z,MI
      type(NEUTRON) :: NEU

      PP=1.D0
C/// calculate values, which are constant in the course of the iteration
      KK=K1(1)**2+K1(2)**2+K1(3)**2
      K0=SQRT(KK)
      Qhkl=CR.Qhkl/10.*CR.DW*CR.Ext1
      MI=CR.MI/10.
      PATH=0.D0
      R2=R1

C ****  Begin of multiple scattering iteration cycle  ****

C// generate random walk step in K direction
C---------------------------------------------------
50    CALL WALKSTEP(CR,1,R2,K1,K0,Qhkl,G,PP,TOUT)
      if (PP.LT.EPS) goto 99   ! don't follow trajectories with low probability
      ! log event
      NEU.R=R2
      NEU.K=K1
      NEU.P=PP
  ! log event
      call AddEventLog(CR%FRAME,NEU)
  ! set K=K+G
      DO I=1,3
        K2(I)=K1(I)+G(I)
      END DO
  ! renormalize
      Z=SQRT(KK/(K2(1)**2+K2(2)**2+K2(3)**2))  ! set |K2| = |K1|
      DO I=1,3
        K2(I)=K2(I)*Z
      END DO

C// generate 2nd random walk step in K+G direction
C-------------------------------------------------
      CALL WALKSTEP(CR,-1,R2,K2,K0,Qhkl,G,PP,TOUT)
      IF (PP.LT.EPS) GOTO 99
      ! log event
      NEU.R=R2
      NEU.K=K2
      NEU.P=PP
  ! log event
      call AddEventLog(CR%FRAME,NEU)
      if (CR%model.eq.0) TOUT=0.D0
      if (TOUT.LE.0) THEN
  ! reached exit
        PP=PP*EXP(-MI*PATH)
      else
  ! continue with random walk
    ! set K=K+G
        DO I=1,3
          K1(I)=K2(I)+G(I)
        END DO
    ! renormalize
        Z=SQRT(KK/(K1(1)**2+K1(2)**2+K1(3)**2))  ! set |K2| = |K1|
        DO I=1,3
          K1(I)=K1(I)*Z
        END DO
    ! repeat cycle
        GOTO 50
      endif

C ****  End of multiple scattering cycle  ****

      RETURN
99    PP=0.D0
      END


C--------------------------------------------------------------
      SUBROUTINE WALKSTEP(CR,DIR,R,K,K0,Q,G,PP,TOUT)
C Move to the reflection point in the crystal - version with random walk
C INPUT:
C   Q  .. Qhkl*DW*Ext1
C   K0 .. |K|
C   R,K  .. position and K-vector
C   DIR  .. 1,-1 for directions K,K+G
C OUTPUT:
C   G(3) .. Local diffraction vector at R
C   PP   .. weight to be attributed to the event
C   TOUT .. Time to reach assembly exit
C NOTE: DIR=-1 AND TOUT=0 stops random walk
C--------------------------------------------------------------
      USE COMPONENTS
      IMPLICIT NONE
      INCLUDE 'crystal.inc'
      TYPE (CRYSTAL) :: CR
      REAL*8 EPS
      PARAMETER(EPS=1D-5)
      INTEGER*4 I,M,DIR,I0(3)
      REAL*8 Z,PHI,P,DT,K0,Q,TOUT,PP
      REAL*8 R(3),K(3),R0(3),G(3)
      REAL*4 RAN1,GASDEV1

C special procedure for mosaic crystals without random walk
      if ((CR%model.eq.0).and.(CR%TYP.ne.ctyp_bent)) then
        call WALKSTEP_SIMPLE(CR,DIR,R,K,K0,Q,G,PP,TOUT)
        return
      endif

C add random vertical mosaic angle
      PHI=CR.VMOS*GASDEV1(0.,3.)
C trace through segments and find times, scatt. probabilities, angular dev. etc..
      call CR_BORDERS(CR,R,K,K0)
      IF (NSEG.LE.0) GOTO 99
      call CR_SCATT(CR,DIR,R,K,Q,PHI)

      TOUT=TSEG(NSEG) ! time to the last segment exit
      P=PSEG(NSEG)   ! scattering probability up to the assembly exit
      IF(DIR.LT.0) P=1. ! for K+G direction, consider the possibility of leaving the assembly
      IF (P*PP.LT.EPS) GOTO 99 ! don't follow trajectories with low probability
C Select the next segment according to scattering probabilities
      Z=P*RAN1()
      M=1
      DO WHILE (M.LE.NSEG.AND.Z.GE.PSEG(M))
          M=M+1
      END DO

C No reflection in any segment
      IF(M.GT.NSEG) THEN
        IF(DIR.LT.0) THEN
            GOTO 90 ! go to the last segment exit and return
        ELSE
            GOTO 99 ! event stopped (no reflection)
        ENDIF
      ENDIF

C Find point of reflection inside the M-th segment
      CALL SEGTRACE(CR,K0,Q,M,DT)

      IF (DT.LE.0) GOTO 99

C accumulate flight-path through the material
      DO I=1,M-1
         PATH=PATH+(TSEG(I)-TSEG0(I))*K0
      ENDDO
      PATH=PATH+DT*K0

C accumulate time-of-flight
      TOF=TOF+DT+TSEG0(M)

C Move to the point of reflection
      DO I=1,3
         R(I)=R(I)+(DT+TSEG0(M))*K(I)
      ENDDO

C Coordinates of the segment centre
      CALL SEGCOORD(CR,R,R0,I0)

C Calculate local G-vector
      CALL LOCALG(CR,DIR,R,R0,I0,ALPHA(M)+GRAD(M)*DT,PHI,G)
      PP=PP*P
      RETURN

C exit when going along K+G (DIR<0)
90    DO I=1,3
         R(I)=R(I)+TSEG(NSEG)*K(I)
      ENDDO
      TOF=TOF+TSEG(NSEG)
      DO I=1,NSEG
          PATH=PATH+(TSEG(I)-TSEG0(I))*K0
      ENDDO
      TOUT=0.D0
      RETURN

C Left without reflection
99    PP=0.D0
      TOUT=0.D0

      END SUBROUTINE WALKSTEP

c---------------------------------------------------------------------C
      SUBROUTINE CR_BORDERS(CR,R,K,K0)
C Traces through all segments along the neutron path
C Calculate the segments entry and exit times
C INPUT:
C   R,K starting position and K-vector
C   K0=|K|
C OUTPUT (in /CRBORDERS/)
C   TSEG    time to the I-th segment exit, started at the assembly entry
C   TSEG0   time to the I-th segment entry, -"-
C   NSEG    number of crossed segments
C ---------------------------------------------------------------------C
      USE COMPONENTS
      implicit none
      INCLUDE 'crystal.inc'
      TYPE (CRYSTAL) :: CR
      REAL(KIND(1.D0)),intent(in) :: R(3),K(3),K0
      INTEGER :: I,J,I0(3)
      REAL(KIND(1.D0)) :: T,DT,TIN1,TIN2,TOUT1,TOUT2,V(3),V0(3)

      REAL(KIND(1.D0)), PARAMETER :: LIMIT=0.5D0

      J=0
      TSEG(0)=0.D0
      TSEG0(0)=0.D0
      T=0.D0 ! measures time-fo-flight
      V=R

! trace while inside the volume of the crystal array
      DO WHILE ((ABS(V(1)/CR.FRAME.SIZE(1)).LT.LIMIT).AND.
     *          (ABS(V(2)/CR.FRAME.SIZE(2)).LT.LIMIT).AND.
     *          (ABS(V(3)/CR.FRAME.SIZE(3)).LT.LIMIT).AND.J.LT.MSEG)


C Get entry (TIN) and exit (TOUT) times of a neutron moving along K and starting at V
C (1) .. crossing segment sector
C (2) .. crossing the segment itself (differs by gaps between the segments, CR.DH ...)
        CALL SEGCOORD(CR,V,V0,I0) ! get segment coordinates
        CALL SEGCROSS(CR,V-V0,K,TIN1,TIN2,TOUT1,TOUT2)

C Count only intersected segments
    ! must start before or inside the segment
        IF(((TOUT2-TIN2)*K0.GT.1.D-3).AND.(TOUT2.GT.1.D-10)) THEN
          J=J+1
C Move to the segment entry, if not already inside
          IF (TIN2.GT.0.) THEN
            DO I=1,3
              V(I)=V(I)+TIN2*K(I)
            ENDDO
            T=T+TIN2
          ELSE
            TIN2=0.D0
          ENDIF
          TSEG0(J)=T     ! J-th segment starting time (with respect to the orig. position R)
          TSEG(J)=T+TOUT2-TIN2  ! J-th segment exit time
        ELSE
          TIN2=0.D0
        ENDIF

C Move to the entry of the next segment sector
    ! move slightly behind, we need to start INSIDE the next sector
        DT=1.000001D0*(TOUT1-TIN2)
        DO I=1,3
           V(I)=V(I)+DT*K(I)
        ENDDO
        T=T+DT
      ENDDO
      NSEG=J
      END SUBROUTINE CR_BORDERS

C ---------------------------------------------------------------------C
      SUBROUTINE CR_SCATT(CR,DIR,R,K,Q,PHI)
C Traces through all segments along the neutron path and calculates
C for entry positions of each segment:
C ALPHA   dThetaB, deviation from Bragg angle
C GRAD    grad(ThetaB), gradient of Bragg angle along K
C PSEG    scattering probability
C INPUT:
C   R,K .. starting position and K-vector
C   DIR .. direction along K (1) or K+G (-1)
C   Q   .. Qhkl*DW*Ext1
C   PHI .. vertical mosaic tilt. angle
C   TSEG,TSEG0 .. cross times through all segments on the path (calculated in CR_BORDERS)
C   NSEG .. number of crossed segments
C OUTPUT:
C   stored in /CRBORDERS/
C ---------------------------------------------------------------------C
      USE COMPONENTS
      implicit none
      INCLUDE 'crystal.inc'
      TYPE (CRYSTAL) :: CR
      integer, intent(in) :: DIR
      REAL(KIND(1.D0)),intent(in) :: R(3),K(3),Q,PHI
      INTEGER :: I,J
      REAL(KIND(1.D0)) :: V(3)

      PSEG(0)=0.D0
      do j=1,NSEG
        DO i=1,3
          V(i)=R(i)+TSEG0(j)*K(i)
        ENDDO
        CALL SEGSCATT(CR,j,DIR,V,K,PHI,Q,TSEG(j)-TSEG0(j))
      enddo
      END SUBROUTINE CR_SCATT

C -------------------------------------------------------------
      SUBROUTINE LOCALG(CR,DIR,R,R0,I0,ETA,PHI,G)
C Calculate local G-vector at R
C I0(3) segment coordinates
C R0(3) segment center physical coordinates
C ETA .. horizontal tilt angle of mosaic domain
C PHI .. vertical tilt angle of mosaic domain
C -------------------------------------------------------------
      USE COMPONENTS
      implicit none
      TYPE (CRYSTAL) :: CR
      INTEGER*4 I,DIR,I1,I0(3)
      REAL*8 R(3),R0(3),G(3),W(3),AT(3),G0(3)
      REAL*8 Z,ETA,PHI,GABS

C Calculate local G-vector (only G-gradient)
      DO I=1,3
          G0(I)=CR.G(I)
          IF (CR.MAPG(I)) THEN
              W(I)=0.D0
            DO I1=1,3
                W(I)=W(I)+CR.DG_DR(I,I1)*(R(I1)-R0(I1))
            ENDDO
            G0(I)=G0(I)+W(I)
          ENDIF
      END DO
      GABS=SQRT(G0(1)**2+G0(2)**2+G0(3)**2)

C Add segment tilt angle and vertical mosaic spread
      CALL SEGTILT(CR,I0,AT)
      G(1)=G0(1)-G0(3)*(AT(1)+AT(3))
      G(3)=G0(3)+G0(1)*(AT(1)+AT(3))
      G(2)=G0(2)-G0(3)*AT(2)+GABS*PHI

C Add the angle of the mosaic block
      DO I=1,3
         G(I)=G(I)+GABS*CR.GAMA(I)*ETA
      END DO
C Renormalize
      Z=GABS/SQRT(G(1)**2+G(2)**2+G(3)**2)
      DO i=1,3
         G(i)=DIR*G(i)*Z
      ENDDO
      END

C ----------------------------------------------------------
      SUBROUTINE SEGTILT(CR,I0,AT)
C     return tilt angles
C ----------------------------------------------------------
      USE COMPONENTS
      implicit none
      TYPE (CRYSTAL) :: CR
      REAL*8 AT(3),da
      INTEGER*4 I0(3)

      IF(I0(1).GT.0.AND.I0(1).LE.CR.NH) THEN
           da=CR.FRAME.SIZE(1)*CR.RH/CR.nh
           AT(1)=(-(CR.nh-1)/2.D0+I0(1)-1.D0)*da
      ELSE
           AT(1)=0.D0
      ENDIF
      IF(I0(2).GT.0.AND.I0(2).LE.CR.NV) THEN
           da=CR.FRAME.SIZE(2)*CR.RV/CR.nv
           AT(2)=(-(CR.nv-1)/2.D0+I0(2)-1.D0)*da
      ELSE
           AT(2)=0.D0
      ENDIF
      IF(I0(3).GT.0.AND.I0(3).LE.CR.NB) THEN
           da=CR.FRAME.SIZE(3)*CR.RB/CR.nb
           AT(3)=(-(CR.nb-1)/2.D0+I0(3)-1.D0)*da
      ELSE
           AT(3)=0.D0
      ENDIF
      END SUBROUTINE SEGTILT


C ---------------------------------------------------------------------------
      SUBROUTINE SEGCOORD(CR,R,R0,I0)
C     return coordinates of the segment R0, in which the particle at R resides
C ---------------------------------------------------------------------------
      USE COMPONENTS
      implicit none
      TYPE (CRYSTAL) :: CR
      REAL*8 R(3),R0(3),HALF
      INTEGER*4 I0(3)
      PARAMETER(HALF=0.5D0)

C ih,iv,ib are the integer-coordinates of the closest segment
      I0(1)=NINT((R(1)/CR.FRAME.SIZE(1)+HALF)*CR.NH+HALF)
      I0(2)=NINT((R(2)/CR.FRAME.SIZE(2)+HALF)*CR.NV+HALF)
      I0(3)=NINT((R(3)/CR.FRAME.SIZE(3)+HALF)*CR.NB+HALF)
C get physical coordinates of the segment center
      R0(1)=CR.FRAME.SIZE(1)*(1.D0*(I0(1)-HALF)/CR.NH-HALF)
      R0(2)=CR.FRAME.SIZE(2)*(1.D0*(I0(2)-HALF)/CR.NV-HALF)
      R0(3)=CR.FRAME.SIZE(3)*(1.D0*(I0(3)-HALF)/CR.NB-HALF)
      END SUBROUTINE SEGCOORD


C ---------------------------------------------------------------------------C
      SUBROUTINE SEGCROSS(CR,R,K,TIN1,TIN2,TOUT1,TOUT2)
C return times when the particle (R,K) crossess the boarders
C of the segment sector (TIN1,TOUT1) and the segment itself (TIN2,TOUT2)
C R is relative to the segment centre !
C ---------------------------------------------------------------------------C
      USE COMPONENTS
      implicit none
      TYPE (CRYSTAL) :: CR
      REAL*8 R(3),R0(3),K(3),TIN1,TIN2,TOUT1,TOUT2
      REAL*8 BSZ(3),T1(3),T2(3),DUM
      INTEGER*4 I

C sector size
      BSZ(1)=CR.FRAME.SIZE(1)/CR.NH
      BSZ(2)=CR.FRAME.SIZE(2)/CR.NV
      BSZ(3)=CR.FRAME.SIZE(3)/CR.NB

C Get entry and exit times for segment sector
      DO I=1,3
      IF (ABS(K(I)).GT.1.0D-10) THEN
               T2(I)=(BSZ(I)/2.-R(I))/K(I)
        T1(I)=(-BSZ(I)/2.-R(I))/K(I)
        IF (T1(I).GT.T2(I)) THEN
          DUM=T1(I)
          T1(I)=T2(I)
          T2(I)=DUM
        ENDIF
      ELSE
          T2(I)=1.0D30
          T1(I)=-1.0D30
      ENDIF
      ENDDO
      TIN1=MAX(T1(1),T1(2),T1(3))
      TOUT1=MIN(T2(1),T2(2),T2(3))
C segment size
      BSZ(1)=CR.FRAME.SIZE(1)/CR.NH-CR.DH
      BSZ(2)=CR.FRAME.SIZE(2)/CR.NV-CR.DV
      BSZ(3)=CR.FRAME.SIZE(3)/CR.NB-CR.DB
C Get entry and exit times for the segment itself
      DO I=1,3
      IF (ABS(K(I)).GT.1.0D-10) THEN
               T2(I)=(BSZ(I)/2.D0-R(I))/K(I)
        T1(I)=(-BSZ(I)/2.D0-R(I))/K(I)
        IF (T1(I).GT.T2(I)) THEN
          DUM=T1(I)
          T1(I)=T2(I)
          T2(I)=DUM
        ENDIF
      ELSE
          T2(I)=1.0D30
          T1(I)=-1.0D30
      ENDIF
      ENDDO

      TIN2=MAX(T1(1),T1(2),T1(3))
      TOUT2=MIN(T2(1),T2(2),T2(3))

      END SUBROUTINE SEGCROSS



C-----------------------------------------------------------------------------
      SUBROUTINE BRAGGDEV(CR,DIR,R,K,PHI,THETA,THETA1)
C Calculate local angular deviation from Bragg condition and its gradient
C INPUT:
C   DIR  .. 1,-1 for directions K,K+G
C   PHI  .. vertical mosaic angle
C   R,K  .. neutron coordinates
C OUTPUT:
C   THETA  .. angular deviation
C   THETA1 .. dTheta/dT, T is the flight path in units of [h_bar/m]
C-----------------------------------------------------------------------------
      USE COMPONENTS
      implicit none
      INCLUDE 'crystal.inc'
      TYPE (CRYSTAL) :: CR
      integer,intent(in) :: DIR
      real(KIND(1.D0)), intent(in) :: R(3),K(3),PHI
      real(KIND(1.D0)), intent(out) :: THETA,THETA1

      INTEGER*4 J,I,I1,I0(3),L
      real(KIND(1.D0)) :: DT,a,b,Z,KK,GABS,DEV
      real(KIND(1.D0)) :: KaG(3),G(3),R0(3),V(3)

      KK=0.
      DO I=1,3
         KK=KK+K(I)**2
      END DO

C get coordinates of the segment centre
      CALL SEGCOORD(CR,R,R0,I0)
      CALL LOCALG(CR,DIR,R,R0,I0,0.D0,PHI,G)

C get angular deviation from the Bragg condition
      GABS=0.
      DO I=1,3
         KaG(I)=K(I)+G(I)
         GABS=GABS+G(I)**2
      END DO
      GABS=SQRT(GABS)
      a=KaG(1)**2+ KaG(2)**2+KaG(3)**2-KK    ! (K+G)^2-K^2
      b=DIR*GABS*(KaG(1)*CR.GAMA(1)+KaG(2)*CR.GAMA(2)+KaG(3)*CR.GAMA(3))

C get angular deviation from the Bragg condition (=alpha)
      THETA=-a/(2.*b)
C get gradient of angular deviation ALPHA
      THETA1=0.D0
      DO I=1,3
        IF (CR.MAPG(I)) THEN
          Z=0.D0
          DO I1=1,3
            Z=Z+DIR*CR.DG_DR(I,I1)*K(I1)
          ENDDO
          THETA1=grad(J)+KaG(I)*Z
        ENDIF
      END DO
      THETA1=-THETA1/b

C 2nd order corrections for large bent crystals
      IF (CR%TYP.EQ.ctyp_bent) THEN
        DEV=THETA
        L=10 ! max. L iterations
        DO WHILE ((ABS(DEV).GT.1D-7).AND.(L.GT.0))
          L=L-1
          DO I=1,3
            V(I)=R(I)-THETA/THETA1*K(I)
          ENDDO
          CALL LOCALG(CR,DIR,V,R0,I0,0.D0,PHI,G)
          GABS=0.
          DO I=1,3
            KaG(I)=K(I)+G(I)
            GABS=GABS+G(I)**2
          END DO
          GABS=SQRT(GABS)
          a=KaG(1)**2+ KaG(2)**2+KaG(3)**2-KK    ! (K+G)^2-K^2
          DEV=-a/(2.*b)
          THETA=THETA+DEV
        ENDDO
      ENDIF

C// add random angle within Darwin box width
c        THETA=THETA+(RAN1(1)-0.5)*Q*CR.DEXT*1D-3/PI


      END SUBROUTINE BRAGGDEV


C-----------------------------------------------------------------------------
      SUBROUTINE SEGSCATT(CR,J,DIR,R,K,PHI,Q,DT)
C Calculate scattering probability for J-th segment on the path DT along K(3)
C INPUT:
C   J    .. segment index
C   DIR  .. 1,-1 for directions K,K+G
C   R,K  .. starting neutron coordinates
C   PHI  .. vertical mosaic angle
C   Q    .. kinematical reflectivity (incl. DW)
C   DT   .. time-of-flight along K(3) through the segment
C OUTPUT:
C   stored in /CRBORDERS/
C-----------------------------------------------------------------------------
      USE COMPONENTS
      implicit none
      INCLUDE 'crystal.inc'
      TYPE (CRYSTAL) :: CR
      integer,intent(in) :: J,DIR
      real(KIND(1.D0)), intent(in) :: R(3),K(3),PHI,Q,DT
      real(KIND(1.D0)) :: Z1,Z2,sigma,K0
      real(KIND(1.D0)) :: HMOS_DIST,ERF

      IF (DT.LE.0) THEN
          PSEG0(J)=0.D0
          PSEG(J)=PSEG(J-1)
          RETURN
      ENDIF

      call BRAGGDEV(CR,DIR,R,K,PHI,alpha(J),grad(J))

      K0=SQRT(K(1)**2+K(2)**2+K(3)**2)

      select case(CR%TYP)
C MOSAIC CRYSTAL
      case(ctyp_mosaic)
        sigma=K0*Q*HMOS_DIST(alpha(J)/CR.HMOS)/CR.HMOS*DT
C MOSAIC CRYSTAL WITH G-GRADIENT
      case(ctyp_gradient)
        Z1=ERF((alpha(J)+grad(J)*DT)/CR.HMOS,0)
        Z2=ERF(alpha(J)/CR.HMOS,0)
        sigma=K0*Q*(Z1-Z2)/grad(J)
        IF ((J.EQ.1).AND.(ABS(ALPHA(J)).LT.SEC)) sigma=0.D0 ! no back refl. in the same segment
C ONLY G-GRADIENT
      case(ctyp_bent)
        Z1=-alpha(J)/grad(J)
        IF ((Z1.GT.0).AND.(Z1.LT.DT)) THEN
          sigma=K0*Q/abs(grad(J))
        ELSE
          sigma=0.D0
        ENDIF
        IF ((J.EQ.1).AND.(ABS(ALPHA(J)).LT.SEC)) sigma=0.D0 ! no back refl. in the same segment
      END select

C// get scattering probability PSEG0 (fr this segment) and PSEG (cumulative)
      IF(sigma.gt.14) then
        PSEG0(J)=1.D0
        PSEG(J)=1.D0
      else
        PSEG0(J)=1.D0-exp(-sigma)
        PSEG(J)=1.D0-(1.D0-PSEG(J-1))*(1.D0-PSEG0(J))
      endif

      END SUBROUTINE SEGSCATT


C-------------------------------------------------------------------------
      SUBROUTINE SEGTRACE(CR,K0,Q,M,DT)
C return random step size to a point of reflection inside the M-th segment
C-------------------------------------------------------------------------
      USE COMPONENTS
      IMPLICIT NONE
      INCLUDE 'crystal.inc'
      TYPE (CRYSTAL) :: CR
      INTEGER*4 M
      REAL*8 SEC,ksi,Z,Z1,sigma,P0,DT,AA
      PARAMETER(SEC=4.85D-6)
      REAL*8 K0,Q
      REAL*4 RAN1
      REAL*8 ERF,HMOS_DIST

      DT=0.D0

C MOSAIC CRYSTAL
      IF(ABS(GRAD(M)).LT.1.D-6*K0*Q) THEN
          ksi=RAN1()
          P0=PSEG0(M)
          AA=ALPHA(M)/CR.HMOS
          sigma=K0*Q*HMOS_DIST(AA)/CR.HMOS
          DT=-Log(1-ksi*P0)/sigma
      ELSE
C MOSAIC CRYSTAL WITH D-GRADIENT
        IF(CR.HMOS.GT.SEC) THEN
          ksi=RAN1()
          P0=PSEG0(M)
          Z=GRAD(M)*Log(1-ksi*P0)/(K0*Q)
          AA=ALPHA(M)/CR.HMOS
C  ERF shoud be used only on (-inf.;0) because of num. precision
          IF(AA.GT.0) THEN
            Z1=ERF(-AA,0)
            IF(1.-Z1-Z.GT.0.5D0) THEN
              DT=-ERF(Z1+Z,1)-AA
            ELSE
              DT=ERF(1.-Z1-Z,1)-AA
            ENDIF
          ELSE
            Z1=ERF(AA,0)
            IF(Z1-Z.GT.0.5D0) THEN
              DT=-ERF(1.-Z1+Z,1)-AA
            ELSE
              DT=ERF(Z1-Z,1)-AA
            ENDIF
          ENDIF
           DT=DT*CR.HMOS/GRAD(M)
C ONLY D-GRADIENT
        ELSE
          DT=-ALPHA(M)/GRAD(M)
        ENDIF
      ENDIF

      END SUBROUTINE SEGTRACE

C-----------------------------------------------------
      SUBROUTINE ERF_INIT(F,AMIN,AMAX)
C Calculate lookup tables for ERF function
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
            pause
          ENDIF
          Z1=XJ+DX*(-B+SQRT(DET))/2/A
          Z2=XJ+DX*(-B-SQRT(DET))/2/A
          IF (ABS(Z2-XJ).LT.ABS(Z1-XJ)) Z1=Z2
          Y1(I)=Z1
        ENDIF
!20    continue
      END DO
      END

C----------------------------------------------------------------
      REAL*8 FUNCTION ERF(ARG,INV)
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
        ERF=ERF_INTERP(XMIN,DX,Y,DIM,ARG)
C else Iverse cumul. function
      ELSE
        IF(ARG.LE.XMIN1) THEN
           ERF=Y1(1) ! left limit
        ELSE IF(ARG.GE.XMIN1+(DIM-1)*DX1) THEN
           ERF=Y1(DIM)  ! right limit
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
             ERF=XJ
          ELSE IF (ABS(A).LT.1D-30) THEN
            ERF=XJ+Y(J+1)
            J=J-1
            GOTO 10
          ELSE
            DET=B**2-4*A*(Y(J)-ARG1)
            DET=SQRT(DET)
            Z1=XJ+DX*(-B+DET)/2/A
            Z2=XJ+DX*(-B-DET)/2/A
            IF (ABS(Z2-XJ).LT.ABS(Z1-XJ)) Z1=Z2
            IF(ARG.GT.9.D-1) Z1=-Z1
            ERF=Z1
          ENDIF
        ENDIF
      ENDIF
      END

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
      end


C--------------------------------------------------
      LOGICAL*4 FUNCTION CR_INSIDE(CRYST,R)
C INSIDE function for CRSYTAL object ... takes into
C account curved surface of a bent crystal plate.
C NOT USED IN THE CURRENT VERSION !
C--------------------------------------------------
      USE COMPONENTS
      implicit none

      TYPE (CRYSTAL) :: CRYST
      REAL*8 R(3),R0(3)
      LOGICAL*4 INSIDE

      R0(3)=R(3)-R(1)**2*CRYST.RH-R(2)**2*CRYST.RV
      R0(1)=R(1)
      R0(2)=R(2)
      CR_INSIDE=INSIDE(CRYST.FRAME,R0)

      END

