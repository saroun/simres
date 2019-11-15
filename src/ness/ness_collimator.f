C $Id: ness_collimator.f,v 1.54 2019/08/16 17:16:26 saroun Exp $
C//////////////////////////////////////////////////////////////////////
C////                                                              ////
C////  NEutron Scattering Simulation - v.2.0, (c) J.Saroun, 2000   ////
C////                                                              ////
C//////////////////////////////////////////////////////////////////////
C////
C////  Subroutines describing objects - ALL COLLIMATOR TYPES
C////  Envelope for collimator segments of any kind:
C////  TUBE, SOLLER , GUIDE/BENDER, PARABOLIC GUIDE
C////
C//////////////////////////////////////////////////////////////////////


C---------------------------------------------------------------------
      REAL(kind(1.D0))FUNCTION BENDER_LAM(OBJ,ID,IL,Z,IO)
C function describing lamella profile
C ID   ... derivative
C IL   ... lamella index
C Z    ... distance from guide entry
C IO   ... horizontal (1) or vertical (2)
C----------------------------------------------------------------------
      USE GUIDES
      IMPLICIT NONE
      TYPE(GUIDE) :: OBJ
      INTEGER :: ID,IL,IO
      REAL(kind(1.D0)) :: PGUIDE_LAM,EGUIDE_LAM,GUIDE_LAM,Z

      select case (OBJ%TYP)
      case(2,3)
         BENDER_LAM=PGUIDE_LAM(OBJ,ID,IL,Z,IO)
      case(4)
         BENDER_LAM=EGUIDE_LAM(OBJ,ID,IL,Z,IO)
      case default
         BENDER_LAM=GUIDE_LAM(OBJ,ID,IL,Z,IO)
      END select
      END FUNCTION BENDER_LAM


C-------------------------------------------------------------------
      SUBROUTINE QUADREQ(A,B,C,X)
C Solve quadratic equation A*X^2 + B*X + C = 0
C Try to find a solution > EPS=1E-10
C 1) no solution .. return 10^30
C 2) 1 solution  .. return this
C 3) 2 solutions .. return the smaller one
C-------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 EPS
      PARAMETER (EPS=1.D-8)
      REAL*8 A,B,C,X,X1,X2,DET

      IF (ABS(A).LT.EPS**2) THEN
        IF (ABS(B).LT.EPS) THEN
          GOTO 20
        ELSE
          X=-C/B
          GOTO 30
        ENDIF
      ELSE
        DET=B**2-4*A*C
        IF (DET.EQ.0.) THEN
           X=-B/2./A
           GOTO 30
        ELSE IF (DET.LT.0.) THEN
           GOTO 20
        ELSE
           DET=SQRT(DET)
           X1=(-B+DET)/2./A
           X2=(-B-DET)/2./A
           IF (X1.GT.EPS.AND.X2.GT.EPS) THEN
              X=MIN(X1,X2)
           ELSE IF (X1.GT.EPS) THEN
              X=X1
           ELSE IF (X2.GT.EPS) THEN
              X=X2
           ELSE
              GOTO 20
           ENDIF
        ENDIF
      ENDIF

30    IF (X.LT.EPS) GOTO 20
      RETURN

20    X=1.D30

      END SUBROUTINE QUADREQ

C--------------------------------------------------------
      LOGICAL*4 FUNCTION BENDER_PASS(OBJ,R,EXCL,IH,IV)
c Checks, whether neutron fits inside any slit and
c returns slit indices
! if EXCL, exclude outer lamella thickness
C--------------------------------------------------------
      USE GUIDES
      IMPLICIT NONE
      TYPE(GUIDE) :: OBJ
      logical :: EXCL
      INTEGER*4 IN_LAM,IH,IV,BENDER_INLAM
      REAL*8 R(3)

      IN_LAM=BENDER_INLAM(OBJ,R,EXCL,IH,IV)
      BENDER_PASS=(IN_LAM.EQ.0)

      END FUNCTION BENDER_PASS


C--------------------------------------------------------
      INTEGER*4 FUNCTION BENDER_INLAM(OBJ,R,EXCL,IH,IV)
c Checks, whether neutron is inside an inner lamella
C Returns the result: no (0) horizontal (1) vertical (2) both (3)
C Returns -1 if neutron is outside the collimator area
c Returns the index of corresponding lamella in IH or IV or both
! if EXCL, exclude outer lamella thickness
C--------------------------------------------------------
      USE GUIDES
      IMPLICIT NONE
      logical :: EXCL
      INTEGER*4 IH,IV
      INTEGER*4 I
      LOGICAL*4 IN_H,IN_V,IN_AREA
      REAL*8 BENDER_LAM
      REAL*8 R(3),jh,jv,W,H
      TYPE(GUIDE) :: OBJ
      REAL*8 OSCD,OSCA
      REAL(KIND(1.D0)) :: WMI,WMA,HMI,HMA,HC,WC
      COMMON /OSCBEND/ OSCD,OSCA
      logical :: dbg
      common /guide_dbg/ dbg

!      dbg=.true.
1     format('BENDER_INLAM: ',a,' ',3(G10.4,2x))
      if (dbg) write(*,1) 'R',R(1:3)

! get width and height of the collimator at the neutron position
  !    W=2.D0*ABS(BENDER_LAM(OBJ,0,0,R(3),1))
  !    H=2.D0*ABS(BENDER_LAM(OBJ,0,0,R(3),2))

  ! limits of the guide horiz. (W) and vert. (H)
      WMI=BENDER_LAM(OBJ,0,0,R(3),1)
      WMA=BENDER_LAM(OBJ,0,OBJ%NLH,R(3),1)
      HMI=BENDER_LAM(OBJ,0,0,R(3),2)
      HMA=BENDER_LAM(OBJ,0,OBJ%NLV,R(3),2)
  ! corresponding widths and centers
      W=ABS(WMA-WMI)
      H=ABS(HMA-HMI)
      WC=(WMA+WMI)/2.D0
      HC=(HMA+HMI)/2.D0

      if (dbg) write(*,1) 'W,H,NLH,NLV',W,H,OBJ%NLH,OBJ%NLV
      jh=((R(1)-OSCD-R(3)*OSCA-WC)/W+0.5)*OBJ%NLH
      jv=((R(2)-HC)/H+0.5)*OBJ%NLV
      IH=NINT(jh-0.5)
      IV=NINT(jv-0.5)

      if (dbg) write(*,1) 'jh,jv',jh,jv
      if (dbg) write(*,1) 'IH,IV',IH,IV
      if (dbg) write(*,1) 'dx,dlx',ABS(jh-NINT(jh))*W/OBJ%NLH,OBJ%DLH/2.
      if (dbg) write(*,1) 'dy,dly',ABS(jv-NINT(jv))*H/OBJ%NLV,OBJ%DLV/2.

      if (OBJ%LAMONLY) then
        IN_AREA=.true.
      else
        IN_AREA=((jh.GT.0).AND.(jv.GT.0).AND.(jh.LT.OBJ%NLH).AND.(jv.LT.OBJ%NLV))
      endif

      if (dbg) write(*,*) 'IN_AREA',IN_AREA

! inside a lamella ?
      IN_H=(ABS(jh-NINT(jh))*W/OBJ%NLH.LT.OBJ%DLH/2.)
! if yes, exclude side wall and set IH to lamella index instead of channel index
      if (IN_H) then
        I=NINT(jh) ! lamella index
      ! side wall ? => no lamella, IH is channel index
        if (EXCL.and.(I.LE.0).or.(I.GE.OBJ%NLH)) then
          IN_H=.FALSE.
      ! otherwise set IH to lamella index
        else
          IH=I
        endif
      endif
      if (dbg) write(*,*) 'IN_H, IH',IN_H,IH
! and the same for vertical ....
      IN_V=(ABS(jv-NINT(jv))*H/OBJ%NLV.LT.OBJ%DLV/2.)
      if (IN_V) then
        I=NINT(jv)
        if (EXCL.and.(I.LE.0).or.(I.GE.OBJ%NLV)) then
          IN_V=.FALSE.
        else
          IV=I
        endif
      endif
      if (dbg) write(*,*) 'IN_V, IV',IN_V,IV

      if (.NOT.IN_AREA) then
        I=-1
      else IF (IN_H.AND.IN_V) THEN
        I=3
      ELSE IF (IN_V) THEN
        I=2
      ELSE IF (IN_H) THEN
        I=1
      ELSE
        I=0
      ENDIF

      BENDER_INLAM=I
      if (dbg) write(*,1) 'result',I
      if (dbg) write(*,*)
!       read(*,*)
      END FUNCTION BENDER_INLAM

C---------------------------------------------------
      REAL*8 FUNCTION BENDER_REF_OLD(ID,OBJ,Q,S)
c returns reflectivity for given momentum transfer Q
C ID identifies which surface is touched
C ID=0  left
C ID=1  right
C ID=2  top
C ID=3  bottom
C---------------------------------------------------
      USE GUIDES
      IMPLICIT NONE
  !    INCLUDE 'ness_common.inc'
      integer, parameter :: MIRROR_DIM=25
      INTEGER*4 m_nrow
      PARAMETER(m_nrow=256)
      INTEGER*4 ID
      TYPE(GUIDE) :: OBJ
      REAL*8 Q,S
      INTEGER*4 iz,NR
      REAL*8 z,dQ,Q1,gamma,R
      INTEGER*4 m_n(MIRROR_DIM)
      CHARACTER*3 m_name(MIRROR_DIM)
      REAL*8 m_alpha(m_nrow,MIRROR_DIM), m_ref1(m_nrow,MIRROR_DIM),m_ref2(m_nrow,MIRROR_DIM)
      COMMON /MIRROR/ m_alpha,m_ref1,m_ref2,m_n,m_name

      logical*4 dbg
      NR=0

c      dbg=(gdea.FRAME%count.eq.90076)
      dbg=.false.
C get critical angle, reflectivity and index to lookup table
      select case(ID)
      case(0)
        IF (S.GE.0) THEN
          gamma=OBJ%GHLU
          R=OBJ%RHLU
          NR=OBJ%NHLU
        ELSE
          gamma=OBJ%GHLD
          R=OBJ%RHLD
          NR=OBJ%NHLD
        endif
      case(1)
        IF (S.GE.0) THEN
          gamma=OBJ%GHRU
          R=OBJ%RHRU
          NR=OBJ%NHRU
        ELSE
          gamma=OBJ%GHRD
          R=OBJ%RHRD
          NR=OBJ%NHRD
        ENDIF
      case(2)
          gamma=OBJ%GVT
          R=OBJ%RVT
          NR=OBJ%NVT
      case(3)
          gamma=OBJ%GVB
          R=OBJ%RVB
          NR=OBJ%NVB
      case DEFAULT
          gamma=0.
          R=0.
      END select

      if (dbg) write(*,*) 'gamma, Q, NR ',gamma, Q, NR

C no reflection for Q<0 or gamma=0
      IF (gamma.LE.0.OR.Q.LT.0) THEN
          BENDER_REF_OLD=0
          RETURN
      ENDIF

C no lookup table, just step function for 0 < Q < 2*pi*gamma
      IF (NR.LE.0.OR.NR.GT.MIRROR_DIM) THEN
        IF (Q.GE.0.AND.Q.LT.2*PI*gamma) THEN
          BENDER_REF_OLD=R
        ELSE
          BENDER_REF_OLD=0
        ENDIF
      ELSE
        Q1=Q/2/PI/GammaNi
        dQ=(m_alpha(m_n(NR),NR)-m_alpha(1,NR))/(m_n(NR)-1)
        z=(Q1-m_alpha(1,NR))/dQ
        iz=INT(z)+1
        if (dbg) write(*,*) 'Q1, dQ, z, iz ',Q1, dQ, z, iz
        IF (z.LT.0.OR.z.GE.m_n(NR).OR.iz.GE.m_n(NR)) THEN
          BENDER_REF_OLD=0
        ELSE IF (S.GE.0) THEN
          BENDER_REF_OLD=m_ref1(iz,NR)+(z-iz+1)*(m_ref1(iz+1,NR)-m_ref1(iz,NR))
        ELSE IF (S.LT.0) THEN
          BENDER_REF_OLD=m_ref2(iz,NR)+(z-iz+1)*(m_ref2(iz+1,NR)-m_ref2(iz,NR))
        ELSE
          BENDER_REF_OLD=0
        ENDIF
      ENDIF
      END FUNCTION BENDER_REF_OLD


C---------------------------------------------------
      REAL*8 FUNCTION BENDER_REF(ID,OBJ,Q,S)
c returns reflectivity for given momentum transfer Q
C ID identifies which surface is touched
C ID=0  left
C ID=1  right
C ID=2  top
C ID=3  bottom
C---------------------------------------------------
      USE GUIDES
      use MIRROR_TABLE
      IMPLICIT NONE
      INTEGER,intent(in) :: ID
      TYPE(GUIDE) :: OBJ
      REAL(KIND(1.D0)),intent(in) :: Q,S
      INTEGER :: NR
      REAL(KIND(1.D0)) :: PR,MC

      call GUIDE_GET_MTAB(ID,OBJ,NR,MC)
      PR=MIRROR_REF(NR,MC,Q,S)
      BENDER_REF=PR
      END FUNCTION BENDER_REF


C--------------------------------------------------------------
      REAL*8 function BENDER_DIV(OBJ,LAMBDA,IDIR)
C Calculate divergence allowed by the bender
C LAMBDA ... wavelength
C IDIR   ... direction: horizontal (1) or vertical (2)
C--------------------------------------------------------------
      USE GUIDES
      IMPLICIT NONE

      TYPE(GUIDE) :: OBJ
      REAL(KIND(1.D0)),intent(in) :: LAMBDA
      integer, intent(in) :: IDIR
      REAL(KIND(1.D0)) :: AMAX,AREF,W1,W2,AL1,AL2,G
      REAL(KIND(1.D0)) :: BENDER_LAM
      INTEGER :: NL

! no limit for a transmitting object
      if ((OBJ%TYP.LT.0).or.(OBJ%FRAME%SIZE(3).LE.0.D0)) then
        BENDER_DIV=1.D30
        RETURN
      endif

! get bender parameters for given direction
! AMAX ... maximum deflection by curved walls
! G    ... critical angle
      select case (IDIR)
      case(2)
        W1=OBJ%FRAME%SIZE(2)
        W2=OBJ%H2
        NL=OBJ%NLV
        G=MAX(OBJ%GVT,OBJ%GVB)
        AL1=ABS(BENDER_LAM(OBJ,1,0,0.D0,2))
        AL2=ABS(BENDER_LAM(OBJ,1,0,OBJ%FRAME%SIZE(3),2))
        AMAX=MAX(AL1,AL2)
      case default
        W1=OBJ%FRAME%SIZE(1)
        W2=OBJ%W2
        NL=OBJ%NLH
        G=MAX(OBJ%GHLU,OBJ%GHRU,OBJ%GHLD,OBJ%GHRD)
        AL1=ABS(BENDER_LAM(OBJ,1,0,0.D0,1))
        AL2=ABS(BENDER_LAM(OBJ,1,0,OBJ%FRAME%SIZE(3),1))
        AMAX=MAX(AL1,AL2)
      end select

      IF (NL.LE.0) then
        AMAX=(W1+W2)/OBJ%FRAME%SIZE(3)
      else
        select case (OBJ%TYP)
! Soller
        case(0)
  ! soller slits
          AMAX=(W1+W2)/NL/OBJ%FRAME%SIZE(3)
! guide/bender
        case DEFAULT
          AREF=2*AMAX+3*G*LAMBDA
          AMAX=MAX(AREF,(W1+W2)/NL/OBJ%FRAME%SIZE(3))
        end select
      endif

      BENDER_DIV=AMAX
      end function BENDER_DIV


C--------------------------------------------------------------
      SUBROUTINE BENDER_SETM(OBJ,M,ISIDE)
C Set critical angle for horizontal or vertical walls
C Adjust dependent parameters
C M ... m-value in [Ni_nat]
C ISIDE ... (1) horizontal (2) vertical
C--------------------------------------------------------------
      USE GUIDES
      use MIRROR_TABLE
      IMPLICIT NONE
      TYPE(GUIDE) :: OBJ
      real(KIND(1.D0)),intent(in) :: M
      integer,intent(in) :: ISIDE
      !INTEGER :: READ_MIRROR
      select case(ISIDE)
      case(1)
        OBJ%GHLU=M*GammaNi
        OBJ%GHRU=OBJ%GHLU
        OBJ%GHLD=OBJ%GHLU
        OBJ%GHRD=OBJ%GHLU
        !OBJ%NHLU=READ_MIRROR(M)
        OBJ%NHLU=MIRROR_GETTABLE(M)
        OBJ%NHLD=OBJ%NHLU
        OBJ%NHRU=OBJ%NHLU
        OBJ%NHRD=OBJ%NHLU
        if ((OBJ%GHLU.gt.0.D0).and.(OBJ%TYP.EQ.0)) OBJ%TYP=1
      case(2)
        OBJ%GVT=M*GammaNi
        OBJ%GVB=OBJ%GVT
        !OBJ%NVT=READ_MIRROR(M)
        OBJ%NVT=MIRROR_GETTABLE(M)
        OBJ%NVB=OBJ%NVT
        if ((OBJ%GVT.gt.0.D0).and.(OBJ%TYP.EQ.0)) OBJ%TYP=1
      end select
      end SUBROUTINE BENDER_SETM

C--------------------------------------------------------------
      SUBROUTINE BENDER_SETREF(OBJ,REF,ISIDE)
C Set reflectivity for horizontal or vertical walls
C Adjust dependent parameters
C REF ... reflectivity value
C ISIDE ... (1) horizontal (2) vertical
C--------------------------------------------------------------
      USE GUIDES
      IMPLICIT NONE
      TYPE(GUIDE) :: OBJ
      real(KIND(1.D0)),intent(in) :: REF
      integer,intent(in) :: ISIDE

      select case(ISIDE)
      case(1)
        OBJ%RHLU=REF
        OBJ%RHRU=REF
        OBJ%RHLD=REF
        OBJ%RHRD=REF
      case(2)
        OBJ%RVT=REF
        OBJ%RVB=REF
      end select
      end SUBROUTINE BENDER_SETREF


