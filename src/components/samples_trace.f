!//////////////////////////////////////////////////////////////////////
!////  $Id: samples_trace.f,v 1.18 2019/08/15 15:02:07 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2011, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.18 $
!////     $Date: 2019/08/15 15:02:07 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Ray-tracing subroutines for SAMPLE
!////
!////////////////////////////////////////////////////////////////////////
      MODULE SAMPLES_TRACE
      use TRACELIB
      USE TRACINGDATA
      USE SAMPLES
      use EVENTLOGS
      use FRAMES_TRACE
      use RNDGEN
      use SAMPLE_POLY_TABLE
      implicit none
      private

      public SAMPLE_GO_INST,SAMPLE_ADJUST

      logical,parameter :: dbg=.false.

      contains


!-------------------------------------------------------------
      logical function SAMPLE_GO_INST(OBJ,PATH,SIGMA)
! call GO on given instance
!-------------------------------------------------------------
      type(TSAMPLE) :: OBJ
      REAL(kind(1.D0)),intent(in) :: PATH
      REAL(kind(1.D0)),intent(out) :: SIGMA
      logical :: LOG1
      REAL(kind(1.D0)) :: XR,lambda,P0

      LOG1=.true.
      lambda=TWOPI/NEUT%K0
      SIGMA=OBJ%SIGA*lambda+OBJ%SIGI+OBJ%SIGSC
  ! nothing to do
      if (.not.(OBJ%TRANS.or.OBJ%SCATT)) goto 99
  ! if only transmission is allowed:
  ! Simres does not integrate over sample thickness and phi if SCATT=false
  ! therefore no weight is necessary
      if (.not.OBJ%SCATT) goto 110

      P0=1.D0
  ! both scattering and transmission:
  ! choose one of them with equal probability
      if (OBJ%TRANS.and.OBJ%SCATT) then
        P0=2.D0
        XR=RAN1()
        ! transmission
        if (XR<0.5) then
          NEUT%P=NEUT%P*P0/PATH
          goto 110
        endif
        ! else continue with scattering
      endif
      SELECT case (OBJ%TYP)
      case(isam_VANAD)
        LOG1=SAM_VANAD(OBJ)
      case default
        LOG1=SAM_GENERIC(OBJ)
      end SELECT
      if (.not.LOG1) goto 99
      NEUT%P=NEUT%P*P0*OBJ%SIGSC
110   SAMPLE_GO_INST=.true.
      RETURN
99    SAMPLE_GO_INST=.false.
      NEUT%P=0.D0
      end function SAMPLE_GO_INST


!-------------------------------------------------
      SUBROUTINE SAMPLE_ADJUST(OBJ,IERR)
! Update position of the component using nominal trajectory represented by NEUT
! Call FRAME_INIT_MAT first.
! Then use SLIT_PRE & SLIT_POST procedures to propagate NEUT
!-------------------------------------------------
      integer,intent(out) :: IERR
      TYPE(TSAMPLE) :: OBJ
      REAL(kind(1.D0)) :: Z
      integer :: i
      IERR=0
      call FRAME_INIT_MAT(OBJ%FRAME)
      call SLIT_PRE(OBJ%FRAME)
      call SLIT_TURN(OBJ%FRAME)
      ! change energy
      if (NEUT%K0.gt.0.D0) then
        Z=OBJ%KF0/NEUT%K0
        do i=1,3
          NEUT%K(i)=NEUT%K(i)*Z
        enddo
        NEUT%K0=OBJ%KF0
      endif
      call SLIT_POST(OBJ%FRAME)
      ! call REPORT_SIGTOT(6, 0.5D0, 5.1D0, 461)
      end SUBROUTINE SAMPLE_ADJUST

!--------------------------------------------------------------------
      LOGICAL FUNCTION SAM_GENERIC(OBJ)
! generic sample for resolution function, R(Q)
! Spreads neutrons to all angles
! Works in Q-coordinates
! integration over theta_H,theta_V (= horizontal and vertical angles))
!--------------------------------------------------------------------
      TYPE(TSAMPLE) :: OBJ
      REAL(kind(1.D0)) :: cos2,a1,a2
      REAL(kind(1.D0)) :: EN,KK
1     format(a12,': ',6(G10.4,1x))
  ! total removal cross-section (capture and incoherent)

! get current nominal (!) angles of kf
      !a1=acos(OBJ%KF(3)/OBJ%KF0)
      if (dbg)  write(*,1) 'GENERIC,IN ',NEUT%R,NEUT%K
      if (dbg)  write(*,1) 'GENERIC,NOM ',OBJ%KF,OBJ%KF0
      a1=atan2(OBJ%KF(1),OBJ%KF(3))
      ! NOTE: KF is in QMAT coordinates, therefore KF(2)=0 always !!
      a2=asin(OBJ%KF(2)/OBJ%KF0)
      if (dbg)  write(*,1) 'nom. a1,a2=',a1/deg,a2/deg
      a1=a1+XRND(IXRND(OBJ%IRND(1)))
      a2=a2+XRND(IXRND(OBJ%IRND(2)))
      if (dbg)  write(*,1) 'act. a1,a2=',a1/deg,a2/deg
    ! integrate over the interval a1 = -PI..+PI, a2=-PI/2 .. +PI/2
      if (abs(a1)>PI) goto 99
      if (abs(a2)>PI/2.D0) goto 99
      cos2=cos(a2)
    ! add energy transfer
      if (OBJ%IRND(3).gt.0) then
        EN=XRND(IXRND(OBJ%IRND(3)))
        if (dbg)  write(*,1) 'act. EN=',EN
        KK=NEUT%K0**2-(EN+OBJ%ESC)/HSQOV2M
        if (KK.LE.0.D0) goto 99
        NEUT%K0=SQRT(KK)
      endif
! construct new kf-vector
      NEUT%K(1)=NEUT%K0*sin(a1)*cos2
      NEUT%K(2)=NEUT%K0*sin(a2)
      NEUT%K(3)=NEUT%K0*cos(a1)*cos2
      if (dbg)  write(*,1) 'GENERIC,OUT ',NEUT%R,NEUT%K
      NEUT%P=NEUT%P*cos2
      SAM_GENERIC=.true.
      return
99    NEUT%P=0.D0
      SAM_GENERIC=.false.
      END function SAM_GENERIC


!--------------------------------------------------------------------
      LOGICAL FUNCTION SAM_GENERIC_OLD(OBJ)
! generic sample for resolution function, R(Q)
! Spreads neutrons to all angles
! Works in Q-coordinates
!--------------------------------------------------------------------
      TYPE(TSAMPLE) :: OBJ
      REAL(kind(1.D0)) :: sin1,a1,a2
      REAL(kind(1.D0)) :: theta,phi,EN,KK,SG
1     format(a12,': ',6(G10.4,1x))
      THETA=XRND(IXRND(OBJ%IRND(1)))
      PHI=XRND(IXRND(OBJ%IRND(2)))
! check limits
      if (dbg) then
          write(*,1) 'SAM_GENERIC ',NEUT%R,NEUT%K
          write(*,1) 'THETA,PHI=',THETA/deg,PHI/deg
          write(*,1) 'QSC=',OBJ%QSC
          write(*,1) 'QN= ',OBJ%QN
          write(*,1) 'KI0,KF0 =',OBJ%KI0,OBJ%KF0
          write(*,1) 'NOM KF =',OBJ%KF
      endif
! get current nominal (!) angles of kf
      a1=acos(OBJ%KF(3)/OBJ%KF0)
      !a1=atan2(OBJ%KF(1),OBJ%KF(3))
      a2=atan2(OBJ%KF(2),OBJ%KF(1))
      SG=sign(1.D0,a1)
      if (dbg)  write(*,1) 'nom. a1,a2=',a1/deg,a2/deg
      a1=a1+THETA
      a2=a2+PHI
    ! integrate over the interval a1 = 0 .. PI, a2=-PI .. +PI
      if ((SG*a1<0.D0).or.(SG*a1>PI)) goto 99
      if (ABS(a2).gt.PI) goto 99
      if (dbg)  write(*,1) '     a1,a2=',a1/deg,a2/deg
      sin1=sin(a1)
    ! add energy transfer
      if (OBJ%IRND(3).gt.0) then
        EN=XRND(IXRND(OBJ%IRND(3)))
        KK=NEUT%K0**2-(EN+OBJ%ESC)/HSQOV2M
        if (KK.LE.0.D0) goto 99
        if (dbg)  write(*,1) 'dE,KI,KF=',EN,NEUT%K0,SQRT(KK)
        NEUT%K0=SQRT(KK)
      endif
! construct new kf-vector
      NEUT%K(1)=NEUT%K0*sin1*cos(a2)
      NEUT%K(2)=NEUT%K0*sin1*sin(a2)
      NEUT%K(3)=NEUT%K0*cos(a1)
      if (dbg)  write(*,1) 'GENERIC,OUT ',NEUT%R,NEUT%K
      NEUT%P=NEUT%P*sin1
      SAM_GENERIC_OLD=.true.
      return
99    NEUT%P=0.D0
      SAM_GENERIC_OLD=.false.
      END function SAM_GENERIC_OLD

!--------------------------------------------------------------------
      LOGICAL FUNCTION SAM_VANAD(OBJ)
! Vanad
! Spreads neutrons to elastically to all angles
! Works in Q-coordinates
! integration over theta,phi (= azimuthal angle)
!--------------------------------------------------------------------
      TYPE(TSAMPLE) :: OBJ
      REAL(kind(1.D0)) :: sin1,a1,a2
      REAL(kind(1.D0)) :: theta,phi
1     format(a12,': ',6(G10.4,1x))
      THETA=XRND(IXRND(OBJ%IRND(1)))
      PHI=XRND(IXRND(OBJ%IRND(2)))
! check limits
      if (ABS(PHI).ge.PI) goto 99
      if (ABS(THETA).gt.PI/2) goto 99
! get current nominal (!) angles of kf
      a1=acos(OBJ%KF(3)/OBJ%KF0)
      a2=atan2(OBJ%KF(2),OBJ%KF(1))
      a1=a1+THETA
      a2=a2+PHI
      sin1=sin(a1)
! construct new kf-vector
      NEUT%K(1)=NEUT%K0*sin1*cos(a2)
      NEUT%K(2)=NEUT%K0*sin1*sin(a2)
      NEUT%K(3)=NEUT%K0*cos(a1)
      NEUT%P=NEUT%P*sin1
      SAM_VANAD=.true.
      return
99    NEUT%P=0.D0
      SAM_VANAD=.false.
      END function SAM_VANAD

!--------------------------------------------------------------------
      LOGICAL FUNCTION SAM_TRANS(OBJ)
! transmission through sample, handles beam attenuation in a container with
! given scattering and absorption cross-section
!--------------------------------------------------------------------
      TYPE(TSAMPLE) :: OBJ
      LOGICAL :: LOG1
      REAL(kind(1.D0)) :: T1,T2,DELTA,P1

      LOG1=.FALSE.
      P1=1.D0
! transform to local coordinates:
      CALL SLIT_PRE(OBJ%FRAME)
      IF (NEUT%P.LE.0.D0) goto 99
! check for the entry and exit times
      IF (FRAME_BORDER(OBJ%FRAME,NEUT%R,NEUT%K,T1,T2)) THEN
  ! TOF inside the sample
        DELTA=T2-T1
        SELECT case(OBJ%TYP)
        case(isam_VANAD)
          P1=exp(-1.D-1*DELTA*NEUT%K0*(OBJ%SIGSC+OBJ%SIGA))
        end SELECT
        NEUT%P=NEUT%P*P1
! transform back to axis coordinates
        CALL SLIT_POST(OBJ%FRAME)
        LOG1=(P1.GT.1.D-10)
      else
        NEUT%P=0.D0
      endif
      if (.NOT.LOG1) goto 99
      call incCounter(OBJ%FRAME%COUNT)
      SAM_TRANS=.TRUE.
      RETURN

99    NEUT%P=0.D0
      SAM_TRANS=.FALSE.
      END FUNCTION SAM_TRANS

      end MODULE SAMPLES_TRACE

