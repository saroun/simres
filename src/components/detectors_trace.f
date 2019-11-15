!//////////////////////////////////////////////////////////////////////
!////  $Id: detectors_trace.f,v 1.29 2019/08/15 15:02:06 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2011, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.29 $
!////     $Date: 2019/08/15 15:02:06 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Ray-tracing subroutines for DETECTOR
!////
!////////////////////////////////////////////////////////////////////////
      MODULE DETECTORS_TRACE
      use TRACELIB
      USE TRACINGDATA
      USE DETECTORS
      use EVENTLOGS
      use FRAMES_TRACE
      use RNDGEN
      use NSTORE
      implicit none
      private

      type(NEUTRON) :: TMPN(DETECTORS_DIM)

      public DETECTOR_GO,DETECTOR_ADJUST,DETECTOR_INIT,DETECTOR_CLR
      public DETECTOR_ENDTRACE,DETECTOR_ANGLE

      logical :: dbg=.false.
      integer :: dbgcnt=0
      contains

!-------------------------------------------------------------
      logical function DETECTOR_GO(INST)
! call GO on given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
        DETECTOR_GO=DETECTOR_GO_INST(ADETECTORS(INST)%X,INST)
      end function DETECTOR_GO

!---------------------------------------------------------------------
      SUBROUTINE DETECTOR_INIT(INST)
!----------------------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(DETECTOR),POINTER :: OBJ
      if (.not.DETECTOR_isValid(INST)) return
      OBJ => ADETECTORS(INST)%X
      CALL FRAME_INIT(OBJ%FRAME)
      dbgcnt=0
      END SUBROUTINE DETECTOR_INIT

!-------------------------------------------------------------
      SUBROUTINE DETECTOR_CLR(INST)
! return reference to given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(DETECTOR),pointer :: OBJ
        if (.not.DETECTOR_isValid(INST)) return
        OBJ => ADETECTORS(INST)%X
        OBJ%FRAME%COUNT=0
        call NSTORE_CLEAR(OBJ%FRAME%REGISTRY)
      end SUBROUTINE DETECTOR_CLR

!--------------------------------------------------------------------------
      subroutine DETECTOR_ENDTRACE(INST,CNT,P)
! Add previously registered event to the monitor storage, using given weight
!--------------------------------------------------------------------------
      integer,intent(in) :: INST,CNT
      REAL(kind(1.D0)),intent(in)  :: P
      TYPE(PDETECTOR) :: OBJ
        OBJ%X => ADETECTORS(INST)%X
        if (OBJ%X%FRAME%REGISTRY>0) then
          ! transform to incident axis coordinates
          call FRAME_LOCAL(-1,OBJ%X%FRAME,TMPN(INST))
          TMPN(INST)%P=P     ! set the final weight to the event
          TMPN(INST)%CNT=CNT ! must set counter to neutron before calling NSTORE_SETEVENT !!
          call NSTORE_SETEVENT(OBJ%X%FRAME%REGISTRY,TMPN(INST))
        endif
      end subroutine DETECTOR_ENDTRACE

!-------------------------------------------------
      SUBROUTINE DETECTOR_ADJUST(INST,IERR)
! Update position of the component using nominal trajectory represented by NEUT
! Call FRAME_INIT_MAT first.
! Then use SLIT_PRE & SLIT_POST procedures to propagate NEUT
!-------------------------------------------------
      integer,intent(in) :: INST
      integer,intent(out) :: IERR
      TYPE(DETECTOR),POINTER :: OBJ
      REAL(kind(1.D0)) :: a
1     format(a,': ',7(G13.6,1x))
      if (.not.DETECTOR_isValid(INST)) return
      OBJ => ADETECTORS(INST)%X
      IERR=0
      call FRAME_INIT_MAT(OBJ%FRAME)
      !write(*,1) 'Move from '//trim(OBJ%FRAME%ID),NEUT%K0*NEUT%T
      !write(*,1) 'K in',NEUT%K
      !write(*,1) 'DETECTOR_ADJUST in  ',(NEUT%T-TOF_SAMPLE)*NEUT%K0, NEUT%R(3), NEUT%K(3)/NEUT%K0
      call SLIT_PRE(OBJ%FRAME)
      !write(*,1) 'K pre',NEUT%K
      !write(*,1) 'DETECTOR_ADJUST pre ',(NEUT%T-TOF_SAMPLE)*NEUT%K0, NEUT%R(3), NEUT%K(3)/NEUT%K0
      call SLIT_POST(OBJ%FRAME)
      !write(*,1) 'K post',NEUT%K
      !write(*,1) 'DETECTOR_ADJUST post',(NEUT%T-TOF_SAMPLE)*NEUT%K0
      !write(*,1) 'Move to '//trim(OBJ%FRAME%ID),NEUT%K0*NEUT%T,OBJ%FRAME%DIST
      ! if rad=0, set rad to the distance from the sample
      if (OBJ%RAD<=0.D0) then
        ! ensure non-zero radius !
        OBJ%RAD=max(1.D0,(NEUT%T-TOF_SAMPLE)*NEUT%K0)
      endif
      ! set cos an sin of the clip angles
      OBJ%COSX=1.D0
      OBJ%SINX=0.D0
      OBJ%COSY=1.D0
      OBJ%SINY=0.D0
      select case (OBJ%TYP)
      case(det_type_cyl,det_type_sph)
        a=0.5D0*min(TWOPI,OBJ%FRAME%SIZE(1)/OBJ%RAD)
        OBJ%COSX=cos(a)
        OBJ%SINX=sin(a)
        if (OBJ%TYP==det_type_sph) then
          a=0.5D0*min(TWOPI,OBJ%FRAME%SIZE(2)/OBJ%RAD)
          OBJ%COSY=cos(a)
          OBJ%SINY=sin(a)
        endif
      end select
      ! set KI in local coordinates

      call FRAME_LOCAL(1,OBJ%FRAME,NEUI)
      ! write(*,1) 'DETECTOR_ADJUST KI, K0, TOF=',NEUI%K, NEUI%K0, NEUT%T/HOVM
      OBJ%KI=NEUI%K
      OBJ%KI0=NEUI%K0
      OBJ%TOF=NEUT%T
      OBJ%LAMBDA=TWOPI/NEUT%K0
      !write(*,1) 'DETECTOR_ADJUST RAD, TOF',OBJ%RAD,NEUT%T
      end SUBROUTINE DETECTOR_ADJUST

!--------------------------------------------------
      LOGICAL FUNCTION DETECTOR_ARRAY(DET,T1,T2)
! tracing for curved detectors
!--------------------------------------------------
      TYPE(DETECTOR) :: DET
      REAL(kind(1.D0)),intent(out) :: T1,T2
      REAL(kind(1.D0)) :: tau1,tau2
      REAL(kind(1.D0)) :: P(0:64),KSI(0:64),T(0:64),R(3),Z
      REAL :: RN
      integer :: J

      !RN=RAN1()
      P(0)=0.D0
      KSI(0)=0.D0
      T1=-1.D-31
      T2=1.D31
      R=NEUT%R
      DO J=1,DET%ND
        R(1)=NEUT%R(1)+(DET%FRAME%SIZE(1)+DET%SPACE)*(J-(DET%ND+1.D0)/2.D0)
        IF (FRAME_BORDER(DET%FRAME,R,NEUT%K,tau1,tau2)) THEN
          KSI(J)=tau2-tau1
          P(J)=P(J-1)+(1.D0-EXP(-DET%ALPHA*TWOPI*KSI(J)))
          T(J)=tau1
        ELSE
          P(J)=P(J-1)
          KSI(J)=0.
          T(J)=0.D0
        ENDIF
      ENDDO
      IF (P(DET%ND).LT.1.E-10) then
        DETECTOR_ARRAY=.false.
        return
      endif

  ! play rus. rulette to decide which tube to enter
      RN=RAN1()
      Z=P(DET%ND)*RN
      J=0
      DO WHILE ((P(J).LE.Z).AND.(J.LT.DET%ND))
        J=J+1
      ENDDO
      T1=T(J)
      T2=T(J)+KSI(J)
      NEUT%P=NEUT%P*(1.D0-P(J-1))
      DETECTOR_ARRAY=.true.
      end FUNCTION DETECTOR_ARRAY


!--------------------------------------------------
      SUBROUTINE DETECTOR_ANGLE(DET,R,T,THETA,LAMBDA)
! Calculate scattering angle and wavelength as if measured
! by the detector for given detection position and time, R,T.
! Takes into account various detector geometries.
! R = recorded neutron position
! T = recorded time from TOF_ZERO (not from source !!)
! THETA = scattering angle calculated from detection position
! LAMBDA = wavelength calculated from the detection time
!--------------------------------------------------
      TYPE(DETECTOR) :: DET
      REAL(kind(1.D0)),intent(in)  :: R(3),T
      REAL(kind(1.D0)),intent(out) :: THETA,LAMBDA
      REAL(kind(1.D0)) :: R1(3),T1,tau,cosa,RN,a,sinthb,t0,dt,tof
      integer :: ndt
1     format(a,6(1x,G13.6))
      !dbg=.false.
      !dbg=(dbgcnt<20)
      !if (dbg) dbgcnt=dbgcnt+1
    ! position relative to the rotation center
      R1=R
      R1(3)=R1(3)+DET%RAD
    ! change position to the detector plane by moving along R1
    ! (in reality, we don't know the detection depth)
      select case(DET%TYP)
      case(det_type_cyl)
        !K=(/R1(1),0.D0,R1(3)/)
        call SURFACE_CYL(R1,R1,DET%RAD,tau,T1)
        R1(1)=R1(1)*(1.D0+T1)
        R1(2)=R1(2)*(1.D0+T1)
        R1(3)=R1(3)*(1.D0+T1)
      case(det_type_sph)
        !K=R1
        call SURFACE_SPH(R1,R1,DET%RAD,tau,T1)
        R1(1)=R1(1)*(1.D0+T1)
        R1(2)=R1(2)*(1.D0+T1)
        R1(3)=R1(3)*(1.D0+T1)
      case default
        !K=R1
        ! T1=-(R1(3)-DET%RAD)/R1(3)
        !R1(1)=R1(1)*(1.D0+T1)
        !R1(2)=R1(2)*(1.D0+T1)
        !R1(3)=R1(3)*(1.D0+T1)
        ! consider aberation for flat PSD => keep x,y as registered,
        ! assume det. centre as the detection plane
        R1(3)=DET%RAD
      end select
    ! calculate scattering angle from measured position
      cosa=DET%KI(1)*R1(1)+DET%KI(2)*R1(2)+DET%KI(3)*R1(3)
      RN=SQRT(R1(1)**2+R1(2)**2+R1(3)**2)
      a=acos(cosa/(RN*DET%KI0))
    ! calculate correction time
      dt=(RN-DET%RAD)/DET%KI0
      if (modulation_d0>0.D0) then
        sinthb=sin(a/2.D0)
        t0=2.D0*modulation_d0*sinthb*(DET%TOF-TOF_ZERO)/DET%LAMBDA
        ndt=NINT((T-t0)/modulation_dt/HOVM)
        dt=dt+ndt*modulation_dt*HOVM
      endif
    ! time of flight
      tof=T-dt
    ! measured wavelengh (only for elastic scattering !!!)
      lambda=tof/(DET%TOF-TOF_ZERO)*DET%LAMBDA
    ! scattering angle
      theta=a
      end subroutine DETECTOR_ANGLE


!--------------------------------------------------
      LOGICAL FUNCTION DETECTOR_PSD(DET,T1,T2)
! in local coordinates: get entry and exit times
! return true if neutron passes through the detection volume
!--------------------------------------------------
      TYPE(DETECTOR) :: DET
      REAL(kind(1.D0)),intent(out) :: T1,T2
! add resolution effect
      if (DET%RESX>0.D0) NEUT%R(1)=NEUT%R(1)+GASDEV()*DET%RESX/R8LN2
      if (DET%RESY>0.D0) NEUT%R(2)=NEUT%R(2)+GASDEV()*DET%RESY/R8LN2
      DETECTOR_PSD=FRAME_BORDER(DET%FRAME,NEUT%R,NEUT%K,T1,T2)
      end FUNCTION DETECTOR_PSD

!------------------------------------------------------
      LOGICAL FUNCTION DETECTOR_CYL(DET,tin,tout)
! tracing for cylindrically curved detectors
! return entry and exit times for the detection volume
! in the begining, NEUT is in local coordinates close to the volume center
!------------------------------------------------------
      TYPE(DETECTOR) :: DET
      REAL(kind(1.D0)),intent(out) :: tin,tout
      LOGICAL :: LOG1
      REAL(kind(1.D0)) :: T1,T2,T3,T4,T5,T6,tau
      REAL(kind(1.D0)) :: R1(3),DX,DY

      REAL(kind(1.D0)) :: RN,a,cosa,R(3)
      LOG1=.FALSE.
! add resolution effect
      if (abs(DET%RESX)+abs(DET%RESY)>0.D0) then
        DX=GASDEV()*DET%RESX/R8LN2
        DY=GASDEV()*DET%RESY/R8LN2
        if (DET%RAD>1.D-3) then
          NEUT%R(1)=NEUT%R(1)+DX*(1.D0+NEUT%R(3)/DET%RAD)
          NEUT%R(3)=NEUT%R(3)+DX*NEUT%R(1)/DET%RAD
        else
          NEUT%R(1)=NEUT%R(1)+DX
        endif
        NEUT%R(2)=NEUT%R(2)+DY
      endif
    ! position relative to the rotation center
      R1=NEUT%R
      R1(3)=R1(3)+DET%RAD
1     format(a,6(1x,G13.6))
      if (NEUT_DBG) then
        write(*,1) 'DETECTOR_CYL R1=',R1
      endif
! get the length of flight path through the detector
      call SURFACE_CYL(R1,NEUT%K,DET%RAD-0.5D0*DET%FRAME%SIZE(3),tau,T1)
      call SURFACE_CYL(R1,NEUT%K,DET%RAD+0.5D0*DET%FRAME%SIZE(3),tau,T2)
      call BORDER_WALL(NEUT%R(2),NEUT%K(2),DET%FRAME%SIZE(2),T3,T4)
      call SURFACE_SECTOR((/R1(1),R1(3)/),(/NEUT%K(1),NEUT%K(3)/),DET%COSX,DET%SINX,T5,T6)
      tin=max(T1,T3,T5)
      tout=min(T2,T4,T6)
      if (NEUT_DBG) then
        write(*,1) 'DETECTOR_CYL  T1..T6=',T1,T2,T3,T4,T5,T6
        write(*,1) 'DETECTOR_CYL  tin,tout=',tin,tout
      endif
      LOG1=(tout>tin)
      ! clip angular range
      if (LOG1) then
        R=NEUT%R-(/0.D0,0.D0,-DET%RAD/)
        cosa=DET%KI(1)*R(1)+DET%KI(2)*R(2)+DET%KI(3)*R(3)
        RN=SQRT(R(1)**2+R(2)**2+R(3)**2)
        a=acos(cosa/(RN*DET%KI0))
        if ((abs(a)<DET%THMIN).or.(abs(a)>DET%THMAX)) then
          LOG1=.false.
        endif
        if (NEUT_DBG) then
          write(*,1) 'DETECTOR_CYL a=',abs(a)*180/PI
          write(*,1) 'DETECTOR_CYL  THMIN, THMAX=',DET%THMIN*180/PI,DET%THMAX*180/PI
          write(*,1) 'DETECTOR_CYL  DET%KI, R=',DET%KI, R
        endif
      endif
      DETECTOR_CYL=LOG1

      !lambda=(NEUT%T-TOF_ZERO)/(DET%TOF-TOF_ZERO)*DET%LAMBDA
      !XVAL=0.5D0*lambda/sin(a/2.D0)

      !write(*,1) 'DETECTOR_CYL L=',(DET%TOF-TOF_ZERO)*DET%KI0,(NEUT%T-TOF_ZERO)*NEUT%K0
      !write(*,1) 'DETECTOR_CYL theta=',TWOPI/NEUT%K0,lambda,XVAL,a/deg

      END FUNCTION DETECTOR_CYL

!------------------------------------------------------
      LOGICAL FUNCTION DETECTOR_SPH(DET,tin,tout)
! tracing for spherically curved detectors
! return entry and exit times for the detection volume
! in the begining, NEUT is in local coordinates close to the volume center
!------------------------------------------------------
      TYPE(DETECTOR) :: DET
      REAL(kind(1.D0)),intent(out) :: tin,tout
      LOGICAL :: LOG1
      REAL(kind(1.D0)) :: T1,T2,T3,T4,T5,T6,tau
      REAL(kind(1.D0)) :: R1(3),DX,DY,sina,sint,cosa,cost
      LOG1=.FALSE.
! add resolution effect
      if (abs(DET%RESX)+abs(DET%RESY)>0.D0) then
        DX=GASDEV()*DET%RESX/R8LN2
        DY=GASDEV()*DET%RESY/R8LN2
        if (DET%RAD>1.D-3) then
          sint=NEUT%R(2)/DET%RAD
          cost=sqrt(NEUT%R(1)**2+NEUT%R(3)**2)
          sina=NEUT%R(1)/cost
          cosa=(DET%RAD+NEUT%R(3))/cost
          NEUT%R(1)=NEUT%R(1)+(DX*cosa-DY*sint*sina)
          NEUT%R(3)=NEUT%R(3)-(DX*sina+DY*sint*cosa)
          NEUT%R(2)=NEUT%R(2)+DY*cost/DET%RAD
        else
          NEUT%R(1)=NEUT%R(1)+DX
          NEUT%R(2)=NEUT%R(2)+DY
        endif
      endif
    ! position relative to the rotation center
      R1=NEUT%R
      R1(3)=R1(3)+DET%RAD
1     format(a,6(1x,G13.5))
      if (dbg) then
        write(*,1) 'DETECTOR_SPH R1=',R1
      endif
! get the length of flight path through the detector
      call SURFACE_SPH(R1,NEUT%K,DET%RAD-0.5D0*DET%FRAME%SIZE(3),tau,T1)
      call SURFACE_SPH(R1,NEUT%K,DET%RAD+0.5D0*DET%FRAME%SIZE(3),tau,T2)
      call SURFACE_SECTOR((/R1(1),R1(3)/),(/NEUT%K(1),NEUT%K(3)/),DET%COSX,DET%SINX,T3,T4)
      call SURFACE_SECTOR((/R1(2),R1(3)/),(/NEUT%K(2),NEUT%K(3)/),DET%COSY,DET%SINY,T5,T6)
      tin=max(T1,T3,T5)
      tout=min(T2,T4,T6)
      if (dbg) then
        write(*,1) 'DETECTOR_SPH T1..T6=',T1,T2,T3,T4,T5,T6
        write(*,1) 'DETECTOR_SPH tin,tout=',tin,tout
      endif
      LOG1=(tout>tin)
      DETECTOR_SPH=LOG1
      END FUNCTION DETECTOR_SPH


!--------------------------------------------------
      LOGICAL FUNCTION DETECTOR_GO_INST(DET,INST)
! tracing for curved detectors
!--------------------------------------------------
      TYPE(DETECTOR) :: DET
      integer,intent(in) :: INST
      LOGICAL :: LOG1
      REAL(kind(1.D0)) :: T1,T2,tau,P0,RN
      dbg=.false.
      !dbg=(dbgcnt<20)
      !dbgcnt=dbgcnt+1

      LOG1=.FALSE.
! move to detector centre and convert to local coordinates
      if (dbg) write(*,1) 'DETECTOR_GO_INST PRE',NEUT%R,NEUT%K
      call DET_PRE(DET)
      if (dbg) write(*,1) 'DETECTOR_GO_INST LOC',NEUT%R,NEUT%K
      call AddEventLog(DET%FRAME,NEUT)
      select case(DET%TYP)
      case(det_type_cyl)
        LOG1=DETECTOR_CYL(DET,T1,T2)
      case(det_type_sph)
        LOG1=DETECTOR_SPH(DET,T1,T2)
      case(det_type_array)
        LOG1=DETECTOR_ARRAY(DET,T1,T2)
      case(det_type_psd)
        LOG1=DETECTOR_PSD(DET,T1,T2)
        if (dbg) then
          write(*,1) 'DETECTOR_GO_INST T1,T2,lambda',LOG1,T1,T2,TWOPI/NEUT%K0
          write(*,1) 'DETECTOR_GO_INST PSD',NEUT%R,NEUT%K
        endif
      case default
        LOG1=TLIB_INSIDE(NEUT%R(1:2),DET%FRAME%SIZE(1:2),DET%FRAME%SHAPE)
        T1=0.D0
        T2=0.D0
      end select
1     format(a,6(1x,G13.5))
      !if (dbg) then
      !  write(*,1) 'DETECTOR_GO_INST LOG1=',LOG1
      !endif

      if (.not.LOG1) goto 99
! sample random peneration depth
      if (DET%TYP.ne.det_type_monit) then
        RN=RAN1()
        P0=1.D0-EXP(-DET%ALPHA*twopi*(T2-T1))
        tau=-log(1.D0-P0*RN)/(DET%ALPHA*twopi)
        if (dbg) then
          write(*,1) 'DETECTOR_GO_INST P0,RN,tau,T1+tau=',P0,RN,tau,T1+tau
        endif
        if (TRACING_UP) then
          call TransportG(T2-tau)
        else
          call TransportG(T1+tau)
        endif
        NEUT%P=NEUT%P*P0
      endif
      dbg=.false.
      call AddEventLog(DET%FRAME,NEUT)
! register event
! IMPORTANT: store NEUT in local coordinates, i.e. before the call to DET_POST
      TMPN(INST)=NEUT
      if (TRACING_UP) TMPN(INST)%K=-TMPN(INST)%K
      call incCounter(DET%FRAME%COUNT)
      call DET_POST(DET)
      DETECTOR_GO_INST=.TRUE.
      RETURN

99    NEUT%P=0.D0
      DETECTOR_GO_INST=.FALSE.
      END FUNCTION DETECTOR_GO_INST

!------------------------------------------------------------
      SUBROUTINE DET_PRE(DET)
! move to the center of detection volume and convert NEUT
! to the local system.
!------------------------------------------------------------
      TYPE(DETECTOR) :: DET
      REAL(kind(1.D0)) :: tin,dt,R(3)
  ! up-stream: use default SLIT procedure
      if (TRACING_UP) then
        call SLIT_PRE(DET%FRAME)
      else
        select case (DET%TYP)
        case(det_type_cyl,det_type_sph)
        ! convert to local coordinates
          NEUT%R(3)=NEUT%R(3)-DET%FRAME%DIST
          call FRAME_LOCAL(1,DET%FRAME,NEUT)
          call GRAVITY_LOCAL(1,DET%FRAME)
        ! get R=position rel. to rotation center
          R=NEUT%R
          R(3)=NEUT%R(3)+DET%RAD
        ! get time to sphere or cylinder surface
          if (DET%TYP==det_type_cyl) then
            call SURFACE_CYL(R,NEUT%K,DET%RAD,tin,dt)
          else
            call SURFACE_SPH(R,NEUT%K,DET%RAD,tin,dt)
          endif
1         format(a,6(1x,G13.5))
          if (dbg) then
            write(*,*)
            write(*,1) 'DET_PRE R,K=',R,NEUT%K
            write(*,1) 'DET_PRE tin,dt=',tin,dt
          endif
        ! transport to the middle of the detection layer (sphere of cylinder)
          call TransportG(DT)
          if (dbg) then
            write(*,1) 'DET_PRE LOC=',NEUT%R
          endif
          !call FRAME_LOCAL(1,DET%FRAME,NEUT)
        case default
          call SLIT_PRE(DET%FRAME)
        end select
      endif
      END SUBROUTINE DET_PRE

!------------------------------------------------------------
      SUBROUTINE DET_POST(DET)
! transform to entry coordinates for detector (up-stream tracing only)
!------------------------------------------------------------
      TYPE(DETECTOR) :: DET
      call SLIT_POST(DET%FRAME)
      END SUBROUTINE DET_POST


      end MODULE DETECTORS_TRACE
