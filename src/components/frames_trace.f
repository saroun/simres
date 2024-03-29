!//////////////////////////////////////////////////////////////////////
!////  $Id: frames_trace.f,v 1.25 2019/06/23 14:51:45 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2011, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.25 $
!////     $Date: 2019/06/23 14:51:45 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Tracing procedures for FRAMES
!////
!////////////////////////////////////////////////////////////////////////
      MODULE FRAMES_TRACE
      use CLASSES
      use SIMATH
      use TRACELIB
      USE TRACINGDATA
      USE FRAMES
      use EVENTLOGS
      implicit none


      contains

!-------------------------------------------------------------
      logical function FRAME_GO(INST)
! call GO on given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
        FRAME_GO=FRAME_GO_INST(AFRAMES(INST)%X)
      end function FRAME_GO

!------------------------------------------------------------
      LOGICAL FUNCTION FRAME_GO_INST(OBJ)
! transport through a slit = entry OBJ window
!------------------------------------------------------------
      TYPE(TFRAME) :: OBJ
      LOGICAL :: LOG1
      LOGICAL :: DBG
      TYPE(NEUTRON) :: NGLOB
      if (OBJ%TRANSPARENT.EQ.1) then
        FRAME_GO_INST=NONE_GO(OBJ)
        return
      endif
      LOG1=.FALSE.
! move to the OBJ entry surface and transform to local coordinates
      CALL SLIT_PRE(OBJ)
! move neutron to z=0
      call SLIT_MIDDLE(.false.)

! log event at the entry of each object
      call AddEventLog(OBJ,NEUT)

2     format(a,2x,10(1x,g10.4))
      !DBG=(OBJ%COUNT<10)
      DBG=.false.
      if (DBG) then
         call Loc2Glob(OBJ,NEUT,NGLOB)
        WRITE(*,*)
        write(*,2) trim(OBJ%name)//' NEUT.R = ',NEUT%R(1:3)
        write(*,2) trim(OBJ%name)//' TLAB = ',OBJ%TLAB(1:3)
        write(*,2) trim(OBJ%name)//' RLAB = ',OBJ%RLAB(1:3,1)
        write(*,2) trim(OBJ%name)//' RLAB = ',OBJ%RLAB(1:3,2)
        write(*,2) trim(OBJ%name)//' RLAB = ',OBJ%RLAB(1:3,3)
        write(*,2) trim(OBJ%name)//' NGLOB.R = ',NGLOB%R(1:3)
        write(*,2) trim(OBJ%name)//' NGLOB.K = ',NGLOB%K(1:3)
      endif

! check that neutron fits within the object
      LOG1=TLIB_INFRAME(NEUT%R,OBJ%SIZE,OBJ%SHAPE)
      IF (LOG1) THEN
! transform to exit axis reference frame
        CALL SLIT_POST(OBJ)
! increment counter
        call incCounter(OBJ%COUNT)
      ELSE
        NEUT%P=0
      ENDIF
      FRAME_GO_INST=LOG1
      END FUNCTION FRAME_GO_INST


!------------------------------------------------------------
      SUBROUTINE FRAME_TRACE0(OBJ)
! Tracing of nominal trajectory (beam axis)
! to be used by guide-like components for alignment by *_ADJUST procedures.
! Starts and ends with NEUT in local coordinates.
! Assume forward tracing.
!------------------------------------------------------------
      TYPE(TFRAME),intent(out) :: OBJ
      REAL(KIND(1.D0)) :: R2(3), K2(3), DX, DY, DT
!1     format('TRACE0, ',a,': ',7(G13.6,1x))
      ! tracing beam axis: we ignore any local stage or gonio movement
      ! revert back to incident axis coord.
      call FRAME_LOCAL(-1,OBJ,NEUT)
      ! transport to the guide exit
      DT = (OBJ%SIZE(3)-NEUT%R(3))/NEUT%K(3)
      if (OBJ%ISDEFL) then
        NEUT%R = OBJ%TDEFL
        K2 = NEUT%K
        !write(*,1) 'deflection in ',NEUT%R,NEUT%K,OBJ%RDEFL(1,3)
        call M3xV3(-1, OBJ%RDEFL, K2, NEUT%K)
        !write(*,1) 'deflection out',NEUT%R,NEUT%K,OBJ%RDEFL(1,3)
        NEUT%T = NEUT%T + DT
      else
        call TransportG(DT)
      endif
!      write(*,1) 'DT, R',DT, NEUT%R
      call FRAME_LOCAL(1,OBJ,NEUT)
      !write(*,1) ' local',NEUT%R, NEUT%K
! JS> OK, it's not very eficient since SLIT_POST does inverse operation just after,
! but this is not a CPU critical code. It is only used for adjustment before simulation.
! For code consistency, I'd prefere to stay in the same reference frame.
      END SUBROUTINE FRAME_TRACE0

!-------------------------------------------------
      SUBROUTINE FRAME_ADJUST(INST)
! Update position of the component using nominal trajectory represented by NEUT
! Call FRAME_INIT_MAT first.
! Then use SLIT_PRE & FRAME_TRACE0 & SLIT_POST procedures to propagate NEUT
!-------------------------------------------------
      integer,intent(in) :: INST
      TYPE(TFRAME),POINTER :: OBJ
      if (.not.FRAME_isValid(INST)) return
      OBJ => AFRAMES(INST)%X
      call SLIT_PRE(OBJ)
      call FRAME_TRACE0(OBJ)
      call SLIT_POST(OBJ)
      end SUBROUTINE FRAME_ADJUST


!---------------------------------------------------------
      integer FUNCTION CONTAINER_BORDER(OBJ,r,k,t,dt)
! Get neutron cross-section times with a container.
! The container can be a hollow object.
! Arguments:
! ----------
! r: current position in local coordinates
! k: wave vector
! Returns:
! --------
! 0: now intersection
! 1: one layer, then t(1), t(2) denote input and output times
! 2: two layers, then in/out times are in t(1),t(3) / t(2),t(4)
! dt then contains sum of thicknesses passed along k before and after the current position
! Assume r,k in local coordinates
! Time is in units [m/h_bar] everywhere   i.e. length=time*K
!---------------------------------------------------------
      TYPE(TFRAME) :: OBJ
      REAL(kind(1.D0)),intent(in)  :: r(3),k(3)
      REAL(kind(1.D0)),intent(out)  :: t(4), dt(2)
      integer :: ires
      logical :: L1, L2, valid, db
      REAL(kind(1.D0)) :: tin, tout, tin2, tout2, dt1, dt2
1     format(a,': ',6(G10.4,1x))
      ires = 0
      L1 = TLIB_BORDER(OBJ%SHAPE,OBJ%SIZE,r,k,tin,tout)
      t = (/tin, tout,0.D0,0.D0/)
      dt = (/-tin, tout/)

      db=.false.
      !db=((r(3)>0).and.(OBJ%COUNT<30))
      if (db) then
        write(*,*) 'CONTAINER_BORDER',L1,OBJ%ISHOLLOW
        write(*,1) '  dt',dt
        write(*,1) '  r,k',r,k
      endif
      ! there is intersection with envelope
      if (L1) then
        ! full objects
        if (.not.OBJ%ISHOLLOW) then
          ires=1
        ! hollow objects
        else
          L2 = TLIB_BORDER(OBJ%SHAPE,OBJ%INSIZE,r,k,tin2,tout2)
          if (L2) then
            dt1 = tin2-tin
            dt2 = tout-tout2
            if ((dt1>0.D0).and.(dt2>0.D0)) then
              ires = 2
              t = (/tin, tin2, tout2, tout/)
              if (t(2)>0.D0) then
                dt(1) = -t(1)
                dt(2) = t(2) + t(4) - t(3)
              else
                dt(1) = t(2) - t(1) - t(3)
                dt(2) = t(4)
              endif
            else if (dt1>0.D0) then
              ires = 1
              t = (/tin, tin2, 0.D0, 0.D0/)
              dt = (/-tin, tin2/)
            else if (dt2>0.D0) then
              ires = 1
              t = (/tout2, tout, 0.D0, 0.D0/)
              dt = (/-tout2, tout/)
            endif
            if (db) then
              write(*,1) '  ires,t',ires,t
              write(*,1) '  dt1,dt2',dt1,dt2
            endif
          else
            ires = 1
          endif
        endif
      endif
      if (ires==1) then
        valid = ((t(1).lt.0.D0).and.(t(2).gt.0.D0))
      else if (ires==2) then
        valid = ((t(1).lt.0.D0).and.(t(4).gt.0.D0))
        if (t(2).ge.0.D0) then
          valid = (valid.and.(t(3).gt.0.D0))
        else
          valid = (valid.and.(t(3).le.0.D0))
        endif
      else
        valid = .true.
      endif
      ! debug condition, this should never happen
      if (.not.valid) then
          write(*,1) 'CONTAINER_BORDER unexpected position, t=',t
          ires=0
      endif
      if (db) then
        write(*,1) '  ires,dt',ires, dt
        write(*,1) '  t',t
      endif
      CONTAINER_BORDER = ires
      END FUNCTION CONTAINER_BORDER

!---------------------------------------------------------
      LOGICAL FUNCTION FRAME_BORDER(OBJ,r,k,tin,tout)
! test neutron cross-section with a frame
! return cross-section times, if any
! assume r,k in local coordinates
! Time is in units [m/h_bar] everywhere   i.e. length=time*K
!---------------------------------------------------------
      TYPE(TFRAME) :: OBJ
      REAL(kind(1.D0)),intent(in)  :: r(3),k(3)
      REAL(kind(1.D0)),intent(out)  :: tin,tout
      FRAME_BORDER = TLIB_BORDER(OBJ%SHAPE,OBJ%SIZE,r,k,tin,tout)
      END FUNCTION FRAME_BORDER


!------------------------------------------------------------
      SUBROUTINE FRAME_AXIS(IT,OBJ,NEU)
! transform neutron coordinates between incident axis and exit axis reference frames
! IT=-1 ... exit -> incident
! else    ... incident -> exit
!------------------------------------------------------------
      TYPE(NEUTRON) :: NEU
      TYPE(TFRAME) :: OBJ
      INTEGER,intent(in) :: IT
      REAL(KIND(1.D0)) :: V(3),VK(3)
!1     format('FRAME_AXIS, ',a,': ',7(G13.6,1x))

      IF (OBJ%MEXI.EQ.0) RETURN
      !write(*,1) trim(OBJ%ID), IT, OBJ%MEXI,OBJ%REXI(1,3)
      select case(IT)
      case(-1)
        CALL M3XV3(-OBJ%MEXI,OBJ%REXI,NEU%R,V)
        NEU%R = V + OBJ%TEXI
        VK = NEU%K
        CALL M3XV3(-OBJ%MEXI,OBJ%REXI,VK, NEU%K)
      case default
!        write(*,1) '  in R, K',NEU%R, NEU%K
        V = NEU%R - OBJ%TEXI
        CALL M3XV3(OBJ%MEXI,OBJ%REXI,V,NEU%R)
        VK = NEU%K
        CALL M3XV3(OBJ%MEXI,OBJ%REXI,VK, NEU%K)
        !write(*,1) ' out VK, K', VK, NEU%K
      end select
      END SUBROUTINE FRAME_AXIS

!------------------------------------------------------------
      SUBROUTINE SLIT_TURN(OBJ)
! turnes neutron by OBJ%AX
! works with NEUT in local coordinates
!------------------------------------------------------------
      TYPE(TFRAME) :: OBJ
      REAL(KIND(1.D0)) :: VKI(3), VKF(3)
      if (OBJ%MEXI.NE.0) then
        ! convert from local to incident axis coord.
        CALL M3XV3(-OBJ%MLOC,OBJ%RLOC,NEUT%K,VKI)
        ! turn neutron by AX
        CALL M3XV3(-OBJ%MEXI,OBJ%REXI,VKI,VKF)
        ! convert back to local coors.
        CALL M3XV3(OBJ%MLOC,OBJ%RLOC,VKF,NEUT%K)
      endif
      END SUBROUTINE SLIT_TURN

!------------------------------------------------------------
      SUBROUTINE SLIT_ENTRY(OBJ)
! move to the object entry surface = the closest cross-section
! in the actual direction of neutron flight
! assume local reference frame
!------------------------------------------------------------
      TYPE(TFRAME) :: OBJ
      REAL(KIND(1.D0)) :: tin,tout
      if (FRAME_BORDER(OBJ,NEUT%R,NEUT%K,tin,tout)) then
  ! move slightly behind the entry, to be really INSIDE the object
  ! this is required by some components (CRYSTAL)
        call TransportG(tin+1.D-7*(tout-tin))
      else
        NEUT%P=0.D0 ! no cross-section with the object
      endif
      END SUBROUTINE SLIT_ENTRY

!---------------------------------------------------------------------------
      SUBROUTINE FRAME_ENTRY(OBJ)
! move to the entry point
! unlike SLIT_ENTRY, this subroutine moves events to the z=0 or z=LENGTH plane
! depending on the tracing direction
!-----------------------------------------------------------------------------
      TYPE(TFRAME) :: OBJ
      REAL(KIND(1.D0)) :: Z0,DT
! trace up-stream: move to z=LENGTH
      if (TRACING_UP) then
        Z0=OBJ%SIZE(3)
! trace down-stream: move to z=0
      else
        Z0=0.D0
      endif
      if (abs(NEUT%K(3)).gt.1.D-20) then
        DT=(Z0-NEUT%R(3))/NEUT%K(3)
      else
        DT=1.D30
        NEUT%P=0.D0
      endif
      call TransportG(DT)
      END subroutine FRAME_ENTRY

!------------------------------------------------------------
      SUBROUTINE SLIT_PRE(OBJ)
! in up-stream mode:
!    rotate position, momentum and gravity to incident coordinates
!    invert momentum
! in down stream mode:
!    shift coordinates by DIST and transport neutron to the component center
!    (usually z=0)
!------------------------------------------------------------
      TYPE(TFRAME) :: OBJ
      INTEGER :: I
!1     format('SLIT_PRE, ',a,': ',7(G13.6,1x))
! transform from previous object ref. frame
! DT is time to R(3)=0 or R(3)=SIZE(3) for down or up stream tracing, resp.
  ! up-stream: convert to incident axis system and invert K
      if (TRACING_UP) then
        call FRAME_AXIS(-1,OBJ,NEUT)
        call GRAVITY_AXIS(-1,OBJ)
        do i=1,3
          NEUT%K(i)=-NEUT%K(i)
        enddo
  ! down-stream: shift origin by OBJ%DIST along z
      else
    ! convert to the incident axis system = shift by +OBJ%DIST
        !write(*,1) trim(OBJ%ID),NEUT%R, NEUT%K
        NEUT%R(3)=NEUT%R(3)-OBJ%DIST
    ! move to R(3)=0
        call SLIT_MIDDLE(.false.)
      endif
! transform to local coordinates
      call FRAME_LOCAL(1,OBJ,NEUT)
      call GRAVITY_LOCAL(1,OBJ)
      !write(*,1) '    local',NEUT%R, NEUT%K
      END SUBROUTINE SLIT_PRE

!------------------------------------------------------------
      SUBROUTINE SLIT_POST(OBJ)
! in down-stream mode:
!    rotate position, momentum and gravity to exit coordinates
!    invert momentum
! in up-stream mode:
!    shift coordinates by -DIST and transport neutron to the component center
!    (usually z=0)
!------------------------------------------------------------
      TYPE(TFRAME) :: OBJ
      integer :: i
!1     format('SLIT_POST, ',a,': ',7(G13.6,1x))
! transform from local coordinates to incident axis reference frame
!      write(*,1) trim(OBJ%ID),NEUT%R, NEUT%K
      call FRAME_LOCAL(-1,OBJ,NEUT)
      call GRAVITY_LOCAL(-1,OBJ)
      !write(*,1) '   incid',NEUT%R, NEUT%K
! up-stream: shift origin by OBJ%DIST along z and invert K
      if (TRACING_UP) then
    ! convert to the exit axis system of the preceding component = shift by -OBJ%DIST
        NEUT%R(3)=NEUT%R(3)+OBJ%DIST
    ! and move back so that R(3)=0
        call SLIT_MIDDLE(.false.)
        do i=1,3
          NEUT%K(i)=-NEUT%K(i)
        enddo
! down-stream: transform from incident to exit axis reference frames
      else
        call FRAME_AXIS(1,OBJ,NEUT)
        call GRAVITY_AXIS(1,OBJ)
!        write(*,1) '   exit',NEUT%R, NEUT%K
      endif
      END SUBROUTINE SLIT_POST

!------------------------------------------------------------
      SUBROUTINE SLIT_MIDDLE(FORW)
! move to the object middle plane (R(3)=0)
! assume local reference frame
! if K(3)~0, then move to R(1)=0 or R(2)=0
! if FORW, allow only DT>0 steps
!------------------------------------------------------------
      logical :: FORW
      INTEGER :: iref
      REAL(KIND(1.D0)) :: DT
      REAL(KIND(1.D0)),parameter :: TOL=1.D-9
!1     format('SLIT_MIDDLE, ',a,': ',7(G12.5,2x))
      !write(*,1) 'before', NEUT%R,NEUT%K,NEUT%K0
      if (abs(NEUT%R(3)).gt.TOL) then
        if (abs(NEUT%K(3)).gt.1.D-6*NEUT%K0) then
          iref=3
        else if (abs(NEUT%K(1)).gt.1.D-6*NEUT%K0) then
          iref=1
        else
          iref=2
        endif
        DT=-NEUT%R(iref)/NEUT%K(iref)
      !write(*,1) 'iref, R,K,DT',iref,NEUT%R(iref),NEUT%K(iref),DT
        if (FORW) STOPTIME=(DT.LT.0.D0)
        call TransportG(DT)
!      write(*,*)
      !write(*,1) 'post',NEUT%R,NEUT%K
      else
        NEUT%R(3)=0.D0
      endif
      END subroutine SLIT_MIDDLE

!------------------------------------------------------------
      SUBROUTINE SLIT_EXIT(OBJ)
! move to the object exit surface
! assume local reference frame
!------------------------------------------------------------
      TYPE(TFRAME) :: OBJ
      REAL(KIND(1.D0)) :: tin,tout
      if (FRAME_BORDER(OBJ,NEUT%R,NEUT%K,tin,tout)) then
        call TransportG(tout)
      else
        NEUT%P=0.D0 ! no cross-section with the object
      endif
      END SUBROUTINE SLIT_EXIT

!------------------------------------------------------------
      LOGICAL FUNCTION CONTAINER_GO(OBJ)
! Transport inside a 3D container
! Leaves neutron in local coordinates !
!------------------------------------------------------------
      TYPE(TFRAME) :: OBJ
      LOGICAL :: L1, L2
      real(kind(1.D0)) :: DT
1     format(a12,': ',6(G10.4,1x))
      L1=.FALSE.
  ! move to the OBJ entry surface and transform to local coordinates
      !write(*,1) 'CONTAINER_GO PRE',NEUT%R(1),NEUT%R(3),NEUT%K(1),NEUT%K(3),NEUT%T
      !write(*,1) 'CONTAINER_GO DT ',XRND(IXRND(OBJ%IRNDP))
      CALL SLIT_PRE(OBJ)
  ! random flight path through the object
      if (OBJ%IRNDP>0) then
        DT=XRND(IXRND(OBJ%IRNDP))/NEUT%K0
        call TransportG(DT)
      endif
  ! log event at the entry of each object
      call AddEventLog(OBJ,NEUT)

  ! check that neutron fits within the object
      L1=TLIB_INFRAME(NEUT%R,OBJ%SIZE,OBJ%SHAPE)
! update for hollow objects:
      if(OBJ%ISHOLLOW) then
        L2=TLIB_INFRAME(NEUT%R,OBJ%INSIZE,OBJ%SHAPE)
        L1 = (L1.and.(.not.L2))
      endif
      !write(*,1) 'CONTAINER_GO SCA',NEUT%R(1),NEUT%R(3),NEUT%K(1),NEUT%K(3),NEUT%T,LOG1
  ! don't increment counter !!
        ! call incCounter(OBJ%COUNT)
      IF (.not.L1) NEUT%P=0
      CONTAINER_GO=L1
      END FUNCTION CONTAINER_GO


!------------------------------------------------------------
      LOGICAL FUNCTION NONE_GO(OBJ)
! transport through a transparent object
!------------------------------------------------------------
      TYPE(TFRAME) :: OBJ
      real(kind(1.D0)) :: DT
      if (TRACING_UP) then
        NEUT%K=-NEUT%K
        NEUT%R(3)=NEUT%R(3)+OBJ%DIST
      else
        NEUT%R(3)=NEUT%R(3)-OBJ%DIST
      endif
      DT=-NEUT%R(3)/NEUT%K(3)
      call TransportG(DT)
      if (TRACING_UP) NEUT%K=-NEUT%K
      call incCounter(OBJ%COUNT)
      NONE_GO=.TRUE.
! log event at the entry of each object
      call AddEventLog(OBJ,NEUT)
      END FUNCTION NONE_GO



!------------------------------------------------------------
      SUBROUTINE GRAVITY_LOCAL(IT,OBJ)
! transform gravity coordinates between incident axis and local reference frames
! IT=-1 ... local -> axis
! else    ... axis -> local
! MUST remain in this module due to dependences
!------------------------------------------------------------
      TYPE(TFRAME) :: OBJ
      INTEGER,intent(in) :: IT
      REAL(KIND(1.D0)) :: GR(3)
      if (.not.use_gravity) return
      GR=NEUT_GR
      CALL M3XV3(IT*OBJ%MLOC,OBJ%RLOC,GR,NEUT_GR)
      END SUBROUTINE GRAVITY_LOCAL


!------------------------------------------------------------
      SUBROUTINE GRAVITY_AXIS(IT,OBJ)
! transform gravity coordinates  between incident axis and exit axis reference frames
! IT=-1 ... exit -> incident
! else    ... incident -> exit
!------------------------------------------------------------
      TYPE(TFRAME) :: OBJ
      INTEGER,intent(in) :: IT
      REAL(KIND(1.D0)) :: GR(3)
      if (.not.use_gravity) return
      IF (OBJ%MEXI.EQ.0) RETURN
      GR=NEUT_GR
      CALL M3XV3(IT*OBJ%MEXI,OBJ%REXI,GR,NEUT_GR)
      END SUBROUTINE GRAVITY_AXIS

!------------------------------------------------------------
      REAL(KIND(1.D0)) FUNCTION CONTAINER_DEPTH(OBJ, R)
! Calculate depth under container surface.
! Depth depends on the container symmetry. It should
! measure distance to the nearest surface (or front z-surface for a box).
!------------------------------------------------------------
      TYPE(TFRAME) :: OBJ
      real(kind(1.D0)), intent(in) :: R(3)
      REAL(kind(1.D0)) :: x1, x2, x3, r2, ksi, sz(3)
      REAL(kind(1.D0)) :: depth
      sz = OBJ%SIZE(1:3)
      select case(OBJ%SHAPE)
        case(FRAME_SHAPE_BOX)
          depth = 0.5*sz(3)-R(3)
        case(FRAME_SHAPE_DISC)
          x1 = R(1)*2.D0/sz(1)
          x2 = R(2)*2.D0/sz(2)
          r2 = R(1)**2 + R(2)**2
          ksi = sqrt(x1**2 + x2**2)
          if (ksi.eq.0.D0) then
            depth = 0.5*min(sz(1), sz(2))
          else
            depth = sqrt(r2)*(1.D0/ksi - 1.D0)
          endif
        case(FRAME_SHAPE_CYLINDER)
          x1 = R(1)*2.D0/sz(1)
          x2 = R(3)*2.D0/sz(3)
          r2 = R(1)**2 + R(3)**2
          ksi = sqrt(x1**2 + x2**2)
          if (ksi.eq.0.D0) then
            depth = 0.5*min(sz(1), sz(3))
          else
            depth = sqrt(r2)*(1.D0/ksi - 1.D0)
          endif
        case default
          x1 = R(1)*2.D0/sz(1)
          x2 = R(2)*2.D0/sz(2)
          x3 = R(3)*2.D0/sz(3)
          r2 = R(1)**2 + R(2)**2 +R(3)**2
          ksi = sqrt(x1**2 + x2**2)
          if (ksi.eq.0.D0) then
            depth = 0.5*min(sz(1), sz(2), sz(3))
          else
            depth = sqrt(r2)*(1.D0/ksi - 1.D0)
          endif
      end select
      CONTAINER_DEPTH = depth
      end FUNCTION CONTAINER_DEPTH


!------------------------------------------------------------
      subroutine CONTAINER_COORD(OBJ, R, depth, TM)
! Calculate depth under container surface.
! Return also corresponding R-coordinate system as a column matrix.
! The R-coordinates depend on the container symmetry:
! BOX: depth is  measured from front z-surface inwards, TM=identity.
! CYLINDER, DISC: x // hoop direction, z // radial direction
! ELLIPSOID: z // R, y // (R x x_loc) or (z_loc x R)
! for alll round shapes, depth is the radial distance to the surface
! NOTE: depth may not be the shortes distance to the surface! 
! On input, R is in local coordinates.
!------------------------------------------------------------
      TYPE(TFRAME) :: OBJ
      real(kind(1.D0)), intent(in) :: R(3)
      real(kind(1.D0)), intent(out) :: depth, TM(3,3)
      REAL(kind(1.D0)) :: sz(3), ksi, rn, x1, x2, x3      
      sz = OBJ%SIZE(1:3)
      TM(:,1) = (/1.D0, 0.D0, 0.D0/)
      TM(:,2) = (/0.D0, 1.D0, 0.D0/)
      TM(:,3) = (/0.D0, 0.D0, 1.D0/)
      select case(OBJ%SHAPE)
        case(FRAME_SHAPE_BOX)
          depth = 0.5*sz(3)-R(3)
        case(FRAME_SHAPE_DISC)
          x1 = R(1)*2.D0/sz(1)
          x2 = R(2)*2.D0/sz(2)
          rn = sqrt(R(1)**2 + R(2)**2)
          ksi = sqrt(x1**2 + x2**2)
          if (ksi.eq.0.D0) then
            depth = 0.5*sz(1)
            TM(:,1) = (/0.D0, 1.D0, 0.D0/)
            TM(:,2) = (/0.D0, 0.D0, 1.D0/)
            TM(:,3) = (/1.D0, 0.D0, 0.D0/)
          else
            depth = rn*(1.D0/ksi - 1.D0)
            TM(:,1) = (/-R(2), R(1), 0.D0/)/rn ! hoop
            TM(:,2) = (/0.D0, 0.D0, 1.D0/) ! axial
            TM(:,3) = (/R(1), R(2), 0.D0/)/rn ! radial
          endif
        case(FRAME_SHAPE_CYLINDER)
          x1 = R(1)*2.D0/sz(1)
          x2 = R(3)*2.D0/sz(3)
          rn = sqrt(R(1)**2 + R(3)**2)
          ksi = sqrt(x1**2 + x2**2)
          if (ksi.eq.0.D0) then
            depth = 0.5*sz(3)
          else
            depth = rn*(1.D0/ksi - 1.D0)
            TM(:,1) = (/-R(3), 0.D0, R(1)/)/rn ! hoop
            TM(:,2) = (/0.D0, 1.D0, 0.D0/) ! axial
            TM(:,3) = (/R(1), 0.D0, R(3)/)/rn ! radial
          endif
        case default ! FRAME_SHAPE_ELLIPSOID
          x1 = R(1)*2.D0/sz(1)
          x2 = R(2)*2.D0/sz(2)
          x3 = R(3)*2.D0/sz(3)
          rn = sqrt(R(1)**2 + R(2)**2 +R(3)**2)
          ksi = sqrt(x1**2 + x2**2 + x3**2)
          if (ksi.eq.0.D0) then
            depth = 0.5*sz(3)
          else
            depth = rn*(1.D0/ksi - 1.D0)
            TM(:,1) = (/-R(1), R(3), 0.D0/)/rn ! hoop
            TM(:,2) = (/0.D0, 1.D0, 0.D0/) ! axial
            TM(:,3) = (/R(1), R(2), R(3)/)/rn ! radial
            ! try to find hoop direction as normal to y_loc and radial
            if (abs(TM(2,3)-1.D0) < 1.D-3) then
               ! R is near y_loc => hoop is // x_loc
               TM(:,1) = (/1.D0, 0.D0, 0.D0/)
            else
               call V3cV3((/0.D0, 1.D0, 0.D0/), TM(:,3), TM(:,1))
            endif
            call V3cV3(TM(:,3), TM(:,1), TM(:,2)) ! axial
          endif
      end select
      
      end subroutine CONTAINER_COORD


      end MODULE FRAMES_TRACE

