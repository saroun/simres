!//////////////////////////////////////////////////////////////////////
!////  $Id: monitors_trace.f,v 1.18 2019/08/15 15:02:07 saroun Exp $
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
!////  Ray-tracing subroutines for MONITOR
!////
!////////////////////////////////////////////////////////////////////////
      MODULE MONITORS_TRACE
      use TRACELIB
      USE TRACINGDATA
      USE MONITORS
      use EVENTLOGS
      use FRAMES_TRACE
      use NSTORE
      implicit none
      private

      logical :: dbg=.false.
      integer :: idbg
      type(NEUTRON) :: TMPN(MONITORS_DIM)
      logical :: TMPIN(MONITORS_DIM)

      public MONITOR_GO,MONITOR_ADJUST,MONITOR_INIT
      public MONITOR_ENDTRACE,MONITOR_CLR
      public MONITOR_DISPOSE_ALL,AddMONITOR,MONITOR_INSIDE,MONITOR_CLIP

      contains


!-------------------------------------------------------------
      SUBROUTINE MONITOR_CLR(INST)
! return reference to given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(MONITOR),pointer :: OBJ
        if (.not.MONITOR_isValid(INST)) return
        OBJ => AMONITORS(INST)%X
        OBJ%IIN=0.D0
        OBJ%IOUT=0.D0
        OBJ%FRAME%COUNT=0
        call NSTORE_CLEAR(OBJ%FRAME%REGISTRY)
      end SUBROUTINE MONITOR_CLR


!-------------------------------------------------------------
      logical function MONITOR_GO(INST)
! call GO on given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
        MONITOR_GO=MONITOR_GO_INST(AMONITORS(INST)%X,INST)
      end function MONITOR_GO


!--------------------------------------------------------------------------
      logical function MONITOR_INSIDE(INST)
! true if event is inside clip area
! limits are applied in LOCAL coordinates !
!--------------------------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(NEUTRON) :: NEU
      TYPE(MONITOR),POINTER :: OBJ
      logical :: LOG1
      REAL(kind(1.D0)) :: X,Y
1     format(a,6(2x,G12.5))
      OBJ => AMONITORS(INST)%X
      NEU=NEUT
      if (TRACING_UP) NEU%K=-NEUT%K
      X=GETVAR(OBJ%IX,NEU,OBJ%K0)
      Y=GETVAR(OBJ%IY,NEU,OBJ%K0)
      LOG1=(abs(X-OBJ%X0)<0.5D0*OBJ%DX)
      if (LOG1) LOG1=(abs(Y-OBJ%Y0)<0.5D0*OBJ%DY)
      !LGX=((X>3.D0).and.(X<4.D0))
      !LGY=((Y>0.005).and.(Y<0.01))
      !dbg=((idbg<100).and.LGX.and.LGY)
      if (dbg) then
        write(*,1) 'MONITOR_INSIDE cnt=',NEU%CNT
        write(*,1) 'MONITOR_INSIDE x,kx=',NEU%R(1),NEU%K(1)/NEU%K0,NEU%P
        write(*,1) 'MONITOR_INSIDE ',OBJ%IX,OBJ%IY,X,Y,LOG1
        idbg=idbg+1
      endif
      MONITOR_INSIDE=LOG1
      end function MONITOR_INSIDE

!--------------------------------------------------------------------------
      subroutine MONITOR_CLIP(INST,NEU,dbg)
! set appropriate weight to the event
!--------------------------------------------------------------------------
      integer,intent(in) :: INST
      logical,intent(in) :: dbg
      TYPE(NEUTRON) :: NEU
      TYPE(MONITOR),POINTER :: OBJ
      logical :: LOG1
      REAL(kind(1.D0)) :: S
1     format(a,6(2x,G12.5))
      OBJ => AMONITORS(INST)%X
      if (OBJ%CNTMODE==monit_mode_add) return
      NEUT=NEU
      ! CALL SLIT_PRE(OBJ%FRAME)
      if (dbg) then
        write(*,1) 'MONITOR_CLIP cnt=',NEU%CNT
      endif
      LOG1=MONITOR_INSIDE(INST)
      select case(OBJ%CNTMODE)
      case(monit_mode_inner)
        if (.not.LOG1) NEU%P=0.D0
      case(monit_mode_outer)
        if (LOG1) NEU%P=0.D0
      case(monit_mode_sub)
        if (.not.LOG1) then
          NEU%P=0.D0
        else
          !if (NEU%CNT<100) write(*,1) 'MONITOR_CLIP ',OBJ%IIN,OBJ%IOUT
          S=(OBJ%IIN+OBJ%IOUT)
          if (S>0.D0) NEU%P=NEU%P*OBJ%IIN/S
        endif
      end select
      if (dbg) then
        write(*,1) 'MONITOR_CLIP x,kx=',NEU%R(1),NEU%K(1)/NEU%K0,NEU%P,LOG1
      endif

      end subroutine MONITOR_CLIP

!--------------------------------------------------------------------------
      subroutine MONITOR_ENDTRACE(INST,CNT,P)
! Add previously registered event to the monitor storage, using given weight
!--------------------------------------------------------------------------
      integer,intent(in) :: INST,CNT
      REAL(kind(1.D0)),intent(in)  :: P
      TYPE(PMONITOR) :: OBJ
        OBJ%X => AMONITORS(INST)%X
        if (OBJ%X%FRAME%REGISTRY>0) then
          ! transform to incident axis coordinates
          call FRAME_LOCAL(-1,OBJ%X%FRAME,TMPN(INST))
          TMPN(INST)%P=P     ! set the final weight to the event
          TMPN(INST)%CNT=CNT ! must set counter to neutron before calling NSTORE_SETEVENT !!
          if (TMPIN(INST)) then
            OBJ%X%IIN=OBJ%X%IIN+P
          else
            OBJ%X%IOUT=OBJ%X%IOUT+P
          endif
          call NSTORE_SETEVENT(OBJ%X%FRAME%REGISTRY,TMPN(INST))
        endif
      end subroutine MONITOR_ENDTRACE

!---------------------------------------------------------------------
      SUBROUTINE MONITOR_INIT(INST)
!----------------------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(MONITOR),POINTER :: OBJ
      if (.not.MONITOR_isValid(INST)) return
      OBJ => AMONITORS(INST)%X
      CALL FRAME_INIT(OBJ%FRAME)
      dbg=.false.
      idbg=0
      END SUBROUTINE MONITOR_INIT

!-------------------------------------------------
      SUBROUTINE MONITOR_ADJUST(INST,IERR)
! Update position of the component using nominal trajectory represented by NEUT
! Call FRAME_INIT_MAT first.
! Then use SLIT_PRE & SLIT_POST procedures to propagate NEUT
!-------------------------------------------------
      integer,intent(in) :: INST
      integer,intent(out) :: IERR
      TYPE(MONITOR),POINTER :: OBJ
1     format(a,': ',7(G13.6,1x))
      if (.not.MONITOR_isValid(INST)) return
      OBJ => AMONITORS(INST)%X
      IERR=0
      !call FRAME_INIT_MAT(OBJ%FRAME)
      !write(*,1) 'Move from '//trim(OBJ%FRAME%ID),NEUT%K0*NEUT%T
      call SLIT_PRE(OBJ%FRAME)
      !write(*,1) 'Inside '//trim(OBJ%FRAME%ID),NEUT%K0*NEUT%T, NEUT%R
      call SLIT_POST(OBJ%FRAME)
      OBJ%K0=NEUT%K0
      !write(*,1) 'Move to '//trim(OBJ%FRAME%ID),NEUT%K0*NEUT%T, OBJ%FRAME%DIST
      end SUBROUTINE MONITOR_ADJUST

!--------------------------------------------------
      LOGICAL FUNCTION MONITOR_GO_INST(OBJ,INST)
! tracing for curved MONITORs
!--------------------------------------------------
      !include 'const.inc'
      TYPE(MONITOR) :: OBJ
      LOGICAL :: LOG1
      integer,intent(in) :: INST

      LOG1=.FALSE.
1     format(a,': ',6(1x,G12.5))

      !dbg=(NINT(trace_tot)==6)
      if (dbg) write(*,1) '  NEUI',NEUI%R(1),NEUI%K(1),NEUI%R(3)
      if (dbg) write(*,1) 'MONITOR_GO',NEUT%R(1),NEUT%K(1),OBJ%FRAME%STA(1)

! move to the OBJ entry surface and transform to local coordinates
      CALL SLIT_PRE(OBJ%FRAME)
      if (dbg) write(*,1) 'MONITOR_GO',NEUT%R(1),NEUT%K(1)
! move neutron to z=0
      call SLIT_MIDDLE(.false.)

! log event at the entry of each object
      call AddEventLog(OBJ%FRAME,NEUT)

! check that neutron fits within the object
      LOG1=TLIB_INSIDE(NEUT%R(1:2),OBJ%FRAME%SIZE(1:2),OBJ%FRAME%SHAPE)
      TMPIN(INST)=LOG1
! handle blocking area
      if (LOG1.and.OBJ%BLOCKAREA) then
        TMPIN(INST)=MONITOR_INSIDE(INST)
        if (OBJ%CNTMODE==monit_mode_inner) then
          LOG1=TMPIN(INST)
          if (dbg) write(*,1) '   inside',LOG1
        else if (OBJ%CNTMODE==monit_mode_outer) then
          LOG1=(.not.TMPIN(INST))
        endif
      endif
      IF (LOG1) THEN
! register event
        TMPN(INST)=NEUT
        if (TRACING_UP) TMPN(INST)%K=-TMPN(INST)%K
! transform to exit axis reference frame
        CALL SLIT_POST(OBJ%FRAME)
        if (dbg) write(*,1) '   post',NEUT%R(1),NEUT%K(1)
        if (dbg) write(*,1) '   registered',LOG1
! increment counter
        call incCounter(OBJ%FRAME%COUNT)
      ELSE
        NEUT%P=0
      ENDIF

      MONITOR_GO_INST=LOG1
      dbg=.false.
      END FUNCTION MONITOR_GO_INST

      end MODULE MONITORS_TRACE
