!//////////////////////////////////////////////////////////////////////
!////  $Id: eventmonitor.f,v 1.16 2019/08/15 15:42:13 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2015, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.16 $
!////     $Date: 2019/08/15 15:42:13 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Handles a list of monitors of events flow
!////  BMONITORS store references to the storage in NSTORE, where
!////  the status of neutrons passing given component is saved for further use
!////
!//////////////////////////////////////////////////////////////////////
      MODULE EVENTMONITOR
      use FIELDDEF
      use FRAMES
      use COMPONENTS
      use VCTABLE_INSTRUMENTS
      use VCTABLE_COMPONENTS
      use VCTABLE_SAMPLES
      implicit none
      private
      SAVE

      ! view coordinates (local, incident or world reference frames)
      INTEGER, parameter :: coord_axis=0 ! default
      INTEGER, parameter :: coord_local=1
      INTEGER, parameter :: coord_world=2
      INTEGER :: coord=coord_axis ! actual frame for plotting and event dump

! references to beam monitors
      integer, parameter :: BMONITORS_DIM=MONITORS_DIM+2
      integer :: BMONITORS_N
      integer :: BMONITOR_PRI,BMONITOR_SEC
      TYPE(TCOMPID) :: BMONITORS(BMONITORS_DIM)
      integer :: BMONITOR_REF=1 ! reference to current data storage for BEAM_1D plots


      public BMONITORS,BMONITOR_REF,REGISTER_BMONITORS,SET_BMONITOR
      public COMPONENTS_SETMCNORM,COMPONENTS_ALLOCATE,BMONITOR_PRI,BMONITOR_SEC
      public coord_axis, coord_local, coord_world
      public BMONITORS_SETCOORD,BMONITORS_GETEVENT,BMONITORS_DUMPMCPL,BMONITORS_READMCPL
      contains



!---------------------------------------------------
      subroutine BMONITORS_READMCPL(FNAME, STORAGE)
! Read events from a MCPL file
! This procedure makes use of "mcplio.dll", which encapsulates
! MCPL library from https://mctools.github.io/mcpl
! mcplio.dll is loaded dynamically on each call to this subroutine
! binding to this dll is implemented in mcplhandle.c
! FNAME = full path name for the output file
! C = coordinate system
! FILT = filter flag
!---------------------------------------------------
      character*(*) :: FNAME
      integer, intent(in) :: STORAGE
      character(256) :: fn
      integer :: dll,i,NCNT,ISTORE, REF
      type(NEUTRON) :: NEU
      REAL(KIND(1.D0)) ::  SUMA, NMAX
      REAL(KIND(1.D0)) ::  r(3), k(3), s(3), t, p, K0(3)
1     format('Loaded from ',a,': records=',G12.6,' counts=',G12.6, 'istore=',G12.4)
2     format('File ',a,' contains ',G12.6,' records.')
3     format('r=',3(1x,G12.6),'k=',3(1x,G12.6),'s=',3(1x,G12.6),'t=',G12.6,'p=',G12.6)

      INTERFACE
        FUNCTION LOADMCPL()
        END FUNCTION LOADMCPL
      END INTERFACE
      select case(STORAGE)
      case(1)
        REF = BMONITOR_PRI
      case(2)
        REF = BMONITOR_SEC
      case default
        REF = BMONITOR_REF
      end select
      !write(*,*) 'BMONITORS_READMCPL REF=',REF
      if (REF>0) then
        ISTORE=BMONITORS(REF)%ISTORE
        if (ISTORE>0) then
          dll=LOADMCPL()
          if (dll>=0) then
            fn=trim(FNAME)//char(0)
            call mcplopenread(trim(fn))
            call mcplnmax(NMAX)
            NCNT = NINT(NMAX)
            write(*,2) trim(FNAME),NCNT
            K0 = (/0.D0, 0.D0, INSOBJ%P_SPEC%KF/)
            call NREC_ALLOCATE(ISTORE,K0,NINT(NMAX),.false.)
            call NSTORE_SETNORM(ISTORE,1.D0)
            SUMA=0.D0
            DO i=1,NCNT
              call MCPLGET(r, k, s, t, p)
              NEU%R = r
              NEU%K = k
              NEU%S = s
              NEU%T = t*HOVM
              NEU%P = p
              NEU%K0=SQRT(NEUT%K(1)**2+NEUT%K(2)**2+NEUT%K(3)**2)
              NEU%REF = 0
              NEU%PHI = 0.D0
              NEU%CNT = i
              NEU%T0 = TOF_ZERO
              call NSTORE_SETEVENT(ISTORE,NEU)
              SUMA=SUMA+NEU%P
            enddo
            PRIMARY_INT = SUMA
            call mcplcloseread
            write(*,1) trim(FNAME), NCNT, SUMA, ISTORE
          else
                write(*,*) 'Cann''t load libmcplio: ',dll
              endif
        ENDIF
      endif
      end subroutine BMONITORS_READMCPL

!---------------------------------------------------
      subroutine BMONITORS_DUMPMCPL(FNAME,C,FILT,STORAGE)
! Dump content of the refernce monitor in a MCPL file
! This procedure makes use of "mcplio.dll", which encapsulates
! MCPL library from https://mctools.github.io/mcpl
! mcplio.dll is loaded dynamically on each call to this subroutine
! binding to this dll is implemented in mcplhandle.c
! FNAME = full path name for the output file
! C = coordinate system
! FILT = filter flag
! STORAGE = storage: (1) BMONITOR_PRI, (2) BMONITOR_SEC, else BMONITOR_REF
!---------------------------------------------------
      character*(*) :: FNAME
      integer, intent(in) :: C,FILT,STORAGE
      character(256) :: fn
      integer :: dll,i,NCNT,ISTORE, REF
      type(NEUTRON) :: NEU
      REAL(KIND(1.D0)) ::  CNORM, SUMA
1     format('Exported to ',a,': records=',G12.6,' counts=',G12.6,' istore=',G12.6)
2     format('Warning: requested monitor ',G10.4,' not defined. No MCPL export.')
3     format('Warning: requested monitor ',G10.4,' is empty. No MCPL export.')
10    format(a,': ',10(1x,G12.5))
      INTERFACE
        FUNCTION LOADMCPL()
        END FUNCTION LOADMCPL
      END INTERFACE
      select case(STORAGE)
      case(1)
        REF = BMONITOR_PRI
      case(2)
        REF = BMONITOR_SEC
      case default
        REF = BMONITOR_REF
      end select
      !write(*,*) 'BMONITORS_DUMPMCPL REF=',REF
      if (REF>0) then
        ISTORE=BMONITORS(REF)%ISTORE
        NCNT=NSTORE_MAXN(ISTORE)
        if (NCNT>0) then
          if (C>=0) then
            call BMONITORS_SETCOORD(C)
          endif
          CNORM=NSTORE_GETNORM(ISTORE)
          dll=LOADMCPL()
          if (dll>=0) then
            fn=trim(FNAME)//char(0)
            call mcplopenwrite(trim(fn))
            SUMA=0.D0
            DO i=1,NCNT
              call BMONITORS_GETEVENT(REF,i,FILT,NEU)
              !if (i<5) then
              !   write(*,10) 'DUMPMCPL',NEU%R,NEU%T/HOVM,NEU%P,CNORM
              !endif
              call mcpladd(NEU%R,NEU%K,NEU%T/HOVM,NEU%P*CNORM)
              SUMA=SUMA+NEU%P*CNORM
            enddo
            call mcplclosewrite
            write(*,1) trim(FNAME),NCNT,SUMA,ISTORE
          else
                write(*,*) 'Cann''t load libmcplio: ',dll
              endif
        else
          write(*,2) STORAGE
        ENDIF
      else
        write(*,2) STORAGE
      endif
      end subroutine BMONITORS_DUMPMCPL


!---------------------------------------------------
      subroutine BMONITORS_SETCOORD(icor)
! define coordinates (local or incident frame)
!---------------------------------------------------
      integer, intent(in) :: icor
      select case (icor)
      case(coord_local)
        coord=coord_local
      case(coord_axis)
        coord=coord_axis
      case(coord_world)
        coord=coord_world
      case default
        coord=coord_axis
      end select
      end subroutine BMONITORS_SETCOORD


!---------------------------------------------------
      subroutine BMONITORS_GETEVENT(IMON,NEV,FILT,NEU)
! get neutron from given monitor
! IMON = monitor index in the BMONITORS list
! NEV = event index in the associated storage
! FILT = filter flag
!    <0 ... all events
!    >=0 ... only events with corresponding LABEL value
!---------------------------------------------------
      integer, intent(in) :: NEV,IMON,FILT
      type(NEUTRON), intent(out) :: NEU
      TYPE(PFRAME):: PF
      TYPE(TFRAME):: F
      type(NEUTRON) :: N
      logical :: trans
      !dbg=(NEV<5)
      ! get neutron data from storage
      call NSTORE_GETEVENT(BMONITORS(IMON)%ISTORE,NEV,NEU)
      ! apply clipping for monitors
      if (BMONITORS(IMON)%OBJ%CCLS==CCLS_MONITOR) then
        call MONITOR_CLIP(BMONITORS(IMON)%OBJ%INST,NEU,.false.)
      endif
      ! transform to local coordinates when required
      select case(coord)
      case(coord_local,coord_world)
        trans=.false.
        if (BMONITORS(IMON)%OBJ%CCLS==0) then ! sample ?
          F=SAMOBJ%P_FRAME
          trans=(SAMOBJ%ICLS.gt.0)
        else
          call PCOMOBJ_GET(BMONITORS(IMON)%OBJ,PF)
          if (associated(PF%X)) then
            F=PF%X
            trans=.true.
          endif
        endif
        if (trans) then
          call FRAME_LOCAL(1,F,NEU)
          if (coord==coord_world) then
            N=NEU
            call Loc2Glob(F,N,NEU)
          endif
        endif
      end select
      ! apply filter
      if (FILT>0) then
        if (iand(NEU%LABEL,FILT).eq.0) NEU%P=0.D0
        ! if (NEU%LABEL.ne.FILT) NEU%P=0.D0
      else if (FILT==0) then
        if (NEU%LABEL.ne.0) NEU%P=0.D0
      endif
      end subroutine BMONITORS_GETEVENT


!--------------------------------------------------------
      logical function is_BMONITOR(IORD,IBEAM)
! true if the CCLS class can be an event monitor
! i.e. it can register its own event storage
!--------------------------------------------------------
      integer,intent(in) :: IORD,IBEAM
      logical :: LOG1
      integer :: CCLS
      CCLS=BEAMLINE(IORD,IBEAM)%OBJ%CCLS
      LOG1=.false.
      select case (CCLS)
    ! all MONITOR components
      case (CCLS_MONITOR)
        LOG1=.true.
    ! other components only if they are at the end
      case default
        select case (IBEAM)
        ! primary beam: only if there is no sample and no secondary beamline
        case (1)
          LOG1=(IORD==BEAMLINE_NC(1))
          LOG1=(LOG1.and.(BEAMLINE_NC(2)==0))
          LOG1=(LOG1.and.(SAMOBJ%ICLS==0))
        ! secondary beam: only the last component
        case (2)
          LOG1=(IORD==BEAMLINE_NC(2))
          LOG1=(LOG1.and.(BEAMLINE_NC(2)>0))
        end select
      end select
      is_BMONITOR=LOG1
      end function is_BMONITOR

!--------------------------------------------------------
      subroutine REGISTER_BMONITOR(IORD,IBEAM)
! register component as a monitor if required
!--------------------------------------------------------
      integer,intent(in) :: IORD,IBEAM
      if (is_BMONITOR(IORD,IBEAM)) then
        BMONITORS_N=BMONITORS_N+1
        BEAMLINE(IORD,IBEAM)%IREG=BMONITORS_N
        BMONITORS(BMONITORS_N)=BEAMLINE(IORD,IBEAM)
      else
        BEAMLINE(IORD,IBEAM)%IREG=0
      endif
      end subroutine REGISTER_BMONITOR

!--------------------------------------------------------
      subroutine CLEAR_BMONITORS
! clear all registered monitors
!--------------------------------------------------------
      integer :: i
      DO i=1,BMONITORS_DIM
        call pcom_clear(BMONITORS(i)%OBJ)
        BMONITORS(i)%IBEAM=0
        BMONITORS(i)%IORD=0
        BMONITORS(i)%ISTORE=0
        BMONITORS(i)%IREG=0
      enddo
      BMONITORS_N=0
      BMONITOR_PRI=0
      BMONITOR_SEC=0
      end subroutine CLEAR_BMONITORS

!--------------------------------------------------------
      subroutine REGISTER_BMONITORS
! register all monitors
! monitors can be only:
! a) all MONITOR components
! b) the last component on the beamline
! c) SAMPLE, if present (registeres two monitors for the initial and final states)
!--------------------------------------------------------
      integer :: i
      TYPE(PFRAME) :: PF
      character(LEN_ID) :: ID
! clear BMONITORS registry
      call CLEAR_BMONITORS
! register monitors with event storage capability from primary beamline
      do i=1,BEAMLINE_NC(1)
        call REGISTER_BMONITOR(i,1)
      enddo
      BMONITOR_PRI=BMONITORS_N
      i=BEAMLINE_NC(1)
! add sample if exists
      if (SAMOBJ%ICLS.gt.0) then
        BMONITORS_N=BMONITORS_N+1
        BMONITORS(BMONITORS_N)%OBJ%INST=1
        BMONITORS(BMONITORS_N)%OBJ%ID='SAMPLE_I'
        BMONITOR_PRI=BMONITORS_N
        BMONITORS_N=BMONITORS_N+1
        BMONITORS(BMONITORS_N)%OBJ%INST=1
        BMONITORS(BMONITORS_N)%OBJ%ID='SAMPLE_F'
        BMONITOR_SEC=BMONITORS_N
      endif
! register monitors with event storage capability from secondary beamline
      do i=1,BEAMLINE_NC(2)
        call REGISTER_BMONITOR(i,2)
      enddo
! register storage for all monitors
      do i=1,BMONITORS_N
        ID=trim(BMONITORS(i)%OBJ%ID)
      ! for beamline components with event storage capability, set pointer to the storage
        if (BMONITORS(i)%IBEAM>0) then
          BMONITORS(i)%ISTORE=NSTORE_REGISTER(ID)
          call PCOMOBJ_GET(BMONITORS(i)%OBJ,PF)
          if (PF%IDX>0) PF%X%REGISTRY=BMONITORS(i)%ISTORE
      ! for sample: register two storages for initial and final states
        else if (trim(ID).eq.'SAMPLE_I') then
          SAMOBJ%ISTORE(1)=NSTORE_REGISTER(ID)
          BMONITORS(i)%ISTORE=SAMOBJ%ISTORE(1)
        else if (trim(ID).eq.'SAMPLE_F') then
          SAMOBJ%ISTORE(2)=NSTORE_REGISTER(ID)
          BMONITORS(i)%ISTORE=SAMOBJ%ISTORE(2)
        endif
      enddo
! ensure that initial and final event references are assigned
! without sample, set both to the last component register
      if (BMONITOR_PRI==0) BMONITOR_PRI=BMONITORS_N
      if (BMONITOR_SEC==0) BMONITOR_SEC=BMONITORS_N
! set the last monitor as the current one
      BMONITOR_REF=BMONITORS_N

      end subroutine REGISTER_BMONITORS

!---------------------------------------------------
      subroutine SET_BMONITOR(ID)
! set reference to current beam monitor
!---------------------------------------------------
      character(*) :: ID
      integer :: i
        i=1
        do while (i<=BMONITORS_N)
          if (trim(BMONITORS(i)%OBJ%ID).eq.trim(ID)) then
            BMONITOR_REF=BMONITORS(i)%ISTORE
          endif
          i=i+1
        enddo
      end subroutine SET_BMONITOR


!--------------------------------------------------------
      SUBROUTINE COMPONENTS_SETMCNORM(IBEAM,NORM)
! set MCNORM to components event storage
!--------------------------------------------------------
      integer, intent(in) :: IBEAM
      real(KIND(1.D0)), intent(in) :: NORM
      integer:: i,iref
      do i=1,BMONITORS_N
        if (IBEAM==BMONITORS(i)%IBEAM) then
          IREF=BMONITORS(i)%ISTORE
          call NSTORE_SETNORM(IREF,NORM)
        endif
      enddo
      end SUBROUTINE COMPONENTS_SETMCNORM

!--------------------------------------------------------
      SUBROUTINE COMPONENTS_ALLOCATE(IBEAM,K0,CNT)
! allocate space for components event storage
!--------------------------------------------------------
      integer, intent(in) :: IBEAM,CNT
      real(KIND(1.D0)), intent(in) :: K0(3)
      integer:: i,iref
      do i=1,BMONITORS_N
        if (IBEAM==BMONITORS(i)%IBEAM) then
          IREF=BMONITORS(i)%ISTORE
          call NREC_ALLOCATE(IREF,K0,CNT,.false.)
        endif
      enddo
      end SUBROUTINE COMPONENTS_ALLOCATE

      end module EVENTMONITOR