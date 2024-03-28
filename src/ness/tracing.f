!//////////////////////////////////////////////////////////////////////
!////  $Id: tracing.f,v 1.138 2019/08/13 15:41:18 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.138 $
!////     $Date: 2019/08/13 15:41:18 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Ray-tracing control
!////
!//////////////////////////////////////////////////////////////////////
      MODULE TRACING
      USE TRACINGDATA
      use TRACINGOPT
      USE NSTORE
      USE COMPONENTS
      use EVENTMONITOR
      USE EVENTLOGS
      use RESMAT
      use REPORTS
      use INSTCONTROL
      use XMLINFO
      use RESULTS
      use FRAMES_TRACE
      use BEAM1D
      use mtmod

      implicit none

c locally used variables
  ! intermediate sums
      INTEGER, parameter :: MPARC=50
      integer :: IPARC ! index of actually filled partial sums
!      integer :: NPARC(0:MPARC) ! cummulative number of trials for each partial cycle
      integer :: CPARC(0:MPARC) ! cummulative number of counts for each partial cycle
      real(KIND(1.D0)) :: SPARC(0:MPARC) ! sum of event weights for each partial cycle
      real(KIND(1.D0)) :: SPARCV(0:MPARC) ! sum of event weights for each partial cycle
      real(KIND(1.D0)) :: PARCINT(MPARC)  ! intensity
      real(KIND(1.D0)) :: PARCVOL(MPARC)  ! volumes
      real(KIND(1.D0)) :: PARCDE(MPARC)   ! energy
      real(KIND(1.D0)) :: PARCDDE(MPARC)  ! energy**2
      real(KIND(1.D0)) :: PARCTOF(MPARC)  ! TOF
      real(KIND(1.D0)) :: PARCDTOF(MPARC)  ! TOF**2


! internal variables
      integer,private :: icycle   ! cycle index for variance reduction mode
      logical,private :: TWOSTEPS
      integer,private :: trigger_count


      CONTAINS

!--------------------------------------------------------
      SUBROUTINE TRACING_INI
! General initialization procedure
! called just before starting the tracing loop after ADJUST calls
!--------------------------------------------------------
      logical,parameter :: dbg=.false.
      integer :: IUDBG = 6
      REAL(KIND(1.D0)) :: EIVALS(CRND),EIVECT(CRND,CRND)
      integer :: N
      type(PDATASET) :: DTS
      NEUT_DBG=.false.
      TOF_ADJ_T0=TROPT%ADJT0

      if (DBG) then
      !  IUDBG=OPENFILEUNIT('TRACING_INI.log',.false.)
      !  if (IUDBG.le.0) IUDBG=6
        IUDBG=6
      endif
      call ResetEventLog
! set gravity parameters
      gravity=TROPT%GRAVITY*G_EARTH
      use_gravity=(TROPT%GRAVITY.ne.0.D0)

! call INIT procedure on all components on the beamline
      if (TROPT%IMODE.ne.tr_secondary) call COMPONENTS_INIT(1)
      if (SAMOBJ%ICLS.gt.0) call PSAMOBJ_INIT(SAMOBJ)
      if (TROPT%IMODE.ge.tr_secondary) call COMPONENTS_INIT(2)

! Initialize generator
      CALL GENERATOR_INIT
      if (DBG)  write(IUDBG,*) 'TRACING_INI GENERATOR_INIT OK'

! Use RESMAT to optimize sampling
  ! creates matrix description of the whole instrument
      call RESMAT_ADDCOMPONENTS

      if (DBG)  write(IUDBG,*) 'TRACING_INI RESMAT_ADDCOMPONENTS OK'

  ! get transformation matrices from the matrix model
  ! N is the number of free parameters = dimension of generated random vectors
      call RESMAT_EIGEN(EIVECT,EIVALS,CRND,N)

      if (DBG)  write(IUDBG,*) 'TRACING_INI RESMAT_EIGEN OK'

! generator uses EIVECT,EIVALS to create transformation from random numbers
! to physical variables (phase-space coordinates etc.)

      call GENERATOR_SETTMAT(EIVECT,EIVALS,CRND,N)
      call GENERATOR_SETMEAN(XY_MEAN,KXY_MEAN,KZ_MEAN)
      if (dbg) then
        write(IUDBG,*)
        write(IUDBG,*) 'RNDIST=      ',RNDIST
        write(IUDBG,*) 'RNDIST_NORM= ',RNDIST_NORM
        write(IUDBG,*) 'RNDIST_RANGE=',RNDIST_RANGE
        write(IUDBG,*) 'RNDIST_VOL=  ',RNDIST_VOL
        write(IUDBG,*) 'TROPT%IMODE=  ',TROPT%IMODE
      endif


! filter options: some optimization methods make sense only for unbiased sampling
      if (RNDIST.ne.rndist_uni) then
        TROPT%MAXV=.FALSE.
      ! TROPT%VARI=.FALSE.
        TROPT%SWPOOL=.FALSE.
      endIF

! set count targets: optimization events are triggered by reaching counting thresholds
      if (TROPT%IMODE.eq.tr_secondary) then
        TCOUNT_REQ = EVENTGEN%IMAX*TROPT%SPLIT
        write(*,*) 'Secondary beam: requested counts = ',TCOUNT_REQ
      else
        TCOUNT_REQ = TROPT%CNT
      endif
  ! variance reduction
      if (TCOUNT_REQ.LT.400) then
        TROPT%VARI=.FALSE.  ! no variance reduction if TROPT%CNT is small
      else if (TCOUNT_REQ.LE.10000) then
        CNT_VARI=500
      else
        CNT_VARI=1000
      endif
  ! max. value optimization
      MAXV_UPD_FLAG=TROPT%MAXV
      if (TCOUNT_REQ.LT.200) then
        MAXV_UPD_FLAG=.FALSE.  ! no max. val. optimization if TCOUNT_REQ is small
      else if (TCOUNT_REQ.LE.10000) then
        MAXV_UPD_COUNT=100
      else
        MAXV_UPD_COUNT=200
      endif
  ! timeout
      CNT_TOUT=MAX(1.D7,1.D5*TCOUNT_REQ)
      TRACING_PRECISION=1.D0
      trigger_count=max(100,min(nint(TROPT%ELIMIT),5000))

! set safety pool flag
      SAFETY_POOL_FLAG=TROPT%SWPOOL

! allocate memory for events storage
      DTS=DSET_DATA(DSET_ISEL)
      ! write(*,*) 'TRACING_INI ',DSET_ISEL,DTS%P%REG1,DTS%P%REG2

    ! primary beamline
      if (TROPT%IMODE.eq.tr_primary) then
        call NREC_ALLOCATE(DTS%P%REG1,(/0.D0,0.D0,DTS%P%KI/),TCOUNT_REQ,.false.)
        call COMPONENTS_ALLOCATE(1,(/0.D0,0.D0,DTS%P%KI/),TCOUNT_REQ)
      endif
    ! both beamlines
      if (TROPT%IMODE.eq.tr_all) then
        call NREC_ALLOCATE(DTS%P%REG1,(/0.D0,0.D0,DTS%P%KI/),TCOUNT_REQ,.false.)
        call COMPONENTS_ALLOCATE(1,(/0.D0,0.D0,DTS%P%KI/),TCOUNT_REQ)
        call NREC_ALLOCATE(DTS%P%REG2,DTS%P%KF,TCOUNT_REQ,.false.)
        call COMPONENTS_ALLOCATE(2,DTS%P%KF,TCOUNT_REQ)
      endif
    ! secondary beamline
      if (TROPT%IMODE.eq.tr_secondary) then
        call NREC_ALLOCATE(DTS%P%REG2,DTS%P%KF,TCOUNT_REQ,.false.)
        call COMPONENTS_ALLOCATE(2,DTS%P%KF,TCOUNT_REQ)
      endif

      call REGISTER_LOGGERS

      IF (.not.REPOPT%MUTE) THEN
        WRITE(SOUT,*) 'Monte-Carlo variables initialized by TRACING_INI.'
!        CALL WRITE_BEAMLINE(20)
c        CALL WRITE_XML(20)
      ENDIF

    ! report transformation matrix
      if(dbg) then
        call XML_VALUE(IUDBG,'DIMENSION','int',' ',RNDLIST%DIM,0.D0,' ')
        call XML_MATRIX(IUDBG,N,N,'TMAT',MATLEG,MATLEG,TMAT(1,1),CRND)
        call XML_ARRAY(IUDBG,N,'LIMIT',MATLEG,RNDLIST%LIMITS(1))
      endif

      if (DBG) then
        close(IUDBG)
      endif

      END SUBROUTINE TRACING_INI



C------------------------------------------------------------------------
      SUBROUTINE TRACING_PROC
C Main procedure for ray-tracing job
C------------------------------------------------------------------------
      integer :: IMODE,NS1,NS2,ND
    !  real(KIND(1.D0)) :: FLUX
    !  type(TREPOPT) :: TMPREP

    !  TWOSTEPS = .false.
      !write(*,*) 'TRACING_PROC '

      if (TROPT%IMODE>tr_all) then
        CALL MSG_WARN('Cannot run simulation for required traing option mode.',1)
        return
      endif
      call XML_STATUS(SMES, "RUNNING", " ", " ")
      call XML_SWITCH(SMES,'MCRUN',1)
      NS1=CSUBSET(1)
      NS2=CSUBSET(2)
      ND=DSET_ISEL
      MCNORM(1:2,NS2,ND)=0.D0
      MCVOLI(1:2,NS1,ND)=0.D0
      MCVOLF(1:2,NS2,ND)=0.D0
      TWOSTEPS = ((TROPT%DBC).and.(TROPT%IMODE.eq.tr_all).and.(SAMOBJ%ICLS.gt.0))
      if (TWOSTEPS) then
    ! save mode
      !  TMPREP=REPOPT
      !  call REPOPT_MUTE(REPOPT)
        IMODE=TROPT%IMODE
    ! trace primary beam
        TROPT%IMODE=tr_primary
        FLUX_NORM=FLXNORM
        call TRACING_PROC_PARC
        PRIMARY_INT=TRACING_INT
        !write(*,*) 'TRACING_PROC ' ,TRACING_INT,TIMEOUT
        if ((TRACING_INT>0.D0).and.(.not.TIMEOUT)) then
    ! trace secondary beam
      !  REPOPT=TMPREP
          TROPT%IMODE=tr_secondary
          FLUX_NORM=PRIMARY_INT   !*MCTRIALS(NS,ND)
          call TRACING_PROC_PARC
        else
          call REPORT_PROG(2)
        endif
        TROPT%IMODE=IMODE
      else
        if (TROPT%IMODE.eq.tr_secondary) then
          ! write(*,*) 'TRACING_PROC secondary', PRIMARY_INT
          if (PRIMARY_INT>0.D0) then
            FLUX_NORM=PRIMARY_INT
            call TRACING_PROC_PARC
          else
            call MSG_WARN('Can''t trace the secondary part, there are no primary events',1)
          endif
        else
        !  IMODE=TROPT%IMODE
    ! trace primary beam
        !  TROPT%IMODE=tr_primary
          FLUX_NORM=FLXNORM
          call TRACING_PROC_PARC
        !  TROPT%IMODE=IMODE
        endif
        if (TROPT%IMODE.eq.tr_primary) PRIMARY_INT=TRACING_INT
      endif
      call NSTORE_XML_LIST(SMES,BMONITOR_REF)
      call XML_SWITCH(SMES,'MCRUN',0)
      call XML_STATUS(SMES, "READY", " ", " ")
      end SUBROUTINE TRACING_PROC

C------------------------------------------------------------------------
      SUBROUTINE TRACING_PROC_PARC
C Main procedure for ray-tracing job:
C controls events flow, counting, logging events etc.
C------------------------------------------------------------------------
      integer, parameter :: LERR = 10  ! loops for error estimation
      LOGICAL :: END_COUNT,END_PRECISION
      logical :: DBG
100   format('.',$)

      DBG=.false.
      if (DBG)  write(*,*) 'TRACING_PROC start IMODE=',TROPT%IMODE
  ! update dataset definition in DATASETS module
      call DEFINE_DSET(DSET_ISEL)
      if (DBG)  write(*,*) 'TRACING_PROC DSET ',DSET_ISEL
  ! initialize beamline
      call TRACING_INI
      if (DBG)  write(*,*) 'TRACING_PROC INI OK'

  ! set optional 2nd cycle for variance reduction technique
      if (TROPT%VARI) then
        icycle=1
        TCOUNT=CNT_VARI
        TCOUNT_ESTIM=CNT_VARI+TCOUNT_REQ
c        write(*,*) 'cycle 1 TCOUNT ',TCOUNT
      else
        icycle=2
        TCOUNT=TCOUNT_REQ
        TCOUNT_ESTIM=TCOUNT_REQ
c        write(*,*) 'cycle 2 TCOUNT ',TCOUNT
      ENDIF
      if (TWOSTEPS.and.(TROPT%IMODE.eq.tr_primary)) TCOUNT_ESTIM=2*TCOUNT_ESTIM


      if (DBG) then
        write(*,*) 'MAXV_UPD_FLAG  ',MAXV_UPD_FLAG
        write(*,*) 'MAXV_UPD_COUNT   ',MAXV_UPD_COUNT
        write(*,*) 'CNT_VARI ',CNT_VARI
        write(*,*) 'CNT_TOUT ',CNT_TOUT
        write(*,*) 'TCOUNT           ',TCOUNT
      !  read(*,*)
      endif

  ! stop event logging if it is still running
      call StopEventLog

  ! report sampling volume
      call REPORT_VOL('INITIAL')
    ! report initial transformation matrix
      call REPORT_VARI


  ! start progress bar
      if ((.not.TWOSTEPS).or.(TROPT%IMODE.ne.tr_secondary)) call REPORT_PROG(0)

  ! reset counters
10    call InitEventCounters
      if (DBG) then
        write(*,*) 'CYCLE=  ',icycle
      endif
      HIT=0
  ! update progess bar
      call REPORT_PROG(1)
  ! clear event counters in components

      if (DBG)  write(*,*) 'TRACING_PROC goto clear'

      if (TROPT%IMODE.ne.tr_secondary) CALL COMPONENTS_CLR(1)
      if (SAMOBJ%ICLS.gt.0) SAMOBJ%P_FRAME%COUNT=0
      if (TROPT%IMODE.ge.tr_secondary) CALL COMPONENTS_CLR(2)
      if (DBG)  write(*,*) 'TRACING_PROC goto clear OK'

  ! reset range checker
      CALL MAXV_UPD(0)
  ! reset covariance counter
      if (TROPT%VARI) CALL COV_CLR

  ! reset tracing end triggers
      TIMEOUT=.FALSE.
      END_COUNT=.FALSE.
      END_PRECISION=.FALSE.
      STOPTRACING=.FALSE.

      if (DBG)  write(*,*) 'TRACING_PROC STARTING CYCLE, dbg=',NEUT_DBG
! *****  TRACING CYCLE  *****
      do while ((.not.END_COUNT).and.(.not.END_PRECISION).and.(.not.TIMEOUT).and.(.not.STOPTRACING))
    ! log events if REPOPT%ELOG=true
        ! NEUT_DBG=(icycle==1)
        call StartEventLog
    ! TRACE EVENT
        !if (NEUT_DBG)  write(*,*) 'TRACING_PROC INI OK'
        call TRACING_RUN
    ! call actions triggered by trace_clock
        call TRACING_OnTrigger
    ! check end conditions
        END_COUNT=(trace_cnt.ge.TCOUNT)
        END_PRECISION=(TRACING_PRECISION.lt.TROPT%PLIMIT)
        if (END_PRECISION) write(*,*) 'precission limit: ',TRACING_PRECISION,IPARC
      !  if (mod(trace_tot,100000.D0).eq.0) then
      !    write(*,*) 'TRACING_PROC_PARC ',trace_tot,trace_cnt
      !  endif
      enddo
  ! tracing stopped for other reason than count target => don't run 2nd cycle
      if (.not.END_COUNT) icycle=2
      if (DBG) then
        write(*,*) 'END CONDITIONS: ',END_COUNT,END_PRECISION,TIMEOUT,STOPTRACING
        write(*,*) 'trace_cnt: ',trace_cnt,', tcount: ',TCOUNT, ', trace_tot: ',trace_tot
      endif

! *****  END OF TRACING CYCLE  *****

    ! restore tracing to DOWNSTREAM
      TRACING_UP=.FALSE.


! Variance reduction optimization - prepare the 2nd cycle:
      if (icycle.ne.2) then
        icycle=2
    ! calculate new covariances from the previous run
        CALL COV_NEW
    ! report new volume
        call REPORT_VOL('NEW')
    ! report new transformation matrix
        call REPORT_VARI
    ! set new count target
        TCOUNT=TCOUNT_REQ
    ! reset MAXV optimization flag
        MAXV_UPD_FLAG=TROPT%MAXV
        goto 10
      endif
! close progress bar
      if ((.not.TWOSTEPS).or.(TROPT%IMODE.ne.tr_primary)) call REPORT_PROG(2)
    ! report final volume
        call REPORT_VOL('FINAL')
    ! report final transformation matrix
        call REPORT_VARI
! report tracing statistics
      call REPORT_STAT
! finalization tasks
      call TRACING_FINAL
      END SUBROUTINE TRACING_PROC_PARC

C------------------------------------------------------------------------
      SUBROUTINE TRACING_RUN
C run a single event
C------------------------------------------------------------------------
! locals
      LOGICAL :: LOG1
      INTEGER :: ihit
      REAL(KIND(1.D0)) :: Z
1     format('MAXV optimization: factor ',G10.4)
2     format(a,': ',G14.8,2x,8(G10.4,1x))
3     format(a,': ',8(G10.4,1x))

! generate the initial neutron state in NEUT
      call GENERATOR_GO

!      if (isProbe) then
!        write(*,3) 'TRACING_RUN probe hit: ',NEUT%R(1),NEUT%R(3),NEUT%K(1),NEUT%K(3),NEUT%P
!        READ(*,*)
!      endif

! count started event
      call CountInitEvent
      call StartNewRayLog

      !if (NINT(trace_tot)==6) then
      !  write(*,3) 'TRACING_RUN: ',NEUT%R(1),NEUT%R(3),NEUT%K(1),NEUT%K(3),NEUT%P
      !endif

! do tracing:
      LOG1=.true.
      if (TROPT%IMODE.ne.tr_secondary) then
        if (LOG1) LOG1=TRACING_GO(1)
      endif
      if (TROPT%IMODE.ge.tr_secondary) then
        if (LOG1) LOG1=TRACING_GO(2)
      endif

! if successful:
      if (LOG1) then

    !  if (isProbe) then
    !    write(*,3) 'TRACING_RUN probe hit: ',NEUT%R(1),NEUT%R(3),NEUT%K(1),NEUT%K(3),NEUT%P
    !    READ(*,*)
    !  endif
        !if (NINT(trace_tot)==6) then
        !  write(*,3) 'TRACING_RUN OK',NEUT%R,NEUT%K
        !endif

  ! update range checker
        CALL MAXV_UPD(1)
  ! check pool: if hit, reset counters and covariances
        IF (SAFETY_POOL(ihit)) THEN
          call TRACING_CLEAR
          RETURN
        endif

  ! check max. value log and update limits when required
        if ((MAXV_UPD_FLAG).and.(MAXV_UPD_COUNT.eq.trace_cnt)) then
          Z=RNDIST_VOL
          CALL MAXV_UPD(2)
          call TRACING_CLEAR
          if (REPOPT%MAXV) write(*,1) Z/RNDIST_VOL
          MAXV_UPD_FLAG=.FALSE.
          ! report sampling volume
          call REPORT_VOL('UPDATED')
          RETURN
        endif
  ! count valid event
        call CountPassedEvent()
! if not successful:
      else
        call ResetRayLog
      endif

      END SUBROUTINE TRACING_RUN

!------------------------------------------------------------------------
      SUBROUTINE TRACING_REP(MSG, NEU)
! write neutron record
!------------------------------------------------------------------------
      TYPE(NEUTRON) :: NEU
      CHARACTER*(*) :: MSG
1     format(a25,': ',8(G10.4,1x))
      write(*,1) MSG,NEU%R, NEU%K,NEU%P
      END SUBROUTINE TRACING_REP


C---------------------------------------------------------------
      LOGICAL FUNCTION TRACING_GO(IBEAM)
C tracing single event through the given beamline.
C Called by TRACING_RUN, which prepares initial NEUT
C IBEAM ... beamline to trace through:
C primary(1), secondary(2) or both(3)
C---------------------------------------------------------------
      integer,intent(in) :: IBEAM
      LOGICAL LOG1
      integer :: i
      TYPE(NEUTRON) :: NEU0,NEU1,NEU2,NEU3,NEU4
1     format(a12,': ',8(G10.4,1x))
      LOG1=(NEUT%P.GT.0)
      if (.not.log1) then
        TRACING_GO=.false.
        return
      endif
      STOPTIME=.false.

      !if (NINT(trace_tot)<10) then
      !  write(*,1) 'TRACING_GO', IBEAM, NINT(trace_tot), isProbe
      !endif

      !NEUT_DBG = ((SAMOBJ%P_FRAME%COUNT>1000 ).and.(SAMOBJ%P_FRAME%COUNT<1100))

! PRIMARY BEAMLINE
      IF (IBEAM.eq.1) then
        NEUT_GR=(/0.D0,-gravity,0.D0/)
        select case(TROPT%DIR)
! down-stream
        case(tr_downstream)
          TRACING_UP=.false.
          NEUI1=NEUT
          i=0
          do while(LOG1.and.(i.lt.BEAMLINE_NC(IBEAM)))
            i=i+1
            LOG1=PCOM_GO(BEAMLINE(i,IBEAM)%OBJ)
            LOG1=(LOG1.and.(.not.STOPTIME))
          enddo
          NEUI=NEUT
          ! trace through sample container
          if (LOG1.and.(SAMOBJ%ICLS.gt.0)) then
            LOG1=PSAMOBJ_CONTAINER_GO()
            LOG1=(LOG1.and.(.not.STOPTIME))
            NEUI=NEUT
            ! NEUT is in sample local coordinates, but we need NEUI in axis (entry) coodinates
            call FRAME_LOCAL(-1,SAMOBJ%P_FRAME,NEUI)
          endif
! up-stream
        case(tr_upstream)
          NEUI=NEUT
        ! trace through sample container first
          if (LOG1.and.(SAMOBJ%ICLS.gt.0)) then
            TRACING_UP=.false. ! trace container forward
            NEUI1 = NEUT ! remember generated position
            ! container expects neutron in the previous component exit frame
            ! it will do oposite conversion at the entry => no move
            NEUT%R(3)=NEUT%R(3)+SAMOBJ%P_FRAME%DIST
            LOG1=PSAMOBJ_CONTAINER_GO()
            LOG1=(LOG1.and.(.not.STOPTIME))
            ! remember the state inside container, it will be used in 2nd part
            NEU0=NEUT
            ! NEUT is in sample local coordinates, but we need it in entry axis coodinates
            call FRAME_LOCAL(-1,SAMOBJ%P_FRAME,NEUT)
            NEUI=NEUT ! this is saved as incident neutron
            ! restore generated neutron and add sample distance
            NEUT = NEUI1
            ! on backtrace, the 1st component before sample
            ! expects neutron in exit axis coordinates ...
            NEUT%R(3)=NEUT%R(3)+SAMOBJ%P_FRAME%DIST

!debug
            if (NEUT_DBG) NEU1 = NEUT

          else
            NEU0=NEUT
          endif
          ! switch temporarily tracing direction in further components to UPSTREAM
          TRACING_UP=.true.
          i=BEAMLINE_NC(IBEAM)
          do while(LOG1.and.(i.gt.0))
          !  if (NINT(trace_tot)==17) then
          !     GUIDES_dbg = (trim(BEAMLINE(i,IBEAM)%OBJ%ID).eq.'GSW')
          !     write(*,1) trim(BEAMLINE(i,IBEAM)%OBJ%ID), mtmod_ngen, LOG1, NEUT%R, NEUT%K
          !  else
          !1    GUIDES_dbg = .false.
          !  endif
          !  mirror_dbg = GUIDES_dbg
            LOG1=PCOM_GO(BEAMLINE(i,IBEAM)%OBJ)
            LOG1=(LOG1.and.(.not.STOPTIME))
            i=i-1
          enddo
        ! correct event records so that they are equivalent to the down-stream case
          ! NEUI is neutron before scattering. Set the right weight and T0
          NEUI%P=NEUT%P
          NEUI%T0=NEUT%T0
          ! NEUI1 is neutron at the source (actually not used ...)
          NEUI1=NEUT
          NEUI1%P=NEUI%P
          ! jump to the initial state with corrected weight
          NEU0%P=NEUT%P
          NEU0%T0=NEUT%T0
          NEUT=NEU0
        END select
! restore tracing direction flag
        TRACING_UP=.false.
!debug
        if (NEUT_DBG) NEU2 = NEUT

! SECONDARY BEAMLINE
      else IF (IBEAM.eq.2) then
    ! trace through sample (scattering) if exists
        if (SAMOBJ%ICLS.gt.0) then
        ! CONTAINER_GO has already prepared NEUT in local coordinates.
        ! PSAMOBJ_GO assumes local coordinates at the entry
        ! for tr_secondary, the events are generated from KSTACK in entry coordinates !!
        ! we need to convert them to sample local coord. as expected by PSAMOBJ_GO
          if (TROPT%IMODE.eq.tr_secondary) then
            ! set gravity to sample entry coordinats
            call M3XV3(-1,SAMOBJ%P_FRAME%RLAB,(/0.D0,-gravity,0.D0/),NEUT_GR)
            call FRAME_LOCAL(1,SAMOBJ%P_FRAME,NEUT)
            call GRAVITY_LOCAL(1,SAMOBJ%P_FRAME)
          endif
!debug
        if (NEUT_DBG) NEU3 = NEUT

          LOG1=PSAMOBJ_GO()
!debug
        if (NEUT_DBG) NEU4 = NEUT

        if (NEUT_DBG.and.LOG1) then
          call TRACING_REP('TRACING_GO after container',NEU1)
          call TRACING_REP('TRACING_GO after   primary',NEU2)
          call TRACING_REP('TRACING_GO before   sample',NEU3)
          call TRACING_REP('TRACING_GO after    sample',NEU4)
          write(*,*)
        endif

          NEUF=NEUT
      ! convert NEUF to sample incident coordinates (the same as for NEUI)
      ! NOTE: PSAMOBJ_GO leaves neutron in exit axis coordinates
          call FRAME_AXIS(-1,SAMOBJ%P_FRAME,NEUF)
        endif
        i=0
        do while(LOG1.and.(i.lt.BEAMLINE_NC(IBEAM)))
          i=i+1
          LOG1=PCOM_GO(BEAMLINE(i,IBEAM)%OBJ)
          LOG1=(LOG1.and.(.not.STOPTIME))
        !  if (NEUT_DBG) then
        !    write(*,1) trim(BEAMLINE(i,IBEAM)%OBJ%ID),NEUT%R(1:3),NEUT%T/HOVM, NEUT%P, LOG1
        !  endif
        enddo
        NEUF1=NEUT
        NEUF%P=NEUT%P
      else
        LOG1=.false.
      endIF

      TRACING_GO=LOG1
      END FUNCTION TRACING_GO

C----------------------------------------------------------------------
      SUBROUTINE TRACING_OnTrigger
C tasks to be performed when trace_clock changes value
C-----------------------------------------------------------------------
! logical :: CHECK_INTERRUPT
      if (trace_clock.gt.trace_clock_old) then
        trace_clock_old=trace_clock
  ! evaluate preliminary results
        call TRACING_EVAL
  ! update progress bar
        call REPORT_PROG(1)
  ! check timeout
        TIMEOUT=(trace_tot.GT.CNT_TOUT)
c        if (TIMEOUT) write(*,*) 'TIMEOUT on count limit'
        if (.not.TIMEOUT) then
          TIMEOUT=(trace_tot/(trace_ctot+1.D0).gt.TROPT%ELIMIT)
          if (TIMEOUT) write(*,*) 'TIMEOUT on efficiency ',int(trace_ctot),int(trace_tot)
        endif
      !  STOPTRACING=CHECK_INTERRUPT()
c        write(*,*) trace_clock,TRACING_PRECISION
      endif
      end SUBROUTINE TRACING_OnTrigger

C----------------------------------------------------------------------
      SUBROUTINE TRACING_CLEAR
C tasks needed to restart tracing loop, e.g. after updating sampling limits
C-----------------------------------------------------------------------
      if (TROPT%IMODE.ne.tr_secondary) then
        CALL COMPONENTS_CLR(1)
      endif
      if (TROPT%IMODE.ge.tr_secondary) then
        CALL COMPONENTS_CLR(2)
      endif
      if (SAMOBJ%ICLS.gt.0) then
        SAMOBJ%P_FRAME%COUNT=0
      endif
      call ResetEventCounters
      if (TROPT%VARI) CALL COV_CLR
      end SUBROUTINE TRACING_CLEAR

C---------------------------------------------------------------
      SUBROUTINE TRACING_FINAL
C final tasks after tracing has finished
C---------------------------------------------------------------
! call FINALIZE procedure on all components on the beamline
      call COMPONENTS_FINALIZE(1)
      if (SAMOBJ%ICLS.gt.0) call PSAMOBJ_FINALIZE(SAMOBJ)
      if (TROPT%IMODE.ge.tr_secondary) call COMPONENTS_FINALIZE(2)
! report event logs
      call StopEventLog
      call ReportEventLog('TRACING.dat')
      call ReportFilteredLog('TRACING_zx.dat',(/4,2/),100)
      call ReportFilteredLog('TRACING_zy.dat',(/4,3/),100)
      if (TOF_MODE) call ReportFilteredLog('TRACING_zt.dat',(/4,9/),100)

! get results
      call TRACING_EVAL

! switch safety pool back on
      if (TROPT%SWPOOL) call SWPOOL(1)
! report volumes
      call REPORT_VOL('FINAL')
  ! compress by rejection
c      call EVARRAY(4,0,I,AUX,DUM)
! report results
      call REPORT_RES

      end SUBROUTINE TRACING_FINAL

!---------------------------------------------------------------------
      subroutine TRACING_EVAL
! evaluate basic results with stat. errors:
! INT ... intensity
! DE  ... energy spread
! TOF ... mean time-of-flight
!---------------------------------------------------------------------
      integer :: i,N
      REAL(KIND(1.D0)) :: Z,Z2,AUX
      integer :: NS,ND,CT
      type(PDATASET) :: DTS
      TRACING_INT=0.D0
      TRACING_DINT=0.D0
      TRACING_VOL=0.D0
      TRACING_DVOL=0.D0
      TRACING_DE=0.D0
      TRACING_DDE=0.D0
      TRACING_TOF=0.D0
      TRACING_DTOF=0.D0
      TRACING_PRECISION=1.D0

      if ((EVENTGEN%SCOUNT.gt.0.D0).and.(trace_sum.gt.0)) then
  ! intensity
        N=0
        AUX=0.D0
        do i=1,IPARC
        ! include only batches with enough events
          CT=CPARC(i)-CPARC(i-1)
          if (CT.gt.TCOUNT/MPARC/2) then
            N=N+1
            Z2=PARCVOL(i)/(SPARCV(i)-SPARCV(i-1))
            TRACING_VOL=TRACING_VOL+Z2*CT/CPARC(IPARC)
            TRACING_DVOL=TRACING_DVOL+(Z2*CT/CPARC(IPARC))**2
          !  Z=PARCINT(i)/(SPARC(i)-SPARC(i-1))
            Z=PARCINT(i)*Z2
            TRACING_INT=TRACING_INT+Z
            TRACING_DINT=TRACING_DINT+Z**2
            AUX=AUX+PARCINT(i)
        !    if (trace_cnt.eq.TCOUNT) then
        !      write(*,"(a,I4,1x,3(2x,G10.3))") '    PARC VOL ',N,PARCVOL(i),Z2,PARCINT(i)*Z2
        !    endif
          endif
        enddo
        if (N.GT.0) TRACING_INT=TRACING_INT/N
        if (N.GT.0) TRACING_VOL=TRACING_VOL/N
        if (N.GT.5) then
          TRACING_DINT=SQRT((TRACING_DINT/N-TRACING_INT**2)/(N-1))
          TRACING_DVOL=SQRT((TRACING_DVOL/N-TRACING_VOL**2)/(N-1))
		  ! disable precission limit if not set
		  if (TROPT%PLIMIT>1.D-10) then
			 TRACING_PRECISION=TRACING_DINT/TRACING_INT
		  endif
c          write(*,*) IPARC,CPARC(IPARC-1),CPARC(IPARC)-CPARC(IPARC-1),trace_cnt
        endif
    ! renormalize to the flux 10^14/s/cm^2 or flux at the sample
        TRACING_INT=TRACING_INT*FLUX_NORM
        TRACING_DINT=TRACING_DINT*FLUX_NORM

  ! energy spread
        N=0
        do i=1,IPARC
          if (PARCINT(i).gt.0) then
            N=N+1
            Z=PARCDE(i)/PARCINT(i)
            Z2=PARCDDE(i)/PARCINT(i)
            Z=SQRT(Z2-Z**2)
            TRACING_DE=TRACING_DE+Z
            TRACING_DDE=TRACING_DDE+Z**2
          endif
        enddo
        if (N.GT.1) TRACING_DDE=R8LN2*SQRT((TRACING_DDE/N-(TRACING_DE/N)**2)/(N-1))
        if (N.GT.0) TRACING_DE=R8LN2*TRACING_DE/N
        TRACING_DE=TRACING_DE*HSQOV2M
        TRACING_DDE=TRACING_DDE*HSQOV2M
  !  TOF spread
        N=0
        do i=1,IPARC
          if (PARCINT(i).gt.0) then
            N=N+1
            Z=PARCTOF(i)/PARCINT(i)
            Z2=PARCDTOF(i)/PARCINT(i)
            Z=SQRT(Z2-Z**2)
            TRACING_TOF=TRACING_TOF+Z
            TRACING_DTOF=TRACING_DTOF+Z**2
          endif
        enddo
        if (N.GT.1) TRACING_DTOF=R8LN2*SQRT((TRACING_DTOF/N-(TRACING_TOF/N)**2)/(N-1))
        if (N.GT.0) TRACING_TOF=R8LN2*TRACING_TOF/N
  ! store results and trials counter
        select case(TROPT%IMODE)
        case(tr_primary)
          NS=CSUBSET(1)
          ND=DSET_ISEL
        case default
          NS=CSUBSET(2)
          ND=DSET_ISEL
        end select
        MCRESULT(1:2,NS,ND)=(/TRACING_INT,TRACING_DINT/)
        MCRESULT(3:4,NS,ND)=(/TRACING_DE,TRACING_DDE/)
        MCRESULT(5:6,NS,ND)=(/TRACING_TOF,TRACING_DTOF/)
        MCTRIALS(NS,ND)=SPARC(IPARC)
        DTS=DSET_DATA(DSET_ISEL)
        if (TROPT%IMODE.eq.tr_primary) then
          MCVOLI(1:2,NS,ND)=(/TRACING_VOL,TRACING_DVOL/)
          MCNORM(1:2,NS,ND)=FLUX_NORM*MCVOLI(1:2,NS,ND)
          call NSTORE_SETNORM(DTS%P%REG1,MCNORM(1,NS,ND))
          call COMPONENTS_SETMCNORM(1,MCNORM(1,NS,ND))
        else
          MCVOLF(1:2,NS,ND)=(/TRACING_VOL,TRACING_DVOL/)
          MCNORM(1:2,NS,ND)=FLUX_NORM*MCVOLF(1:2,NS,ND)
          if (TROPT%IMODE.eq.tr_all) call NSTORE_SETNORM(DTS%P%REG1,MCNORM(1,NS,ND))
          call NSTORE_SETNORM(DTS%P%REG2,MCNORM(1,NS,ND))
          call COMPONENTS_SETMCNORM(1,MCNORM(1,NS,ND))
          call COMPONENTS_SETMCNORM(2,MCNORM(1,NS,ND))
        endif
      !    if (trace_cnt.eq.TCOUNT) then
      !      write(*,"(a,3(2x,G10.3))") 'TRACING_EVAL ',AUX*TRACING_VOL,MCNORM(1,NS,ND),TRACING_VOL
      !    endif
      endif
      end subroutine TRACING_EVAL

!---------------------------------------------------------------------
      subroutine CountInitEvent
! increment counters of started events
!---------------------------------------------------------------------
      if (isProbe) return
      trace_tot=trace_tot+1
      trace_clock=INT(trace_tot/trigger_count)
      end subroutine CountInitEvent

!---------------------------------------------------------------------
      subroutine CountPassedEvent
! increment counters of successful events
!---------------------------------------------------------------------
      REAL(KIND(1.D0)) :: E,T
      logical :: db=.false.
      type(PDATASET) :: DTS
      logical :: dbg=.false.
1     format(a,': ',6(G10.3,1x),I5)
    ! don't count probe events
       if (isProbe) return
       ! write(*,*) 'CountPassedEvent ',icycle,trace_cnt
    ! counters
        trace_ctot=trace_ctot+1
        trace_crun=trace_crun+1
        trace_cnt=trace_cnt+1
        trace_sum=trace_sum+NEUT%P
        if (icycle.eq.1) then
          CALL COV_UPD(1.D0)
          return
        endif
        DTS=DSET_DATA(DSET_ISEL)
    ! must set cnt to neutron before calling NSTORE_SETEVENT !!
    ! event storage
! PRIMARY
        if (TROPT%IMODE.eq.tr_primary) then
          call COMPONENTS_ENDTRACE(1,NEUI)
          if (SAMOBJ%ICLS>0) call PSAMOBJ_ENDTRACE(1,NEUI)
        endif
! SECONDARY
        if (TROPT%IMODE.eq.tr_secondary) then
          call COMPONENTS_ENDTRACE(2,NEUF1)
          if (SAMOBJ%ICLS>0) call PSAMOBJ_ENDTRACE(2,NEUF)
        endif
! ALL
        if (TROPT%IMODE.eq.tr_all) then
          call COMPONENTS_ENDTRACE(1,NEUI)
          call NSTORE_SETEVENT(DTS%P%REG1,NEUI)
          NEUF%REF=(/DTS%P%REG1,NEUI%CNT/)
          NEUF1%REF=NEUF%REF
          call COMPONENTS_ENDTRACE(2,NEUF1)
          if (SAMOBJ%ICLS>0) then
            call PSAMOBJ_ENDTRACE(1,NEUI)
            call PSAMOBJ_ENDTRACE(2,NEUF)
          endif
        endif
    ! partial sums
        IPARC=MIN(MPARC,INT(1.0*MPARC*trace_cnt/float(TCOUNT))+1)
    !switch safety-pool off for the 2nd half of tracing job
        if (IPARC.gt.MPARC/2) then
          call SWPOOL(0)
        endif
        CPARC(IPARC)=trace_cnt
  ! NOTE: NEUT%P is weighted by EVENTGEN%VOL
        SPARC(IPARC)=EVENTGEN%OLDSCOUNT + (EVENTGEN%SCOUNT-EVENTGEN%OLDSCOUNT)*EVENTGEN%VOL/RNDIST_VOL
        SPARCV(IPARC)=EVENTGEN%SCOUNT
        PARCINT(IPARC)=PARCINT(IPARC)+NEUT%P
        PARCVOL(IPARC)=RNDIST_VOL
        E=NEUT%K0**2
        PARCDE(IPARC)=PARCDE(IPARC)+NEUT%P*E
        PARCDDE(IPARC)=PARCDDE(IPARC)+NEUT%P*E**2
      ! PARCTOF evaluates true time spread, i.e. time resolution
        T=NEUT%T-NEUI1%T
        PARCTOF(IPARC)=PARCTOF(IPARC)+NEUT%P*T
        PARCDTOF(IPARC)=PARCDTOF(IPARC)+NEUT%P*T**2
      end subroutine CountPassedEvent

!---------------------------------------------------------------------
! reset counters incl. total events counters
!---------------------------------------------------------------------
      subroutine InitEventCounters
        call ResetEventCounters
        trace_tot=0.D0
        trace_ctot=0.D0
        if ((.not.TWOSTEPS).or.(TROPT%IMODE.ne.tr_secondary)) trace_crun=0.D0
      !  write(*,*) 'InitEventCounters ctot=',trace_ctot,' CRUN=',CRUN
      !  trace_ctot=0.D0
        trace_clock=0
        trace_clock_old=0
      end subroutine InitEventCounters

!---------------------------------------------------------------------
! reset counters (except total events count)
!---------------------------------------------------------------------
      subroutine ResetEventCounters
      integer :: i
        call GENERATOR_CLEAR
        trace_cnt=0
        trace_sum=0.D0
        IPARC=0
      !  NPARC(0)=0
        SPARC(0)=0.D0
        SPARCV(0)=0.D0
        CPARC(0)=0
        do i=1,MPARC
          PARCINT(i)=0.D0
          PARCVOL(i)=0.D0
          PARCDE(i)=0.D0
          PARCDDE(i)=0.D0
          PARCTOF(i)=0.D0
          PARCDTOF(i)=0.D0
      !    NPARC(i)=0
          SPARC(i)=0.D0
          SPARCV(i)=0.D0
          CPARC(i)=0
        enddo
        TRACING_PRECISION=1.D0
    !    call ResetEventLog
      end subroutine ResetEventCounters

      end MODULE TRACING