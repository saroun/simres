!//////////////////////////////////////////////////////////////////////
!////  $Id: reports.f,v 1.36 2019/08/15 15:02:08 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2008, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.36 $
!////     $Date: 2019/08/15 15:02:08 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Subrouting for reporting on ray-tracing progress
!////  (statistics, event flow, progress bar etc.)
!////
!//////////////////////////////////////////////////////////////////////
      MODULE REPORTS
      USE TRACINGDATA
      use REPORTSOPT
      use COMPONENTS
      use GENERATOR
      USE XMLINFO
      use RESULTS
      use FILETOOLS
      use EVENTMONITOR
      USE RESMAT,only: MATLEG
      use MIRRLOG
      implicit none
      SAVE
      private

! CPU times (start, run, estimated)
      REAL :: TSTART            ! start
      REAL :: TRUN              ! current
      REAL :: TRUN0,CRUN,CRUN0  ! auxilliary time and count records needed for estimations
      REAL :: TESTIM            ! estimation
      REAL :: TCOUNT_ESTIM      ! target count used for run time estimation

      integer, parameter :: TCOUNT_SPEED_MAX=10 ! maximum number of values involved in average speed measurement
      integer :: TCOUNT_SPEED_NUM               ! actual number of values involved in average speed measurement
      real :: TCOUNT_SPEED(TCOUNT_SPEED_MAX)    ! measured average count speed

      LOGICAL :: PROG_ON=.false.

      PUBLIC REPORT_PROG,REPORT_VARI,REPORT_VOL,REPORT_STAT,REPORT_RES,TCOUNT_ESTIM,REPORT_GUIDE_PROF
      public REPORT_GUIDE_REC,REPORT_MCPL,LOAD_MCPL,SAVE_MCPL
      public REPORT_MIRROR_REF, REPORT_BEAMLINE, REPORT_SAMPLE


      contains


!-------------------------------------------------------------
      subroutine REPORT_MCPL(SARG)
! Dump content of the curretly selected monitor to MCPL file
! without a path sepcified, the output is in the project output directory
! optional arguments are
! int, coordinate system
! int, filter flag
!-------------------------------------------------------------
      character*(*) :: SARG
      integer :: is,il,c,FILT
      character*256 :: fout,fdir,fn,fname
      character*32 :: fe
      real(kind(1.D0)) :: ARG(3)
      integer :: NARG
      if (len_trim(SARG)>0) then
        !call FINDSTRPAR(SARG,' ',1,IS,IL)
        IS=1
        call FINDPAR(SARG,1,IS,IL)
        fname=SARG(IS:IS+IL-1)
        CALL GETLINPARG(SARG(IS+IL+1:),ARG,3,NARG)
        call SplitPathName(fname,fdir,fn,fe)
        if (len_trim(fdir).eq.0) fdir=trim(OUTPATH)
        if (len_trim(fe).eq.0) fe='mcpl'
        fout=trim(fdir)//'/'//trim(fn)//'.'//trim(fe)
        if (NARG>0) then
          c=NINT(ARG(1))
        else
          c=-1
        endif
        if (NARG>1) then
          FILT=NINT(ARG(2))
        else
          FILT=-1
        endif
        call BMONITORS_DUMPMCPL(fout,c,FILT,0)
      endif
      end subroutine REPORT_MCPL

!-------------------------------------------------------------
      subroutine SAVE_MCPL(SARG)
! Dump content of the specified monitor to MCPL file
! without a path sepcified, the output is in the project output directory
! optional arguments are
! int, storage definition: (1) BMONITOR_PRI, (2) BMONITOR_SEC, else BMONITOR_REF
!-------------------------------------------------------------
      character*(*) :: SARG
      integer :: is,il,s
      character*256 :: fout,fdir,fn,fname
      character*32 :: fe
      real(kind(1.D0)) :: ARG(3)
      integer :: NARG
      if (len_trim(SARG)>0) then
        !call FINDSTRPAR(SARG,' ',1,IS,IL)
        IS=1
        call FINDPAR(SARG,1,IS,IL)
        fname=SARG(IS:IS+IL-1)
        CALL GETLINPARG(SARG(IS+IL+1:),ARG,3,NARG)
        call SplitPathName(fname,fdir,fn,fe)
        if (len_trim(fdir).eq.0) fdir=trim(OUTPATH)
        if (len_trim(fe).eq.0) fe='mcpl'
        fout=trim(fdir)//'/'//trim(fn)//'.'//trim(fe)
        if (NARG>0) then
          s=NINT(ARG(1))
        else
          s=0
        endif
        call BMONITORS_DUMPMCPL(fout,coord_axis,-1, s)
      endif
      end subroutine SAVE_MCPL

!-------------------------------------------------------------
      subroutine LOAD_MCPL(SARG)
! Dump content of the curretly selected monitor to MCPL file
! without a path sepcified, the output is in the project output directory
! optional arguments are
! int, storage definition: (1) BMONITOR_PRI, (2) BMONITOR_SEC, else BMONITOR_REF
!-------------------------------------------------------------
      character*(*) :: SARG
      integer :: is,il,s
      character*256 :: fout,fdir,fn,fname
      character*32 :: fe
      real(kind(1.D0)) :: ARG(3)
      integer :: NARG
      if (len_trim(SARG)>0) then
        !call FINDSTRPAR(SARG,' ',1,IS,IL)
        IS=1
        call FINDPAR(SARG,1,IS,IL)
        fname=SARG(IS:IS+IL-1)
        CALL GETLINPARG(SARG(IS+IL+1:),ARG,3,NARG)
        call SplitPathName(fname,fdir,fn,fe)
        if (len_trim(fdir).eq.0) fdir=trim(OUTPATH)
        if (len_trim(fe).eq.0) fe='mcpl'
        fout=trim(fdir)//'/'//trim(fn)//'.'//trim(fe)
        if (NARG>0) then
          s=NINT(ARG(1))
        else
          s=0
        endif
        call BMONITORS_READMCPL(fout, s)
      endif
      end subroutine LOAD_MCPL

C------------------------------------------------------------------------
      SUBROUTINE REPORT_PROG(ICOM)
C report on progress
C update count time, progress bars etc.
C------------------------------------------------------------------------
      integer,intent(in) :: ICOM
      REAL :: speed,DT,DC,T1
      integer :: i,j
    !  write(*,*) 'REPORT_PROG ',REPOPT%PROG,ICOM
1     format(a,': ',6(G10.4,1x))
  ! initial and final time stamps
      select case(ICOM)
      case(0)
        TSTART=SECNDS(0.0)
      case(2)
        TRUN=SECNDS(TSTART) ! timestamp 2
      end select
  ! report, only when required
      if (REPOPT%PROG) then
        select case(ICOM)
      ! initialize progress bar
        case(0)
          call XML_NESS(smes,0)
          PROG_ON=.true.
          TRUN=0.0
          CRUN=0.0
          TCOUNT_SPEED_NUM=0
      !    write(*,*) 'REPORT_PROG init'
      ! update progress bar
        case(1)
          if (PROG_ON) then
            T1=SECNDS(TSTART)
            if ((T1-TRUN.gt.REPOPT%TIME).and.(trace_crun-CRUN.gt.10)) then
              TRUN0=TRUN
              TRUN=T1
              CRUN0=CRUN
              CRUN=trace_crun
              TESTIM=0.D0
              DT=TRUN-TRUN0
              DC=CRUN-CRUN0
              if ((DT.gt.0.0).and.(DC.GT.10.0)) then
      ! measure count rate
                if (TCOUNT_SPEED_NUM.ge.TCOUNT_SPEED_MAX) then
                  j=TCOUNT_SPEED_MAX
                  do i=1,TCOUNT_SPEED_MAX-1
                    TCOUNT_SPEED(i)=TCOUNT_SPEED(i+1)
                  enddo
                else
                  TCOUNT_SPEED_NUM=TCOUNT_SPEED_NUM+1
                  j=TCOUNT_SPEED_NUM
                endif
                TCOUNT_SPEED(j)=DC/DT
          ! calculate average count rate
                speed=0.0
                do i=1,TCOUNT_SPEED_NUM
                  speed=speed+TCOUNT_SPEED(i)
                enddo
                speed=speed/TCOUNT_SPEED_NUM
                TESTIM=TRUN+(TCOUNT_ESTIM-trace_cnt)/speed
                call XML_NESS(smes,2)
              endif
c              write(*,1) 'REPORT_PROG',DT,int(DC),speed,PROGRES_UPDATE_CNT,trace_crun,trace_tot
            endif
          endif
      ! finalize progress bar - write a report
        case(2)
          call XML_NESS(smes,3)
          PROG_ON=.false.
        end select
      endif
      end SUBROUTINE REPORT_PROG

C------------------------------------------------------------------------
      SUBROUTINE REPORT_VARI
C report on new transformation matrix
C------------------------------------------------------------------------
      if (REPOPT%MUTE) return
      if (REPOPT%VARI) then
          call XML_RSXDUMP(sout,' ',1)
          call XML_MATRIX(sout,RNDLIST%DIM,RNDLIST%DIM,'TMAT','x|y|kx|ky|kz|t','x|y|kx|ky|kz|t',TMAT(1,1),CRND)
          call XML_ARRAY(sout,RNDLIST%DIM,'LIMIT','x|y|kx|ky|kz|t',RNDLIST%LIMITS(1))
          call XML_RSXDUMP(sout,' ',0)
      endif
      end SUBROUTINE REPORT_VARI

C------------------------------------------------------------------------
      SUBROUTINE REPORT_VOL(CAP)
C report on sampling volume
C------------------------------------------------------------------------
      character*(*),intent(in) :: CAP
      if (REPOPT%MUTE) return
      if (REPOPT%VOL) then
          call XML_RSXDUMP(sout,' ',1)
          call XML_VALUE(sout,trim(CAP)//' VOLUME','float',' ',0,RNDIST_VOL,' ')
          call XML_RSXDUMP(sout,' ',0)
      endif
      end SUBROUTINE REPORT_VOL

C------------------------------------------------------------------------
      SUBROUTINE REPORT_STAT
C report on tracing statistics
C------------------------------------------------------------------------
      if (REPOPT%MUTE) return
      if (REPOPT%STATS) then
        write(sout,*)
        CALL EVENTSTAT
1       format('counts: ',I8,'  total: ',G11.5,'  efficiency: ',G11.5,'  time: ',G10.4)
        if (trace_tot.gt.0.D0) write(sout,1) trace_cnt,trace_tot,trace_cnt/trace_tot,TRUN
        if (TROPT%SWPOOL) then
          call XML_RSXDUMP(sout,' ',1)
          call XML_IARRAY(sout,RNDLIST%DIM,'SAFETY POOL',MATLEG,HIT(1))
          call XML_RSXDUMP(sout,' ',0)
        endif
      !  write(*,*) 'REPORT_STAT RNDLIST%DIM=',RNDLIST%DIM
      endif
      end SUBROUTINE REPORT_STAT

C------------------------------------------------------------------------
      SUBROUTINE REPORT_RES
C report on tracing basic results
C------------------------------------------------------------------------
      if (REPOPT%MUTE) return
      if (REPOPT%RES) then
        call XML_RSXDUMP(sout,' ',1)
        call XML_FVALUE(sout,'INTENSITY',' ',TRACING_INT,TRACING_DINT)
        call XML_FVALUE(sout,'ENERGY_WIDTH','meV',TRACING_DE,TRACING_DDE)
        call XML_FVALUE(sout,'TOF_WIDTH','us',TRACING_TOF,TRACING_DTOF)
        call XML_RSXDUMP(sout,' ',0)
      endif
      end SUBROUTINE REPORT_RES

C--------------------------------------------
      SUBROUTINE EVENTSTAT
C print statistics of events through beamline
C---------------------------------------------
      TYPE(TCOMOBJ) :: OBJ
      integer :: NC,IP,IS,IL,i,di
      real(KIND(1.D0)) :: COUNTS(BEAMLINE_DIM*2+1,3)
      character*2048 HDR
      character*128 STR
      TYPE(PFRAME) :: F
      if (REPOPT%MUTE) return
! no events => return
      if (EVENTGEN%COUNT.le.0) RETURN

! Collect components and their counter values (only if count>0)
      select case(TROPT%DIR)
      case(tr_upstream)
        i=BEAMLINE_NC(1)
        di=-1
      case default
        i=1
        di=1
      end select
  ! PRIMARY
      IP=0
      HDR = 'generator'
      NC=1
      COUNTS(1,1)=trace_tot
      COUNTS(1,2)=1
      COUNTS(1,3)=1
      if (TROPT%IMODE.ne.tr_secondary) IP=1
      DO WHILE ((i.gt.0).and.(i.le.BEAMLINE_NC(1)).and.(IP.gt.0))
        OBJ=BEAMLINE(i,1)%OBJ
        call PCOMOBJ_GET(OBJ,F)
        if (.not.associated(F%X)) EXIT
        STR=trim(F%X%NAME)
        call BOUNDS(STR,IS,IL)
        !write(*,*) 'EVENTSTAT ',trim(STR),' ',F%X%COUNT
        if (F%X%COUNT.GT.0) then
      ! add a component to the list
          NC=NC+1
          HDR=trim(HDR)//'|'//STR(IS:IS+IL-1)
          COUNTS(NC,1)=F%X%COUNT
        else
          IP=0 ! stop on OBJ%COUNT=0
        endif
        i=i+di
      enddo
  ! SPECIMEN
      if (SAMOBJ%ICLS.gt.0) then
        HDR=trim(HDR)//'|SAMPLE'
        NC=NC+1
        COUNTS(NC,1)=SAMOBJ%P_FRAME%COUNT
      endif
  ! SECONDARY
      IP=0
      if (TROPT%IMODE.gt.tr_primary) IP=1
      do i=1,BEAMLINE_NC(2)
        if (IP.GT.0) then
          OBJ=BEAMLINE(i,2)%OBJ
          call PCOMOBJ_GET(OBJ,F)
          if (.not.associated(F%X)) EXIT
          STR=trim(F%X%NAME)
          call BOUNDS(STR,IS,IL)
          if (F%X%COUNT.GT.0) then
      ! add a component to the list
            NC=NC+1
            HDR=trim(HDR)//'|'//STR(IS:IS+IL-1)
            COUNTS(NC,1)=F%X%COUNT
          else
            IP=0 ! stop on OBJ%COUNT=0
          endif
        endif
      enddo

! calculate efficiences
      do i=2,NC
        COUNTS(i,2)=COUNTS(i,1)/COUNTS(i-1,1) ! relative
        COUNTS(i,3)=COUNTS(i,1)/COUNTS(1,1)   ! absolute
      enddo


! print table
      call XML_RSXDUMP(sout,' ',1)
      if (NC.gt.0) call XML_MATRIX(sout,3,NC,'COUNTS','count|eff_rel|eff_abs',trim(HDR),COUNTS(1,1),BEAMLINE_DIM*2+1)
      call XML_RSXDUMP(sout,' ',0)

      END SUBROUTINE EVENTSTAT


C--------------------------------------------------------------------
      SUBROUTINE XML_NESS(IU,STATE)
C Print XML tag: information about MC simulation progress
C INPUT:
C   IU          ... output unit
C   STATE   ... (0) start (1) time estimation (2) progress (3) end
C   TESTIM  ... estimeted time
C   TPASSED ... passed time
C   NCMAX   ... requested events
C   NC      ... passed events
C   NEVENT  ... total attempts
C------------------------------------------------------------------
      INTEGER*4 IU,STATE
      real*8 Z
      CHARACTER*32 CTESTIM,CTPASSED,CEFF,CNC,CNCMAX,CNEVENT,CPREC

10    FORMAT('<SIMULATION nc="',a,'" ncmax="',a,'" >')
11    format('<TIME units="s">',a,'</TIME>')
12    format('<EFFICIENCY units="%">',a,'</EFFICIENCY>')
13    format('<TIMEGUESS units="s">',a,'</TIMEGUESS>')
14    format(a)

1     FORMAT('.',$)
3     FORMAT(' Wait please ',$)
4     FORMAT(' Time to wait [s]: ',a)
5     FORMAT(' Timeout reached, simulation not completed.')
6     FORMAT(' Reached required precision.')
8     FORMAT(' CPU: ',a8,' events: ',a10,' counts: ',a7,' err: ',a6)
9     FORMAT(' CPU: ',a8,' estim: ',a6,' events: ',a10,' counts: ',a7,' err: ',a6)

c        write(*,*) 'XML_NESS ',trace_cnt,TRACING_PRECISION

c format numbers
        call FLOAT2STR(TRACING_PRECISION,CPREC)
        call FLOAT2STR(1.D0*trace_cnt,CNC)
        call FLOAT2STR(1.D0*TCOUNT,CNCMAX)
        call FLOAT2STR(1.D0*EVENTGEN%COUNT,CNEVENT)
        call FLOAT2STR(1.D-1*int(TESTIM*10),CTESTIM)
        call FLOAT2STR(1.D0*TRUN,CTPASSED)
        if (EVENTGEN%COUNT.GT.0) then
          Z=1.D0*trace_cnt/EVENTGEN%COUNT
        else
          Z=0.D0
        endif
        call FLOAT2STR(Z*100,CEFF)

c XML output
      if (XMLOUT.GT.0) then
c initial tag
        call XML_RSXDUMP(IU,'NESS',1)
        write(IU,10) trim(CNC),trim(CNCMAX)
c inner tags
        select case(state)
        case(0) ! start
          call XML_TAG(IU,'STARTSIM',' ',2)
        case(1) ! time estimation
          write(IU,13) trim(CTESTIM)
        case(2) ! progress
          call XML_TAG(IU,'PROGRESS',' ',1)
          write(IU,11) trim(CTPASSED)
          write(IU,12) trim(CEFF)
          write(IU,13) trim(CTESTIM)
          call XML_TAG(IU,'PROGRESS',' ',0)
        case(3) ! end
          call XML_TAG(IU,'ENDSIM',' ',1)
          if (REPOPT%MUTE) then
            write(IU,14) 'mute'
          endif
          write(IU,11) trim(CTPASSED)
          write(IU,12) trim(CEFF)
          call XML_TAG(IU,'ENDSIM',' ',0)
        end select
c end tag
        call XML_TAG(IU,'SIMULATION',' ',0)
        call XML_RSXDUMP(IU,'NESS',0)
c console output
      else
        select case(state)
        case(0) ! start
          WRITE(IU,3)
          WRITE(IU,*) ' ... '
        case(1) ! time estimation
          WRITE(IU,4) trim(CTESTIM)
        case(2) ! progress
          IF (SILENT.LT.1) WRITE(IU,9) trim(CTPASSED),trim(CTESTIM),trim(CNEVENT),trim(CNC),trim(CPREC)
        case(3) ! end
          WRITE(IU,*)
          IF (SILENT.LT.1) WRITE(IU,8) trim(CTPASSED),trim(CNEVENT),trim(CNC),trim(CPREC)
          if (TIMEOUT) WRITE(IU,5)
          if (TRACING_PRECISION.lt.TROPT%PLIMIT) WRITE(IU,6)
        end select
      endif
      END SUBROUTINE XML_NESS



!-------------------------------------------------------------
      subroutine REPORT_GUIDE_REC(FNAME,CCLS,INST)
! write record of the mirror bounces in a file
! of given name (each component adds suffix to the file name
! such as FNAME_ID.dat
!-------------------------------------------------------------
      character*(*) :: FNAME
      integer, intent(in) :: CCLS,INST
      integer :: irec,IO

      irec=MLOGS_GET_REC(CCLS,INST)
      write(*,*) 'REPORT_GUIDE_REC ',trim(FNAME), ' inst=',INST, 'irec=',irec

      if (irec>0) then
        IO=6
        if (len_trim(FNAME)>0) IO=OPENFILEOUT_SFX(fname,MLOGS(irec)%ID)
        if (IO>0) then
          call MLOGS_REC_OUT(IREC,IO)
          if (IO.ne.6) close(IO)
        endif
      endif
      end subroutine REPORT_GUIDE_REC

!-------------------------------------------------------------
      subroutine REPORT_GUIDE_PROF(FNAME,CCLS,INST)
! write guide profiles of the primary part of the instrument
! into files of given name (each component adds suffix to the file name
! such as FNAME_ID.dat
!-------------------------------------------------------------
      character*(*) :: FNAME
      integer, intent(in) :: CCLS,INST
      TYPE(TCOMOBJ) :: OBJ
      integer :: i,IO
      ! write(*,*) 'REPORT_GUIDE_PROF ',trim(FNAME), ' inst=',INST
      do i=1,BEAMLINE_NC(1)
        OBJ=BEAMLINE(i,1)%OBJ
        if ((INST<1).or.(INST==OBJ%INST)) then
        if (CCLS==OBJ%CCLS) then
          ! write(*,*) trim(OBJ%ID), ' obj.inst=',OBJ%INST
          select case(OBJ%CCLS)
          case(CCLS_GUIDE,CCLS_SGUIDE)
            IO=6
            if (len_trim(FNAME)>0) IO=OPENFILEOUT_SFX(fname,OBJ%ID)
            if (IO>0) then
              select case(OBJ%CCLS)
              case(CCLS_GUIDE)
                call GUIDE_INIT(OBJ%INST)
                call GUIDE_PROFILE(OBJ%INST,IO)
              case(CCLS_SGUIDE)
                call SGUIDE_INIT(OBJ%INST)
                call SGUIDE_CALC_PARAM(OBJ%INST)
                call SGUIDE_PROFILE(OBJ%INST,IO)
              end select
              if (IO.ne.6) close(IO)
            endif
          end select
        endif
        endif
      enddo
      end subroutine REPORT_GUIDE_PROF

!-------------------------------------------------------------
      subroutine REPORT_MIRROR_REF(CMD)
! Write reflectivity data for given m-value and angle or wavelength
! Input:
! FNAME = filename without extension. Suffix _mx.x will be added.
! CMD format is "filename mval angle lambda"
! mval = m-value of the supermirror
! angle = incidence angle [rad]
! lambda = wavelength
!-------------------------------------------------------------
      character*(*) :: CMD
      integer :: IO, IS, IL, NA
      character(32) :: CNUM
      REAL(KIND(1.D0)) :: A(3)
1     format(a,': ',7(1x,G12.5))
      call FINDSTRPAR(CMD,' ',1,IS,IL)
      call READ_ARRAY8(CMD(IS+IL:),3,A,NA)
      !write(*,1) 'REPORT_MIRROR_REF: '//CMD(IS:IS+IL-1), NA, A
      if (NA>2) then
        IO=6
        call FLOAT2STR(A(1),CNUM)
        IO=OPENFILEOUT_SFX(CMD(IS:IS+IL-1),'m'//trim(CNUM))
        if (IO>0) then
          call WriteReflectivity(A(1), A(2), A(3), IO)
          if (IO.ne.6) close(IO)
        endif
      endif
      end subroutine REPORT_MIRROR_REF

!-------------------------------------------------------------
      subroutine REPORT_BEAMLINE(FNAME)
! Write a table with all components and their coordinates.
! The table columns:
! ID, x,y,z (world coordinates, mm), distance (mm), tof (us)
! Input:
! FNAME = filename for output.
!-------------------------------------------------------------
      character*(*) :: FNAME
      integer :: IO, ierr
      ! open output file
      IO=6
      if (len_trim(FNAME)>0) IO=OPENFILEOUT_SFX(fname,' ')
      call AdjustBeamline(ierr, IO)
      if (IO.ne.6) then
        close(IO)
      ! info to GUI that the file was saved.
        call MSG_INFO('Beamline log file saved in '//trim(fname),1)
      endif
      end subroutine REPORT_BEAMLINE

!-------------------------------------------------------------
      subroutine REPORT_SAMPLE(CMD)
! Write a report with sample scattering cross-sections in 1/cm
! The table columns:
! wavelength, sig_tot, sig_abs, sig_diff, sig_inc, sig_sph, sig_mph
! Input:
! FNAME = filename for output.
!-------------------------------------------------------------
      character*(*) :: CMD
      integer :: IO, IS, IL, NA
      character*(MAX_FNAME_LENGTH) :: FNAME
      REAL(KIND(1.D0)) :: A(3), lmin, lmax
      integer :: nlam

      call FINDSTRPAR(CMD,' ',1,IS,IL)
      FNAME=CMD(IS:IS+IL-1)
      lmin = 0.5D0
      lmax = 5.5D0
      nlam = 501
      if (IL<len_trim(CMD)) then
        call READ_ARRAY8(CMD(IS+IL:),3,A,NA)
        if (NA>2) then
          lmin = A(1)
          lmax = A(2)
          nlam = nint(A(3))
        endif
      endif
      ! open output file
      IO=6
      if (len_trim(FNAME)>0) IO=OPENFILEOUT_SFX(CMD(IS:IS+IL-1),' ')
      call REPORT_SIGTOT(IO, lmin, lmax, nlam)
      if (IO.ne.6) then
        close(IO)
      ! info to GUI that the file was saved.
        call MSG_INFO('Sample report saved in '//trim(fname),1)
      endif
      end subroutine REPORT_SAMPLE

      end module REPORTS