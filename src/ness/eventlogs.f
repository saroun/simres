!//////////////////////////////////////////////////////////////////////
!////  $Id: eventlogs.f,v 1.12 2012/04/26 16:40:40 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2008, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.12 $
!////     $Date: 2012/04/26 16:40:40 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Module for logging neutron histories during ray-tracing
!////
!//////////////////////////////////////////////////////////////////////
      MODULE EVENTLOGS
      USE TRACINGDATA
      use REPORTSOPT
      use FRAMES
      USE FILETOOLS
      implicit none

c arrays for event logs
      integer, parameter :: MAXLOG=1024
      integer, parameter :: MAXRAY=512
      integer :: LOGINDX                   ! current record index
      integer :: NLOGS                    ! maximum allocated record
      integer :: RAYINDX                  ! current ray index
      integer :: NRAYS                    ! maximum allocated ray
      integer :: IRAY(MAXRAY)             ! indices of records where new rays start
      LOGICAL :: LOGGING                  ! current flag for logging
      real(KIND(1.D0)) :: EVLOG(10,MAXLOG) ! log records
      character*10 :: EVLOGOBJ(MAXLOG)    ! source names for each record
      logical :: DBGLOG=.false.            ! flag for debug output


      contains


!---------------------------------------------------------------------
! Initialize event logging
!---------------------------------------------------------------------
      subroutine ResetEventLog
1     format(a,': ',6(G10.4,1x))
        NLOGS=0
        LOGINDX=0
        NRAYS=0
        RAYINDX=0
        if (DBGLOG) then
          write(*,1) 'ResetEventLog ',RAYINDX,NRAYS,LOGINDX,NLOGS,trace_cnt
         ! read(*,*)
        endif
      end subroutine ResetEventLog

!---------------------------------------------------------------------
! Dispose first 50% of rays and continue logging
!---------------------------------------------------------------------
      subroutine ShiftEventLog
      integer :: ir,ilog,ilog0,i
1     format(a,': ',6(G10.4,1x))
      if (LOGGING) then
c        DBGLOG=.true.
        if (RAYINDX.GT.10) then
          ir=RAYINDX/2
          ilog=IRAY(ir) ! the first log to be preserved
        if (DBGLOG) write(*,1) 'ShiftEventLog ir,ilog',ir,NRAYS,ilog,NLOGS
          do i=1,NLOGS-ilog+1
            EVLOG(1:10,i)=EVLOG(1:10,i+ilog-1)
            EVLOGOBJ(i)=EVLOGOBJ(i+ilog-1)
          enddo
          ilog0=IRAY(ir)-1
          do i=1,NRAYS-ir+1
            IRAY(i)=IRAY(i+ir-1)-ilog0
          enddo
          NLOGS=NLOGS-ilog+1
          NRAYS=NRAYS-ir+1
          LOGINDX=NLOGS
          RAYINDX=NRAYS
        if (DBGLOG) then
          write(*,1) 'ShiftEventLog ',RAYINDX,NRAYS,LOGINDX,NLOGS,trace_cnt
          read(*,*)
        endif
        else
          LOGGING=.FALSE.
        endif
      endif
      end subroutine ShiftEventLog

!---------------------------------------------------------------------
! Start logging of a new ray
!---------------------------------------------------------------------
      subroutine StartNewRayLog
1     format(a,': ',6(G10.4,1x))
      if (LOGGING) then
!      DBGLOG=(trace_tot.ge.31350)
        if (RAYINDX.ge.REPOPT%NRAYS) then
!          DBGLOG=.true.
!          write(*,1) 'shift on REPOPT%NRAYS: ',RAYINDX,REPOPT%NRAYS
          call ShiftEventLog
        endif
        if (NLOGS.ge.MAXLOG-50) then
!          DBGLOG=.true.
!          write(*,1) 'shift on MAXLOG: ',NLOGS
          call ShiftEventLog
        endif
        if (LOGGING) then
          RAYINDX=RAYINDX+1
          NRAYS=NRAYS+1
          IRAY(RAYINDX)=LOGINDX+1
        endif
        if (DBGLOG) then
          write(*,1) 'StartNewRayLog',RAYINDX,NRAYS,LOGINDX,NLOGS,trace_cnt
          read(*,*)
        endif
      endif
      end subroutine StartNewRayLog

!---------------------------------------------------------------------
! Start a new ray, discard the old one
!---------------------------------------------------------------------
      subroutine ResetRayLog
1     format(a,': ',6(G10.4,1x))
      if (LOGGING) then
        LOGINDX=max(IRAY(RAYINDX)-1,0)
        NRAYS=max(NRAYS-1,0)
        RAYINDX=max(RAYINDX-1,0)
        NLOGS=LOGINDX
        if (DBGLOG) then
          write(*,1) 'ResetRayLog',RAYINDX,NRAYS,LOGINDX,NLOGS,trace_cnt
c          read(*,*)
        endif
      endif
      end subroutine ResetRayLog

!---------------------------------------------------------------------
! Start event logging
!---------------------------------------------------------------------
      subroutine StartEventLog
      logical :: COND
      if (REPOPT%ELOG.and.(.NOT.LOGGING)) then
    ! log the last 100 events
c        COND=(trace_cnt.eq.TRACING_CNT-99)
        COND=.true.
        LOGGING=COND
c        if (COND) write(*,*) 'StartEventLog '
      endif
      end subroutine StartEventLog

!---------------------------------------------------------------------
! Start event logging
!---------------------------------------------------------------------
      subroutine StopEventLog
      logical :: COND
    ! stop anyway
        COND=.TRUE.
c        if (LOGGING) then
c          write(*,*) 'StopEventLog ', COND, trace_tot,trace_cnt
c          if (COND) call ReportEventLog('StopEventLog.dat')
c          READ(*,*)
c        endif
        LOGGING=(.NOT.COND)
      end subroutine StopEventLog

!---------------------------------------------------------------------
! Add event log
!---------------------------------------------------------------------
      subroutine AddEventLog(OBJ,NEU)
      TYPE(TFRAME) :: OBJ
      TYPE(NEUTRON) :: NEU,NGLOB
      integer :: i
1     format(a10,': ',6(G10.4,1x))
      if (LOGGING) then
        LOGINDX=MIN(LOGINDX+1,MAXLOG)
        call Loc2Glob(OBJ,NEU,NGLOB)
    !      write(*,1) trim(OBJ%NAME),NEU.R(1),NEU.R(3),NEU.K(1),NEU.K(3)
    !      write(*,1) trim(OBJ%NAME),NGLOB.R(1),NGLOB.R(3),NGLOB.K(1),NGLOB.K(3)
    !      read(*,*)
        EVLOG(1,LOGINDX)=trace_tot
        do i=1,3
          EVLOG(i+1,LOGINDX)=NGLOB%R(i)
          EVLOG(i+4,LOGINDX)=NGLOB%K(i)
        enddo
        EVLOG(8,LOGINDX)=NGLOB%P
        EVLOG(9,LOGINDX)=NGLOB%T
        EVLOG(10,LOGINDX)=OBJ%COUNT
        EVLOGOBJ(LOGINDX)=trim(OBJ%ID)
        NLOGS=MAX(NLOGS,LOGINDX)
        if (DBGLOG) then
          write(*,1) trim(EVLOGOBJ(LOGINDX)),RAYINDX,LOGINDX,trace_cnt,int(trace_tot)
      !    read(*,*)
        endif
      endif
      end subroutine AddEventLog

!---------------------------------------------------------------------
      subroutine PrintEvent(IU,OBJ,LABEL)
! Print a single row with event data
! = global neutron coordinates
!---------------------------------------------------------------------
      integer,intent(in) :: IU
      TYPE(TFRAME),intent(in) :: OBJ
      character(*) :: LABEL
      TYPE(NEUTRON) :: NEU
      character(64) :: SR,SK,SC,ST
1     format(a,' r=',a,' k=',a,' t=',a,' p=',G15.9)
      call Loc2Glob(OBJ,NEUT,NEU)
      call ARRAY2STR(NEU%R,3,SR)
      call ARRAY2STR(NEU%K,3,SK)
      call FLOAT2STR(trace_ctot,SC)
      call FLOAT2STR(NEUT%T/HOVM,ST)
      SC=trim(SC)//' '//trim(LABEL)
      write(IU,1) trim(SC),trim(SR),trim(SK),trim(ST),NEUT%P
      end subroutine PrintEvent

!---------------------------------------------------------------------
      subroutine ReportEventLog(FNAME)
! Write log-table to a file
!---------------------------------------------------------------------
      character(*),intent(in) :: FNAME
      integer :: i,IO
1     format(G15.9,1x,8(G12.5,1x),' ',a10,' ',G12.5)
      if (.NOT.REPOPT%ELOG.or.(NLOGS.LE.0)) goto 999
      !Open(Unit=22,File=FNAME,err=999,Status='Unknown')
      IO=OPENFILEOUT_SFX(trim(FNAME),' ')
      if (IO>0) then
        write(IO,*) 'id x y z kx ky kz p t obj cnt'
        do i=1,NLOGS
          write(IO,1) EVLOG(1:9,i),EVLOGOBJ(i),EVLOG(10,i)
        enddo
        if (IO.ne.6) CLOSE(IO)
      endif
999   return
      end subroutine ReportEventLog

!---------------------------------------------------------------------
      subroutine ReportFilteredLog(FNAME,INDX,MAXRAYS)
! Write filtered table with selected projection of trajectories .
! INDX(2) ... selected coordinates (1st is the independent variable)
! MAXRAYS: maximum number of rays reported
! Each ray is represented by a pair of columns (X,Y)
!---------------------------------------------------------------------
      character(*),intent(in) :: FNAME
      integer,intent(in) :: MAXRAYS,INDX(2)
      integer :: i,j,ir,iray0,inode,is,il,is1,il1,maxnode
      character*2048 :: HDR,TABLE(128)
      integer :: LTAB,IO
      character*13 :: CNUM,CHDR(2)
      integer,parameter :: LITEM=26
      character*(LITEM) :: CNUM2,HDRITEM
      character*(LITEM),parameter :: CEMPTY=' -           -            '
      character*20,parameter :: VARNAMES='i:x:y:z:kx:ky:kz:p:t'

2     format(1x,F12.5,1x,F12.5)
4     format(1x,a12,1x,a12)
  ! check input
      if (.NOT.REPOPT%ELOG.or.(NLOGS.LE.0)) goto 999
      if ((INDX(1).LT.1).or.(INDX(1).GT.9).or.(INDX(2).LT.1).or.(INDX(2).GT.9)) then
        goto 999
      endif
  ! init. variables
      ir=0
      inode=0
      iray0=0
      maxnode=0
      LTAB=1
      do i=1,128
        TABLE(i)=' '
      enddo
      HDR=' '
  ! collect output table
      do i=1,NLOGS
      if ((ir.le.NRAYS).AND.(inode.lt.128)) then
      ! new ray
        if ((ir.lt.NRAYS).and.(i.eq.IRAY(ir+1))) then
          ir=ir+1
        ! add empty cells to the rest of rows in the previous ray
          if (ir.gt.1) then
            do j=inode+1,128
              TABLE(j)=TABLE(j)(1:LTAB)//CEMPTY
            enddo
        ! new table row length
            LTAB=LTAB+LITEM
          endif
          inode=0
      ! append header item
          write(CNUM,*) ir
          call BOUNDS(CNUM,IS1,IL1)
          do j=1,2
            call FINDSTRPAR(VARNAMES,':',INDX(j),is,il)
            CHDR(j)=VARNAMES(IS:IS+IL-1)//'_'//CNUM(is1:is1+il1-1)
          enddo
          write(HDRITEM,4) (CHDR(j),j=1,2)
          HDR=HDR(1:LTAB)//HDRITEM
        endif
      ! new node
        inode=inode+1
        maxnode=MAX(inode,maxnode)
        write(CNUM2,2) EVLOG(INDX(1),i),EVLOG(INDX(2),i)
        if (LTAB+LITEM.LE.2048) then
          TABLE(inode)=TABLE(inode)(1:LTAB)//CNUM2
        endif
      endif
      enddo
  ! write table to a file
      IO=OPENFILEOUT_SFX(trim(FNAME),' ')
      if (IO.gt.0) then
        write(IO,*) trim(HDR)
        do i=1,maxnode
          write(IO,*) trim(TABLE(i))
        enddo
        if (IO.ne.6) CLOSE(IO)
      endif
999   return
      end subroutine ReportFilteredLog



      end module EVENTLOGS