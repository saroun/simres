!//////////////////////////////////////////////////////////////////////
!////  $Id: io.f,v 1.9 2019/08/15 17:24:08 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2008, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.9 $
!////     $Date: 2019/08/15 17:24:08 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Console input handling. Implements command buffer.
!////
!//////////////////////////////////////////////////////////////////////
      MODULE IO
      implicit none
      SAVE
      PRIVATE

      LOGICAL :: isPaused=.false.       ! pause status
      LOGICAL :: BatchMode=.false.      ! batch mode status
      LOGICAL :: ConsoleEnabled=.true. ! input from STDIN is enabled

  ! input modes
      integer, parameter :: imCONSOLE=0
      integer, parameter :: imBUFFER=1
      integer :: InputMode=imCONSOLE
  ! handler modes
      integer, parameter :: hmCOMMAND=0
      integer, parameter :: hmECHO=1
      integer, parameter :: hmFEED=2
      integer :: HandleMode=hmCOMMAND


  ! command buffer
      integer,PARAMETER :: LEN_LINE=1024    ! default length of input line
      integer,PARAMETER :: MAXSIZE=128   ! maximum number of pages
      integer,PARAMETER :: MAXLINES=16384 ! maximum number of lines
      integer,PARAMETER :: MPAGE=4096 ! size of one memory page

      character, ALLOCATABLE :: SBUFF(:) ! string buffer
      integer :: NALLOC=0 ! number of allocated pages
      integer :: LCOUNT=0 ! number of lines
      integer :: TOP=0 ! top of assigned strings
      integer :: NP(0:MAXLINES) ! partioning of the SBUFF array (line ends)
      integer :: IERR=0


      INTEGER :: SILENT=1     ! silence flag
      INTEGER :: SINP=5       ! input unit
      INTEGER :: SOUT=6       ! output unit
      INTEGER :: SMES=6       ! message unit

      LOGICAL:: GOEND=.false.

  ! published items
      public SMES,SINP,SOUT,SILENT,GOEND,LEN_LINE
      public isPaused,BatchMode,ConsoleEnabled
      public imCONSOLE,imBUFFER,InputMode
      public hmCOMMAND,hmECHO,hmFEED,HandleMode
      public RSXREAD,RSXSTOP,INIT_IO,FREE_IO
      public CHANGE_INPUT,CHANGE_OUTPUT,HLPOUT
      public SBUFF_DELETEALL,SBUFF_ADD,SBUFF_LCOUNT
      public SHOWINPUT, SHOWHANDLE
      contains

!---------------------------------------------------------------
      SUBROUTINE INIT_IO
! initialization of console (buffer allocation etc.)
!---------------------------------------------------------------
      call SBUFF_INIT
      end SUBROUTINE INIT_IO

!---------------------------------------------------------------
      SUBROUTINE FREE_IO
! dispose console (deallocate buffer etc.)
!---------------------------------------------------------------
      call SBUFF_FREE
      end SUBROUTINE FREE_IO

!---------------------------------------------------------------
      SUBROUTINE RSXSTOP
! set flag for program stop
!---------------------------------------------------------------
      GOEND=.true.
      end SUBROUTINE RSXSTOP

!---------------------------------------------------------------
      SUBROUTINE RSXREAD(LINE)
! Read a line from string buffer (SBUFF) or console
! Handle input redirection, PAUSEMODE
!---------------------------------------------------------------
      character*(*) :: LINE
      character(LEN_LINE) :: S
      integer :: LCOUNT,IS,IL

1     format(a)
!2     format('RSXREAD ',a,'|',a,'| count=',I3)

      call SBUFF_LCOUNT(LCOUNT)
      S=' '
c// take a line from buffer . If empty, switch to CONSOLE
      if (InputMode.EQ.imBUFFER) then
        if (LCOUNT.GT.0) then
          call SBUFF_GET(1,S)
          call SBUFF_DELETE(1)
    !      write(*,2) 'buffer ',trim(S),LCOUNT
        !  read(*,*)
        else
          InputMode=imCONSOLE
        endif
C// read a line from SINP or stdin
      else if (InputMode.EQ.imCONSOLE) then
  !        write(*,*) 'console enabled? ',ConsoleEnabled
c// in pause mode, read from stdinp only
        if (ConsoleEnabled) then
          if (isPaused) then
            read(*,1) S
c// batch file => read from SINP
          else if (SINP.NE.5) then
            read(sinp,1,end=20) S
          else if (ConsoleEnabled) then
c// read a line from console
            read(*,1,end=20) S
          endif
        endif
      endif

c// format output and exit
      call BOUNDS(S,IS,IL)
      IL=MIN(IL,LEN(LINE))
      if (IL.GT.0) then
        LINE=S(IS:IS+IL-1)
      else
        LINE=' '
      endif
    !  READ(*,*)
      return
c// on END file, reset IO to standard units
20    CALL CHANGE_INPUT(5)
      CALL LINPSETIO(SINP,SOUT,SMES)
      LINE=' '
      end SUBROUTINE RSXREAD

!--------------------------------------------------------
      SUBROUTINE SBUFF_INIT
!-------------------------------------------------------
      integer :: I
      IERR=0
      IF (ALLOCATED(SBUFF)) then
        DEALLOCATE(SBUFF,STAT=IERR)
      endif
      IF (IERR.NE.0) GOTO 98
      ALLOCATE (SBUFF(1:MPAGE),STAT=IERR) ! allocate one page
      IF (IERR.NE.0) GOTO 99
      NALLOC=1
  ! initialize partion
      do i=0,MAXLINES
        NP(i)=0
      enddo
      TOP=0
      LCOUNT=0
      return
98    write(*,*) 'SBUFF_INIT: deallocation of string buffer failed'
      write(*,*) 'ierr=',ierr
      RETURN
99    write(*,*) 'SBUFF_INIT: allocation of string buffer failed'
      write(*,*) 'ierr=',ierr
      RETURN
      end SUBROUTINE SBUFF_INIT

!--------------------------------------------------------
      SUBROUTINE SBUFF_FREE
!-------------------------------------------------------
      integer :: I
      IERR=0
      IF (ALLOCATED(SBUFF)) DEALLOCATE(SBUFF,STAT=IERR)
      IF (IERR.NE.0) GOTO 98
      do i=0,MAXLINES
        NP(i)=0
      enddo
      TOP=0
      LCOUNT=0
      return
98    WRITE(*,*) 'SBUFF_FREE: deallocation of string buffer failed'
      write(*,*) 'ierr=',ierr
      RETURN
      end SUBROUTINE SBUFF_FREE

!--------------------------------------------------------
      SUBROUTINE SBUFF_LCOUNT(N)
!-------------------------------------------------------
      integer,intent(out) :: N
      IERR=0
      N=LCOUNT
      end SUBROUTINE SBUFF_LCOUNT

!--------------------------------------------------------
      SUBROUTINE SBUFF_ADD(LINE)
!-------------------------------------------------------
      CHARACTER*(*),intent(in) :: LINE
      integer :: IL,I,IS,NAL
      character, ALLOCATABLE :: STMP(:) ! temporary string buffer
      IERR=0
      call bounds(LINE,IS,IL)
      if ((LCOUNT.LT.MAXLINES).and.(MAXSIZE*MPAGE.GT.TOP+IL)) then
c// we need to allocate another page
        IF (NALLOC*MPAGE.LE.TOP+IL) then
          NAL=INT((TOP+IL-NALLOC*MPAGE)/MPAGE)+1 ! number of required pages
          ALLOCATE (STMP(1:NALLOC*MPAGE),STAT=IERR) ! allocate temporary buffer
          IF (IERR.NE.0) GOTO 99
          STMP(1:NALLOC*MPAGE)=SBUFF(1:NALLOC*MPAGE) ! copy SBUFF to STMP
          IF (ALLOCATED(SBUFF)) DEALLOCATE(SBUFF,STAT=IERR)
          IF (IERR.NE.0) GOTO 98
          ALLOCATE (SBUFF(1:(NALLOC+NAL)*MPAGE),STAT=IERR) ! allocate new SBUFF
          IF (IERR.NE.0) GOTO 99
          SBUFF(1:NALLOC*MPAGE)=STMP(1:NALLOC*MPAGE) ! copy STMP back to SBUFF
          IF (ALLOCATED(STMP)) DEALLOCATE(STMP,STAT=IERR) ! deallocate STMP
          IF (IERR.NE.0) GOTO 98
          NALLOC=NALLOC+NAL ! update NALLOC
        endif
c// add additional line
        DO I=1,IL
          SBUFF(TOP+I)=LINE(IS+I-1:IS+I-1)
        enddo
        TOP=TOP+IL
        LCOUNT=LCOUNT+1
        NP(LCOUNT)=NP(LCOUNT-1)+IL
      else
        goto 97 ! buffer overflow
      endif
      return
97    write(*,*) 'SBUFF_ADD: input buffer overflow'
      RETURN
98    write(*,*) 'SBUFF_ADD: deallocation of string buffer failed'
      write(*,*) 'ierr=',ierr
      RETURN
99    write(*,*) 'SBUFF_ADD: allocation of string buffer failed'
      write(*,*) 'ierr=',ierr

      end SUBROUTINE SBUFF_ADD

!--------------------------------------------------------
      SUBROUTINE SBUFF_DELETE(INDX)
!-------------------------------------------------------
      integer,intent(in) :: INDX
      integer :: IL,I
      IERR=0
      if ((INDX.GT.0).AND.(INDX.LE.LCOUNT)) then
        IL=NP(INDX)-NP(INDX-1)
        DO I=NP(INDX-1)+1,NP(LCOUNT)-IL
          SBUFF(I)=SBUFF(I+IL)
        enddo
        do I=INDX,LCOUNT-1
          NP(I)=NP(I+1)-IL
        enddo
        LCOUNT=LCOUNT-1
        TOP=TOP-IL
      endif
      end SUBROUTINE SBUFF_DELETE

!--------------------------------------------------------
      SUBROUTINE SBUFF_DELETEALL
!-------------------------------------------------------
      IERR=0
      DO WHILE (LCOUNT.GT.0)
        call SBUFF_DELETE(LCOUNT)
      enddo
      end SUBROUTINE SBUFF_DELETEALL

!--------------------------------------------------------
      SUBROUTINE SBUFF_GET(INDX,LINE)
!-------------------------------------------------------
      CHARACTER*(*),intent(out) :: LINE
      integer,intent(in) :: INDX
      integer :: IL,IS,I
      IERR=0
      IF ((INDX.GT.0).AND.(INDX.LE.LCOUNT)) then
        IL=NP(INDX)-NP(INDX-1)
        IL=MIN(IL,LEN(LINE))
        IS=NP(INDX-1)+1
        IF (IL.GT.0) then
          DO I=1,IL
            LINE(I:I)=SBUFF(IS+I-1)
          enddo
        else
          LINE=' '
        endif
      else
        LINE=' '
      endif
      end SUBROUTINE SBUFF_GET

!-----------------------------------------------------------------
      SUBROUTINE FILL_BUFFER(LINE)
! Add LINE to the commands buffer (SBUFF, see strings.f)
!----------------------------------------------------------------
      CHARACTER(*) :: LINE
      INTEGER :: IS,IL
      call BOUNDS(LINE,IS,IL)
      IL=MIN(IL,LEN(LINE))
    !  write(*,*) 'FILL_BUFFER ',trim(LINE)
      if (IL.GT.0) then
        call SBUFF_ADD(LINE(IS:IS+IL-1))
      else
        call SBUFF_ADD(' ')
      endif
      call SBUFF_LCOUNT(IS)
      END SUBROUTINE FILL_BUFFER


C-----------------------------------------------------------------
      logical function CHECK_INTERRUPT()
C Processes all input in the command buffer (SBUFF, see strings.f)
C-----------------------------------------------------------------
      logical :: RES
      CHARACTER(LEN_LINE) :: LINE
      INTEGER :: I,N
      RES=.false.
      call SBUFF_LCOUNT(N)
      DO I=1,N
        call SBUFF_GET(I,LINE)
        if (trim(LINE).eq.'INTERRUPT') then
          RES=.true.
          EXIT
        endif
      ENDDO
      CHECK_INTERRUPT=.false.
      END function CHECK_INTERRUPT

C***************************************************************
      SUBROUTINE SHOWHANDLE
C***************************************************************
1     format('Handle=' ,a)
      select case(HandleMode)
      case(hmCOMMAND)
        write(*,1) 'COMMAND'
      case(hmECHO)
        write(*,1) 'ECHO'
      case(hmFEED)
        write(*,1) 'FEED'
      end select
      END SUBROUTINE SHOWHANDLE


C***************************************************************
      SUBROUTINE SHOWINPUT
C***************************************************************
1     format(' Input=' ,a)
      select case(InputMode)
      case(imCONSOLE)
        write(*,1) 'CONSOLE'
      case(imBUFFER)
        write(*,1) 'BUFFER'
      end select
      END SUBROUTINE SHOWINPUT

!-------------------------------------------------------------
      SUBROUTINE CHANGE_INPUT(IU)
! redirection of input
!-------------------------------------------------------------
      integer,intent(in) :: IU
      IF (SINP.NE.5) CLOSE(SINP)
      if (IU.gt.0) THEN
        SINP=IU
        BatchMode=.TRUE.
      else
        SINP=5
        HandleMode=hmCOMMAND
        InputMode=imCONSOLE
        BatchMode=.FALSE.
        call CHANGE_OUTPUT(6)
      endif
      END SUBROUTINE CHANGE_INPUT

!-------------------------------------------------------------
      SUBROUTINE CHANGE_OUTPUT(IU)
! redirection of output
!-------------------------------------------------------------
      integer,intent(in) :: IU
      IF (SOUT.NE.6) then
        call flush(SOUT)
        CLOSE(SOUT)
      endif
      if (IU.le.0) then
        SOUT=6
      else
        SOUT=IU
      endif
      END SUBROUTINE CHANGE_OUTPUT

!-------------------------------------------------------------------------
      subroutine HLPOUT
! print command help
!-------------------------------------------------------------------------
      integer :: J
      CHARACTER(LEN=60) :: HLPOPT(19)  ! hints to command line options
      HLPOPT(1)='RESTRAX options:'
      HLPOPT(2)='-------------------- '
      HLPOPT(3)=' '
      HLPOPT(4)='-d        debug mode (no sampling optimization)'
      HLPOPT(5)='-sx       randomize with x, Example: -s10001'
      HLPOPT(6)='-tx       test random generator with dim=x, e.g. -t5 '
      HLPOPT(7)='-flxfile  read flux distribution dPhi/dLambda from file'
      HLPOPT(8)='            Example: -flxdist.dat '
      HLPOPT(9)='-flhx     set horizontal source divergence limit to x [deg] '
      HLPOPT(10)='            Example: -flh1.5 '
      HLPOPT(11)='-flvx     set vertical source divergence limit to x [deg]'
      HLPOPT(12)='            Example: -flv2.3'
      HLPOPT(13)='-Voigt    use pseudo-Voigt mosaic distribution'
      HLPOPT(14)='-ran1     use Numerical Recipes RAN1 generator'
      HLPOPT(15)='-rand     use system random number generator (RAN)'
      HLPOPT(16)='-noopt    no automatic sampling optimization'
      HLPOPT(17)='-nmon     weight by monitor efficiency ~ 1/ki'
      HLPOPT(18)='-help     show this help'
      HLPOPT(19)='  '
      DO J=1,18
          write(*,*) HLPOPT(J)
      ENDDO
      call RSXSTOP
      end subroutine HLPOUT





      end MODULE IO
