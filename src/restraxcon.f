C---------------------------------------------
C   SIMRES console interface for DLL export
C $Id: restraxcon.f,v 1.39 2019/06/12 17:55:13 saroun Exp $

C***************************************************************
      INTEGER FUNCTION RESTRAX_BATCH(CMD)
C In batch mode, return one line from a batch file in CMD
C Except of PAUSEMODE => must wait for console input
C for DLL export only
C return number of entered characters (ignore empty lines)
C***************************************************************
      use FILETOOLS
      IMPLICIT NONE
      INCLUDE 'linp.inc'
      character*(*) CMD
      character*1024 LINE
      integer*4 IS,IL

1     format(a)
      if (BatchMode.AND.(.NOT.isPaused)) then
        call RSXREAD(LINE)
c        read(sinp,1,end=20) LINE
        call BOUNDS(LINE,IS,IL)
        IL=MIN(IL,LEN(CMD)-1)
        if (IL.GT.0) then
          CMD=LINE(IS:IS+IL-1)//char(0)
          RESTRAX_BATCH=IL
          return
        endif
      endif
      CMD=char(0)
      RESTRAX_BATCH=0
      return

c// on END file, reset IO to standrd units
20    CALL REINP(' ')
      CALL LINPSETIO(SINP,SOUT,SMES)
      CMD=char(0)
      RESTRAX_BATCH=0
      end FUNCTION RESTRAX_BATCH

C***************************************************************
      SUBROUTINE RESTRAX_PROMPT(PSTR)
C  Return prompt string in PSTR including char(0) for DLL export
C***************************************************************
      use XMLINFO
      IMPLICIT NONE
      INCLUDE 'linp.inc'
      character*(*) PSTR

      if ((XMLOUT.LE.0).AND.(linp_in.EQ.5)) then
        PSTR=linp_p(1:linp_np)//'> '//char(0)
      else
        PSTR=' '//char(0)
      endif

      end SUBROUTINE RESTRAX_PROMPT


!***************************************************************
      INTEGER FUNCTION RESTRAX_HANDLE(LINE)
! Dispatches a single command line
! Handles ECHO and PAUSE modes and sends commands to CMD_HANDLE
! Returns:
! 0 .. continue
! 1 .. end session
!***************************************************************
      USE FILETOOLS
      USE CMDHANDLER
      use XMLINFO
      USE IO
      IMPLICIT NONE
      INCLUDE 'linp.inc'

      CHARACTER*(*) LINE
      CHARACTER*1024 S,UPCASE
      character*64 key
      INTEGER*4 IS,IL

1     FORMAT(a)
3     FORMAT(/,'press ENTER to continue or Q to quit...')
      !write(*,*) 'RESTRAX_HANDLE: |',trim(LINE),'|',InputMode,HandleMode, GOEND, BatchMode
C run
c      dbgref=.TRUE.
      IF (.not.GOEND) THEN
        CALL BOUNDS(LINE,IS,IL)
        S=UPCASE(LINE(IS:IS+IL-1))
C RSXFEED, switch handler to FEED mode
        if (S(1:IL).EQ.'<RSXFEED>') THEN
          HandleMode=hmFEED
c          InputMode=imCONSOLE
c          isPaused=.false.
          goto 90
        endif

C Handle input according to the InputMode value
        select case (HandleMode)
C ECHO mode => copy input to output until S='END'
          case(hmECHO)
            IF (S(1:IL).EQ.'END') THEN ! END ECHO
      !        write(*,*) 'ECHO OFF'
              HandleMode=hmCOMMAND
              if (XMLOUT.GT.0) call XML_RSXDUMP(SOUT,'ECHO',0)
            ELSE
              if (XMLOUT.GT.0) then ! parse entities
                call PARSE_ENTITY(LINE)
                call BOUNDS(LINE,IS,IL)
              endif
              WRITE(linp_out,1) LINE(IS:IS+IL-1)
            ENDIF
            goto 90
C COMMAND mode => handle input as commands
          case(hmCOMMAND)
C handle pause
            IF (isPaused) then
              IF (S(1:IL).EQ.'Q') then
                linp_eof=1 ! skip a rest of input file other than SINP
              endif
              isPaused=.FALSE.! end PAUSEMODE on any enter
! wait signal
            ELSE IF (index(S(1:IL),'WAIT=').eq.1) THEN
              if (IL>5) then
                key =trim(S(6:IL))
                call XML_STATUS(SMES, "WAIT", trim(key), " ")
              else
                call XML_STATUS(SMES, "WAIT", " ", " ")
              endif
! notify signal
            ELSE IF (index(S(1:IL),'NOTIFY=').eq.1) THEN
              if (IL>7) then
                key =trim(S(8:IL))
                call XML_STATUS(SMES, "NOTIFY", trim(key), " ")
              else
                call XML_STATUS(SMES, "NOTIFY", " ", " ")
              endif
C echo
            ELSE IF (S(1:IL).EQ.'ECHO') THEN
              HandleMode=hmECHO
      !        write(*,*) 'ECHO ON'
              if (XMLOUT.GT.0) call XML_RSXDUMP(SOUT,'ECHO',1)
C pause
            ELSE IF (S(1:IL).EQ.'PAUSE') THEN
              isPaused=.TRUE.
              if (XMLOUT.GT.0) then
                call XML_RSXDUMP(SOUT,' ',1)
                call XML_TAG(SOUT,'PAUSE',' ',2)
                call XML_RSXDUMP(SOUT,' ',0)
              else
                write(*,3)
              endif
C handle input except empty lines, # comments and </RSXFEED>
            ELSE IF ((IL.GT.0).AND.(LINE(1:1).NE.'#').AND.(S(1:IL).NE.'</RSXFEED>')) THEN
              !if (index(line,'TIMING')>0) write(*,*) 'RESTRAX_HANDLE: [',LINE(IS:IS+IL-1),']'
              CALL CMD_HANDLE(LINE(IS:IS+IL-1))
              IF (SOUT.NE.6) call flush(SOUT)
            ENDIF
C FEED mode => add input to the command buffer until </RSXFEED>
          case(hmFEED)
            if (S(1:IL).EQ.'</RSXFEED>') THEN
              InputMode=imBUFFER
              HandleMode=hmCOMMAND
            else
              call SBUFF_ADD(LINE(IS:IS+IL-1))
              !if (index(line,'TIMING')>0) write(*,*) '   buffered: [',LINE(IS:IS+IL-1),']'
            endif
        end select
      endif
90    if (GOEND) then
        RESTRAX_HANDLE=1
      else
        RESTRAX_HANDLE=0
      endif

C// debug info
      if (.FALSE.) then
        call SHOWINPUT
        call SHOWHANDLE
        !write(*,*) 'GOEND=',GOEND, 'linp_eof=',linp_eof
      endif

! handle requests on I/O reset to STDIN/STDOUT
      if (GOEND.OR.(linp_eof.NE.0)) then
        CALL REINP(' ')
        CALL LINPSETIO(SINP,SOUT,SMES)
      endif
      END FUNCTION RESTRAX_HANDLE

!-----------------------------------------------------------------
      INTEGER FUNCTION PROCESS_BUFFER()
! Processes all input in the command buffer (SBUFF, see strings.f)
! Returns the result of RESTRAX_HANDLE , which is either
! 0 .. continue
! 1 .. end session
!-----------------------------------------------------------------
      USE XMLSIMRES
      use IO
      USE MESSAGES
      use XMLINFO
      IMPLICIT NONE

      CHARACTER*1024 LINE
      INTEGER*4 IRES,RESTRAX_HANDLE,I
1     format(/,'ProcessBuffer ',a,2x,I3)
C run
      IRES=0
      call SBUFF_LCOUNT(I)
      !write(*,*) 'ProcessBuffer count: ',I
      DO WHILE ((I.GT.0).AND.(IRES.EQ.0).AND.(InputMode.EQ.imBUFFER))
    !    write(*,1) 'START:',I
        call RSXREAD(LINE)
        !write(*,*) 'ProcessBuffer [',trim(LINE),']',I
    ! process XML command here: read instr. definition from buffer
        !if (trim(LINE).eq.'XML') then
        !  CALL SETXML
        !  call SBUFF_DELETEALL
        !  call FLUSH_MESSAGES
        !else
          IRES=RESTRAX_HANDLE(LINE)
        !endif
        call SBUFF_LCOUNT(I)
        !if (I.EQ.0) call XML_STATUS(SMES, "READY", " ")
        if (IRES.EQ.1) call XML_STATUS(SMES, "CLOSING", " ", "received GOEND signal")
      ENDDO
      InputMode=imCONSOLE
      PROCESS_BUFFER=IRES
      END FUNCTION PROCESS_BUFFER


