!//////////////////////////////////////////////////////////////////////
!////  $Id: cmdhandler.f,v 1.40 2019/08/15 15:02:06 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2009, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.40 $
!////     $Date: 2019/08/15 15:02:06 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Main command dispatcher
!////
!////////////////////////////////////////////////////////////////////////
      MODULE CMDHANDLER
      USE CONSTANTS
      USE XMLSIMRES
      use XMLWRITER
      USE COMPONENTS
      use TRACINGOPT
      use EVENTMONITOR
      USE COMPONENTS_IO
      USE FILETOOLS
      use DIALOGS
      use COMMANDS
      use IO
      USE MESSAGES
      use GRFEXEC
      use SIMEXEC
      use MIRROR_TABLE
      use STRAIN_TABLE
      use REPORTS
      implicit none
      SAVE
      PRIVATE

      INTEGER :: NOS=0            ! number of command arguments
      REAL(KIND(1.D0)) :: RET(40)             ! numerical command arguments
      CHARACTER(LEN_LINE) :: RETSTR     ! string command argument

      logical :: isCmdBuffered=.false.
      CHARACTER(LEN_LINE) :: CmdLine
      character(4) :: Space='    '

      integer, PARAMETER :: RES_NCMD=25   	! Number of commands
      CHARACTER*5 RES_NAM(RES_NCMD)         ! command names
      CHARACTER*60 RES_HLP(RES_NCMD)        ! command hints

      public CMD_EXECUTE,CMD_HANDLE,INIT_COMMANDS
      CONTAINS

!---------------------------------------------------------------
      SUBROUTINE CMD_EXEC_COMP(SCOMM)
! component specific commands
! SCOMM should contain: OBJ CMD ARG
! OBJ ... object ID
! CMD ... command string
! ARG ... command arguments
!---------------------------------------------------------------
      use SGUIDES
      CHARACTER*(*) SCOMM
      type(TCLASS) :: OBJ
      integer :: IS1,IS2,IL1,IL2,IS3,IL3
      character(256) :: CMD,ARG
      IL1=0
      IL2=0
      IL3=0
      ARG=''
      call FINDSTRPAR(SCOMM,' ',1,IS1,IL1)
      if (IL1>0) call FINDSTRPAR(SCOMM,' ',2,IS2,IL2)
      if (IL2>0) call FINDSTRPAR(SCOMM,' ',3,IS3,IL3)
      if (IL2<=0) return

      ! write(*,*) 'CMD_EXEC_COMP ',SCOMM(IS1:IS1+IL1-1)
      ! write(*,*) '         CMD= ',SCOMM(IS2:IS2+IL2-1)
      ! write(*,*) '          ID= ',SCOMM(IS3:IS3+IL3-1)

      OBJ=getClassByID(SCOMM(IS1:IS1+IL1-1))
! only components are accepted
      if ((OBJ%ICLS>0).and.(OBJ%TCLS==TCLS_COM)) then
        CMD=SCOMM(IS2:IS2+IL2-1)
        call MKUPCASE(CMD)
        select case(OBJ%P_COM%CCLS)
        case(CCLS_GUIDE,CCLS_SGUIDE)
          select case(trim(CMD))
          case('PROF')
            if (IL3>0) ARG=SCOMM(IS3:IS3+IL3-1)
            call REPORT_GUIDE_PROF(ARG,OBJ%P_COM%CCLS,OBJ%P_COM%INST)
          case('MREC')
            write(*,*) 'CMD_EXEC_COMP REPORT_GUIDE_REC ',trim(SCOMM)
            if (IL3>0) ARG=SCOMM(IS3:IS3+IL3-1)
            call REPORT_GUIDE_REC(ARG,OBJ%P_COM%CCLS,OBJ%P_COM%INST)
          end select
        end select
      endif
      end SUBROUTINE CMD_EXEC_COMP

!---------------------------------------------------------------
      SUBROUTINE CMD_HANDLE(SCOMM)
!  Main command handler for RESTRAX
!  call CMD_HANDLE(' ') as the first command of RESTRAX and then
!  any time you need to reset default menu
!---------------------------------------------------------------
      INCLUDE 'linp.inc'
      CHARACTER*(*) SCOMM
      INTEGER*4 ICOM,NPAR,L,I,LC, IS, IL
      CHARACTER(LEN_LINE) :: LINE,LINPEXECSTR
      LOGICAL :: LOG1,IsInteger
2     FORMAT(a5,' = ',G12.6)
100   FORMAT(1X,'Unknown command')
200   FORMAT(1X,'RESTRAX Error: ',I4,/,1x,A)

      L=len_trim(SCOMM)
! ignore comments and empty lines
      IF ((L.LE.0).OR.(SCOMM(1:1).EQ.'#')) RETURN
      IF (SCOMM(1:L).EQ.'BUFFER') THEN
        isCmdBuffered=.true.
        CmdLine=''
        Space=''
        return
      endif

! flush any pending messages
      call FLUSH_MESSAGES

      if (isCmdBuffered) then
        !write(*,*) 'CMD_HANDLE BUFFER [',trim(CmdLine),']'
        IF (SCOMM(1:L).EQ.'ENDBUFFER') then
          isCmdBuffered=.false.
          !write(*,*) 'CMD_HANDLE ENDBUFFER '
          !CmdLine=''
        else IF (SCOMM(1:L).EQ.'[SP]') then
          Space='[SP]'
          return
        else
          LC=len_trim(CmdLine)
          if (LC<LEN_LINE) then
            if (len_trim(Space)>0) then
              CmdLine=trim(CmdLine)//' '//SCOMM(1:L)
              Space=''
            else
              CmdLine=trim(CmdLine)//SCOMM(1:L)
            endif
            return
          else
            isCmdBuffered=.false.
            !write(*,*) 'CMD_HANDLE ENDBUFFER 2 '
            !CmdLine=''
          endif
        endif
      else
        CmdLine=SCOMM(1:L)
      endif

      ! handle handle handshake by sending response tag to GUI...
      if (trim(CmdLine(1:5)).eq.'HNDSK') then
        !write(*,*) 'CMD_HANDLE handshake'
        call FINDSTRPAR(trim(CmdLine),' ',1,IS,IL)
        call XML_RSXDUMP(smes,' ',1)
        if (IL.gt.0) then
          call XML_TAG(smes,'HNDSK',CmdLine(IS:IS+IL-1),2)
        else
          call XML_TAG(smes,'HNDSK',' ',2)
        endif
        call XML_RSXDUMP(smes,' ',0)
        return
      endif

C send commands to linp
      IF (trim(CmdLine).EQ.'SETLINP') THEN
        CALL LINPSET(RES_NCMD,'SimRes',RES_NAM,RES_HLP)
        CALL LINPSETIO(SINP,SOUT,SMES)
        RETURN
      ENDIF

      RETSTR=' '
      LINE=' '

C// process command string through LINP
      LINE=LINPEXECSTR(trim(CmdLine),ICOM,NPAR)
      !write(*,*) 'CMD_HANDLE ',ICOM, ' ', trim(LINE)
C// get numeric arguments
      CALL GETLINPARG(LINE,RET,40,NOS)
      !write(*,*) 'CMD_HANDLE ',NOS
C// get the whole line as a string argument
      IF (NPAR.GT.0) RETSTR=LINE
      IF (ICOM.LT.0) then ! not a LINP command
      ! try the new handler
        LOG1=CMD_EXECUTE(trim(RETSTR))
        RETURN
      endif

C// process I/O commands
      IF (LINE(1:4).EQ.'LIST') THEN
          I=DLG_PARAM(' ',.false.)
      ELSE IF (LINE(1:4).EQ.'QUIT') THEN
          call RSXSTOP
C// process execution commands
      ELSE IF (ICOM.GT.0) THEN
        !write(*,*) 'CMD_HANDLE: ',trim(RES_NAM(ICOM)), trim(RETSTR)
        SELECT CASE (RES_NAM(ICOM))
          CASE('SET'); I=DLG_PARAM(trim(RETSTR),.true.)
          CASE('CMD'); I=DLG_PARAM('CMD '//trim(RETSTR),.true.)
          CASE('DO'); LOG1=CMD_EXECUTE(RETSTR)
          CASE('EXEC'); call CMD_EXEC_COMP(RETSTR)
          CASE('MRTAB'); call SET_MTABLE(RETSTR)
          CASE('STTAB'); call SET_STABLE(RETSTR)
          CASE('CFG'); ! CALL SETCFG(RETSTR)
      ! XML sends information to GUI via console IO
          CASE('XML')
            if (XMLOUT.gt.0) then
              call SETXML(RETSTR)
            endif
          CASE('CPATH'); call SETRESPATH(RETSTR) !
          CASE('OPATH'); call SETOUTPATH(RETSTR) !
          CASE('IMP'); ! CALL OPEN_RESCAL(RETSTR,IERR)
          CASE('BAT')
            CALL REINP(RETSTR)
            CALL LINPSETIO(SINP,SOUT,SMES)
          CASE('REPLOT'); call DOPLOT(.true.)
          CASE('GRFDE'); CALL SELGRFDEV(RETSTR)
          CASE('SEED')
            if (isInteger(RETSTR,I)) then
              ISEED=I
              CALL RAN1SEED(ISEED)
            endif
          CASE('CNT')
            if (isInteger(RETSTR,I)) then
              TROPT%CNT=min(10000000,max(50,I))
            endif
          ! adjust crystal
          CASE('ADJCR'); call MKUPCASE(RETSTR); call DO_ADJCRYST(RETSTR)
          CASE('REPCR'); call MKUPCASE(RETSTR); call REPCRYST(RETSTR)
          CASE('REPRE'); call REPORT_MIRROR_REF(RETSTR)
          CASE('REPXT'); call MKUPCASE(RETSTR); call REPXTAL(RETSTR)
          CASE('REPBM'); call MKUPCASE(RETSTR); call REPORT_BEAMLINE(RETSTR)
          case('BREF'); call SET_BMONITOR(RETSTR); call NSTORE_XML_LIST(SMES,BMONITOR_REF)
          case('MCPL'); call REPORT_MCPL(RETSTR)
          case('MCPIN'); call LOAD_MCPL(RETSTR)
          case('MCPEX'); call SAVE_MCPL(RETSTR)
          case('NDUMP'); call nstore_dump(RETSTR)
          CASE('EXIT','EXFF'); call RSXSTOP
          CASE DEFAULT
            write(sout,100)
        END SELECT
      ENDIF

! flush any pending messages
      call FLUSH_MESSAGES

      if (GOEND) then
        !write(*,*) 'CMD_HANDLE: GOEND signal'
        call RESEND
      endif

      END SUBROUTINE CMD_HANDLE


!--------------------------------------------------------------
      logical function CMD_EXECUTE(CMDSTR)
!  dispatch to command handlers
!--------------------------------------------------------------
      TYPE(TCMDOBJ) :: OBJ
      TYPE(TCLASS) :: CLS
      character(*),intent(in) :: CMDSTR
      character(MAX_FNAME_LENGTH) :: FN
      character(LEN_ID) :: id,cmd
      character(LEN_LINE) :: arg
      integer :: ierr,IDCOMPARE,i,IS,IL
      logical :: LOG1
      real(kind(1.D0)) :: Z1(3),Z2(3)
      CMD_EXECUTE=.false.
      LOG1=.true.
! empty command => print help
      if (len_trim(CMDSTR).le.0) then
        call COMMANDS_PRN(SOUT)
        return
      endif


      i=1
      cmd=''
      arg=''
! DO is stripped off from the command
      if ((len_trim(CMDSTR).gt.2).and.(IDCOMPARE('DO',CMDSTR(1:3)).eq.0)) i=3
      call FINDSTRPAR(trim(CMDSTR(i:)),' ',1,IS,IL)
      cmd=CMDSTR(i+IS-1:i+IS+IL-2)
      OBJ=getCmdObjByID(cmd)
      call FINDSTRPAR(trim(CMDSTR(i:)),' ',2,IS,IL)
      if (IL>0) arg=CMDSTR(i+IS-1:i+IS+IL-2)
      CALL GETLINPARG(arg,RET,40,NOS)
! dispatch to handlers
      if (OBJ%ICLS.gt.0) then
        select case(trim(OBJ%P_CMD%ID))
        CASE('LOAD')
          call getCmdParamS('LOAD',1,'FILE',FN)
          call load_config(trim(FN),ierr)
        CASE('SAVE')
          call getCmdParamS('SAVE',1,'FILE',FN)
          checkOverwrite=(getCmdParam('SAVE','OVER').eq.0)
          call save_config(trim(FN))
        CASE('CMDSAVE')
          call getCmdParamS('CMDSAVE',1,'FILE',FN)
          checkOverwrite=(getCmdParam('CMDSAVE','OVER').eq.0)
          call save_commands(trim(FN))
        CASE('CSAVE')
          ID=' '
          call getCmdParamS('CSAVE',1,'FILE',FN)
          call getCmdParamS('CSAVE',1,'ID',ID)
          checkOverwrite=(getCmdParam('CSAVE','OVER').eq.0)
          if (len_trim(ID).gt.0) then
            call save_component(trim(FN),trim(id))
          else
            call MSG_WARN('Unrecognized component ID ['//trim(ID)//']',1)
          endif
        CASE('CLOAD')
          call getCmdParamS('CLOAD',1,'FILE',FN)
          call load_config(trim(FN),ierr)
        CASE('GRSAVE')
          call getCmdParamS('GRSAVE',1,'FILE',FN)
          checkOverwrite=(getCmdParam('GRSAVE','OVER').eq.0)
          call SaveAllPorts(trim(FN))
          call SaveAuxData(trim(FN))
        CASE('MC')
          call TRACING_PROC
        case('FMERIT')
          if (NOS>0) then
            call DO_FMERIT(NINT(RET(1)),Z1,Z2)
          else
            call DO_FMERIT(-1,Z1,Z2)
          endif
      !  case('SWARM')
      !    CLS=PCLASS(OBJ)
      !    info%lun=SOUT
      !    call xml_open(info,'STDIN', .false.)
      !    call OBJ2XML(info,CLS,.true.)
      !    call xml_close(info)
        CASE('PLOT')
          call DOPLOT(.false.)
        CASE('REPLOT')
          call DOPLOT(.true.)
        case('BEAM2D')
          call setCmdParam('PLOT','TYPE',1.D0*pl_BEAM2D)
          call DOPLOT(.false.)
        case('DET2D')
          call setCmdParam('PLOT','TYPE',1.D0*pl_DET2D)
          call DOPLOT(.false.)
        case('BEAM1D')
          call setCmdParam('PLOT','TYPE',1.D0*pl_BEAM1D)
          call DOPLOT(.false.)
        case('DET1D')
          call setCmdParam('PLOT','TYPE',1.D0*pl_DET1D)
          call DOPLOT(.false.)
        CASE('SCAN1D')
          call DOSCAN1D
          call setCmdParam('PLOT','TYPE',1.D0*pl_SCAN1D)
          call DOPLOT(.false.)
        CASE('SCAN2D')
          call DOSCAN2D
          call setCmdParam('PLOT','TYPE',1.D0*pl_SCAN2D)
          call DOPLOT(.false.)
        CASE('RES1D')
          call setCmdParam('PLOT','TYPE',1.D0*pl_RES1D)
          call DOPLOT(.false.)
        case('RES2D')
          call setCmdParam('PLOT','TYPE',1.D0*pl_RES2D)
          call DOPLOT(.false.)
        case('GAUGE2D')
          call setCmdParam('PLOT','TYPE',1.D0*pl_GAUGE2D)
          call DOPLOT(.false.)
        case default
          LOG1=.false.
        end select
      else
        call MSG_WARN('Undefined command: '//trim(CMDSTR),1)
        LOG1=.false.
      endif
      call FLUSH_MESSAGES
      CMD_EXECUTE=LOG1
      end function CMD_EXECUTE


!-------------------------------------------------------------------------
      subroutine INIT_COMMANDS
!-------------------------------------------------------------------------
      integer :: i
      i = 0
      i = i+1
      RES_NAM(i)='CPATH'
      RES_HLP(i)='[path] set search path for configuration files'
      i = i+1
      RES_NAM(i)='OPATH'
      RES_HLP(i)='[path] set path for output'
      i = i+1
      RES_NAM(i)= 'SET'
      RES_HLP(i)='[component ID] [variable ID] [value] set parameter'
      i = i+1
      RES_NAM(i)='CMD'
      RES_HLP(i)='[command ID] set command parameters'
      i = i+1
      RES_NAM(i)='DO'
      RES_HLP(i)='[command] execute command'
      i = i+1
      RES_NAM(i)='ADJCR'
      RES_HLP(i)='[ID] [theta] adjust crystal for given wavelength'
      i = i+1
      RES_NAM(i)='REPCR'
      RES_HLP(i)='[ID] report CRYSTAL reflecting properties'
      i = i+1
      RES_NAM(i)='REPXT'
      RES_HLP(i)='[ID] report XTAL reflecting properties'
      i = i+1
      RES_NAM(i)='REPRE'
      RES_HLP(i)='fname mval angle lambda; report mirror reflectivity'
      i = i+1
      RES_NAM(i)='REPBM'
      RES_HLP(i)='fname write a table with all components coordinates '
      i = i+1
      RES_NAM(i)='XML'
      RES_HLP(i)='load XML definition from console (GUI only)'
      i = i+1
      RES_NAM(i)='SEED'
      RES_HLP(i)='[int] set seed for random number generator'
      i = i+1
      RES_NAM(i)='CNT'
      RES_HLP(i)='[int] set number of requested counts'
      i = i+1
      RES_NAM(i)='REPLOT'
      RES_HLP(i)='replot the latest graph'
      i = i+1
      RES_NAM(i)='GRFDE'
      RES_HLP(i)='set output device name for PGPLOT'
      i = i+1
      RES_NAM(i)='BAT'
      RES_HLP(i)='[jobname] execute a batch file'
      i = i+1
      RES_NAM(i)='BREF'
      RES_HLP(i)='[index] select repository for BEAM1D and 2D plots'
      i = i+1
      RES_NAM(i)='MCPL'
      RES_HLP(i)='[fname][coord, filt] dump current registry to MCPL file'
      i = i+1
      RES_NAM(i)='MCPIN'
      RES_HLP(i)='[fname][0..2] load MCPL file to given registry'
      i = i+1
      RES_NAM(i)='MCPEX'
      RES_HLP(i)='[fname][0..2] export MCPL file from given registry'
      i = i+1
      RES_NAM(i)='NDUMP'
      RES_HLP(i)='[nr] dump given neutron registry to a file'
      i = i+1
      RES_NAM(i)='MRTAB'
      RES_HLP(i)='[m file] read mirror reflectivity table'
      i = i+1
      RES_NAM(i)='STTAB'
      RES_HLP(i)='[id file] read table with strain profile'
      i = i+1
      RES_NAM(i)='EXEC'
      RES_HLP(i)='[CID] [ARG] execute component specific command'
      i = i+1
      RES_NAM(i)='EXFF'
      RES_HLP(i)='exit'

      CALL LINPSET(RES_NCMD,'SimRes',RES_NAM,RES_HLP)
      CALL LINPSETIO(SINP,SOUT,SMES)

      end subroutine INIT_COMMANDS

      end MODULE CMDHANDLER