!//////////////////////////////////////////////////////////////////////
!////  $Id: commands_2.f,v 1.1 2009/02/16 00:56:31 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.1 $
!////     $Date: 2009/02/16 00:56:31 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Defines classes for commands definition
!////
!//////////////////////////////////////////////////////////////////////
      MODULE COMMANDS
      use XMLINFO
      USE XMLUTILS
      implicit none


! COMMAND definition
      integer,parameter :: CMD_MAXARG=16
      integer,parameter :: CMD_MAX=32
      integer,parameter :: CMD_SLEN=256 ! max. length of string argument
! string field dimensions: keep private
      integer,parameter,PRIVATE :: LEN_ID=16
      integer,parameter,PRIVATE :: LEN_NAME=32
      type TCOMMAND;sequence
        INTEGER :: ICLS ! class index
        INTEGER :: ICMD ! command index
        character(LEN_ID)     :: ID       ! ID string
        character(LEN_NAME)   :: NAME   ! name string
        CHARACTER(64) :: HINT           ! description
        integer :: NARG ! number of numerical arguments
        integer :: NSARG ! number of string arguments
        character(CMD_MAXARG*(LEN_ID+1)) :: IDSTR  ! ID strings for numerical arguments
        character(8*(LEN_ID+1)) :: SIDSTR  ! ID strings for string arguments
        real(kind(1.D0)) :: ARG(CMD_MAXARG)
        INTEGER :: ARG_NUM(CMD_MAXARG) ! positions of the last element of each numeric parameter (for arrays)
        character(CMD_SLEN) :: SARG
      end TYPE TCOMMAND

  ! number of commands defined
      integer :: COMMANDS_NUM
  ! instances of command classes
      type(TCOMMAND),target :: myCOMMANDS(CMD_MAX)

      interface COMMANDS_INP
        module procedure COMMANDS_INP_R
        module procedure COMMANDS_INP_S
      end interface

      interface COMMANDS_OUT
        module procedure COMMANDS_OUT_R
        module procedure COMMANDS_OUT_S
      end interface


! instance of TCOMMAND used at loaing the commands definitions
!      type(TCOMMAND), private :: myCMD

      private COMMANDS_INP_R
      private COMMANDS_INP_S
      private COMMANDS_OUT_R
      private COMMANDS_OUT_S
      contains

C---------------------------------------------------------
      SUBROUTINE COMMANDS_DEFAULT(OBJ)
C clear buffers of the object
C---------------------------------------------------------
      TYPE(TCOMMAND),intent(out) :: OBJ
      OBJ%ICLS=0
      OBJ%ICMD=0
      OBJ%ID=' '
      OBJ%HINT='undefined'
      OBJ%NAME=' '
      OBJ%NARG=0
      OBJ%NSARG=0
      OBJ%ARG=0.D0
      OBJ%ARG_NUM=0
      OBJ%SARG=' '
      OBJ%IDSTR=' '
      OBJ%SIDSTR=' '
      end SUBROUTINE COMMANDS_DEFAULT


!---------------------------------------------------------
      SUBROUTINE COMMANDS_PARINFO(OBJ,PNAME,IDX,TID,NARG)
! return actual class, of which PNAME is a member.
! RETURNS
!   IDX   ... index of the parameter in IDSTR or SIDSTR list
!   TID   ... ID string of variable type (S..string, N..numeric)
!   NARG  ... size of the parameter (>1 for array types)
!---------------------------------------------------------
      TYPE(TCOMMAND),intent(in) :: OBJ
      character(*),intent(in) :: PNAME
      integer,intent(out) :: IDX,NARG
      character(*),intent(out) :: TID
      integer :: I
      call GETPARINDX(OBJ%IDSTR,':',trim(PNAME),IDX)
      if (IDX.gt.0) then
        call GETPARINDX(OBJ%SIDSTR,':',trim(PNAME),I)
        if (I.gt.0) then
          IDX=I
          NARG=1
          TID='S'
        else
          TID='N'
          NARG=OBJ%ARG_NUM(IDX)-OBJ%ARG_NUM(IDX-1)
        endif
      else
        IDX=-1
        NARG=0
        TID=' '
      endif
      end SUBROUTINE COMMANDS_PARINFO

C---------------------------------------------------------
      SUBROUTINE COMMANDS_PRN(IO)
C Print the list of commands
C INPUT
C    IO      .... output unit
C---------------------------------------------------------
      INTEGER,intent(in) :: IO
      integer :: i
1     FORMAT('   ',a10,' ',a)
2     FORMAT(a)
      write(IO,2) 'COMMANDS:'
      do i=1,COMMANDS_NUM
        write(IO,1) adjustl(myCOMMANDS(i)%ID),trim(adjustl(myCOMMANDS(i)%HINT))
      enddo
      END SUBROUTINE COMMANDS_PRN


C---------------------------------------------------------
      subroutine AssignCommand(OBJ)
! assign contents of OBJ to actual command instance
! OBJ%ICMD determines command index
C---------------------------------------------------------
      type(TCOMMAND) :: OBJ
      if ((OBJ%ICMD.gt.0).and.(OBJ%ICMD.le.COMMANDS_NUM)) then
        myCOMMANDS(OBJ%ICMD)=OBJ
      endif
      end subroutine AssignCommand

C---------------------------------------------------------
      type(TCOMMAND) function getCmdByID(IDSTR)
! return command structure for given ID string
! unique substring is handled as an exact match => abbreviations are allowed
! NOTE: only a copy of the command data is returned. To receive a pointer
! to actual command instance, use getCmdObjByID from vctable_commands.mod
C---------------------------------------------------------
      character(*),intent(in) :: IDSTR
      type(TCOMMAND) :: CMD
      integer :: ic
      call COMMANDS_DEFAULT(CMD)
      ic=getCmdIndex(IDSTR)
      if (ic.gt.0) CMD=myCOMMANDS(ic)
      getCmdByID=CMD
      end function getCmdByID

C---------------------------------------------------------
      integer function getCmdIndex(IDSTR)
! return command index for given ID string
! unique substring is handled as an exact match => abbreviations are allowed
C---------------------------------------------------------
      character(*),intent(in) :: IDSTR
      integer :: i,ic,IDCOMPARE,is,icomp
      ic=0
      is=0
      do i=1,COMMANDS_NUM
        icomp=IDCOMPARE(trim(IDSTR),trim(myCOMMANDS(i)%ID))
!        write(*,*) 'getCmdIndex ',trim(IDSTR),' ',trim(myCOMMANDS(i)%ID),icomp
        if (icomp.eq.0) then
          ic=i
          is=1
          EXIT
        else if (icomp.eq.1) then
          is=is+1 ! count substring matches
          ic=i
        endif
      enddo
      if (is.ne.1) ic=0
      getCmdIndex=ic
      end function getCmdIndex

C---------------------------------------------------------
      REAL(KIND(1.D0)) function getCmdParam(CMD,PARID)
! return numerical value of command parameter
C---------------------------------------------------------
      character(*),intent(in) :: CMD
      character(*),intent(in) :: PARID
      integer :: ID,IOBJ
      getCmdParam=-1.D0
      IOBJ=getCmdIndex(CMD)
      if (IOBJ.gt.0) then
        call GETPARINDX(myCOMMANDS(IOBJ)%IDSTR,':',trim(PARID),ID)
        if ((ID.gt.0).and.(ID.le.myCOMMANDS(IOBJ)%NARG)) then
          getCmdParam=myCOMMANDS(IOBJ)%ARG(ID)
        else
          call MSG_WARN('unknown command parameter: '//trim(PARID),1)
        endif
      else
        call MSG_WARN('unknown command ID: '//trim(CMD),1)
      endif
      end function getCmdParam

C---------------------------------------------------------
      subroutine setCmdParam(CMD,PARID,ARG)
! set numerical value of command parameter
C---------------------------------------------------------
      character(*),intent(in) :: CMD
      character(*),intent(in) :: PARID
      REAL(KIND(1.D0)),intent(in) :: ARG
      integer :: ID,IOBJ
      IOBJ=getCmdIndex(CMD)
      if (IOBJ.gt.0) then
        call GETPARINDX(myCOMMANDS(IOBJ)%IDSTR,':',trim(PARID),ID)
        if ((ID.gt.0).and.(ID.le.myCOMMANDS(IOBJ)%NARG)) then
          myCOMMANDS(IOBJ)%ARG(ID)=ARG
        else
          call MSG_WARN('unknown command parameter: '//trim(PARID),1)
        endif
      else
        call MSG_WARN('unknown command ID: '//trim(CMD),1)
      endif
      end subroutine setCmdParam

C---------------------------------------------------------
      subroutine getCmdParamS(CMD,PARID,SARG)
! return string value of command parameter
C---------------------------------------------------------
      character(*),intent(in) :: CMD
      character(*),intent(in) :: PARID
      character(*),intent(out) :: SARG
      integer :: ID,IOBJ
      SARG=' '
      IOBJ=getCmdIndex(CMD)
      if (IOBJ.gt.0) then
        call GETPARINDX(myCOMMANDS(IOBJ)%SIDSTR,':',trim(PARID),ID)
        if ((ID.gt.0).and.(ID.le.myCOMMANDS(IOBJ)%NSARG)) then
          call GETSTRPAR(myCOMMANDS(IOBJ)%SARG,':',ID,SARG)
        else
          call MSG_WARN('unknown command parameter: '//trim(PARID),1)
        endif
      else
        call MSG_WARN('unknown command ID: '//trim(CMD),1)
      endif
      end subroutine getCmdParamS

C---------------------------------------------------------
      subroutine setCmdParamS(CMD,PARID,SARG)
! set string value of command parameter
C---------------------------------------------------------
      character(*),intent(in) :: CMD
      character(*),intent(in) :: PARID
      character(*),intent(in) :: SARG
      integer :: ID,IOBJ
      IOBJ=getCmdIndex(CMD)
      if (IOBJ.gt.0) then
        call GETPARINDX(myCOMMANDS(IOBJ)%SIDSTR,':',trim(PARID),ID)
        if ((ID.gt.0).and.(ID.le.myCOMMANDS(IOBJ)%NSARG)) then
!         write(*,*) 'setCmdParamS ',trim(PARID)//'['//trim(SARG)//']'
          call SETSTRPAR(myCOMMANDS(IOBJ)%SARG,':',ID,trim(PARID))
        else
          call MSG_WARN('unknown command parameter: '//trim(PARID),1)
        endif
      else
        call MSG_WARN('unknown command ID: '//trim(CMD),1)
      endif
      end subroutine setCmdParamS

C---------------------------------------------------------
      SUBROUTINE COMMANDS_INP_R(OBJ,PNAME,ARG,NARG)
C---------------------------------------------------------
      TYPE(TCOMMAND),intent(out) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(in) :: ARG(*)
      integer,intent(out) :: NARG
      character(LEN_NAME) :: TID
      integer :: CLS,IDX,NA
      NARG=0
      call COMMANDS_PARINFO(OBJ,PNAME,IDX,TID,NA)
      if ((IDX.gt.0).and.(IDX+NA.le.OBJ%NARG+1).and.(NA.gt.0)) then
        OBJ%ARG(IDX:IDX+NA-1)=ARG(1:NA)
        NARG=NA
      endif
      END SUBROUTINE COMMANDS_INP_R


C---------------------------------------------------------
      SUBROUTINE COMMANDS_OUT_R(OBJ,PNAME,ARG,NARG)
C---------------------------------------------------------
      TYPE(TCOMMAND),intent(out) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(out) :: ARG(*)
      integer,intent(out) :: NARG
      character(LEN_NAME) :: TID
      integer :: CLS,IDX,NA
      NARG=0
      call COMMANDS_PARINFO(OBJ,PNAME,IDX,TID,NA)
!      write(*,*) 'COMMANDS_OUT_R ',trim(PNAME),CLS,IDX,trim(TID),NA,OBJ%NARG
      if ((IDX.gt.0).and.(IDX+NA.le.OBJ%NARG+1).and.(NA.gt.0)) then
        ARG(1:NA)=OBJ%ARG(IDX:IDX+NA-1)
!        write(*,*) 'COMMANDS_OUT_R ',ARG(1:NA)
        NARG=NA
      endif
      END SUBROUTINE COMMANDS_OUT_R


C---------------------------------------------------------
      SUBROUTINE COMMANDS_INP_S(OBJ,PNAME,SARG,NARG)
C---------------------------------------------------------
      TYPE(TCOMMAND),intent(out) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      CHARACTER(*),intent(in) :: SARG
      integer,intent(out) :: NARG
      character(LEN_NAME) :: TID
      integer :: CLS,IDX,NA
      NARG=0
      call COMMANDS_PARINFO(OBJ,PNAME,IDX,TID,NA)
      if (trim(TID).eq.'S') then
        if ((IDX.gt.0).and.(IDX.le.OBJ%NSARG)) then
!        write(*,*) 'COMMANDS_INP ',trim(PNAME)//'['//trim(SARG)//']'
          call SETSTRPAR(OBJ%SARG,':',IDX,trim(SARG))
          NARG=NA
        endif
      endif
      END SUBROUTINE COMMANDS_INP_S

C---------------------------------------------------------
      SUBROUTINE COMMANDS_OUT_S(OBJ,PNAME,SARG,NARG)
C---------------------------------------------------------
      TYPE(TCOMMAND),intent(out) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      CHARACTER(*),intent(out) :: SARG
      integer,intent(out) :: NARG
      character(LEN_NAME) :: TID
      integer :: CLS,IDX,NA
      SARG=' '
      NARG=0
      call COMMANDS_PARINFO(OBJ,PNAME,IDX,TID,NA)
      if (trim(TID).eq.'S') then
        if ((IDX.gt.0).and.(IDX.le.OBJ%NSARG)) then
          call GETSTRPAR(OBJ%SARG,':',IDX,SARG)
!          write(*,*) 'COMMANDS_OUT ',trim(PNAME)//'['//trim(SARG)//']',CLS,IDX,trim(TID)
          NARG=NA
        endif
      endif
      END SUBROUTINE COMMANDS_OUT_S

!--------------------------------------------------------
      subroutine CLEAR_COMMANDS
! clear definitions of all commands
!--------------------------------------------------------
      integer :: i
      COMMANDS_NUM=0
      do i=1,CMD_MAX
        call COMMANDS_DEFAULT(myCOMMANDS(i))
      enddo
      end subroutine CLEAR_COMMANDS

!--------------------------------------------------------
      subroutine ADD_COMMAND(CMD,error)
! defines new command instance from a template !
!--------------------------------------------------------
      type(TCOMMAND) :: CMD
      logical,intent(out)           :: error
        if (COMMANDS_NUM.lt.CMD_MAX) then
          COMMANDS_NUM=COMMANDS_NUM+1
          CMD%ICMD=COMMANDS_NUM
          myCOMMANDS(COMMANDS_NUM)=CMD
        else
          error=.true.
      !    call PARSEINFO
          call MSG_ERROR('COMMANDS','too many commands defined',0,1)
        endif
      end subroutine ADD_COMMAND

      end module COMMANDS
