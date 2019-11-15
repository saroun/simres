!//////////////////////////////////////////////////////////////////////
!////  $Id: commands.f,v 1.31 2018/10/29 17:50:48 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.31 $
!////     $Date: 2018/10/29 17:50:48 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Defines classes for commands
!////
!//////////////////////////////////////////////////////////////////////
      MODULE COMMANDS
      use XMLINFO
      USE XMLUTILS
      use FIELDDEF
      use CLASSDEF
      USE FIELDDATA
      implicit none

      private
      save

! COMMAND definition
      integer,parameter :: CMD_MAXARG=16
      integer,parameter :: CMD_MAX=32

! common type for all commands
      type TCOMMAND
        INTEGER :: ICLS=0 ! class index (see global definitions in CLASSES module)
        INTEGER :: ICMD=0 ! command index (pointer to myCOMMANDS array in this module)
        character(LEN_ID)     :: ID=''       ! ID string
        character(LEN_NAME)   :: NAME=''   ! name string
        integer :: NARG=0 ! number of numerical arguments
        TYPE(TFIELD) :: FIELDS(CMD_MAXARG)
        type(TCLSDEF),pointer :: CDEF
      end TYPE TCOMMAND

  ! number of commands defined
      integer :: COMMANDS_NUM
  ! instances of command classes
      type(TCOMMAND),target :: myCOMMANDS(CMD_MAX)

      interface getCmdParam
        module procedure getCmdParam1
        module procedure getCmdParam2
      end interface

      public getCmdParam,getCmdParamS,setCmdParam,TCOMMAND,getCmdIndex,COMMANDS_NUM,myCOMMANDS
      public ALLOC_COMMAND,ASSIGN_COMMAND,COMMANDS_PRN, DISPOSE_COMMAND
      public COMMANDS_INP,COMMANDS_OUT,ALLOCATE_COMMANDS, DISPOSE_ALL_COMMANDS

      contains

!--------------------------------------------------------
      subroutine DISPOSE_ALL_COMMANDS
! clear definitions of all commands
!--------------------------------------------------------
      integer :: i
      do i=1,COMMANDS_NUM
        call DISPOSE_COMMAND(myCOMMANDS(i))
      enddo
      COMMANDS_NUM=0
      end subroutine DISPOSE_ALL_COMMANDS

!--------------------------------------------------------
      subroutine ALLOCATE_COMMANDS
! allocates and adds all commands defined in classes.xml
!--------------------------------------------------------
      type(PCLSDEF) :: CDEF
      integer :: i
      logical :: error
      call DISPOSE_ALL_COMMANDS
      do i=1,CLS_NUM
        call GetClassDef(i,CDEF)
        if (CDEF%C%isCMD) then
          call ADD_COMMAND(CDEF%C,error)
        endif
      enddo
      end subroutine ALLOCATE_COMMANDS

!--------------------------------------------------------
      subroutine ADD_COMMAND(CDEF,error)
! defines new command instance from a template !
!--------------------------------------------------------
      type(TCLSDEF),pointer :: CDEF
      logical,intent(out)           :: error
      integer :: IC
      if (associated(CDEF)) then
        if (COMMANDS_NUM.lt.CMD_MAX) then
          IC=COMMANDS_NUM+1
          myCOMMANDS(IC)%ICMD=IC
          myCOMMANDS(IC)%ICLS=CDEF%IDX
          call ALLOC_COMMAND(myCOMMANDS(IC),CDEF,error)
          if (.not.error) COMMANDS_NUM=IC
        else
          error=.true.
          call MSG_ERROR('COMMANDS','too many commands defined',0,1)
        endif
      endif
      end subroutine ADD_COMMAND

!--------------------------------------------------------
      subroutine ALLOC_COMMAND(OBJ,CDEF,error)
! allocate command fields and initiate from definition object TCLSDEF !
!--------------------------------------------------------
      type(TCOMMAND) :: OBJ
      type(TCLSDEF),pointer :: CDEF
      logical,intent(out)           :: error
      integer :: i
        OBJ%CDEF=>CDEF
        OBJ%ID=trim(CDEF%ID)
        OBJ%NAME=trim(CDEF%NAME)
        OBJ%NARG=0
        do i=1,CDEF%NF
          OBJ%FIELDS(i)%DEF=>CDEF%FIELDS(i)%F
          if (ALLOC_FIELD(OBJ%FIELDS(i),CDEF%FIELDS(i)%F%DIM)) then
            OBJ%NARG=OBJ%NARG+1
          else
            error=.true.
            write(*,*) 'COMMANDS','cannot allocate field '//trim(CDEF%ID),' obj=',trim(OBJ%ID),' cls=',trim(OBJ%CDEF%ID)
            write(*,*) '       field=',i,trim(CDEF%FIELDS(i)%F%ID)
            call MSG_ERROR('COMMANDS','cannot allocate field '//trim(CDEF%ID),0,1)
            exit
          endif
        enddo
        error=(CDEF%NF.ne.OBJ%NARG)
        if (error) call MSG_WARN('Some fields of '//trim(CDEF%ID)//' where not allocated.',1)
      end subroutine ALLOC_COMMAND

!---------------------------------------------------------

      SUBROUTINE DISPOSE_COMMAND(OBJ)
! deallocate field buffers of the object
! preserves CDEF = pointer to CLSDEF
!--------------------------------------------------------

      TYPE(TCOMMAND) :: OBJ
      integer :: i
      do i=1,OBJ%NARG
        call DISPOSE_FIELD(OBJ%FIELDS(i))
      enddo
      OBJ%NARG=0









      OBJ%ICLS=0
      OBJ%ICMD=0
      OBJ%ID=' '
      OBJ%NAME=' '      
      end SUBROUTINE DISPOSE_COMMAND



C---------------------------------------------------------
      SUBROUTINE ASSIGN_COMMAND(CMD)
C---------------------------------------------------------
      TYPE(TCOMMAND) :: CMD
      integer :: i,j
      logical :: error
      i=getCmdIndex(trim(CMD%ID))
      if (i.gt.0) then
  !      write(*,*) 'ASSIGN_COMMAND i=',i,' ID=',trim(CMD%ID)
        do j=1,CMD%NARG
          call FIELD_COPY(CMD%FIELDS(j),myCOMMANDS(i)%FIELDS(j),error)
          if (error) then
            call send_message('ADD_FIELD',CMD%ID,CMD%FIELDS(j)%ID)
            exit
          endif
        enddo

      endif
      end SUBROUTINE ASSIGN_COMMAND

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
        write(IO,1) adjustl(myCOMMANDS(i)%ID),trim(adjustl(myCOMMANDS(i)%CDEF%HINT))
      enddo
      END SUBROUTINE COMMANDS_PRN

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
        if (icomp.eq.0) then
          ic=i
          is=1
          EXIT
        else if (icomp.eq.1) then
          is=is+1
          ic=i
        endif
      enddo
  !    write(*,*) 'getCmdIndex ',trim(IDSTR),' is=',is,' ic=',ic
      if (is.ne.1) ic=0
      getCmdIndex=ic
      end function getCmdIndex

C---------------------------------------------------------
      integer function getParIndex(CMD,IDSTR)
! return command index for given ID string
! unique substring is handled as an exact match => abbreviations are allowed
C---------------------------------------------------------
      type(TCOMMAND) :: CMD
      character(*) :: IDSTR
      integer i,IDCOMPARE
      logical :: found
      found=.false.
      i=0
      do WHILE((.not.found).and.(i.lt.CMD%NARG))
        i=i+1
        found=(IDCOMPARE(trim(IDSTR),trim(CMD%FIELDS(i)%ID)).eq.0)
      enddo
      if (found) then
        getParIndex=i
      else
        getParIndex=0
      endif
      end function getParIndex

C---------------------------------------------------------
      REAL(KIND(1.D0)) function getCmdParam1(CMD,PARID)
! return numerical value of command parameter
! for arrays, return the first field defined by IDX
C---------------------------------------------------------
      character(*),intent(in) :: CMD
      character(*),intent(in) :: PARID
      getCmdParam1=getCmdParam2(CMD,1,PARID)
      end function getCmdParam1

C---------------------------------------------------------
      REAL(KIND(1.D0)) function getCmdParam2(CMD,IDX,PARID)
! return numerical value of command parameter
! for arrays, return the first field defined by IDX
C---------------------------------------------------------
      character(*),intent(in) :: CMD
      character(*),intent(in) :: PARID
      integer,intent(in) :: IDX
      integer :: IOBJ,IPAR,N,IX
      REAL(KIND(1.D0)) :: AUX(1)
      AUX=0.D0
      if (IDX.gt.0) then
          IOBJ=getCmdIndex(CMD)
          if (IOBJ.gt.0) then
            IPAR=getParIndex(myCOMMANDS(IOBJ),trim(PARID))
            if (IPAR.gt.0) then
              IX=min(IDX,myCOMMANDS(IOBJ)%FIELDS(IPAR)%NP)
		!  write(*,*) 'getCmdParam ',trim(CMD),'.',trim(PARID),' IX=',IX,' IPAR=',IPAR
              call FIELD2ARRAY(myCOMMANDS(IOBJ)%FIELDS(IPAR),IX,AUX,1,N)
		!  if (CMD.eq.'PLOT') write(*,*) 'getCmdParam PARID='//trim(PARID)//' result=',AUX(1:N)
            else

              call send_message('P',PARID,' ')
            endif
          else


            call send_message('C',CMD,' ')
          endif
      endif
      getCmdParam2=AUX(1)
      end function getCmdParam2

C---------------------------------------------------------
      SUBROUTINE setCmdParam(CMD,PARID,A)
! set numerical value of command parameter
C---------------------------------------------------------
      character(*),intent(in) :: CMD
      character(*),intent(in) :: PARID
      REAL(KIND(1.D0)),intent(in) :: A
      integer :: IOBJ,N
      REAL(KIND(1.D0)) :: AUX(1)
      type(TFIELD) :: ARG
      TYPE(PCLSDEF) :: CDEF
      TYPE(TFIELDDEF),target :: FDEF
      AUX=A
      IOBJ=getCmdIndex(CMD)
      if (IOBJ.gt.0) then
        call GetClassDef(trim(CMD),CDEF)
        if (associated(CDEF%C)) then
          call GetFieldDef(CDEF%C%IDX,trim(PARID),FDEF)
          if (FDEF%IDX.gt.0) then
            ARG%DEF=> FDEF
            if (ALLOC_FIELD(ARG,FDEF%DIM)) then
              call ARRAY2FIELD(ARG,1,AUX,1,N)
              call COMMANDS_INP(myCOMMANDS(IOBJ),trim(PARID),ARG,N)
              call DISPOSE_FIELD(ARG)
            endif
          endif
        endif
      endif
      end SUBROUTINE setCmdParam

C---------------------------------------------------------
      subroutine getCmdParamArray(CMD,PARID,A,NARG)
! return numerical values of command parameter as an array
C---------------------------------------------------------
      character(*),intent(in) :: CMD
      character(*),intent(in) :: PARID
      REAL(KIND(1.D0)),intent(out) :: A(:)
      integer,intent(out) :: NARG
      integer :: IOBJ,IPAR,i
      REAL(KIND(1.D0)) :: AUX(FIELD_BUFFER_FLOAT)
      AUX=0.D0
      NARG=0
      IOBJ=getCmdIndex(CMD)
      if (IOBJ.gt.0) then
        IPAR=getParIndex(myCOMMANDS(IOBJ),trim(PARID))
        if (IPAR.gt.0) then
          call FIELD2ARRAY(myCOMMANDS(IOBJ)%FIELDS(IPAR),0,AUX,FIELD_BUFFER_FLOAT,NARG)
          do i=1,min(NARG,SIZE(A))
            A(i)=AUX(i)
          enddo
        else
          call send_message('P',PARID,' ')
        endif
      else
        call send_message('C',CMD,' ')
      endif
      end subroutine getCmdParamArray

C---------------------------------------------------------
      subroutine getCmdParamS(CMD,IDX,PARID,SARG)
! return numerical value of command parameter
! for arrays, return the first field defined by IDX
C---------------------------------------------------------
      character(*),intent(in) :: CMD
      character(*),intent(in) :: PARID
      integer,intent(in) :: IDX
      character(*),intent(out) :: SARG
      integer :: IOBJ,IPAR,N,IX
      SARG=' '
      if (IDX.le.0) return
      IOBJ=getCmdIndex(CMD)
      if (IOBJ.gt.0) then
        IPAR=getParIndex(myCOMMANDS(IOBJ),trim(PARID))
        if (IPAR.gt.0) then
          IX=min(IDX,myCOMMANDS(IOBJ)%FIELDS(IPAR)%NP)
          call FIELD2STR(myCOMMANDS(IOBJ)%FIELDS(IPAR),IX,SARG,N)
        else
          call send_message('P',PARID,' ')
        endif
      else
        call send_message('C',CMD,' ')
      endif
      end subroutine getCmdParamS

!---------------------------------------------------------
      SUBROUTINE COMMANDS_INP(OBJ,PNAME,ARG,NARG)
!---------------------------------------------------------
      TYPE(TCOMMAND) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: IP
      logical :: error
      NARG=0
  !    write(*,*) 'COMMANDS_INP  PNAME=',trim(PNAME)
      if (trim(PNAME).ne.trim(ARG%ID)) then
        write(*,*) 'COMMANDS_INP: unequal argument IDs: ['//trim(PNAME)//','//trim(ARG%ID)//']'
        return
      endif
      IP=getParIndex(OBJ,PNAME)
      if (IP.gt.0) then
        call FIELD_COPY(ARG,OBJ%FIELDS(IP),error)
        if (.not.error) NARG=OBJ%FIELDS(IP)%NP
      else
        write(*,*) 'COMMANDS_INP: undefined parameter : ['//trim(PNAME)//']'
      endif
      END SUBROUTINE COMMANDS_INP

!---------------------------------------------------------
      SUBROUTINE COMMANDS_OUT(OBJ,PNAME,ARG,NARG)
!---------------------------------------------------------
      TYPE(TCOMMAND) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: IP
      logical :: error
      NARG=0
      if (trim(PNAME).ne.trim(ARG%ID)) then
        write(*,*) 'COMMANDS_OUT: unequal argument IDs: ['//trim(PNAME)//','//trim(ARG%ID)//']'
        return
      endif
      IP=getParIndex(OBJ,PNAME)
      if (IP.gt.0) then
        call FIELD_COPY(OBJ%FIELDS(IP),ARG,error)
        if (.not.error) NARG=ARG%NP
      else
        write(*,*) 'COMMANDS_OUT: undefined parameter : ['//trim(PNAME)//']'
      endif
      END SUBROUTINE COMMANDS_OUT

!--------------------------------------------------------
      subroutine send_message(what,ID,TID)
!--------------------------------------------------------
      character(*),intent(in) :: what,ID,TID
      select case(what)
      case('NR')
        call MSG_WARN('cannot read numerical value from: '//trim(ID)//':'//trim(tid),1)
      case('NA')
        call MSG_WARN('cannot assign numerical value to: '//trim(ID)//':'//trim(tid),1)
      case('SR')
        call MSG_WARN('cannot read string value from: '//trim(ID)//':'//trim(tid),1)
      case('SA')
        call MSG_WARN('cannot assign string value to: '//trim(ID)//':'//trim(tid),1)
      case('P')
        call MSG_WARN('unknown command parameter: '//trim(ID),1)
      case('C')
        call MSG_WARN('unknown command ID: '//trim(ID),1)
      case('ASSIGN_ICLS')
        call MSG_WARN('ASSIGN_COMMAND: wrong ICLS for '//trim(ID)//':'//trim(tid),1)
      case('ASSIGN_ICMD')
        call MSG_WARN('ASSIGN_COMMAND: wrong ICMD for '//trim(ID)//':'//trim(tid),1)
      case('ADD_FIELD')
        call MSG_WARN('ASSIGN_COMMAND: cannot add field '//trim(ID)//':'//trim(tid),1)
      case default
        call MSG_WARN('undefined error: '//trim(what),1)
      end select

      !READ(*,*)
      end subroutine send_message

      end module COMMANDS
