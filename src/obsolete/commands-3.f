!//////////////////////////////////////////////////////////////////////
!////  $Id: commands-3.f,v 1.1 2009/02/17 23:35:42 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.1 $
!////     $Date: 2009/02/17 23:35:42 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Defines classes for commands
!////
!//////////////////////////////////////////////////////////////////////
      MODULE COMMANDS
      use XMLINFO
      USE XMLUTILS
      implicit none


! COMMAND definition
      integer,parameter,PRIVATE :: CMD_MAXARG=16
      integer,parameter,PRIVATE :: CMD_MAX=32
      integer,parameter :: CMD_SLEN=256 ! max. length of string argument
! string field dimensions: keep private
      integer,parameter,PRIVATE :: LEN_ID=16
      integer,parameter,PRIVATE :: LEN_NAME=32

! Intrinsic type VRANGE defines a range of a component parameter
      TYPE VRANGE;SEQUENCE
        CHARACTER(2*LEN_ID+1) :: ID ! parameer ID in format COMPONENT.PAR
        REAL(KIND(1.D0)) :: VMI     ! minimum value
        REAL(KIND(1.D0)) :: VMA     ! maximum value
        REAL(KIND(1.D0)) :: VINC    ! increment
      end type VRANGE

! intrinsic field type, can represent an array of integers, floats, strings or value ranges
! fields are dynamically allocated
      TYPE TFIELD
        CHARACTER(LEN_ID) :: ID=''
        CHARACTER(LEN_ID) :: TID=''
        integer :: LENGTH=0
        real(kind(1.D0)), DIMENSION(:),POINTER :: RFIELD => NULL()
        INTEGER, DIMENSION(:),POINTER :: IFIELD => NULL()
        CHARACTER(CMD_SLEN), DIMENSION(:), POINTER :: SFIELD => NULL()
        TYPE(VRANGE), DIMENSION(:), POINTER :: VFIELD => NULL()
      end type TFIELD

! common type for all commands
      type TCOMMAND
        INTEGER :: ICLS=0 ! class index (see global definitions in CLASSES module)
        INTEGER :: ICMD=0 ! command index (pointer to myCOMMANDS array in this module)
        character(LEN_ID)     :: ID=''       ! ID string
        character(LEN_NAME)   :: NAME=''   ! name string
        CHARACTER(64) :: HINT=''           ! description
        integer :: NARG=0 ! number of numerical arguments
        TYPE(TFIELD) :: FIELDS(CMD_MAXARG)
      end TYPE TCOMMAND


  ! number of commands defined
      integer :: COMMANDS_NUM
  ! instances of command classes
      type(TCOMMAND),target :: myCOMMANDS(CMD_MAX)

      interface COMMANDS_INP
        module procedure COMMANDS_INP_R
        module procedure COMMANDS_INP_S
        module procedure COMMANDS_INP_V
      end interface

      interface COMMANDS_OUT
        module procedure COMMANDS_OUT_R
        module procedure COMMANDS_OUT_S
        module procedure COMMANDS_OUT_V
      end interface

      interface COMMANDS_CLONE
        module procedure COMMANDS_CLONE_ID
        module procedure COMMANDS_CLONE_OBJ
      end interface


! instance of TCOMMAND used at loaing the commands definitions
!      type(TCOMMAND), private :: myCMD

      private COMMANDS_INP_R
      private COMMANDS_INP_S
      private COMMANDS_INP_V
      private COMMANDS_OUT_R
      private COMMANDS_OUT_S
      private COMMANDS_OUT_V

      private COMMANDS_CLONE_ID
      private COMMANDS_CLONE_OBJ

      private send_message
      contains

!--------------------------------------------------------
      subroutine CLEAR_COMMANDS
! clear definitions of all commands
!--------------------------------------------------------
      integer :: i
      do i=1,COMMANDS_NUM
        call COMMANDS_DEFAULT(myCOMMANDS(i))
      enddo
      COMMANDS_NUM=0
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
        !  call COMMANDS_DISPOSE(myCOMMANDS(CMD%ICMD))
          call COMMANDS_CLONE_OBJ(CMD,myCOMMANDS(CMD%ICMD))
  !        write(*,*) CMD%ICMD,' ADD_COMMAND  ',trim(myCOMMANDS(CMD%ICMD)%ID),' loc=',POINTER(myCOMMANDS(CMD%ICMD))
        !  read(*,*)
        else
          error=.true.
      !    call PARSEINFO
          call MSG_ERROR('COMMANDS','too many commands defined',0,1)
        endif
      end subroutine ADD_COMMAND

C---------------------------------------------------------
      SUBROUTINE FIELD_COPY(SRC,TAR)
C---------------------------------------------------------
      TYPE(TFIELD),intent(inout) :: SRC
      TYPE(TFIELD),intent(inout) :: TAR
      REAL(KIND(1.D0)),allocatable :: R(:)
      integer,allocatable :: I(:)
      character(CMD_SLEN),allocatable :: S(:)
      type(VRANGE),allocatable :: V(:)
      integer :: L
        if (trim(SRC%TID).eq.trim(TAR%TID)) then
          if (SRC%LENGTH.eq.TAR%LENGTH) then
            L=SRC%LENGTH
            TAR%ID=SRC%ID
            select case (SRC%TID)
            case('F')
              if ((associated(SRC%RFIELD)).and.(associated(TAR%RFIELD))) then
                allocate(R(L))
                R(1:L)=SRC%RFIELD(1:L)
                TAR%RFIELD(1:L)=R(1:L)
                deallocate(R)
              endif
            case('I')
              if ((associated(SRC%IFIELD)).and.(associated(TAR%IFIELD))) then
                allocate(I(L))
                I(1:L)=SRC%IFIELD(1:L)
                TAR%IFIELD(1:L)=I(1:L)
                deallocate(I)
              endif
            case('S')
              if ((associated(SRC%SFIELD)).and.(associated(TAR%SFIELD))) then
                allocate(S(L))
                S(1:L)=SRC%SFIELD(1:L)
                TAR%SFIELD(1:L)=S(1:L)
                deallocate(S)
              endif
            case('V')
              if ((associated(SRC%VFIELD)).and.(associated(TAR%VFIELD))) then
                allocate(V(L))
                V(1:L)=SRC%VFIELD(1:L)
                TAR%VFIELD(1:L)=V(1:L)
                deallocate(V)
              endif
            end select
          else
            write(*,*) trim(SRC%ID),' FIELD_COPY unequal lengths: [',SRC%LENGTH,',',TAR%LENGTH,']'
          endif
        else
          write(*,*) trim(SRC%ID),' FIELD_COPY unequal types: [',trim(SRC%TID),',',trim(TAR%TID),']'
        endif
      end SUBROUTINE FIELD_COPY

C---------------------------------------------------------
      SUBROUTINE ADD_FIELD(OBJ,id,tid,isize,error)
C---------------------------------------------------------
      TYPE(TCOMMAND) :: OBJ
      character(*),intent(in) :: id
      character(*),intent(in) :: tid
      integer,intent(in) :: isize
      logical,intent(out) :: error
      integer ierr
      ierr=0

    !  call ListCommandPar(OBJ)
      if (OBJ%NARG.lt.CMD_MAXARG) then
  !      write(*,*) OBJ%ICMD,'    ADD_FIELD '//trim(OBJ%id)//'.'//trim(id),' ',trim(tid),' ',isize
  !      call ListCommandPar(OBJ)
        OBJ%NARG=OBJ%NARG+1
        OBJ%FIELDS(OBJ%NARG)%ID=ID
        OBJ%FIELDS(OBJ%NARG)%TID=TID
        ierr=ALLOCATE_FIELD(OBJ,OBJ%NARG,isize)
  !      call ListCommandPar(OBJ)
  !      write(*,*) OBJ%ICMD,'    ADD_FIELD done '//trim(OBJ%ID)//'.'//trim(OBJ%FIELDS(OBJ%NARG)%TID)
      else
        ierr=1
        call MSG_ERROR('COMMANDS','too many fields in '//trim(OBJ%ID),0,1)
      endif
      error=(ierr.ne.0)
      if (error) write(*,*) OBJ%ICMD,'    ADD_FIELD  error'
    !  call ListCommandPar(OBJ)
    !  read(*,*)
      end subroutine ADD_FIELD


C---------------------------------------------------------
      integer function DISPOSE_FIELD(OBJ,i)
C clear buffers of the object
C---------------------------------------------------------
      TYPE(TCOMMAND) :: OBJ
      integer i
      integer :: ierr
      logical :: dbg
      dbg=.false.
      ierr=0
      if ((i.gt.0).and.(i.le.OBJ%NARG)) then
        if (dbg) write(*,*) OBJ%ICMD,' DISPOSE_FIELD  ',trim(OBJ%ID),'.',trim(OBJ%FIELDS(i)%ID),' i=',i
        if (dbg) write(*,*) '     fields loc=',pointer(OBJ%FIELDS(i))
        select case(OBJ%FIELDS(i)%TID)
        case('F')
          if (dbg) write(*,*) '     R, loc=',pointer(OBJ%FIELDS(i)%RFIELD)
          if (associated(OBJ%FIELDS(i)%RFIELD)) deallocate(OBJ%FIELDS(i)%RFIELD,STAT=ierr)
        case('I')
          if (dbg)  write(*,*) '     I, loc=',pointer(OBJ%FIELDS(i)%IFIELD)
          if (associated(OBJ%FIELDS(i)%IFIELD)) deallocate(OBJ%FIELDS(i)%IFIELD,STAT=ierr)
        case('S')
          if (dbg) write(*,*) '     S, loc=',pointer(OBJ%FIELDS(i)%SFIELD)
          if (associated(OBJ%FIELDS(i)%SFIELD)) deallocate(OBJ%FIELDS(i)%SFIELD,STAT=ierr)
        case('VRANGE')
          if (dbg) write(*,*) '     V, loc=',pointer(OBJ%FIELDS(i)%VFIELD)
          if (associated(OBJ%FIELDS(i)%VFIELD)) deallocate(OBJ%FIELDS(i)%VFIELD,STAT=ierr)
        end select
        if (ierr.ne.0) then
          write(*,*) 'DISPOSE_FIELD err=',ierr,i,trim(OBJ%FIELDS(i)%ID)
        endif
        OBJ%FIELDS(i)%RFIELD=>NULL()
        OBJ%FIELDS(i)%IFIELD=>NULL()
        OBJ%FIELDS(i)%SFIELD=>NULL()
        OBJ%FIELDS(i)%VFIELD=>NULL()
        OBJ%FIELDS(i)%LENGTH=0
      endIF
      DISPOSE_FIELD=ierr
      if (dbg) write(*,*) OBJ%ICMD,' DISPOSE_FIELD  done ',ierr
      end function DISPOSE_FIELD

C---------------------------------------------------------
      integer function ALLOCATE_FIELD(OBJ,i,isize)
C clear buffers of the object
C---------------------------------------------------------
      TYPE(TCOMMAND) :: OBJ
      integer i,isize
      integer :: ierr
      logical :: dbg
      dbg=.false.
      ierr=0
      if ((i.gt.0).and.(i.le.OBJ%NARG)) then
         if (dbg) write(*,*) OBJ%ICMD,' ALLOCATE_FIELD  ',trim(OBJ%ID),'.',trim(OBJ%FIELDS(i)%ID),' i=',i,' size=',isize
      !  ierr=DISPOSE_FIELD(OBJ,i)
        if ((ierr.eq.0).and.(isize.gt.0)) then
          if (dbg) write(*,*) '     fields loc=',pointer(OBJ%FIELDS(i))
          select case (OBJ%FIELDS(i)%TID)
            case('F')
              allocate(OBJ%FIELDS(i)%RFIELD(1:isize),STAT=ierr)
              if (dbg) write(*,*) '     R, loc=',pointer(OBJ%FIELDS(i)%RFIELD)
            case('I')
              allocate(OBJ%FIELDS(i)%IFIELD(1:isize),STAT=ierr)
              if (dbg) write(*,*) '     I, loc=',pointer(OBJ%FIELDS(i)%IFIELD)
            case('S')
              allocate(OBJ%FIELDS(i)%SFIELD(1:isize),STAT=ierr)
              if (dbg) write(*,*) '     S, loc=',pointer(OBJ%FIELDS(i)%SFIELD)
            case('VRANGE')
              allocate(OBJ%FIELDS(i)%VFIELD(1:isize),STAT=ierr)
              if (dbg) write(*,*) '     V, loc=',pointer(OBJ%FIELDS(i)%VFIELD)
          end select
          if (ierr.ne.0) then
            call MSG_ERROR('ALLOCATE_FIELD','cannot allocate for '//trim(OBJ%ID)//'.'//trim(OBJ%FIELDS(i)%ID),0,1)
            OBJ%FIELDS(i)%LENGTH=0
          ELSE
            OBJ%FIELDS(i)%LENGTH=isize
          endif
        endif
      endIF
      ALLOCATE_FIELD=ierr
      if (dbg) write(*,*) OBJ%ICMD,' ALLOCATE_FIELD  done ',ierr
      end function ALLOCATE_FIELD


C---------------------------------------------------------
      SUBROUTINE COMMANDS_DISPOSE(OBJ)
C clear buffers of the object
C---------------------------------------------------------
      TYPE(TCOMMAND) :: OBJ
      integer :: i,ierr,ip
      ip=POINTER(OBJ)
    !  if (OBJ%NARG.gt.0) write(*,*) 'COMMANDS_DISPOSE ',trim(OBJ%ID),' loc=',ip
      do i=1,OBJ%NARG
        ierr=DISPOSE_FIELD(OBJ,i)
        if (ierr.ne.0) then
          write(*,*) 'COMMANDS_DISPOSE err=',ierr,i,trim(OBJ%FIELDS(i)%ID)
    !      STOP
        endif
        OBJ%FIELDS(i)%RFIELD=>NULL()
        OBJ%FIELDS(i)%IFIELD=>NULL()
        OBJ%FIELDS(i)%SFIELD=>NULL()
        OBJ%FIELDS(i)%VFIELD=>NULL()
        OBJ%FIELDS(i)%LENGTH=0
      enddo
  !    if (OBJ%NARG.gt.0) write(*,*) 'COMMANDS_DISPOSE done for ',trim(OBJ%ID),' loc=',ip
      OBJ%NARG=0
    !  read(*,*)
      end SUBROUTINE COMMANDS_DISPOSE


C---------------------------------------------------------
      SUBROUTINE COMMANDS_DEFAULT(OBJ)
C clear buffers of the object
C---------------------------------------------------------
      TYPE(TCOMMAND) :: OBJ
      integer :: i
    !  write(*,*) 'COMMANDS_DEFAULT ',OBJ%ID,OBJ%NARG
    !  read(*,*)
      call COMMANDS_DISPOSE(OBJ)
      OBJ%ICLS=0
      OBJ%ICMD=0
      OBJ%ID=' '
      OBJ%HINT='undefined'
      OBJ%NAME=' '
      end SUBROUTINE COMMANDS_DEFAULT


!---------------------------------------------------------
      SUBROUTINE COMMANDS_PARINFO(OBJ,PNAME,IDX,TID,NARG)
! return actual class, of which PNAME is a member.
! RETURNS
!   IDX   ... index of the parameter in IDSTR or SIDSTR list
!   TID   ... ID string of variable type (S..string, N..numeric)
!   NARG  ... size of the parameter (>1 for array types)
!---------------------------------------------------------
      TYPE(TCOMMAND) :: OBJ
      character(*),intent(in) :: PNAME
      integer,intent(out) :: IDX,NARG
      character(*),intent(out) :: TID
      integer :: I
  !    write(*,*) OBJ%ICMD,'COMMANDS_PARINFO '//trim(OBJ%ID)//' '//trim(PNAME)
      i=getParIndex(OBJ,PNAME)
  !    write(*,*) OBJ%ICMD,'COMMANDS_PARINFO '//trim(OBJ%ID)//' '//trim(PNAME),' ',i
  !    if ((i.eq.7).and.(trim(OBJ%ID).eq.'SCAN2D')) then
  !      write(*,*) 'COMMANDS_PARINFO   FIELDS(7)=',OBJ%FIELDS(7)%IFIELD(1)
  !      write(*,*) 'COMMANDS_PARINFO             ',OBJ%FIELDS(i)%LENGTH,OBJ%FIELDS(i)%TID
  !    endif

      if (i.gt.0) then
        IDX=I
        NARG=OBJ%FIELDS(i)%LENGTH
        TID=trim(OBJ%FIELDS(i)%TID)
      else
        IDX=-1
        NARG=0
        TID=' '
      ENDIF
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
      SUBROUTINE COMMANDS_CLONE_OBJ(SRC,TAR)
! Copy command object SRC to TAR
! As = on types with pointer fields makes soft copy only,
!  we need a special procedure to make a real clone.
C---------------------------------------------------------
      type(TCOMMAND) :: SRC
      type(TCOMMAND) :: TAR
      integer :: i,j
      logical :: error
      character(LEN_NAME) :: fid,ftyp
  !    write(*,*) 'COMMANDS_CLONE_OBJ ',trim(SRC%ID)
  !    call ListCommandPar(SRC)
      call COMMANDS_DISPOSE(TAR)
      TAR%ID=SRC%ID
      TAR%ICLS=SRC%ICLS
      TAR%ICMD=SRC%ICMD
      TAR%NAME=SRC%NAME
      TAR%HINT=SRC%HINT
  !    call ListCommandPar(TAR)
      do i=1,SRC%NARG
        fid=SRC%FIELDS(i)%id
        ftyp=SRC%FIELDS(i)%tid
        call ADD_FIELD(TAR,trim(fid),trim(ftyp),SRC%FIELDS(i)%length,error)
        call FIELD_COPY(SRC%FIELDS(i),TAR%FIELDS(i))
      enddo
  !    call ListCommandPar(TAR)
  !    write(*,*) 'COMMANDS_CLONE_OBJ done ',trim(TAR%ID),' src=',POINTER(SRC),' tar=',POINTER(TAR)
  !    read(*,*)
      end SUBROUTINE COMMANDS_CLONE_OBJ


C---------------------------------------------------------
      SUBROUTINE COMMANDS_CLONE_ID(IDSTR,TAR)
! return command structure for given ID string
! unique substring is handled as an exact match => abbreviations are allowed
! NOTE: a clone of the command data is returned. To receive a pointer
! to actual command instance, use getCmdObjByID from vctable_commands.mod
C---------------------------------------------------------
      character(*),intent(in) :: IDSTR
      type(TCOMMAND) :: TAR
      integer :: ic
      logical :: error
      ! CMD%NARG=0
  !    write(*,*) 'COMMANDS_CLONE 1 ',trim(IDSTR)

    !  read(*,*)
    !  call COMMANDS_DEFAULT(CMD)
    !  write(*,*) 'COMMANDS_CLONE 2 '
    !  read(*,*)
      ic=getCmdIndex(IDSTR)

    !  write(*,*) 'COMMANDS_CLONE 3 ',ic
    !  read(*,*)
      if (ic.gt.0) then
        call COMMANDS_CLONE_OBJ(myCOMMANDS(ic),TAR)
      !  write(*,*) 'COMMANDS_CLONE_ID ',trim(CMD%ID),' src=',POINTER(myCOMMANDS(ic)),' tar=',POINTER(CMD)
      else
        call MSG_ERROR('COMMANDS_CLONE','undefined command '//trim(IDSTR),0,1)
      endif
      end SUBROUTINE COMMANDS_CLONE_ID

C---------------------------------------------------------
      subroutine COMMANDS_ASSIGN(OBJ)
! assign contents of OBJ to actual command instance
! ID and ICLS must agree
! OBJ%ICMD determines command index
C---------------------------------------------------------
      type(TCOMMAND),POINTER :: OBJ
      integer :: i
      logical :: error
      character(LEN_NAME) :: fid,ftyp
      if ((OBJ%ICMD.gt.0).and.(OBJ%ICMD.le.COMMANDS_NUM)) then
        if (OBJ%ICLS.ne.myCOMMANDS(OBJ%ICMD).ICLS) then
          call send_message('ASSIGN_ICLS',OBJ%ID,myCOMMANDS(OBJ%ICMD).ID)
        else if (OBJ%ID.ne.myCOMMANDS(OBJ%ICMD).ID) then
          call send_message('ASSIGN_ICMDS',OBJ%ID,myCOMMANDS(OBJ%ICMD).ID)
        ELSE
        !  call COMMANDS_DISPOSE(myCOMMANDS(OBJ%ICMD))
          call COMMANDS_CLONE_OBJ(OBJ,myCOMMANDS(OBJ%ICMD))
        endif
      endif
      end subroutine COMMANDS_ASSIGN

C---------------------------------------------------------
      SUBROUTINE ListCommandPar(OBJ)
! return command index for given ID string
! unique substring is handled as an exact match => abbreviations are allowed
C---------------------------------------------------------
      type(TCOMMAND) :: OBJ
      integer :: i
      write(*,*) '     List of ',trim(OBJ%ID),' loc=',POINTER(OBJ)
      do i=1,OBJ%NARG
        write(*,*) '     ',i,'   ',trim(OBJ%ID),'.',OBJ%FIELDS(i)%ID,OBJ%FIELDS(i)%TID,OBJ%FIELDS(i)%LENGTH
      enddo
      end SUBROUTINE ListCommandPar

C---------------------------------------------------------
      integer function getCmdIndex(IDSTR)
! return command index for given ID string
! unique substring is handled as an exact match => abbreviations are allowed
C---------------------------------------------------------
      character(*),intent(in) :: IDSTR
      integer :: i,ic,IDCOMPARE,is,icomp
      ic=0
      is=0
    !  do i=1,COMMANDS_NUM
    !    write(*,*) i,'   getCmdIndex ',trim(myCOMMANDS(i)%ID),' loc=',POINTER(myCOMMANDS(i))
   !   enddo

      do i=1,COMMANDS_NUM
        icomp=IDCOMPARE(trim(IDSTR),trim(myCOMMANDS(i)%ID))
      !  if (trim(IDSTR).eq.'PLOT') write(*,*) '   getCmdIndex ',trim(IDSTR),' ',trim(myCOMMANDS(i)%ID),' ',icomp
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
    !  write(*,*) '   getCmdIndex ',ic
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
      logical :: dbg
      found=.false.
      i=0
      dbg=.true.
      if (dbg) write(*,*) CMD%ICMD,'getParIndex '//trim(CMD%ID)//' '//trim(IDSTR),' NARG=',CMD%NARG

      do WHILE((.not.found).and.(i.lt.CMD%NARG))
        i=i+1
        if (dbg) write(*,*) '    ',i,trim(CMD%FIELDS(i)%ID)
        if (dbg) write(*,*) '    ',i,trim(IDSTR)
        found=(IDCOMPARE(trim(IDSTR),trim(CMD%FIELDS(i)%ID)).eq.0)
      enddo
      if (found) then
        getParIndex=i
        if (dbg) write(*,*) '     found ',i,trim(IDSTR)
      else
        getParIndex=0
      endif
      end function getParIndex


C---------------------------------------------------------
      REAL(KIND(1.D0)) function getCmdParam(CMD,PARID)
! return numerical value of command parameter
! for arrays, return the first field
C---------------------------------------------------------
      character(*),intent(in) :: CMD
      character(*),intent(in) :: PARID
      integer :: ID,IOBJ
      character(LEN_ID) :: tid
      REAL(KIND(1.D0)) :: res
      res=-1.D0
      IOBJ=getCmdIndex(CMD)
      if (IOBJ.gt.0) then
        ID=getParIndex(myCOMMANDS(IOBJ),PARID)
        if (ID.gt.0) then
          tid=myCOMMANDS(IOBJ)%FIELDS(ID)%TID
          select case(tid)
          case('F')
            res=myCOMMANDS(IOBJ)%FIELDS(ID)%RFIELD(1)
          case('I')
            res=myCOMMANDS(IOBJ)%FIELDS(ID)%IFIELD(1)
          case default
            call send_message('NR',PARID,TID)
          end select
        else
          call send_message('P',PARID,' ')
        endif
      else
        call send_message('C',CMD,' ')
      endif
      getCmdParam=res
      end function getCmdParam

C---------------------------------------------------------
      subroutine setCmdParam(CMD,PARID,ARG)
! set numerical value of command parameter
C---------------------------------------------------------
      character(*),intent(in) :: CMD
      character(*),intent(in) :: PARID
      REAL(KIND(1.D0)),intent(in) :: ARG
      integer :: ID,IOBJ
      character(LEN_ID) :: tid
      IOBJ=getCmdIndex(CMD)
      if (IOBJ.gt.0) then
        ID=getParIndex(myCOMMANDS(IOBJ),PARID)
        if (ID.gt.0) then
          tid=myCOMMANDS(IOBJ)%FIELDS(ID)%TID
          select case(tid)
          case('F')
            myCOMMANDS(IOBJ)%FIELDS(ID)%RFIELD(1)=ARG
          case('I')
            myCOMMANDS(IOBJ)%FIELDS(ID)%IFIELD(1)=NINT(ARG)
          case default
            call send_message('NA',PARID,TID)
          end select
        else
          call send_message('P',PARID,' ')
        endif
      else
        call send_message('C',CMD,' ')
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
      character(LEN_ID) :: tid
      character(LEN_NAME) :: SID,SMI,SMA,SINC
      REAL(KIND(1.D0)) :: res
      SARG=' '
      IOBJ=getCmdIndex(CMD)
      if (IOBJ.gt.0) then
        ID=getParIndex(myCOMMANDS(IOBJ),PARID)
        if (ID.gt.0) then
          tid=myCOMMANDS(IOBJ)%FIELDS(ID)%TID
          select case(tid)
          case('S')
            SARG=myCOMMANDS(IOBJ)%FIELDS(ID)%SFIELD(1)
          case('F')
            call FLOAT2STR(myCOMMANDS(IOBJ)%FIELDS(ID)%RFIELD(1),SARG)
          case('I')
            call FLOAT2STR(1.D0*myCOMMANDS(IOBJ)%FIELDS(ID)%IFIELD(1),SARG)
          case('VRANGE')
            SID=myCOMMANDS(IOBJ)%FIELDS(ID)%VFIELD(1)%ID
            call FLOAT2STR(1.D0*myCOMMANDS(IOBJ)%FIELDS(ID)%VFIELD(1)%VMI,SMI)
            call FLOAT2STR(1.D0*myCOMMANDS(IOBJ)%FIELDS(ID)%VFIELD(1)%VMA,SMA)
            call FLOAT2STR(1.D0*myCOMMANDS(IOBJ)%FIELDS(ID)%VFIELD(1)%VINC,SINC)
            SARG=trim(SID)//' '//trim(SMI)//' '//trim(SMA)//' '//trim(SINC)
          case default
            call send_message('SR',PARID,TID)
          end select
        else
          call send_message('P',PARID,' ')
        endif
      else
        call send_message('C',CMD,' ')
      endif
      end subroutine getCmdParamS

C---------------------------------------------------------
      subroutine setCmdParamS(CMD,PARID,SARG)
! set string value of command parameter
C---------------------------------------------------------
      character(*),intent(in) :: CMD
      character(*),intent(in) :: PARID
      character(*),intent(in) :: SARG
      integer :: ID,IOBJ,ires
      logical :: IsFloat,IsInteger
      character(LEN_NAME) :: tid,SID
      character(256) :: AUX
      REAL(KIND(1.D0)) :: res,A(3)
      integer :: IS,IL,N
      IOBJ=getCmdIndex(CMD)
      if (IOBJ.gt.0) then
        ID=getParIndex(myCOMMANDS(IOBJ),PARID)
        if (ID.gt.0) then
          tid=myCOMMANDS(IOBJ)%FIELDS(ID)%TID
          select case(tid)
          case('S')
            myCOMMANDS(IOBJ)%FIELDS(ID)%SFIELD(1)=SARG
          case('F')
            if (IsFloat(SARG,res)) then
              myCOMMANDS(IOBJ)%FIELDS(ID)%RFIELD(1)=res
            endif
          case('I')
            if (IsInteger(SARG,ires)) then
              myCOMMANDS(IOBJ)%FIELDS(ID)%IFIELD(1)=ires
            endif
          case('VRANGE')
            IS=1
            CALL FINDPAR(SARG,1,IS,IL)
            if (IL.gt.0) then
              AUX=trim(SARG(IS+IL:))
              call STRCOMPACT(AUX,IS,IL)
              call STR2ARRAY8(SARG(IS:IS+IL-1),' ',A,3,N)
              if (N.ge.3) then
                myCOMMANDS(IOBJ)%FIELDS(ID)%VFIELD(1)%ID=SARG(IS:IS+IL-1)
                myCOMMANDS(IOBJ)%FIELDS(ID)%VFIELD(1)%VMI=A(1)
                myCOMMANDS(IOBJ)%FIELDS(ID)%VFIELD(1)%VMA=A(2)
                myCOMMANDS(IOBJ)%FIELDS(ID)%VFIELD(1)%VINC=A(3)
              endif
            endif
          case default
            call send_message('SA',PARID,TID)
          end select
        else
          call send_message('P',PARID,' ')
        endif
      else
        call send_message('C',CMD,' ')
      endif
      end subroutine setCmdParamS

C---------------------------------------------------------
      SUBROUTINE COMMANDS_INP_R(OBJ,PNAME,ARG,NARG)
C---------------------------------------------------------
      TYPE(TCOMMAND) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(in) :: ARG(*)
      integer,intent(out) :: NARG
      character(LEN_ID) :: TID
      integer :: IDX,NA,i
      NARG=0
    !  write(*,*) OBJ%ICMD,' COMMANDS_INP_R '//trim(OBJ%ID),' ',trim(PNAME)
      call COMMANDS_PARINFO(OBJ,PNAME,IDX,TID,NA)
      write(*,*) OBJ%ICMD,' COMMANDS_INP_R '//trim(OBJ%ID),'.',trim(PNAME),IDX,trim(TID),NA
      if ((IDX.gt.0).and.(NA.gt.0)) then
        select case(TID)
        case('F')
          DO i=1,NA
            OBJ%FIELDS(IDX)%RFIELD(i)=ARG(i)
          enddo
          NARG=NA
        case('I')
          do i=1,NA
            OBJ%FIELDS(IDX)%IFIELD(i)=NINT(ARG(i))
          enddo
      !    write(*,*) '    ARG=',OBJ%FIELDS(IDX)%IFIELD(1:NA)
          NARG=NA
        end select
      endif
      END SUBROUTINE COMMANDS_INP_R

C---------------------------------------------------------
      SUBROUTINE COMMANDS_INP_V(OBJ,PNAME,ARG,NARG)
C---------------------------------------------------------
      TYPE(TCOMMAND) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(VRANGE), intent(in) :: ARG(*)
      integer,intent(out) :: NARG
      character(LEN_ID) :: TID
      integer :: IDX,NA
      NARG=0
      call COMMANDS_PARINFO(OBJ,PNAME,IDX,TID,NA)
      if ((IDX.gt.0).and.(NA.gt.0)) then
        select case(TID)
        case('VRANGE')
          OBJ%FIELDS(IDX)%VFIELD(1:NA)=ARG(1:NA)
          NARG=NA
        end select
      endif
      END SUBROUTINE COMMANDS_INP_V


C---------------------------------------------------------
      SUBROUTINE COMMANDS_OUT_R(OBJ,PNAME,ARG,NARG)
C---------------------------------------------------------
      TYPE(TCOMMAND) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(out) :: ARG(*)
      integer,intent(out) :: NARG
      character(LEN_NAME) :: TID
      integer :: IDX,NA,i
      NARG=0
  !    write(*,*) 'COMMANDS_OUT_R   FIELDS(7)=',OBJ%FIELDS(7)%IFIELD(1)
      call COMMANDS_PARINFO(OBJ,PNAME,IDX,TID,NA)
  !    write(*,*) 'COMMANDS_OUT_R   FIELDS(7)=',OBJ%FIELDS(7)%IFIELD(1)
  !    write(*,*) OBJ%ICMD,' COMMANDS_OUT_R '//trim(OBJ%ID),' ',trim(PNAME),IDX,trim(TID),NA

      if ((IDX.gt.0).and.(NA.gt.0)) then
        select case(TID)
        case('F')
          do i=1,NA
            ARG(i)=OBJ%FIELDS(IDX)%RFIELD(i)
          enddo
          NARG=NA
        case('I')
          do i=1,NA
            ARG(i)=OBJ%FIELDS(IDX)%IFIELD(i)
          enddo
  !        write(*,*) '    ARG=',OBJ%FIELDS(IDX)%IFIELD(1:NA)
          NARG=NA
        end select
      endif
      END SUBROUTINE COMMANDS_OUT_R

C---------------------------------------------------------
      SUBROUTINE COMMANDS_OUT_V(OBJ,PNAME,ARG,NARG)
C---------------------------------------------------------
      TYPE(TCOMMAND) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(VRANGE), intent(out) :: ARG(*)
      integer,intent(out) :: NARG
      character(LEN_NAME) :: TID
      integer :: IDX,NA,i
      NARG=0
      call COMMANDS_PARINFO(OBJ,PNAME,IDX,TID,NA)
      if ((IDX.gt.0).and.(NA.gt.0)) then
        select case(TID)
        case('VRANGE')
          ARG(1:NA)=OBJ%FIELDS(IDX)%VFIELD(1:NA)
          NARG=NA
        end select
      endif
      END SUBROUTINE COMMANDS_OUT_V

C---------------------------------------------------------
      SUBROUTINE COMMANDS_INP_S(OBJ,PNAME,SARG,NARG)
C---------------------------------------------------------
      TYPE(TCOMMAND) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      CHARACTER(*),intent(in) :: SARG
      integer,intent(out) :: NARG
      character(LEN_NAME) :: TID
      integer :: IDX,NA
      NARG=0
    !  write(*,*) 'COMMANDS_INP_S '//trim(OBJ%ID),' ',OBJ%NARG
      call COMMANDS_PARINFO(OBJ,PNAME,IDX,TID,NA)
    !  write(*,*) '      TID= '//trim(TID),' ',IDX
      if (trim(TID).eq.'S') then
        if (IDX.gt.0) then
          OBJ%FIELDS(IDX)%SFIELD(1)=trim(SARG)
          NARG=NA
        endif
      endif
      END SUBROUTINE COMMANDS_INP_S

C---------------------------------------------------------
      SUBROUTINE COMMANDS_OUT_S(OBJ,PNAME,SARG,NARG)
C---------------------------------------------------------
      TYPE(TCOMMAND) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      CHARACTER(*),intent(out) :: SARG
      integer,intent(out) :: NARG
      character(LEN_NAME) :: TID
      integer :: IDX,NA
      SARG=' '
      NARG=0
      call COMMANDS_PARINFO(OBJ,PNAME,IDX,TID,NA)
      if (trim(TID).eq.'S') then
        if (IDX.gt.0) then
          SARG=OBJ%FIELDS(IDX)%SFIELD(1)
          NARG=NA
        endif
      endif
      END SUBROUTINE COMMANDS_OUT_S

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
        call MSG_WARN('COMMANDS_ASSIGN: wrong ICLS for '//trim(ID)//':'//trim(tid),1)
      case('ASSIGN_ICMD')
        call MSG_WARN('COMMANDS_ASSIGN: wrong ICMD for '//trim(ID)//':'//trim(tid),1)
      case('ADD_FIELD')
        call MSG_WARN('COMMANDS_ASSIGN: cannot add field '//trim(ID)//':'//trim(tid),1)
      case default
        call MSG_WARN('undefined error: '//trim(what),1)
      end select

      !READ(*,*)
      end subroutine send_message

      end module COMMANDS
