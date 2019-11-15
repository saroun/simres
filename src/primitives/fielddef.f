!//////////////////////////////////////////////////////////////////////
!////  $Id: fielddef.f,v 1.6 2012/04/09 22:35:04 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron spectrometers
!////
!////     Copyright (C) 1995-2009, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.6 $
!////     $Date: 2012/04/09 22:35:04 $
!//////////////////////////////////////////////////////////////////////
!////
!////   Defines basic types used by RESTRAX
!////   Includes XML handlers for field definitions in classes.xml
!////
!//////////////////////////////////////////////////////////////////////
      MODULE FIELDDEF
      use XMLUTILS
      use XMLINFO
      use ENUMDEF
!      use FILETOOLS!

      implicit none

      private
      save

      integer, parameter :: LEN_ID=16           ! length of ID strings
      integer, parameter :: LEN_NAME = 32       ! length of name attributes
      !integer, parameter :: LEN_LINE=1024        ! default length of input line

    ! values of TFIELDDEF%TCLS
    !  integer,parameter :: IDX_FIELD_INT=1
    !  integer,parameter :: IDX_FIELD_FLOAT=2
    !  integer,parameter :: IDX_FIELD_STR=3
    !  integer,parameter :: IDX_FIELD_RANGE=4
    !  integer,parameter :: IDX_FIELD_ENUM=5
    !  integer,parameter :: IDX_FIELD_SELECT=6
    !  character(34),parameter :: ID_FIELD_NAMES='INT:FLOAT:STRING:RANGE:ENUM:SELECT'

    ! values of TFIELDDEF%TID
      character(3),parameter :: ID_FIELD_INT='INT'
      character(5),parameter :: ID_FIELD_FLOAT='FLOAT'
      character(6),parameter :: ID_FIELD_STR='STRING'
      character(5),parameter :: ID_FIELD_RANGE='RANGE'
      character(4),parameter :: ID_FIELD_ENUM='ENUM'
      character(6),parameter :: ID_FIELD_SELECT='SELECT'
      character(8),parameter :: ID_FIELD_CLASS='CLASSOBJ'
      character(5),parameter :: ID_FIELD_TABLE='TABLE'

      TYPE TFIELDDEF
        integer :: OWNER            ! ICLS of the owner (TCLSDEF type)
        INTEGER :: IDX=0            ! type index in the owner's fields list
        character(LEN_ID) :: ID     ! field ID
        character(LEN_NAME) :: NAME       ! name string
        character(LEN_NAME) :: UNITS      ! units string
        character(LEN_ID) :: TID    ! type ID string
        character(LEN_ID) :: CID    ! class type ID for class fields
        integer :: DIM              ! dimension (for arrays)
        integer :: NROW             ! number of rows (for tables)
        integer :: NCOL             ! number of columns (for tables)
        integer :: ENUM             ! pointer to enumerated type definition
        LOGICAL :: CAL              ! true for calculated fields, default=false
      end TYPE TFIELDDEF


! Intrinsic type VRANGE defines a range of a component parameter
      TYPE VRANGE
        CHARACTER(2*LEN_ID+1) :: ID ! parameter ID in format COMPONENT.PAR
        REAL(KIND(1.D0)) :: VMI     ! minimum value
        REAL(KIND(1.D0)) :: VMA     ! maximum value
        REAL(KIND(1.D0)) :: VINC    ! increment
      end type VRANGE

! encapsulated pointer to TFIELDDEF
      type PFIELDDEF
        TYPE(TFIELDDEF),pointer :: F
  !      TYPE(TFIELDDEF),allocatable :: FA(:) ! field definitions
      end TYPE PFIELDDEF

! private instances for XML import
      type(PFIELDDEF) :: myFIELDDEF
      TYPE(TFIELDDEF) :: ZERO_FIELD

! flag that enables XML parsing

      public LEN_NAME,LEN_ID,LEN_LINE
      public PFIELDDEF,TFIELDDEF
      public DISPOSE_TMP_TYPES, CLEAR_TYPES,ALLOC_FIELDDEF,DISPOSE_FIELDDEF,startFIELDDEF
    !  public IDX_FIELD_INT,IDX_FIELD_FLOAT,IDX_FIELD_STR,IDX_FIELD_RANGE,IDX_FIELD_ENUM,IDX_FIELD_SELECT
      public ID_FIELD_INT,ID_FIELD_FLOAT,ID_FIELD_STR,ID_FIELD_RANGE,ID_FIELD_ENUM
      public ID_FIELD_SELECT,ID_FIELD_CLASS,iD_FIELD_TABLE
      public myFIELDDEF,GET_EMPTY_FIELD
      contains

!----------------------------------------------
      TYPE(TFIELDDEF) function GET_EMPTY_FIELD()
!----------------------------------------------
      GET_EMPTY_FIELD=ZERO_FIELD
      end function GET_EMPTY_FIELD

!////////////////////////////////////////////////////////////
!/////                                                  /////
!/////                 PRIVATE AREA                     /////
!/////                                                  /////
!////////////////////////////////////////////////////////////

!---------------------------------------------------------
      subroutine ALLOC_FIELDDEF(F,ID,NAME)
!---------------------------------------------------------
        character(*) :: ID
        character(*) :: NAME
        TYPE(PFIELDDEF) :: F
        integer :: ierr
        allocate(F%F,STAT=ierr)
        if(ierr.eq.0) then
          call CLEAR_FIELDDEF(F%F)
          F%F%ID=trim(ID)
          F%F%NAME=trim(NAME)
        endif
      end subroutine ALLOC_FIELDDEF

!---------------------------------------------------------
      subroutine CLEAR_FIELDDEF(F)
!---------------------------------------------------------
        TYPE(TFIELDDEF),intent(out) :: F
          F%ID=' '
          F%NAME=' '
          F%UNITS=' '
          F%DIM=0
          F%TID=' '
          F%CID=' '
          F%IDX=0
          F%ENUM=0
          F%OWNER=0
          F%CAL=.false.
          F%NROW=0
          F%NCOL=0
      end subroutine CLEAR_FIELDDEF

!---------------------------------------------------------
      subroutine DISPOSE_FIELDDEF(F)
!---------------------------------------------------------
        TYPE(PFIELDDEF) :: F
        if (associated(F%F)) then
          DEallocate(F%F)
          nullify(F%F)
        endif
      end subroutine DISPOSE_FIELDDEF

!--------------------------------------------------------
      subroutine CLEAR_TYPES
! clear definitions of all types
!--------------------------------------------------------
        call CLEAR_FIELDDEF(ZERO_FIELD)
        call DISPOSE_FIELDDEF(myFIELDDEF)
      end subroutine CLEAR_TYPES

!--------------------------------------------------------
      subroutine DISPOSE_TMP_TYPES
! dispose temporary instaces of type definitions
!--------------------------------------------------------
        call DISPOSE_FIELDDEF(myFIELDDEF)
      end subroutine DISPOSE_TMP_TYPES

!--------------------------------------------------------
      subroutine startFIELDDEF(CLSID,tag,attribs,error)
!--------------------------------------------------------
        character(len=*)                  :: CLSID
        character(len=*)                  :: tag
        character(len=*), dimension(:,:)  :: attribs
        logical                           :: error
        integer :: isize
        character(LEN_ID)  :: PARAM_ID,cid
        character(LEN_NAME) :: uname,ename,pname
          if (.not.associated(myFIELDDEF%F)) then
            select case(trim(tag))
          ! skip some tags
            case('GROUP','ENUMITEM','TRANSLATE')
              return
            case default
          ! read common attributes to all fields except ENUMITEM
              call readAttrib(attribs,'id',PARAM_ID,error)
              call MKUPCASE(PARAM_ID)
              if (.not.error) call readAttrib(attribs,'name',pname,error)
              if (.not.error) call ALLOC_FIELDDEF(myFIELDDEF,trim(PARAM_ID),trim(pname))
            end select
          ! field should be allocated now
            if ((error).or.(.not.associated(myFIELDDEF%F))) then
              call MSG_ERROR('CLASSES','cannot create new field definition: '//trim(trim(PARAM_ID)),0,1)
              error=.true.
              return
            endif
        ! flag for calculated fields
            myFIELDDEF%F%CAL=readAttribBoolean(attribs,'calc',.false.)
        ! handle particular field types
            select case(trim(tag))
            case(ID_FIELD_FLOAT)
              myFIELDDEF%F%TID=ID_FIELD_FLOAT
              call readAttrib(attribs,'units',uname,error)
              isize=readAttribInt(attribs,'size',1)
              myFIELDDEF%F%DIM=isize
              myFIELDDEF%F%UNITS=trim(uname)
            case(ID_FIELD_INT)
              myFIELDDEF%F%TID=ID_FIELD_INT
              isize=readAttribInt(attribs,'size',1)
              myFIELDDEF%F%DIM=isize
            case(ID_FIELD_STR)
              myFIELDDEF%F%TID=ID_FIELD_STR
              isize=readAttribInt(attribs,'size',1)
              myFIELDDEF%F%DIM=isize
            case(ID_FIELD_RANGE)
              myFIELDDEF%F%TID=ID_FIELD_RANGE
              isize=readAttribInt(attribs,'size',1)
              myFIELDDEF%F%DIM=isize
            case(ID_FIELD_SELECT)
              myFIELDDEF%F%TID=ID_FIELD_SELECT
              myFIELDDEF%F%DIM=1
            case(ID_FIELD_TABLE)
              myFIELDDEF%F%TID=ID_FIELD_TABLE
              isize=readAttribInt(attribs,'cols',1)
              myFIELDDEF%F%NCOL=isize
              isize=readAttribInt(attribs,'rows',1)
              myFIELDDEF%F%NROW=isize
              myFIELDDEF%F%DIM=myFIELDDEF%F%NROW*myFIELDDEF%F%NCOL
            case(ID_FIELD_CLASS)
              myFIELDDEF%F%TID=ID_FIELD_CLASS
              call readAttrib(attribs,'cid',cid,error)
              if (.not.error)  myFIELDDEF%F%CID=trim(cid)
              myFIELDDEF%F%DIM=1
            case(ID_FIELD_ENUM)
          ! enum definition is embeded in the class => read the definition now
              ! derive name from class name and enum. id
              ename=trim(CLSID)//'_'//adjustl(PARAM_ID)
              call startENUMDEF(ename,error)
              myFIELDDEF%F%TID=ID_FIELD_ENUM
              myFIELDDEF%F%DIM=1 ! arrays of enumerators are not possible
            case default
          ! field type can be one of already defined enumerators
              if (isEnumDefined(trim(tag))) then
                myFIELDDEF%F%ENUM=GET_ENUMINDEX(trim(tag))
                myFIELDDEF%F%TID=ID_FIELD_ENUM
                myFIELDDEF%F%DIM=1
              else
                call MSG_ERROR('startFIELDDEF','unknown field type: '//trim(tag),0,1)
                error=.true.
                return
              endif
            end select
            if (error) call MSG_ERROR('startFIELDDEF','wrong or missing attributes for tag '//trim(tag),0,1)
          endif
      end subroutine startFIELDDEF



      end module FIELDDEF
