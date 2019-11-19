!//////////////////////////////////////////////////////////////////////
!////  $Id $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.45 $
!////     $Date: 2019/08/15 15:02:09 $
!//////////////////////////////////////////////////////////////////////
!////
!////  XML handler for classes used in SIMRES
!////  General procedures, which apply to any class defined in CLASSES module
!////
!//////////////////////////////////////////////////////////////////////
      MODULE XMLHANDLER
      use xmlparse
      USE CLASSES
      use ENUMDEF
      use FIELDDEF
      use CLASSDEF
      ! use CLASSDATA
      USE FIELDDATA
      use VCTABLE
      USE XMLUTILS
      implicit none

      private
      save

!  startDEFAULT,dataDEFAULT and endDEFAULT handle all field data (tags inside classes)
!  endDEFAULT should be called also by end hadlers of classes => returns to the previous handler
!  All field ID's must be defined in XHND_NAMES, otherwise they are ignored


! variables that control parsing
      character(LEN_ID*16) :: XHANDLERS        ! handler IDs tree (: separated ID strings)
      character(LEN_NAME) :: XHND             ! actual handler ID
      character(DIM_ID*6) :: XHND_NAMES       ! current namespace for value tags
      integer :: XHND_ICLS                    ! id number of currently handled classs
    !  integer :: XHND_CDATA                   ! instance index of temporary subclass data
      character(LEN_ID+1) :: XHND_PREFIX      ! subclass prefix
      character(LEN_ID+1) :: XHND_SUBCLS      ! subclass prefix

      character(1024) :: XHND_ITEM                ! string for items
    !  integer :: XHND_LIST_LEN                ! declared length of list (number of following ITEM tags)
    !  integer :: XHND_LIST_SEL                ! selection item
      TYPE(TFIELDDEF),target :: XHND_FIELD      ! current parameter definition
      TYPE(TENUMDEF) :: XHND_ENUM               ! current enumerator type definition
      TYPE(TFIELD) :: XHND_ARG                 ! current field data
    !  TYPE(PCLSDEF) :: XHND_CDEF               ! temporary class definition
    !  TYPE(PCLSDEF) :: XHND_CDEF               ! temporary class definition

      integer :: ENUM_VALUE                   ! ord value of enumerated variable (-2 for other types)
      logical :: DOPARSING=.false.            ! flag for start/stop parsing
      logical :: REDEFINE=.false.             ! flag indicating definition of a new instrument instead of parameters update

      public startDEFAULT,dataDEFAULT,endDEFAULT,CLEAR_XHNDARG,CLEAR_XMLHANDLER
      public XHND,XHANDLERS,XHND_NAMES,XHND_ICLS
      public REDEFINE,DOPARSING,XHNDPass,XHNDReturn

      contains

!--------------------------------------------------------
      subroutine CLEAR_XMLHANDLER
! clear all handler data
!--------------------------------------------------------
        DOPARSING=.false.
        XHND_NAMES=' '
        XHND_ICLS=0
        XHANDLERS=' '
        XHND=' '
        XHND_PREFIX=' '
        XHND_SUBCLS=' '
        call CLEAR_XHNDARG
      end subroutine CLEAR_XMLHANDLER

!--------------------------------------------------------
      subroutine CLEAR_XHNDARG
! clear field data for the handler
!--------------------------------------------------------
      call DISPOSE_FIELD(XHND_ARG)
      XHND_ARG%DEF=>null()
      XHND_ITEM =' '
      ENUM_VALUE=-2
      call GetFieldDef(0,' ',XHND_FIELD)
      call GetEnumDef(0,XHND_ENUM)
      end subroutine CLEAR_XHNDARG

!--------------------------------------------------------
      subroutine XHNDPass(newHND)
! pass task to a new handler
!--------------------------------------------------------
      character(*) :: newHND
        call PREPENDPAR(XHANDLERS,':',trim(XHND))
        XHND_NAMES=' '
        XHND_PREFIX=' '
        XHND_SUBCLS=' '
        XHND=adjustl(newHND)
      end subroutine XHNDPass

!--------------------------------------------------------
      subroutine XHNDReturn
! Return tasks to the previous handler
!--------------------------------------------------------
      integer :: IS,IL
        call FINDSTRPAR(XHANDLERS,':',1,IS,IL)
        if (IL.gt.0) then
          XHND=XHANDLERS(IS:IS+IL-1)
          XHANDLERS=trim(XHANDLERS(IS+IL+1:))
        else
          XHND=' '
        endif
        XHND_NAMES=' '
        XHND_PREFIX=' '
        XHND_SUBCLS=' '
      end subroutine XHNDReturn

!--------------------------------------------------------
      subroutine startVALUES(attribs,error)
! handle START of a value tag: handle attributes and allocate XHND_ARG
! subclass dot delimiter is supported
! Assumes XHND_FIELD has been defined in startDEFAULT
!--------------------------------------------------------
      character(len=*), dimension(:,:)  :: attribs
      logical                           :: error
      character(LEN_NAME) :: CNUM
      character(2*LEN_ID+1) :: clstag
      integer :: i
      logical :: IsInteger
      ! get enumerator definition (if XHND_FIELD is enumerated type)
        call GetEnumDef(XHND_FIELD%ENUM,XHND_ENUM)
        ENUM_VALUE=-1
        clstag=trim(XHND_PREFIX)//trim(XHND_FIELD%ID)
      ! associate field definition with XHND_ARG
        XHND_ARG%DEF=>XHND_FIELD
        select case (trim(XHND_FIELD%TID))
          case(ID_FIELD_ENUM)
          ! NOTE: value attribute is not obligatory
            call readAttrib(attribs,'value',CNUM,error)
            if ((.not.error).and.(IsInteger(CNUM,i))) ENUM_VALUE=i
            error=(.not.ALLOC_FIELD(XHND_ARG,XHND_FIELD%DIM))
          case(ID_FIELD_SELECT)
            call readAttrib(attribs,'length',CNUM,error)
            if ((.not.error).and.(IsInteger(CNUM,i))) then
              error=(.not.ALLOC_FIELD(XHND_ARG,i))
            else
              call MSG_ERROR('XMLHANDLER','missing length attribute for '//trim(clstag),0,1)
              error=.true.
              return
            endif
            call readAttrib(attribs,'selected',CNUM,error)
            if ((.not.error).and.(IsInteger(CNUM,i))) XHND_ARG%EORD=i
          case(ID_FIELD_RANGE)
            call readAttrib(attribs,'length',CNUM,error)
            if ((.not.error).and.(IsInteger(CNUM,i))) then
              error=(.not.ALLOC_FIELD(XHND_ARG,i))
            else
              call MSG_ERROR('XMLHANDLER','missing length attribute for '//trim(clstag),0,1)
              error=.true.
              return
            endif
          case(ID_FIELD_TABLE)
            call readAttrib(attribs,'rows',CNUM,error)
            if ((.not.error).and.(IsInteger(CNUM,i))) then
              error=(.not.ALLOC_TABLE(XHND_ARG,i))
              XHND_ARG%NP=0
            else
              call MSG_ERROR('XMLHANDLER','missing rows attribute for '//trim(clstag),0,1)
              error=.true.
              return
            endif
          case default
            error=(.not.ALLOC_FIELD(XHND_ARG,XHND_FIELD%DIM))
        end select
      end subroutine startVALUES

!--------------------------------------------------------
      subroutine endVALUES(OBJ,tag,error)
! handle end of a value tag (already includes XML_PREFIX)
! send data to OBJ and clear handler data
!--------------------------------------------------------
      type(TCLASS) :: OBJ
      character(len=*)                  :: tag
      logical                           :: error
      integer :: LR
        select case (trim(XHND_FIELD%TID))
          ! apply limits to selection index
          case(ID_FIELD_SELECT)
            XHND_ARG%EORD=min(XHND_ARG%EORD,XHND_ARG%NP)
            XHND_ARG%EORD=max(XHND_ARG%EORD,-1)
        end select
        !if (trim(tag).eq.'SEG') then
!1         format(a,6(1x,G12.5))
          !write(*,1) 'endVALUES '//trim(XHND_ARG%DEF%TID),XHND_ARG%NP,XHND_ARG%LENGTH,XHND_ARG%MAXROW
          !ncol=XHND_ARG%DEF%NCOL
          !do i=1,XHND_ARG%NP
          !  write(*,1)  'endVALUES ',i,XHND_ARG%RFIELD((i-1)*ncol+1)
          !enddo
        !endif

        call PCLASS_INP(OBJ,trim(tag),XHND_ARG,LR)
        call CLEAR_XHNDARG
      end subroutine endVALUES

!--------------------------------------------------------
      subroutine dataENUM(xmldata)
! enum data handler
! check consistence of ORD and string values
! If string value is defined => use it and get corresponding ENUM_VALUE
! otherwise use ordinal value given as atribute
!--------------------------------------------------------
      character(len=*)                  :: xmldata
      integer :: IE
      character(LEN_NAME) :: CNUM
      character(2*LEN_ID+1) :: clstag
      character(FIELD_BUFFER_STR) :: strval,LINE
      logical :: dbg = .false.
      !  dbg=(trim(XHND_ENUM%ID).eq.'SIGN')
        strval=trim(adjustl(xmldata))
        clstag=trim(XHND_PREFIX)//trim(XHND_FIELD%ID)
      ! defined string for SIGN is '+1', not '1', but '1' should be legal
        if ((trim(XHND_ENUM%ID).eq.'SIGN').and.(trim(strval).eq.'1')) strval='+1'
        call GetEnumOrd(XHND_FIELD%ENUM,strval,IE)
! string value is not defined for given enumerated type => use ordinal value from attributes
        if (IE.lt.0) then
          call INT2STR(XHND_FIELD%ENUM,CNUM)
          LINE=trim(XHND_ENUM%ID)//' ['//trim(strval)//'] id='//trim(clstag)//' enum='//trim(CNUM)
          call MSG_WARN('XMLHANDLER: undefined enumerated value, '//trim(LINE),1)
          call GetEnumString(XHND_ENUM%ID,ENUM_VALUE,strval)
          call MSG_INFO('Assumed value is '//trim(strval),1)
! string and ordinal values do not agree => prefer string value
        else if (ENUM_VALUE.NE.IE) then
        ! warning if value attribute is given, but does not agree with string value
          if (ENUM_VALUE.ge.0) then
            call INT2STR(ENUM_VALUE,CNUM)
            call STRIM(trim(XHND_ENUM%ID)//' index='//trim(CNUM)//' value='//trim(strval),LINE)
            !LINE=trim(XHND_ENUM%ID)//' index='//trim(CNUM)//' value='//trim(strval)
            call MSG_WARN('XMLHANDLER: incompatible enumerator index and value: '//trim(LINE),1)
            call MSG_INFO('Assumed value is ['//trim(strval)//']',1)
          endif
          ENUM_VALUE=IE
        endif
        if (dbg) write(*,*) 'dataENUM ',trim(XHND_FIELD%ID),' [',trim(strval),']  ord=',IE,' value=',ENUM_VALUE
      end subroutine dataENUM

!--------------------------------------------------------
      subroutine dataVALUES(xmldata,error)
! default data handler
!--------------------------------------------------------
      character(len=*)                  :: xmldata
      logical                           :: error
      !integer,parameter :: MAXA = 128
      integer :: NA
      REAL(KIND(1.D0)) :: BUF(FIELD_BUFFER_FLOAT)
      logical :: dbg = .false.
      !  dbg=(trim(XHND_ENUM%ID).eq.'SIGN')
! get parameter values (ARG) and number (NA) for each value type

        select case (trim(XHND_FIELD%TID))
        case(ID_FIELD_ENUM)
          call dataENUM(xmldata)
          BUF(1)=ENUM_VALUE
          call ARRAY2FIELD(XHND_ARG,0,BUF,FIELD_BUFFER_FLOAT,NA)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_STR)
          call STR2FIELD(XHND_ARG,0,trim(xmldata),NA)
        case(ID_FIELD_RANGE,ID_FIELD_SELECT,ID_FIELD_CLASS,ID_FIELD_TABLE)
      ! nothing to do
        case default
          if (len_trim(XHND_FIELD%TID).gt.0) then
            call MSG_ERROR('XMLHANDLER','unknown type for '//trim(XHND_FIELD%ID)//' ['//trim(XHND_FIELD%TID)//']',0,1)
            error=.true.
            XHND_ARG%DEF=>null()
          endif
        end select
      end subroutine dataVALUES


!--------------------------------------------------------
      subroutine startDEFAULT(tag,attribs,error)
! default handler of start tags (assume value tags)
! 1) handle ITEM
! 2) get field definition from tag
! 3) set XHND_PREFIX for class fields, call startVALUES for normal fields
!--------------------------------------------------------
      character(len=*)                  :: tag
      character(len=*), dimension(:,:)  :: attribs
      logical                           :: error
      integer :: IVAL
      character(2*LEN_ID+1) :: clstag

! do not handle fields if calling object is not defined
      if (XHND_ICLS.le.0) return

      select case(trim(tag))
      case('ITEM')
          XHND_ITEM=' '
      case default
          call CLEAR_XHNDARG
    ! ignor parent class tags
          IVAL=isParent(XHND_ICLS,tag)
          if (IVAL.gt.0) return
        !  if (trim(tag).eq.'SAMPLE') write(*,*) 'startDEFAULT result=',IVAL
    ! get tag including subclass prefix
          clstag=trim(XHND_PREFIX)//trim(tag)
        !  write(*,*) 'startDEFAULT [',trim(clstag),']'
        !  if (INDEX(tag,'VFOC').gt.0) write(*,*) 'startDEFAULT [',trim(clstag),']'
          call GetFieldDef(XHND_ICLS,trim(clstag),XHND_FIELD)
          if (XHND_FIELD%IDX.le.0) then
            write(*,*) 'startDEFAULT ['//trim(clstag)//']',IVAL,XHND_ICLS

            call MSG_ERROR('XMLHANDLER','undefined type for field '//trim(clstag),0,1)
            error=.true.
            return
          endif
          select case (trim(XHND_FIELD%TID))
          case(ID_FIELD_CLASS)
        ! field=class
            XHND_PREFIX=trim(XHND_FIELD%ID)//"."
            XHND_SUBCLS=trim(XHND_FIELD%ID)
          !  write(*,*) 'startDEFAULT ['//trim(XHND_SUBCLS)//']'
        ! other field types
          case default
          ! check that tag is a member of current namespace
            call GETPARINDX(XHND_NAMES,':',trim(clstag),IVAL)
          !  write(*,*) 'startVALUES for ',trim(XHND_FIELD%ID),' OWNER=',XHND_FIELD%OWNER,' ',IVAL
            if (IVAL.gt.0) call startVALUES(attribs,error)
          end select
      end select
      end subroutine startDEFAULT

!--------------------------------------------------------
      subroutine endDEFAULT(OBJ,tag,error)
! defaullt handler for end tags (value tags or return handler to caller)
! handle also ITEM tags (list and range types)
!--------------------------------------------------------
      type(TCLASS) :: OBJ
      character(len=*)                  :: tag
      logical                           :: error
      integer :: LR,IVAL
    ! root tag
    !  if (INDEX(tag,'NSEG').gt.0) write(*,*) 'endDEFAULT [',trim(XHND_PREFIX)//trim(tag),']',trim(XHND_SUBCLS)
    ! ignor parent class tags
      if (isParent(XHND_ICLS,tag).gt.0) return
      !write(*,*) 'endDEFAULT tag=',trim(tag),' ',trim(XHND)
! standard behaviour => only return to handling of previous tag
      if (trim(tag).eq.trim(XHND)) then
        call XHNDReturn
        return
      endif
      !if (trim(tag).eq.'ITEM') write(*,*) 'endDEFAULT 0'
! do not handle end tags if OBJ is not defined
      if (.not.associated(OBJ%DEF)) return
      !if (trim(tag).eq.'ITEM') write(*,*) 'endDEFAULT 1'
! handle end of class fields  => clear prefix
      if ((len_trim(XHND_PREFIX).gt.0).and.(trim(tag).eq.trim(XHND_SUBCLS))) then
        XHND_PREFIX=' '
        XHND_SUBCLS=' '
        call CLEAR_XHNDARG
        return
      endif
      !if (trim(tag).eq.'ITEM') write(*,*) 'endDEFAULT 2'
! do not handle end tags if XHND_ARG is not defined
      if (.not.associated(XHND_ARG%DEF)) return
      !if (trim(tag).eq.'ITEM') write(*,*) 'endDEFAULT 3 '//trim(XHND_FIELD%TID)
      select case(trim(tag))
! handle ITEM
      case('ITEM')
        select case (trim(XHND_FIELD%TID))
        case(ID_FIELD_TABLE)
          !write(*,*) 'endDEFAULT 4 ',XHND_ARG%NP,XHND_ARG%MAXROW
          if (XHND_ARG%NP.lt.XHND_ARG%MAXROW) then
            XHND_ARG%NP=XHND_ARG%NP+1
            !write(*,*) 'endDEFAULT ITEM=',XHND_ARG%NP,' ',trim(XHND_ITEM)
            call STR2ROW(XHND_ARG,XHND_ARG%NP,0,trim(XHND_ITEM),LR)
          endif
        case(ID_FIELD_SELECT,ID_FIELD_RANGE)
          if (XHND_ARG%NP.lt.XHND_ARG%LENGTH) then
            XHND_ARG%NP=XHND_ARG%NP+1
            call ITEM2FIELD(XHND_ARG,XHND_ARG%NP,trim(XHND_ITEM))
          endif
        end select
      case default
! handle field tags
    ! process only defined identifiers
        call GETPARINDX(XHND_NAMES,':',trim(XHND_PREFIX)//trim(tag),IVAL)
      ! write(*,*) 'endDEFAULT ['//trim(XHND_SUBCLS)//']'
        if (IVAL.gt.0) call endVALUES(OBJ,trim(XHND_PREFIX)//trim(tag),error)
      end select
      end subroutine endDEFAULT

!--------------------------------------------------------
      subroutine dataDEFAULT(tag,xmldata,error)
! default data handler
!--------------------------------------------------------
      character(len=*)                  :: tag
      character(len=*)                  :: xmldata
      logical                           :: error
      integer :: IVAL
! do not handle empty data
      if (len_trim(xmldata).le.0) return
! do not handle data if XHND_ARG is not defined
      if (.not.associated(XHND_ARG%DEF)) return

      select case (trim(tag))
      case('ITEM')
        select case (trim(XHND_FIELD%TID))
        case(ID_FIELD_RANGE,ID_FIELD_SELECT,ID_FIELD_TABLE)
        !write(*,*) 'dataDEFAULT ITEM=',trim(xmldata)
          XHND_ITEM=trim(xmldata)
        end select
      case default
      ! process only defined identifiers
        call GETPARINDX(XHND_NAMES,':',trim(XHND_PREFIX)//trim(tag),IVAL)
        if (XHND_ICLS.le.8) then
        !  write(*,*) 'dataVALUES for ',trim(XHND_FIELD%ID),' OWNER=',XHND_FIELD%OWNER,' ',trim(xmldata),IVAL
        endif
        if (IVAL.gt.0) call dataVALUES(xmldata,error)
      end select
      end subroutine dataDEFAULT

      end module XMLHANDLER
