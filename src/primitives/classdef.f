!//////////////////////////////////////////////////////////////////////
!////  $Id: classdef.f,v 1.11 2019/08/16 17:16:26 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron spectrometers
!////
!////     Copyright (C) 1995-2009, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.11 $
!////     $Date: 2019/08/16 17:16:26 $
!//////////////////////////////////////////////////////////////////////
!////
!////   Definition of TCLSDEF ... generic class definition
!////   Includes XML handlers for class definitions in classes.xml
!////
!//////////////////////////////////////////////////////////////////////
      MODULE CLASSDEF
      use CLASSES
      use FIELDDEF
      use ENUMDEF
      use XMLUTILS
      use XMLINFO
      use FILETOOLS

      implicit none

      private
      save

      integer, parameter :: MAX_CLASSES=64      ! max. number of classes
      integer, parameter :: MAX_FIELDS=32       ! max. number of fields in one class
      integer, parameter :: DIM_ID=MAX_FIELDS*(LEN_ID+1)  ! length of the string with field IDs for one class

! common type for all commands
      type TCLSDEF
        INTEGER :: IDX=0                       ! class index (see global definitions in CLASSES module)
        character(LEN_ID)     :: ID            ! ID string
        character(LEN_NAME)   :: NAME          ! name string
        CHARACTER(64) :: HINT                  ! description
        integer :: NF                          ! number of fields
        integer :: PARENT                      ! parent type index
        logical :: isCMD                       ! indicates that parent is COMMAND
        integer :: cType                       ! class type (e.g. TCLS_COM)
        TYPE(PFIELDDEF) :: FIELDS(MAX_FIELDS)  ! field definitions
      end TYPE TCLSDEF

! encapsulated pointer to TCLSDEF
      type PCLSDEF
        TYPE(TCLSDEF),pointer :: C
  !      TYPE(TCLSDEF),TARGET,allocatable :: CA(:) ! field definitions
      end TYPE PCLSDEF

! physical storage of classes definitions

      integer :: CLS_NUM=0
      TYPE(PCLSDEF) :: CLSDEF(MAX_CLASSES)

! private instances for XML import
      TYPE(PCLSDEF) :: myCLSDEF
      TYPE(TCLSDEF) ::  ZERO_CLS

! flag that enables XML parsing
      logical :: isParsing=.false.

      public MAX_FIELDS,DIM_ID

      interface GetClassDef
        module procedure GetClassDef_ID
        module procedure GetClassDef_IX
      end interface GetClassDef

      public startcls,datacls,endcls,CLEAR_CLASSES
      public new_report_classes
      public getICLS
      public GetClassParents,GetClassDef,GetFieldIndex,GetFieldDef,isParent
      public CLASSES_NAMESPACE,CLASSES_LENGTH
      public PCLSDEF,TCLSDEF
      public CLS_NUM

      contains

!--------------------------------------------------------
! get index in CLS_IDS from CLS index in given namespace (CLSNAMES)
!--------------------------------------------------------
      integer function getICLS(CCLS,CLSNAMES)
      integer,intent(in) :: CCLS
      character(*),intent(in) :: CLSNAMES
      integer :: IRES,IS,IL
      IRES=-1
      call FINDSTRPAR(CLSNAMES,':',CCLS,IS,IL)
      if (IL.gt.0) IRES=GET_CLSINDEX(CLSNAMES(IS:IS+IL-1))
      getICLS=IRES
      end function getICLS


!---------------------------------------------------------
      SUBROUTINE CLASSES_NAMESPACE(ICLS,NAMESPACE)
! return parameter namespace as : delimitted string
! including parent(s)
! exclude calculated fields
!---------------------------------------------------------
      integer,intent(in) :: ICLS
      character(*),intent(out) :: NAMESPACE
      integer :: NCLS,PCLS(1:10),i,j,k
      TYPE(TFIELDDEF),pointer :: FD
      character(LEN_ID+1) :: prefix
      TYPE(PCLSDEF) :: CSUB

      call GetClassParents(ICLS,PCLS,NCLS)
      NAMESPACE=' '
      if (NCLS.le.0) return
      i=NCLS
      do while (i.gt.0)
!        write(*,*) 'CLASSES_NAMESPACE ',i,' parent=',PCLS(i),' ',trim(CLSDEF(PCLS(i))%C%ID)
        do j=1,CLSDEF(PCLS(i))%C%NF
          FD => CLSDEF(PCLS(i))%C%FIELDS(j)%F
!          write(*,*) '    field=',trim(FD%ID),' cal=',FD%CAL
          if (trim(FD%TID).eq.ID_FIELD_CLASS) then
            prefix=trim(FD%ID)//"."
            call getClassDef(trim(FD%CID),CSUB)
            if (associated(CSUB%C)) then
              do k=1,CSUB%C%NF
                FD => CSUB%C%FIELDS(k)%F
                if (.not.FD%CAL) call APPENDPAR(NAMESPACE,':',trim(prefix)//trim(FD%ID))
              enddo
            endif
          else
            if (.not.FD%CAL) call APPENDPAR(NAMESPACE,':',trim(FD%ID))
          endif
        enddo
        i=i-1
      enddo
    !  write(*,*) 'CLASSES_NAMESPACE=[',trim(NAMESPACE),']'
      end SUBROUTINE CLASSES_NAMESPACE

!-------------------------------------------------------------
      SUBROUTINE GetFieldDef(ICLS,PNAME,OBJ)
! get field definition for given class ICLS and field name
! Handles dot convention for subclass field names (because GetFieldIndex does)
!-------------------------------------------------------------
      integer, intent(in) :: ICLS
      character(*) :: PNAME
      type(TFIELDDEF),intent(out) :: OBJ
      integer :: ires,OWNER
      ires=0
      if ((ICLS.gt.0).and.(ICLS.le.CLS_NUM)) then
        ires=GetFieldIndex(ICLS,PNAME,OWNER)
      !  if (INDEX(PNAME,'VFOC').gt.0) write(*,*) 'GetFieldDef [',trim(PNAME),']',ires,OWNER
        if (ires.gt.0)  then
          OBJ=CLSDEF(OWNER)%C%FIELDS(ires)%F
        endif
      endif
      if (ires.le.0) OBJ=GET_EMPTY_FIELD()
      end SUBROUTINE GetFieldDef

!-------------------------------------------------------------
      SUBROUTINE GetClassDef_ID(ID,OBJ)
! get class definition for given class ID
!-------------------------------------------------------------
      character(*) :: ID
      type(PCLSDEF) :: OBJ
      integer :: ICLS
      ICLS=GET_CLSINDEX(ID)
      call GetClassDef_IX(ICLS,OBJ)
      end SUBROUTINE GetClassDef_ID

!-------------------------------------------------------------
      SUBROUTINE GetClassDef_IX(ICLS,OBJ)
! get class definition for given class ICLS
!-------------------------------------------------------------
      integer, intent(in) :: ICLS
      type(PCLSDEF) :: OBJ
      if ((ICLS.gt.0).and.(ICLS.le.CLS_NUM)) then
        OBJ=CLSDEF(ICLS)
      else
        ! OBJ=myCLSDEF
        NULLIFY(OBJ%C)
      endif
      end SUBROUTINE GetClassDef_IX

!---------------------------------------------------------
      subroutine GetClassParents(ICLS,PLIST,NLIST)
! return integer array with intexes of all class parents
! the class itself is the first
! NLIST is the number of items
!---------------------------------------------------------
      integer,intent(in) :: ICLS
      integer,intent(inout) :: PLIST(:)
      integer,intent(out) :: NLIST
      integer :: IC
            NLIST=0
            if (ICLS.gt.0) then
                  if ((ICLS.le.CLS_NUM).and.associated(CLSDEF(ICLS)%C)) then
                        NLIST=1
                        PLIST(NLIST)=ICLS
                        IC=CLSDEF(ICLS)%C%PARENT
                        do while ((IC.gt.0).and.(NLIST.lt.SIZE(PLIST)))
                              if (associated(CLSDEF(IC)%C)) then
                                    NLIST=NLIST+1
                                    PLIST(NLIST)=CLSDEF(IC)%C%IDX
                                    IC=CLSDEF(IC)%C%PARENT
                              else
                                    exit
                              endif
                        endDO
                  endIF
            endif
      end subroutine GetClassParents

!---------------------------------------------------------
      integer function isParent(ICLS,clsName)
! return >0 if clsName is a name of any parent of the class ICLS.
! including itself.
! result=ICLS of the identified parent or 0
!---------------------------------------------------------
      integer,intent(in) :: ICLS
      character(*) :: clsName
      integer :: IC,RES
      logical :: LOG1
        RES=0
        IC=ICLS
        LOG1=.false.
        do while ((IC.gt.0).and.(IC.le.CLS_NUM).and.(.not.LOG1))
            if (associated(CLSDEF(IC)%C)) then
                  LOG1=(trim(clsName).eq.trim(CLSDEF(IC)%C%ID))
                  if (LOG1) RES=IC
                  IC=CLSDEF(IC)%C%PARENT
            else
                  exit
            endif
        endDO
        isParent=RES
      end function isParent


!---------------------------------------------------------
      integer function GetFieldIndex(ICLS,PNAME,OWNER)
! return index of a parameter of given class, and index of owner (member of CLSDEF array)
! includes parents' fields, in the order as returned by CLASSES_NAMESPACE
! i.e. PARENT1/PARENT2/.../CLS
! Handles dot convention for subclass field names.
!---------------------------------------------------------
      integer,intent(in) :: ICLS
      character(*),intent(in) :: PNAME
      integer,intent(out) :: OWNER
      integer :: IC,IX,i,k
      integer :: PCLS(1:10),NCLS
      type(PCLSDEF) :: COWN
      TYPE(TFIELDDEF),pointer :: FD
      if ((ICLS.le.0).or.(ICLS.gt.CLS_NUM)) then
        OWNER=0
        GetFieldIndex=0
        return
      endif
! collect the class tree indices in PCLS
      call GetClassParents(ICLS,PCLS,NCLS)
      IC=NCLS
      IX=0
      OWNER=0
    ! handle dot convention for subclass field names
      i=INDEX(PNAME,'.')
    ! scan given class and all parent classes
      DO WHILE ((IC.GT.0).and.(IX.le.0))
    ! has prefix => expect ClassField type
        if (i.gt.1) then
        ! get ClassField type identified by prefix
          k=GET_FIELDINDEX(CLSDEF(PCLS(IC))%C,PNAME(:i-1))
          if (k.gt.0) then
            ! get subclass definition for the ClassField type
            FD => CLSDEF(PCLS(IC))%C%FIELDS(k)%F
            call GetClassDef(FD%CID,COWN)
            ! get field index from the subclass definition
            if (associated(COWN%C)) then
              IX=GET_FIELDINDEX(COWN%C,trim(PNAME(i+1:)))
              if (IX.gt.0) OWNER=COWN%C%IDX
            endif
          endif
      ! no prefix => expect normal field
        ELSE
          IX=GET_FIELDINDEX(CLSDEF(PCLS(IC))%C,trim(PNAME))
          if (IX.gt.0) OWNER=PCLS(IC)
        endIF
        IC=IC-1
      endDO
      GetFieldIndex=IX
      end function GetFieldIndex

!---------------------------------------------------------
      integer function CLASSES_LENGTH(ICLS)
! return number of parameters (array elements!) for given class
! including parent(s) and suclass fields
!---------------------------------------------------------
      integer,intent(in) :: ICLS
      integer :: IRES,i,j,k,NCLS
      integer :: PCLS(1:10)
      TYPE(TFIELDDEF),pointer :: FD
      TYPE(PCLSDEF) :: CSUB
      CLASSES_LENGTH=0
      if ((ICLS.le.0).or.(ICLS.gt.CLS_NUM)) return
      IRES=0
      call GetClassParents(ICLS,PCLS,NCLS)
      do j=1,NCLS
        do i=1,CLSDEF(PCLS(j))%C%NF
          FD => CLSDEF(PCLS(j))%C%FIELDS(i)%F
          if (trim(FD%TID).eq.ID_FIELD_CLASS) then
            call getClassDef(trim(FD%CID),CSUB)
            if (associated(CSUB%C)) then
              do k=1,CSUB%C%NF
                FD => CSUB%C%FIELDS(k)%F
                IRES=IRES+FD%DIM
              enddo
            endif
          else
            IRES=IRES+FD%DIM
          endif
        enddo
      enddo
      CLASSES_LENGTH=IRES
      end function CLASSES_LENGTH

!////////////////////////////////////////////////////////////
!/////                                                  /////
!/////                 PRIVATE AREA                     /////
!/////                                                  /////
!////////////////////////////////////////////////////////////

!---------------------------------------------------------
      logical function isClassDefined(ID)
!---------------------------------------------------------
      character(*) :: ID
      logical :: b
      integer :: i
      b=.false.
      do i=1,CLS_NUM
        IF (ASSOCIATED(CLSDEF(i)%C)) THEN
        !  if (trim(ID).eq.'FOCARRAY') write(*,*) 'isClassDefined '//trim(CLSDEF(i)%C%ID)
          b=(trim(ID).eq.trim(CLSDEF(i)%C%ID))
          if (b) exit
        ENDIF
      enddo
      isClassDefined=b
      end function isClassDefined

!---------------------------------------------------------
      logical function isFieldDefined(CLS,ID)
!---------------------------------------------------------
      type(TCLSDEF),pointer :: CLS
      character(*) :: ID
      logical :: b
      integer :: i
      b=.false.
      do i=1,CLS%NF
        if (associated(CLS%FIELDS(i)%F)) then
          b=(trim(ID).eq.trim(CLS%FIELDS(i)%F%ID))
          if (b) exit
        endif
      enddo
      isFieldDefined=b
      end function isFieldDefined

!---------------------------------------------------------
      subroutine ALLOC_CLSDEF(CLS,ID,NAME)
!---------------------------------------------------------
        character(*) :: ID,NAME
        TYPE(PCLSDEF) :: CLS
        integer :: ierr
        allocate(CLS%C,STAT=ierr)
        if(ierr.eq.0) then
          call CLEAR_CLSDEF(CLS%C)
          CLS%C%ID=trim(ID)
          CLS%C%NAME=trim(NAME)
        endif
      end subroutine ALLOC_CLSDEF

!---------------------------------------------------------
      subroutine CLEAR_CLSDEF(CLS)
!---------------------------------------------------------
        TYPE(TCLSDEF),intent(out) :: CLS
          CLS%IDX=0
          CLS%ID=' '
          CLS%NAME=' '
          CLS%HINT=' '
          CLS%NF=0
          CLS%PARENT=0
          CLS%isCMD=.false.
          CLS%cType=0
      end subroutine CLEAR_CLSDEF


!---------------------------------------------------------
      subroutine DISPOSE_CLSDEF(C)
!---------------------------------------------------------
        TYPE(PCLSDEF) :: C
        integer :: i
        if (associated(C%C)) then
          do i=1,C%C%NF
            call DISPOSE_FIELDDEF(C%C%FIELDS(i))
          enddo
          C%C%PARENT=0
          DEallocate(C%C)
          nullify(C%C)
        endif
      end subroutine DISPOSE_CLSDEF

!---------------------------------------------------------
      integer function GET_CLSINDEX(ID)
! get index of class with given ID
!---------------------------------------------------------
      character(*) :: ID
      integer :: i,ix
      ix=0
      i=0
      do while ((i.lt.CLS_NUM).and.(ix.eq.0))
        i=i+1
        if (trim(ID).eq.trim(CLSDEF(i)%C%ID)) ix=i
      enddo
      GET_CLSINDEX = ix
      end function GET_CLSINDEX

!---------------------------------------------------------
      integer function GET_FIELDINDEX(CLS,ID)
! get index of field with given ID in given class
!---------------------------------------------------------
      character(*) :: ID
      type(TCLSDEF),pointer :: CLS
      integer :: i,ix
      ix=0
      i=0
      do while ((i.lt.CLS%NF).and.(ix.eq.0))
          i=i+1
          if (trim(ID).eq.trim(CLS%FIELDS(i)%F%ID)) ix=i
      enddo
      GET_FIELDINDEX = ix
      end function GET_FIELDINDEX

!---------------------------------------------------------
      subroutine ADD_CLSDEF(CLS,error)
! add new class definition
!---------------------------------------------------------
      type(TCLSDEF),pointer :: CLS
      logical,intent(out) :: error
      integer :: I,J
      error=.true.
      if ((CLS_NUM.lt.MAX_CLASSES).and.(.not.isClassDefined(CLS%ID))) then
        I=CLS_NUM+1
        call ALLOC_CLSDEF(CLSDEF(I),CLS%ID,CLS%NAME)
        if (associated(CLSDEF(I)%C)) then
          error=.false.
          CLS_NUM=I
          CLSDEF(I)%C=CLS
          CLS%NF=0
          do j=1,CLS%NF
            call ADD_FIELD(CLSDEF(I)%C,CLS%FIELDS(j)%F,error)
            if (error) exit
            CLSDEF(I)%C%FIELDS(j)%F%OWNER=I
          enddo
          CLSDEF(I)%C%IDX=I
        !  write(*,*) 'ADD_CLSDEF ',trim(CLSDEF(I)%C%ID)
        !  CLSDEF(I)%C=CLS
        !  CLSDEF(I)%C%HINT=
        endif
      endif
      end subroutine ADD_CLSDEF

!---------------------------------------------------------
      subroutine ADD_FIELD(CLS,FIELD,error)
! add new field definition to a class definition
!---------------------------------------------------------
      type(TCLSDEF),pointer :: CLS
      type(TFIELDDEF),pointer :: FIELD
      logical,intent(out) :: error
      integer :: I
      error=.true.
      if ((CLS%NF.lt.MAX_FIELDS).and.(.not.isFieldDefined(CLS,FIELD%ID))) then
        I=CLS%NF+1
        call ALLOC_FIELDDEF(CLS%FIELDS(I),FIELD%ID,FIELD%NAME)
        if (associated(CLS%FIELDS(I)%F)) then
          CLS%NF=I
          CLS%FIELDS(I)%F=FIELD
          CLS%FIELDS(I)%F%IDX=I
          CLS%FIELDS(I)%F%OWNER=CLS%IDX
          if (FIELD%TID.eq.ID_FIELD_CLASS) then
            error=(.not.isClassDefined(FIELD%CID))
          !  if (error) write(*,*) 'ADD_FIELD','Class field not defined: '//trim(CLS%ID)//'.'//trim(FIELD%CID)
            if (error) call MSG_ERROR('ADD_FIELD','Class field not defined: '//trim(CLS%ID)//'.'//trim(FIELD%CID),0,1)
          endif
          error=.false.
        endif
      endif
      end subroutine ADD_FIELD

!--------------------------------------------------------
      subroutine CLEAR_CLASSES
! clear definitions of all classes and types
!--------------------------------------------------------
        integer :: i
        call CLEAR_CLSDEF(ZERO_CLS)
        do i=1,CLS_NUM
          call DISPOSE_CLSDEF(CLSDEF(i))
        enddo
        CLS_NUM=0
        call DISPOSE_CLSDEF(myCLSDEF)
        call CLEAR_TYPES
        call CLEAR_ENUMS
      end subroutine CLEAR_CLASSES

!--------------------------------------------------------
      subroutine DISPOSE_TMP_CLASSES
! dispose temporary instaces of class definitions
!--------------------------------------------------------
        call DISPOSE_CLSDEF(myCLSDEF)
      end subroutine DISPOSE_TMP_CLASSES

 !--------------------------------------------------------
      integer function GET_CLASS_TYPE(IDSTR)
! dconvert type ID string to integer value
! one of TCLS_xxxx values defined in classes.f
!--------------------------------------------------------
      character(*) :: IDSTR
      integer :: RES
      RES=0
      select case  (trim(IDSTR))
      case('FRAME')
        RES=TCLS_COM
      case('IGROUP')
        RES=TCLS_IGP
      case('OPTION')
        RES=TCLS_OPT
      case('COMMAND')
        RES=TCLS_CMD
      case('SAMPLE')
        RES=TCLS_SAM
      case('SPECTROMETER')
        RES=TCLS_INS
      end select
      GET_CLASS_TYPE=RES
      end function GET_CLASS_TYPE

!--------------------------------------------------------
      subroutine startcls(tag,attribs,error)
!--------------------------------------------------------
        character(len=*)                  :: tag
        character(len=*), dimension(:,:)  :: attribs
        logical                           :: error
        character(LEN_ID)  :: parent,CLASS_ID
        character(LEN_NAME) :: ename,CLASS_NAME
        character(64) :: MSG
        logical :: hasParent
        MSG='invalid attributes in '//trim(tag)
        select case (trim(tag))
    ! new class
        case('CLASSES')
          call CLEAR_CLASSES
          isParsing=.true.
          return
    ! ignored tags
        case('GROUP','ENUMITEM','TRANSLATE')
          RETURN
        end select
        if (isParsing.and.(.not.error)) then
          if (associated(myCLSDEF%C)) then
        ! start of field definition
            call startFIELDDEF(myCLSDEF%C%ID,tag,attribs,error)
          else
            select case(trim(tag))
        ! new class
            case('CLASS')
              if (CLS_NUM.GE.MAX_CLASSES) then
                call MSG_ERROR('CLASSDEF','too many classes defined',0,1)
                error=.true.
                return
              endif
              call readAttrib(attribs,'id',CLASS_ID,error)
              if (error) goto 99
              call readAttrib(attribs,'name',CLASS_NAME,error)
              if (error) goto 99
              call readAttrib(attribs,'inherits',parent,error)
              hasParent=(.not.error)
              call MKUPCASE(parent)
              call MKUPCASE(CLASS_ID)
              call ALLOC_CLSDEF(myCLSDEF,CLASS_ID,CLASS_NAME)
              error=(.not.associated(myCLSDEF%C))
              if (error) then
                MSG='class not allocated: '//trim(tag)
                goto 99
              endif
            ! determine class type
              if (hasParent) then
                if (isClassDefined(parent)) then
                  myCLSDEF%C%PARENT=GET_CLSINDEX(parent)
                  myCLSDEF%C%cType=GET_CLASS_TYPE(parent)
                  myCLSDEF%C%isCMD=(myCLSDEF%C%cType.eq.TCLS_CMD)
                ! exception: SAMPLE inherits from FRAME
                  if (trim(CLASS_ID).eq.'SAMPLE') myCLSDEF%C%cType=TCLS_SAM
                ELSE
                  hasParent=.false.
                  MSG='undefined parent class: '//trim(parent)
                  goto 99
                endif
              else
                ! get root class type
                myCLSDEF%C%cType=GET_CLASS_TYPE(CLASS_ID)
              endif
        ! new enumerator, can be defined outside class
            case(ID_FIELD_ENUM)
              call readAttrib(attribs,'id',ename,error)
              if (.not.error) then
                call MKUPCASE(ename)
                call startENUMDEF(ename,error)
              endif
            case default
              MSG='undefined tag '//trim(tag)
              error=.true.
            end select
          endif
        endif
        return
99      error=.true.
        call MSG_ERROR('CLASSDEF',trim(MSG),0,1)
      end subroutine startcls


!--------------------------------------------------------
      subroutine datacls( tag, xmldata, error )
!--------------------------------------------------------
        character(len=*)                  :: tag
        character(len=*)                  :: xmldata
        logical                           :: error
      ! process only non-empty data inside class entries
        if (.not.isParsing) return
        if (len_trim(xmldata).le.0) return
        select case(trim(tag))
    ! ignored tags
          case('GROUP','TRANSLATE')
            return
    ! enumerator item
          case('ENUMITEM')
            call dataENUMDEF(xmldata)
    ! hint text for class definition
          case('CLASS')
            if (associated(myCLSDEF%C)) then
              myCLSDEF%C%HINT=trim(adjustl(xmldata))
            endif
        end select
      end subroutine datacls

!--------------------------------------------------------
      subroutine endcls( tag, error )
!--------------------------------------------------------
      character(len=*)                  :: tag
      logical                           :: error
      if (.not.isParsing) return
      select case (trim(tag))
    ! ignored tags
      case ('GROUP','ENUMITEM','TRANSLATE')
        return
    ! end of classes definition, stop pasing
      case('CLASSES')
        isParsing=.false.
        call DISPOSE_TMP_CLASSES
        call DISPOSE_TMP_TYPES
        call DISPOSE_TMP_ENUMS
        return
    ! enum definition
      case(ID_FIELD_ENUM)
        call endENUMDEF(error)
        if (associated(myFIELDDEF%F)) myFIELDDEF%F%ENUM=ENUM_NUM
      end select
      if (.not.associated(myCLSDEF%C)) return
    ! do this only if class is associated
      select case(trim(tag))
    ! class definition finished
      case('CLASS')
        call ADD_CLSDEF(myCLSDEF%C,error)
        call DISPOSE_CLSDEF(myCLSDEF)
      case default
        if (associated(myFIELDDEF%F)) then
          call ADD_FIELD(myCLSDEF%C,myFIELDDEF%F,error)
          if (error) call MSG_ERROR('CLASSDEF','error when adding field: '//trim(myFIELDDEF%F%ID),0,1)
          call DISPOSE_FIELDDEF(myFIELDDEF)
        endif
      end select
      end subroutine endcls

!--------------------------------------------------------
      subroutine report_class(CLS,IU)
! report class definition
!--------------------------------------------------------
      type(PCLSDEF) :: CLS
      integer :: IU
      character(16) :: CUNIT
      integer :: NP,j
      integer, parameter :: nt=3
      type(TFIELDDEF),POINTER :: FD
      character*64  :: TABS
      integer :: ierr
      TABS=' '
1     format(a,I5,a,' ',a16,a16,a16,' size=',I5)
2     format(a,I5,a,' ',a16,a16,a16,' size=',I5,' ',a)
3     format(I5,a,a,' fields: ',I5,' parent: ',I5)
11    format(' [',a,'] ')
      if (associated(CLS%C)) then
        NP=CLS%C%NF
        write(IU,3) CLS%C%IDX,TABS(1:nt),trim(CLS%C%ID),NP,CLS%C%PARENT
        do j=1,NP
          FD=>CLS%C%FIELDS(j)%F
          WRITE(CUNIT,11,IOSTAT=ierr) trim(FD%UNITS)
          if (FD%ENUM>0) then
            write(IU,2) TABS(1:nt),FD%IDX,TABS(1:nt),trim(FD%ID),trim(CUNIT),trim(FD%TID),FD%DIM,trim(GetEnumID(FD%ENUM))
          else
            write(IU,1) TABS(1:nt),FD%IDX,TABS(1:nt),trim(FD%ID),trim(CUNIT),trim(FD%TID),FD%DIM
          endif
        enddo
      endif
      end subroutine report_class

!--------------------------------------------------------
      subroutine new_report_classes
! report classes definition as read from classes.xml
!--------------------------------------------------------
      integer :: IU,i

      IU=OPENFILEUNIT('new_report_classes.txt',.false.)
  !   IU=6
      if (IU.le.0) return
      write(IU,*) 'CLASSES:'
      do i=1,CLS_NUM
        if (.not.CLSDEF(i)%C%isCMD) then
          call report_class(CLSDEF(i),IU)
        endif
      enddo
      write(IU,*) 'COMMANDS:'
      do i=1,CLS_NUM
        if (CLSDEF(i)%C%isCMD) then
          call report_class(CLSDEF(i),IU)
        endif
      enddo
      write(IU,*) 'ENUM:'
      do i=1,ENUM_NUM
        call report_enum(i,IU)
      enddo
      close(IU)
      end subroutine new_report_classes


      end module CLASSDEF


