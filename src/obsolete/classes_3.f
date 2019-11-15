!//////////////////////////////////////////////////////////////////////
!////  $Id: classes_3.f,v 1.1 2009/02/22 00:45:16 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.1 $
!////     $Date: 2009/02/22 00:45:16 $
!//////////////////////////////////////////////////////////////////////
!////
!////   Definitions of all classes and parameters parsed from classes.xml
!////  - class IDs, indexes, parameter IDs, names, units, hints, enum. definitions etc.
!////
!//////////////////////////////////////////////////////////////////////
      MODULE CLASSES
      use FILETOOLS,ONLY: OPENFILEUNIT
      use XMLINFO
      USE COMMANDS
      USE XMLPARSE
      implicit none

      private

! ELEMENTS DEFINED BY SIMRES:
! only those with named indices

! class IDs used internally in dispatch procedures
      integer,parameter :: TCLS_COM=1
      integer,parameter :: TCLS_SAM=2
      integer,parameter :: TCLS_INS=3
      integer,parameter :: TCLS_OPT=4
      integer,parameter :: TCLS_CMD=5
      character(45),parameter :: TCLS_NAMES='COMPONENT:SPECIMEN:INTERFACE:OPTIONS:COMMANDS'

! option classes
      integer,parameter :: OCLS_TRACING=1
      integer,parameter :: OCLS_REPORTS=2
      character(15),parameter :: OCLS_NAMES='TRACING:REPORTS'

! device interfaces
      integer, PARAMETER :: DCLS_SPEC = 1
      integer, PARAMETER :: DCLS_PWDIF = 2
      integer, PARAMETER :: DCLS_TAS = 3
      character(20),parameter :: DCLS_NAMES='SPECTROMETER:PWD:TAS'

! component classes IDs
      integer, PARAMETER :: CCLS_FRAME = 1
      integer, PARAMETER :: CCLS_SOURCE = 2
      integer, PARAMETER :: CCLS_DETECTOR = 3
      integer, PARAMETER :: CCLS_GUIDE = 4
      integer, PARAMETER :: CCLS_CRYSTAL = 5
      integer, parameter :: CCLS_NUMBER = CCLS_CRYSTAL
  ! class identificator (keep the same order as class ID !)
      character(35),parameter :: CCLS_NAMES='FRAME:SOURCE:DETECTOR:GUIDE:CRYSTAL'

! sample classes
      integer, PARAMETER :: SCLS_SAMPLE = 1
      integer, PARAMETER :: SCLS_POLY = 2
      integer, PARAMETER :: SCLS_SINGLE = 3
      character(20),parameter :: SCLS_NAMES='SAMPLE:PCRYST:SCRYST'

! ENUMERATED TYPES:
! FRAME_SHAPE
      integer, PARAMETER :: FRAME_SHAPE_ELLIPSOID = 0
      integer, PARAMETER :: FRAME_SHAPE_CYLINDER = 1  ! cylinder, axis // y
      integer, PARAMETER :: FRAME_SHAPE_DISC = 2      ! cylinder, axis // z
      integer, PARAMETER :: FRAME_SHAPE_BOX = 3

! SAMPLE_TYPE
      integer, PARAMETER :: SAMPLE_TYPE_INELASTIC = 0
      integer, PARAMETER :: SAMPLE_TYPE_ELASTIC= 1
      integer, PARAMETER :: SAMPLE_TYPE_POWDER= 2
      integer, PARAMETER :: SAMPLE_TYPE_VANAD= 3
      integer, PARAMETER :: SAMPLE_TYPE_PHONON= 4

! GUIDE_TYPE
      integer, PARAMETER :: GUIDE_TYPE_SOLLER = 0
      integer, PARAMETER :: GUIDE_TYPE_GUIDE= 1
      integer, PARAMETER :: GUIDE_TYPE_PARA= 2
      integer, PARAMETER :: GUIDE_TYPE_PARA2= 3
      integer, PARAMETER :: GUIDE_TYPE_ELL= 4

! CRYSTAL_MODEL
      integer, PARAMETER :: CRYSTAL_MODEL_SIMPLE = 0
      integer, PARAMETER :: CRYSTAL_MODEL_RW = 1

! DETECTOR_TYPE
      integer, PARAMETER :: DETECTOR_TYPE_AREA = 0
      integer, PARAMETER :: DETECTOR_TYPE_ARRAY = 1
      integer, PARAMETER :: DETECTOR_TYPE_PSD = 2

! BOOLEAN
      integer, PARAMETER :: BOOL_FALSE = 0
      integer, PARAMETER :: BOOL_TRUE = 1

! NOYES
      integer, PARAMETER :: NOYES_NO = 0
      integer, PARAMETER :: NOYES_YES= 1

! KFMODE (Kf axis orientation)
      integer,parameter :: KFMODE_FLAT=0
      integer,parameter :: KFMODE_TOF=1


! ELEMENTS DEFINED IN classes.xml:
! DIMENSIONS - should be large enough to fit all parameters in
      integer, parameter :: LEN_ID=16             ! length of ID strings
      integer, parameter :: LEN_NAME=32           ! length of NAME strings
      integer, parameter :: LEN_LINE=255          ! default length of input line
      integer, parameter :: DIM_CLS=64            ! max. number of classes
      integer, parameter :: DIM_FLDS=32           ! max. number of fields (any type) in one class
      integer, parameter :: DIM_FLDSS=8           ! max. number of string fields in one class
      integer, parameter :: DIM_ID=DIM_FLDS*(LEN_ID+1)      ! length of the string with field IDs for one class
      integer, parameter :: DIM_NAMES=DIM_FLDS*(LEN_NAME+1)  ! length of the string with field NAMEs for one class
! BEAMLINE (keep it here, it is neded to define dimensions by many other modules)
      integer, parameter :: BEAMLINE_DIM=100  ! max. number of components

! CLASSES
  ! string with class IDs
      character(DIM_CLS*(LEN_ID+1)) :: CLS_IDS
  ! string with class names
      character(DIM_CLS*(LEN_NAME+1)) :: CLS_NAMES
  ! number of items in CLS_IDS
      integer :: CLS_NUM
  ! parents of each class
      character(DIM_ID) :: CLS_PARENT(DIM_CLS)
      integer :: CLS_PARENTID(DIM_CLS)
  ! flag for command classes
      logical :: CLS_COMMAND(DIM_CLS)
! FIELDS
  ! parameter IDs (including string type)
      character(DIM_ID) :: CPAR_IDS(DIM_CLS)
  ! parameter IDs (ONLY string type)
      character(DIM_FLDSS*(LEN_ID+1)) :: SPAR_IDS(DIM_CLS)
  ! types of parameters (F,I,S or an enum type)
      character(DIM_ID) :: CPAR_TYPES(DIM_CLS)
  ! strings with parameter descriptive names
      character(DIM_NAMES) :: CPAR_NAMES(DIM_CLS)
  ! strings with physical units for each parameter
      character(DIM_NAMES) :: CPAR_UNITS(DIM_CLS)
  ! partition table for parameter arrays
      integer :: CPAR_NPAR(DIM_CLS)            ! number of fields (any type) for each class
      integer :: SPAR_NPAR(DIM_CLS)            ! number of string fieldsfor each class
      integer :: CPAR_INDX(0:DIM_FLDS,DIM_CLS) ! addresses of each field in each class
! ENUMERATED TYPES
  ! definition of enumerated types read from classes.xml
  ! enum names are constructed as (classname)_(variable name), e.g. SAMPLE_TYPE
      integer, parameter :: ENUM_MAX_TYPES=32               ! max. number of enumerated types
      integer, parameter :: ENUM_BUFF_LEN=1024              ! buffer length for storage of enum items IDs
      character(ENUM_MAX_TYPES*(2*LEN_ID+1)) ::  ENUM_TYPES ! list of type names
      character*(ENUM_BUFF_LEN) :: ENUM_TYPE_LIST(ENUM_MAX_TYPES)    ! list of items for each type

! local variables useful during parsing
      integer,private :: CLS_INDEX
      integer,private :: ENUM_INDEX
      CHARACTER(LEN_ID),private  :: CLASS_ID
      CHARACTER(LEN_ID),private  :: PARAM_ID
      CHARACTER(LEN_NAME),private  :: CLASS_NAME
      logical,private :: isCOMMAND
! instance of TCOMMAND used at loaing the commands definitions
      type(TCOMMAND), private, SAVE :: myCMD


  !    interface getEnumTypeOrd
  !      module procedure getEnumTypeOrd_1
  !      module procedure getEnumTypeOrd_2
  !    end interface


      public LEN_ID,LEN_NAME,DIM_ID
      public DCLS_SPEC
      public OCLS_TRACING,OCLS_REPORTS,OCLS_NAMES
      public SCLS_SAMPLE,SCLS_SINGLE,SCLS_POLY,SCLS_NAMES
      public CCLS_FRAME,CCLS_GUIDE,CCLS_DETECTOR,CCLS_SOURCE,CCLS_CRYSTAL,CCLS_NAMES
      public DCLS_TAS,DCLS_PWDIF,DCLS_NAMES
      public TCLS_COM,TCLS_SAM,TCLS_OPT,TCLS_INS,TCLS_CMD,TCLS_NAMES
      public KFMODE_FLAT,KFMODE_TOF
      public FRAME_SHAPE_BOX,FRAME_SHAPE_CYLINDER,FRAME_SHAPE_DISC,FRAME_SHAPE_ELLIPSOID

      public getICLS
      public DIM_FLDS,BEAMLINE_DIM

      contains


!---------------------------------------------------------
      integer function getEnumTypeOrd_1(ICLS,PNAME)
! return type index for a parameter given by object ICLS and ID string
! return 0 if not defined
!---------------------------------------------------------
      integer :: ICLS
      character(*) :: PNAME
      integer :: CLS,IDX,NARG,IT
      character(LEN_NAME) :: ETYPE
      call CLASSES_PARINFO(ICLS,PNAME,CLS,IDX,ETYPE,NARG)
      call GETPARINDX(ENUM_TYPES,':',trim(ETYPE),IT)
      getEnumTypeOrd_1=IT
      end function getEnumTypeOrd_1

C-------------------------------------------------------------
      integer function getEnumTypeOrd_2(ETYPE)
! convert enumerator type to integer
C-------------------------------------------------------------
      character(*) ETYPE
      integer :: IT
      IT=0
      call GETPARINDX(ENUM_TYPES,':',trim(ETYPE),IT)
      getEnumTypeOrd_2=IT
      end function getEnumTypeOrd_2

!---------------------------------------------------------
      SUBROUTINE CLASSES_NAMESPACE(ICLS,NAMESPACE)
! return parameter namespace as : delimitted string
! including parent(s)
!---------------------------------------------------------
      integer,intent(in) :: ICLS
      character(*),intent(out) :: NAMESPACE
      integer :: IPAR,NC
      NAMESPACE=' '
      call COUNTPAR(trim(CLS_IDS),':',NC)
      if ((ICLS.le.0).or.(ICLS.gt.NC)) return
      NAMESPACE=trim(CPAR_IDS(ICLS))
      IPAR=CLS_PARENTID(ICLS)
      do while (IPAR.gt.0)
        call PREPENDPAR(NAMESPACE,' ',CPAR_IDS(IPAR))
        IPAR=CLS_PARENTID(IPAR)
      enddo
      end SUBROUTINE CLASSES_NAMESPACE

!---------------------------------------------------------
      integer function CLASSES_PARINDX(ICLS,PNAME)
! return index of a parameter of given class
! including parent(s), in order as returned by CLASSES_NAMESPACE
! i.e. PARENT1/PARENT2/.../CLS
!---------------------------------------------------------
      integer,intent(in) :: ICLS
      character(*),intent(in) :: PNAME
      integer :: IC,IP,IX
      integer :: PCLS(0:10),NCLS
! collect all ancestors IDs in PCLS
      NCLS=0
      IC=ICLS
      PCLS(0)=IC
      do while ((IC.gt.0).and.(NCLS.lt.10).and.(CLS_PARENTID(IC).gt.0))
        NCLS=NCLS+1
        PCLS(NCLS)=CLS_PARENTID(IC)
        IC=CLS_PARENTID(IC)
      enddo
      IC=NCLS
      IP=0
      IX=0
      DO WHILE ((IC.GE.0).and.(IP.le.0))
        CALL GETPARINDX(CPAR_IDS(PCLS(IC)),':',trim(PNAME),IP)
        if (IP.gt.0) then
          IX=IX+CPAR_INDX(IP-1,PCLS(IC))+1
        else
          IX=IX+CPAR_INDX(CPAR_NPAR(PCLS(IC)),PCLS(IC))
        endif
        IC=IC-1
      enddo
      if (IP.le.0) IX=0
      CLASSES_PARINDX=IX
      end function CLASSES_PARINDX

!---------------------------------------------------------
      SUBROUTINE CLASSES_PARINFO(ICLS,PNAME,CLS,IDX,TID,NARG)
! return actual class, of which PNAME is a member.
! ICLS is the class ID, but it may include parents.
! This procedure searches in the parents namespaces and returns ID
! of the parent or ICLS
! RETURNS
!   CLS   ... class ID of the base parent or ICLS
!   IDX   ... index of the parameter in CPAR_IDS or SPAR_IDS list
!   TID   ... ID string of variable type
!   NARG  ... size of the parameter (>1 for array types)
!---------------------------------------------------------
      integer,intent(in) :: ICLS
      character(*),intent(in) :: PNAME
      integer,intent(out) :: CLS,IDX,NARG
      character(*),intent(out) :: TID
      integer :: IC,IP,IL,IS
      IC=ICLS
      IP=0
      NARG=0
      TID=' '
      do while ((IC.gt.0).and.(IP.le.0))
        CALL GETPARINDX(CPAR_IDS(IC),':',trim(PNAME),IP)
        if (IP.le.0) then
          IC=CLS_PARENTID(IC)
        else
          CALL FINDSTRPAR(CPAR_TYPES(IC),':',IP,IS,IL)
          TID=CPAR_TYPES(IC)(IS:IS+IL-1)
          NARG=CPAR_INDX(IP,IC)-CPAR_INDX(IP-1,IC)
    ! for strings, IDX=index in the SPAR_IDS  list
          if (trim(TID).eq.'S') CALL GETPARINDX(SPAR_IDS(IC),':',trim(PNAME),IP)
        endif
      enddo
      CLS=IC
      IDX=IP
      end SUBROUTINE CLASSES_PARINFO


!---------------------------------------------------------
      integer function CLASSES_NPAR(ICLS)
! return number of parameters (array elements) for given class
! including parent(s)
!---------------------------------------------------------
      integer,intent(in) :: ICLS
      integer :: IPAR,NC,IRES
      CLASSES_NPAR=0
      call COUNTPAR(trim(CLS_IDS),':',NC)
      if ((ICLS.le.0).or.(ICLS.gt.NC)) return
      IRES=CPAR_INDX(CPAR_NPAR(ICLS),ICLS)-CPAR_INDX(0,ICLS)
      IPAR=CLS_PARENTID(ICLS)
      do while (IPAR.gt.0)
        IRES=IRES+CPAR_INDX(CPAR_NPAR(IPAR),IPAR)-CPAR_INDX(0,IPAR)
        IPAR=CLS_PARENTID(IPAR)
      enddo
      CLASSES_NPAR=IRES
      end function CLASSES_NPAR

C-------------------------------------------------------------
      SUBROUTINE GetEnumOrd(ETYPE,STRVAL,ITYPE,IORD)
! convert enumerator type and value given by strings to integers
! note: values are indexed from 0 !
! IORD=-1 for valid enum tyoe, but undefined value
! IORD=-2 for invalid enum type
C-------------------------------------------------------------
      character(*) ETYPE,STRVAL
      integer, intent(out) :: ITYPE,IORD
      integer :: IT,IE
      ITYPE=-1
      IORD=-1
      IT=getEnumTypeOrd(ETYPE)
      if (IT.gt.0) then
        ITYPE=IT
        call GETPARINDX(ENUM_TYPE_LIST(IT),':',trim(adjustl(STRVAL)),IE)
        if (IE.gt.0) IORD=IE-1
      else
        IORD=-2
      endif
      end SUBROUTINE GetEnumOrd


!--------------------------------------------------------
      subroutine GetEnumString(ETYPE,IORD,ESTR)
! return string value for given enumerated type and ord number
!--------------------------------------------------------
      character(*),intent(in) :: ETYPE
      integer,intent(in) :: IORD
      character(*),intent(out) :: ESTR
      integer :: IE,IS,IL
        ESTR='undefined'
        CALL GETPARINDX(ENUM_TYPES,':',trim(ETYPE),IE)
        if ((IE.gt.0).and.(IORD.ge.0)) then
          call FINDSTRPAR(ENUM_TYPE_LIST(IE),':',IORD+1,IS,IL)
          if (IL.GT.0) ESTR=ENUM_TYPE_LIST(IE)(IS:IS+IL-1)
        endif
      END subroutine GetEnumString


C-------------------------------------------------------------
      SUBROUTINE GetParAttrib(ICLS,PID,PTYPE,PNAME,PUNIT)
! Get attributes (name, units, size) of a named parameter
! INPUT
!   ICLS   .. class index in CLS_IDS
!   PID    .. parameter ID string (may belong to parent)
! OUTPUT
!   PTYPE  .. type ID string
!   PNAME  .. descriptive name
!   PUNIT  .. phys. units
C-------------------------------------------------------------
      integer, intent(in) :: ICLS
      character(*),intent(in) :: PID
      character(*),intent(out) :: PTYPE,PNAME,PUNIT
      integer :: IC,IVAR,IS,IL
      IC=0
      IVAR=0
      PTYPE=' '
      PNAME=' '
      PUNIT=' '
      if (ICLS.le.0) return
! try to find the variable first in the specific set
      CALL GETPARINDX(CPAR_IDS(ICLS),':',PID,IVAR)
!      write(*,*) 'GetParAttrib PID=',trim(PID),' ICLS=',ICLS,' IVAR=',IVAR
      IC=ICLS
      if (IVAR.le.0) then
! try parents
        do while ((CLS_PARENTID(IC).gt.0).and.(IVAR.le.0))
          IC=CLS_PARENTID(IC)
          CALL GETPARINDX(CPAR_IDS(IC),':',PID,IVAR)
!          write(*,*) 'GetParAttrib PID=',trim(PID),' IC=',IC,' IVAR=',IVAR
!          if (ICLS.eq.6) read(*,*)
        enddo
      endif
! if IVAR>0, then also IC>0
      if (IVAR.gt.0) then
!        write(*,*) 'GetParAttrib CPAR_TYPES(IC)=',trim(CPAR_TYPES(IC))
        CALL FINDSTRPAR(CPAR_TYPES(IC),':',IVAR,IS,IL)
        PTYPE=CPAR_TYPES(IC)(IS:IS+IL-1)
        CALL FINDSTRPAR(CPAR_NAMES(IC),':',IVAR,IS,IL)
        PNAME=CPAR_NAMES(IC)(IS:IS+IL-1)
        CALL FINDSTRPAR(CPAR_UNITS(IC),':',IVAR,IS,IL)
        PUNIT=CPAR_UNITS(IC)(IS:IS+IL-1)
!        write(*,*) 'GetParAttrib PNAME=',trim(PNAME),' PTYPE=',trim(PTYPE)
      endif
      end SUBROUTINE GetParAttrib


!--------------------------------------------------------
! get index in CLS_IDS from CLS index in given namespace (CLSNAMES)
!--------------------------------------------------------
      integer function getICLS(CCLS,CLSNAMES)
      integer,intent(in) :: CCLS
      character(*),intent(in) :: CLSNAMES
      integer :: IRES,IS,IL
      IRES=-1
      call FINDSTRPAR(CLSNAMES,':',CCLS,IS,IL)
      if (IL.gt.0) then
        call GETPARINDX(CLS_IDS,':',CLSNAMES(IS:IS+IL-1),IRES)
      endif
      getICLS=IRES
      end function getICLS

!--------------------------------------------------------
! get index in CLSNAMES from ICLS (index from CLS_IDS)
!--------------------------------------------------------
      integer function getCCLS(ICLS,CLSNAMES)
      integer,intent(in) :: ICLS
      character(*),intent(in) :: CLSNAMES
      integer :: IRES,IS,IL
      IRES=-1
      call FINDSTRPAR(CLS_IDS,':',ICLS,IS,IL)
      if (IL.gt.0) then
        call GETPARINDX(CLSNAMES,':',CLS_IDS(IS:IS+IL-1),IRES)
      endif
      getCCLS=IRES
      end function getCCLS


!************************************************************************
! Parser handlers - CLASSES
!************************************************************************

!--------------------------------------------------------
      subroutine CLEAR_CLASSES
! clear definitions of all classes and types
!--------------------------------------------------------
        character(LEN_ID) :: ename
        CLS_IDS=' '
        CLS_NAMES=' '
        CLS_NUM=0
        CLS_COMMAND(1:DIM_CLS)=.false.
        CLS_PARENT(1:DIM_CLS)=' '
        CPAR_IDS(1:DIM_CLS)=' '
        SPAR_IDS(1:DIM_CLS)=' '
        CPAR_NAMES(1:DIM_CLS)=' '
        CPAR_UNITS(1:DIM_CLS)=' '
        CPAR_TYPES(1:DIM_CLS)=' '
        CPAR_NPAR(1:DIM_CLS)=0
        SPAR_NPAR(1:DIM_CLS)=0
        CLS_PARENTID(1:DIM_CLS)=0
        CPAR_INDX(:,:)=0
        ENUM_TYPES=' '
        ENUM_TYPE_LIST(1:ENUM_MAX_TYPES)=' '
        ENUM_INDEX=0
        CALL CLEAR_COMMANDS
    ! add default enumerators

        ename='NOYES'
        call APPENDPAR(ENUM_TYPES,':',ename)
        ENUM_INDEX=ENUM_INDEX+1
        ENUM_TYPE_LIST(ENUM_INDEX)=' '
        call APPENDPAR(ENUM_TYPE_LIST(ENUM_INDEX),':','no')
        call APPENDPAR(ENUM_TYPE_LIST(ENUM_INDEX),':','yes')

        ename='SIGN'
        call APPENDPAR(ENUM_TYPES,':',ename)
        ENUM_INDEX=ENUM_INDEX+1
        ENUM_TYPE_LIST(ENUM_INDEX)=' '
        call APPENDPAR(ENUM_TYPE_LIST(ENUM_INDEX),':','-1')
        call APPENDPAR(ENUM_TYPE_LIST(ENUM_INDEX),':','+1')


      end subroutine CLEAR_CLASSES


!--------------------------------------------------------
      subroutine startcls(tag,attribs,error)
!--------------------------------------------------------
        character(len=*)                  :: tag
        character(len=*), dimension(:,:)  :: attribs
        logical                           :: error
        integer :: i,isize,iat,NP
        character(LEN_ID)  :: parent
        character(LEN_NAME) :: uname,ename,pname
    ! new class
        if (CLS_INDEX.le.0) then
          select case(trim(tag))
    ! START: clear all fields
          case('CLASSES')
            call CLEAR_CLASSES
    ! new class or a reference to a class from another one
          case('CLASS')
            CLASS_ID=' '
            CLASS_NAME=' '
            parent=' '
            iat=0
            if (CLS_NUM.GE.DIM_CLS) then
              call MSG_ERROR('CLASSES','too many classes defined',0,1)
              error=.true.
              return
            endif
            do i=1,size(attribs,2)
              select case(trim(attribs(1,i)))
              case('id')
                CLASS_ID=trim(attribs(2,i))
                iat=iat+1
              case('name')
                CLASS_NAME=trim(attribs(2,i))
                iat=iat+1
              case('inherits')
                parent=trim(attribs(2,i))
              end select
            enddo
            if (iat.lt.2) goto 99
            isCOMMAND=(trim(parent).eq.'COMMAND')
        ! we have a new valid class definition:
            call APPENDPAR(CLS_IDS,':',CLASS_ID)
            call APPENDPAR(CLS_NAMES,':',CLASS_NAME)
            CLS_NUM=CLS_NUM+1
            CLS_INDEX=CLS_NUM
            if (len_trim(parent).gt.0) then
              call GETPARINDX(CLS_IDS,':',trim(adjustl(parent)),i)
              if (i.gt.0) then
                CLS_PARENTID(CLS_INDEX)=i
              else
                CLS_PARENTID(CLS_INDEX)=0
                call MSG_ERROR('CLASSES','undefined parent class: '//trim(parent),0,1)
              endif
            else
              CLS_PARENTID(CLS_INDEX)=0
            endif
            CLS_PARENT(CLS_INDEX)=trim(parent)
            CLS_COMMAND(CLS_INDEX)=isCOMMAND
          end select
    ! inside class entry
        else
          iat=0
          isize=1
          uname=" "
        !  if (isCOMMAND) then
        !    write(*,*) 'startcls ',myCMD%NARG,' ',trim(tag)
        !    read(*,*)
        !  endif
          select case(trim(tag))
          case('FLOAT')
            do i=1,size(attribs,2)
              select case(trim(attribs(1,i)))
              case('id')
                PARAM_ID=trim(attribs(2,i))
                iat=iat+1
              case('name')
                pname=trim(attribs(2,i))
                iat=iat+1
              case('units')
                uname=trim(attribs(2,i))
                iat=iat+1
              case('size')
                READ(trim(attribs(2,i)),*,err=99) isize
              end select
            enddo
            if (iat.lt.3) goto 99
            call APPENDPAR(CPAR_IDS(CLS_INDEX),':',PARAM_ID)
            call APPENDPAR(CPAR_NAMES(CLS_INDEX),':',pname)
            call APPENDPAR(CPAR_UNITS(CLS_INDEX),':',uname)
            call APPENDPAR(CPAR_TYPES(CLS_INDEX),':','F')
            NP=CPAR_NPAR(CLS_INDEX)+1
            CPAR_NPAR(CLS_INDEX)=NP
            CPAR_INDX(NP,CLS_INDEX)=CPAR_INDX(NP-1,CLS_INDEX)+isize
          case('INT','NOYES','SIGN','VRANGE')
            do i=1,size(attribs,2)
              select case(trim(attribs(1,i)))
              case('id')
                PARAM_ID=trim(attribs(2,i))
                iat=iat+1
              case('name')
                pname=trim(attribs(2,i))
                iat=iat+1
              case('size')
                READ(trim(attribs(2,i)),*,err=99) isize
              end select
            enddo
            if (iat.lt.2) goto 99
            if (trim(tag).eq.'INT') then
              ename='I'
            else
              ename=trim(tag)
            endif
            call APPENDPAR(CPAR_IDS(CLS_INDEX),':',PARAM_ID)
            call APPENDPAR(CPAR_NAMES(CLS_INDEX),':',pname)
            call APPENDPAR(CPAR_UNITS(CLS_INDEX),':'," ")
            call APPENDPAR(CPAR_TYPES(CLS_INDEX),':',trim(ename))
            NP=CPAR_NPAR(CLS_INDEX)+1
            CPAR_NPAR(CLS_INDEX)=NP
            CPAR_INDX(NP,CLS_INDEX)=CPAR_INDX(NP-1,CLS_INDEX)+isize
          case('STRING')
            do i=1,size(attribs,2)
              select case(trim(attribs(1,i)))
              case('id')
                PARAM_ID=trim(attribs(2,i))
                iat=iat+1
              case('name')
                pname=trim(attribs(2,i))
                iat=iat+1
              end select
            enddo
            if (iat.lt.2) goto 99
            call APPENDPAR(CPAR_IDS(CLS_INDEX),':',PARAM_ID)
            call APPENDPAR(CPAR_NAMES(CLS_INDEX),':',pname)
            call APPENDPAR(CPAR_UNITS(CLS_INDEX),':'," ")
            call APPENDPAR(CPAR_TYPES(CLS_INDEX),':','S')
            call APPENDPAR(SPAR_IDS(CLS_INDEX),':',PARAM_ID)
            SPAR_NPAR(CLS_INDEX)=SPAR_NPAR(CLS_INDEX)+1
            NP=CPAR_NPAR(CLS_INDEX)+1
            CPAR_NPAR(CLS_INDEX)=NP
            CPAR_INDX(NP,CLS_INDEX)=CPAR_INDX(NP-1,CLS_INDEX)+1
          case('ENUM')
            do i=1,size(attribs,2)
              select case(trim(attribs(1,i)))
              case('id')
                PARAM_ID=trim(attribs(2,i))
                iat=iat+1
              case('name')
                pname=trim(attribs(2,i))
                iat=iat+1
              end select
            enddo
            if (iat.lt.2) goto 99
            ename=trim(adjustl(CLASS_ID))//'_'//trim(adjustl(PARAM_ID))
            call APPENDPAR(CPAR_IDS(CLS_INDEX),':',PARAM_ID)
            call APPENDPAR(CPAR_NAMES(CLS_INDEX),':',pname)
            call APPENDPAR(CPAR_UNITS(CLS_INDEX),':'," ")
            call APPENDPAR(CPAR_TYPES(CLS_INDEX),':',ename)
            NP=CPAR_NPAR(CLS_INDEX)+1
            CPAR_NPAR(CLS_INDEX)=NP
            CPAR_INDX(NP,CLS_INDEX)=CPAR_INDX(NP-1,CLS_INDEX)+1
        ! update list of enumerators only if not yet defined
            if (getEnumTypeOrd(ename).le.0) then
              call APPENDPAR(ENUM_TYPES,':',ename)
              ENUM_INDEX=ENUM_INDEX+1
              ENUM_TYPE_LIST(ENUM_INDEX)=' '
            endif
          end select
c    1     format(a16,a16,a16)
c          write(*,1) trim(PARAM_ID),trim(pname),trim(uname)
c          read(*,*)
        endif
        if (isCOMMAND) call startCOMMAND(tag,isize,error)
        return
99      error=.true.
        call MSG_ERROR('CLASSES','invalid attributes in '//trim(tag),0,1)
      end subroutine startcls


!--------------------------------------------------------
      subroutine startCOMMAND(tag,isize,error)
! special handling of start elements for COMMAND
!--------------------------------------------------------
      character(*),intent(in) :: tag
      integer,intent(in) :: isize
      logical,intent(out) :: error
    !  integer :: new_index,old_index
      !  write(*,*) 'startCOMMAND ',trim(tag),' ',trim(PARAM_ID)
      !  read(*,*)
  !      if (trim(tag).ne.'CLASS') then
  !        write(*,*) 'startCOMMAND ',trim(tag),'.',trim(PARAM_ID)
  !        call ListCommandPar(myCMD)
  !      endif

        select case(trim(tag))
        case ('CLASS')
        !  call COMMANDS_DEFAULT(myCMD)
        !  write(*,*) 'startCOMMAND ',trim(CLASS_ID),' loc=',POINTER(myCMD)
          myCMD%ID=trim(CLASS_ID)
        case ('FLOAT')
          call ADD_FIELD(myCMD,trim(PARAM_ID),'F',isize,error)
        case ('INT')
          call ADD_FIELD(myCMD,trim(PARAM_ID),'I',isize,error)
        case ('STRING')
          call ADD_FIELD(myCMD,trim(PARAM_ID),'S',isize,error)
        case ('NOYES')
          call ADD_FIELD(myCMD,trim(PARAM_ID),'I',isize,error)
        case ('SIGN')
          call ADD_FIELD(myCMD,trim(PARAM_ID),'I',isize,error)
        case ('ENUM')
          call ADD_FIELD(myCMD,trim(PARAM_ID),'I',isize,error)
        case ('VRANGE')
          call ADD_FIELD(myCMD,trim(PARAM_ID),'VRANGE',isize,error)
        end select
  !      if (trim(tag).ne.'CLASS') then
  !        write(*,*) 'startCOMMAND done ',trim(tag),'.',trim(PARAM_ID)
  !        call ListCommandPar(myCMD)
  !        read(*,*)
  !      endif

      !  write(*,*) '          n=',myCMD%NARG
      !  read(*,*)
      end subroutine startCOMMAND

!--------------------------------------------------------
      subroutine datacls( tag, xmldata, error )
!--------------------------------------------------------
        character(len=*)                  :: tag
        character(len=*)                  :: xmldata
        logical                           :: error
      ! process only non-empty data inside class entries
        if ((CLS_INDEX.gt.0).and.(len_trim(xmldata).gt.0)) then
          select case(trim(tag))
          case('ENUMITEM')
            if (ENUM_INDEX.gt.0) call APPENDPAR(ENUM_TYPE_LIST(ENUM_INDEX),':',trim(adjustl(xmldata)))
          case('CLASS')
            if (isCOMMAND) myCMD%HINT=trim(adjustl(xmldata))
          end select
        endif
      end subroutine datacls

!--------------------------------------------------------
      subroutine endcls( tag, error )
!--------------------------------------------------------
      character(len=*)                  :: tag
      logical                           :: error
      ! process only non-empty data inside class entries
      if (CLS_INDEX.gt.0) then
        select case(trim(tag))
        case('CLASS')
          if (isCOMMAND) then
            myCMD%ICLS=CLS_INDEX
            myCMD%NAME=trim(CLASS_NAME)
            call ADD_COMMAND(myCMD,error)
            call COMMANDS_DEFAULT(myCMD)
          endif
        !  write(*,*) 'CLASSES [',trim(CLS_IDS),']',CLS_INDEX
          CLS_INDEX=0
        case('CLASSES')
          isCOMMAND=.false.
        end select
      endif

      end subroutine endcls


      end MODULE CLASSES