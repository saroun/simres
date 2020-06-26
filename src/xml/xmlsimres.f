!//////////////////////////////////////////////////////////////////////
!////  $Id: xmlsimres.f,v 1.92 2019/04/01 12:59:33 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.92 $
!////     $Date: 2019/04/01 12:59:33 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Module for parsing instrument definition from XML data
!////
!//////////////////////////////////////////////////////////////////////
      MODULE XMLSIMRES
      use XMLINFO
      USE XMLUTILS
      USE FILETOOLS
   !  use FIELDS
      USE CLASSES
      use COMMANDS
      USE VCTABLE_COMMANDS
      USE VCTABLE_COMPONENTS
      use xmlparse
      USE XMLHANDLER
      USE COMPONENTS
      USE COMPONENTS_IO
      USE INSTCONTROL
      use XTALS

      implicit none

      private
      save

      logical :: DBG=.false.

      type(TCLASS) :: myOBJ
      type(TCOMOBJ) :: myCOMOBJ
! local instances for all component class which might be found in the instrument definition
  ! components
      !TYPE(TFRAME) :: myFRAME
      !TYPE(TFRAME) :: mySLIT
      !TYPE(SOURCE) :: mySOURCE
      !TYPE(DETECTOR) :: myDETECTOR
      !TYPE(BENDER) :: myGUIDE
      !TYPE(CRYSTAL) :: myCRYSTAL
      !TYPE(XTAL) :: myXTAL
  ! samples
      TYPE(TSAMPLE) :: mySAMPLE
      TYPE(TPCRYST) :: myPCRYST
      TYPE(TSCRYST) :: mySCRYST
  ! interfaces
      TYPE(TSPEC)   :: mySPEC
  !    TYPE(TPWDIF)  :: myPWDIF
  !    TYPE(TTAS)    :: myTAS
  ! options
      TYPE(TTROPT)   :: myTROPT
      TYPE(TREPOPT)   :: myREPOPT

      type(TCOMMAND) :: myCMD

  ! commands


      integer :: ibeam   ! current beamline: (1) primary, (2) secondary, (0) sample
! publish only this:
      public parse_xmlclasses,parse_xmlcrystals
      public load_config

      contains


!--------------------------------------------------------
      subroutine parse_xmlcrystals(fname,ierr)
! parse XML file with definition of crystals (crystals.xml)
!--------------------------------------------------------
      use TABLE_CRYSTALS
      use xmlparse
      use FILETOOLS
      character(*),intent(in)     :: fname
      integer,intent(out) :: ierr
      logical          :: error
      character(MAX_FNAME_LENGTH) :: FRES,FN
      integer :: IRES
      ierr=1
  ! prepare prompt and current filename
      if (len_trim(FNAME).le.0) then
        FN=adjustl(XCONFIG_CRYSTALS)
      else
        FN=adjustl(FNAME)
      endif
  ! run dialog for file selection
      CALL DLG_FILEOPEN(FNAME,'|'//trim(RESPATH)//'|'//trim(CFGPATH),'xml',0,0,IRES,FRES)
      ! write(*,*) 'parse_xmlcrystals DLG_FILEOPEN IRES=',IRES,' ',trim(FRES)
      if (IRES.gt.0) then
        call xml_process_file(trim(FRES),startCRTABLE, dataCRTABLE, endCRTABLE,error )
        if (.not.error) then
          ierr=0
          call MSG_INFO('Definition of crystals loaded from '//trim(FRES),1)
          XCONFIG_CRYSTALS=trim(FRES)
        else
          write(*,*) 'parse_xmlclasses processed error=',error
          call CLEAR_TABLE_CRYSTALS
        endif
      endif
      !write(*,*) 'parse_xmlclasses ierr==',ierr
      end subroutine parse_xmlcrystals

!--------------------------------------------------------
      subroutine parse_xmlclasses(fname,ierr)
! parse XML file with classes definition (classes.xml)
!--------------------------------------------------------
      use CLASSDEF,ONLY:startcls,datacls,endcls,CLEAR_CLASSES
      USE COMMANDS,ONLY: ALLOCATE_COMMANDS,COMMANDS_PRN
      use xmlparse
      use FILETOOLS
      character(*)     :: fname
      integer,intent(out) :: ierr
      logical          :: error
      character(MAX_FNAME_LENGTH) :: FRES,FN
      integer :: IRES
      ierr=1
  ! prepare prompt and current filename
      if (len_trim(FNAME).le.0) then
        FN=adjustl(XCONFIG_CLASSES)
      else
        FN=adjustl(FNAME)
      endif
  ! run dialog for file selection
      CALL DLG_FILEOPEN(FNAME,'|'//trim(RESPATH)//'|'//trim(CFGPATH),'xml',0,0,IRES,FRES)
      !write(*,*) 'parse_xmlclasses DLG_FILEOPEN IRES=',IRES
      if (IRES.gt.0) then
        call xml_process_file(trim(FRES),startcls, datacls, endcls,error )
        !write(*,*) 'parse_xmlclasses processed error=',error
        if (.not.error) then
          ierr=0
          call MSG_INFO('Definition of classes loaded from '//trim(FRES),1)
          XCONFIG_CLASSES=trim(FRES)
        else
          call CLEAR_CLASSES
        endif
      endif
      if (.not.error) then
        call ALLOCATE_COMMANDS
      endif
      !write(*,*) 'parse_xmlclasses error =',error

      ! if(error) read(*,*)
      end subroutine parse_xmlclasses

!--------------------------------------------------------
      subroutine load_config(fname,ierr)
! parse XML file with instrument definition
!--------------------------------------------------------
      character(*),intent(in)     :: fname
      integer,intent(out) :: ierr
      logical          :: error
      character(MAX_FNAME_LENGTH) :: FRES,FN
      integer :: IRES,IPROMPT,i
      logical :: std=.false.
      std=(trim(fname).eq.'STDIN')
    !  write(*,*) 'load_config '//trim(fname)
      IERR=1
      if (len_trim(FNAME).le.0) then
        FN=adjustl(XCONFIG_NAME)
        IPROMPT=1
      else if (std) then
        FN='STDIN'
      ELSE
        FN=adjustl(FNAME)
        IPROMPT=0
      endif
      CALL DLG_FILEOPEN(trim(FN),'|'//trim(RESPATH)//'|'//trim(CFGPATH),'xml',IPROMPT,0,IRES,FRES)
      if (IRES.gt.0) then
    ! reset parser !!
        XHND=' '
        DOPARSING=.false.
    ! parse the file
        call xml_process_file(trim(FRES), startMAIN, dataMAIN, endMAIN, error )
        !write(*,*) 'load_config process OK',error
    ! tracing options must be defined
    ! create default if not read from XML file)
        if (.not.associated(TRAOBJ%P_TRA)) then
          myOBJ=CREATE_CLASS('TRACING','TR','tracing options')
          i=AddOptions(myOBJ%P_OPT)
          call PCLASS_CLEAR(myOBJ)
          call MSG_WARN('Tracing options not defined in '//trim(FRES)//': default is used',1)
        endif
    ! reports options must be defined
    ! create default if not read from XML file)
        if (.not.associated(REPOBJ%P_REP)) then
          myOBJ=CREATE_CLASS('REPORTS','REP','reports options')
          i=AddOptions(myOBJ%P_OPT)
          call PCLASS_CLEAR(myOBJ)
          call MSG_WARN('Report options not defined in '//trim(FRES)//': default is used',1)
        endif
    ! interface must be defined
    ! create default if not read from XML file)
        if (.not.associated(INSOBJ%P_SPEC)) then
          myOBJ=CREATE_CLASS('SPECTROMETER','INS','default interface')
          i=AddInterface(myOBJ%P_INS)
          call PCLASS_CLEAR(myOBJ)
          call MSG_WARN('Instrument interface not defined in '//trim(FRES)//': default is used',1)
        endif
    ! there must be a source
        if (SOURCES_NC.le.0) then
          call MSG_ERROR('LOAD_CONFIG','No source is defined in this setup.',1,1)
        endif
        if (.not.error) then
          !write(*,*) 'load_config REDEFINE=',REDEFINE
        !  write(*,*) 'load_config INSOBJ%P_SPEC%KI=',INSOBJ%P_SPEC%KI
          if (REDEFINE) then
      ! new instrument definition => clear all datasets
            call DSET_CLEAR
      ! define the 1st dataset. INST_ADJUST and BEAM_OUT are called inside DEFINE_DSET
            call DEFINE_DSET(1)
            IERR=0
            if (.not.std) then
              call MSG_INFO('Instrument setup loaded from '//trim(FRES),1)
              XCONFIG_NAME=trim(FRES)
            endif
            write(*,*) 'load_config '//trim(fname)//' REDEFINE OK'
          else
            call DEFINE_DSET(DSET_ISEL)
            write(*,*) 'load_config '//trim(fname)//' update OK'
          endif
        endif
      else
        write(*,*) 'load_config, can''t read file '//trim(fname)
      endif
  ! only one dataset and subset is considered in this version:
      CSUBSET(1)=1
      CDATASET(1)=1
      CSUBSET(2)=1
      CDATASET(2)=1
      call CLEAR_XHNDARG
      !call dbgreport
      end subroutine load_config

      subroutine dbgreport
        TYPE(TCLASS) :: OBJ
        character(LEN_ID) :: CLSID,OBJID
        character(LEN_NAME) :: OBJNAME
        OBJ=PCLASS(SAMOBJ)
        call GetClassAttrib(OBJ,CLSID,OBJID,OBJNAME)
  !      write(*,*) ' XMLSIMRES sample=',OBJ%ICLS,OBJ%TCLS,trim(CLSID),' ',trim(OBJID),' ',trim(OBJNAME)
      end subroutine dbgreport

!--------------------------------------------------------
      TYPE(TCLASS) function CREATE_CLASS(CID,IDNAME,CNAME)
! return a local instance of a class for given class ID string
!--------------------------------------------------------
      character(*) :: CID,IDNAME,CNAME
      type(TCLASS) :: OBJ
      call PCLASS_CLEAR(OBJ)
      if (PCOM_IS_COMPONENT(CID)) then
        myCOMOBJ=PCOM_CREATE(CID,IDNAME,CNAME)
        OBJ=PCLASS(myCOMOBJ)
        !write(*,*) 'CREATE_CLASS ',trim(CID),' ',myCOMOBJ%ICLS,myCOMOBJ%CCLS,myCOMOBJ%INST,' ',trim(myCOMOBJ%ID)
      else
      select case(CID)
  ! samples
        case('SAMPLE')
          OBJ=PCLASS(PSAMOBJ(mySAMPLE))
        case('PCRYST')
          OBJ=PCLASS(PSAMOBJ(myPCRYST))
        case('SCRYST')
          OBJ=PCLASS(PSAMOBJ(mySCRYST))
  ! interfaces
        case('SPECTROMETER')
          OBJ=PCLASS(PINSOBJ(mySPEC))
        case('PWD')
  !    TYPE(TPWDIF)  :: myPWDIF
        case('TAS')
  !    TYPE(TTAS)    :: myTAS
  ! options
        case('TRACING')
          OBJ=PCLASS(POPTOBJ(myTROPT))
        case('REPORTS')
          OBJ=PCLASS(POPTOBJ(myREPOPT))
  ! commands (return global instance)
        case('COMMAND')
          OBJ=CREATE_CMDCLASS(IDNAME)
      end select
      endif
      if (OBJ%ICLS.gt.0) then
        call PCLASS_DEFAULT(OBJ)
        call SetClassAttrib(OBJ,IDNAME,CNAME)
      endif
      CREATE_CLASS=OBJ
      end function CREATE_CLASS

!--------------------------------------------------------
      subroutine CLEAR_TMP_OBJECT
! clear all data associated with temporary object myOBJ
! including namespace
!--------------------------------------------------------
        if (REDEFINE) then
          call PCLASS_DISPOSE(myOBJ)
        else
          call PCLASS_CLEAR(myOBJ)
        endif
        XHND_NAMES=' '
        XHND_ICLS=0
      end subroutine CLEAR_TMP_OBJECT

!--------------------------------------------------------
      TYPE(TCLASS) function CREATE_CMDCLASS(CID)
! creating commands is different, requires special soubroutine !
!--------------------------------------------------------
      character(*),intent(in) :: CID
      type(TCLASS) :: OBJ
      TYPE(PCLSDEF) :: CDEF
      TYPE(TCMDOBJ) :: COBJ
      logical :: error
      call PCLASS_CLEAR(OBJ)
      call getClassDef(trim(CID),CDEF)
      if (associated(CDEF%C)) then
        if (CDEF%C%isCMD) then
          myCMD%ICLS=CDEF%C%IDX
          call ALLOC_COMMAND(myCMD,CDEF%C,error)
          if (.not.error) then
             COBJ=PCMDOBJ(myCMD)
             OBJ=PCLASS(COBJ)
          endif

        endif
      endif
      CREATE_CMDCLASS=OBJ
      end function CREATE_CMDCLASS

!************************************************************************
! Parser handlers - INSTRUMENT
!************************************************************************

!--------------------------------------------------------
      subroutine startFRAME(tag,attribs,error)
! handle START tags in the FRAME section
!--------------------------------------------------------
      character(len=*)                  :: tag
      character(len=*), dimension(:,:)  :: attribs
      logical                           :: error
      call startDEFAULT(tag,attribs,error)
      end subroutine startFRAME

!--------------------------------------------------------
      subroutine endFRAME(tag,error)
! handle END tags in the FRAME section
!--------------------------------------------------------
      character(len=*)                  :: tag
      logical  :: error
      integer :: i
      select case(trim(tag))
      case('FRAME')
        if (REDEFINE) then
          select case(myOBJ%TCLS)
          case(TCLS_COM)
            i=AddComponent(myOBJ%P_COM,IBEAM)
          case(TCLS_SAM)
            i=AddSpecimen(myOBJ%P_SAM)
          end select
        endif
        !write(*,*) ' endFRAME ',trim(myOBJ%ID),' ',trim(myOBJ%P_COM%ID),' ',XHND_ICLS
        call CLEAR_TMP_OBJECT
      end select
      call endDEFAULT(myOBJ,tag,error)
      end subroutine endFRAME

!--------------------------------------------------------
      subroutine startCOMMANDS(tag,attribs,error)
! handle START tags in the COMMANDS section
!--------------------------------------------------------
      character(len=*)                  :: tag
      character(len=*), dimension(:,:)  :: attribs
      logical                           :: error
      character(LEN_ID) :: CID
      select case(trim(tag))
      case('COMMAND')
        call readAttrib(attribs,'id',CID,error)
        if (error) return
        myOBJ=CREATE_CLASS('COMMAND',CID,CID)
        XHND_ICLS=myOBJ%ICLS
        if (myOBJ%ICLS.gt.0) then
          call CLASSES_NAMESPACE(myOBJ%ICLS,XHND_NAMES)
        else
          call MSG_ERROR('XMLSIMRES','Unknown command: '//trim(CID),0,1)
          error=.true.
        endif
      case default
        call startDEFAULT(tag,attribs,error)
      end select
      end subroutine startCOMMANDS


!--------------------------------------------------------
      subroutine endCOMMANDS(tag,error)
! handle END tags in the COMMANDS section
!--------------------------------------------------------
      character(len=*)                  :: tag
      logical  :: error
      select case(trim(tag))
      case('COMMAND')
        call ASSIGN_COMMAND(myOBJ%P_CMD%P_CMD)
        call CLEAR_TMP_OBJECT
      end select
      call endDEFAULT(myOBJ,tag,error)
      end subroutine endCOMMANDS

!--------------------------------------------------------
      subroutine startOPTIONS(tag,attribs,error)
! handle START tags in the INTERFACE section
!--------------------------------------------------------
      character(len=*)                  :: tag
      character(len=*), dimension(:,:)  :: attribs
      logical                           :: error
      character(LEN_ID) :: CLSNAME
      character(LEN_ID) :: IDNAME
      character(LEN_NAME) :: CNAME
      select case(trim(tag))
      case('OPTION')
        call readClassAttrib(attribs,CLSNAME,IDNAME,CNAME,error)
        if (error) return
        myOBJ=CREATE_CLASS(CLSNAME,IDNAME,CNAME)
        XHND_ICLS=myOBJ%ICLS
        if (myOBJ%ICLS.gt.0) then
          call CLASSES_NAMESPACE(myOBJ%ICLS,XHND_NAMES) ! define namespace
        else
          call MSG_ERROR('XMLSIMRES','Unknown options class: '//trim(CLSNAME),0,1)
          error=.true.
        endif
      case default
        call startDEFAULT(tag,attribs,error)
      end select
      end subroutine startOPTIONS


!--------------------------------------------------------
      subroutine endOPTIONS(tag,error)
! handle END tags in the FRAME section
!--------------------------------------------------------
      character(len=*)                  :: tag
      logical  :: error
      integer :: i
      select case(trim(tag))
      case('OPTION')
        i=AddOptions(myOBJ%P_OPT)
        call CLEAR_TMP_OBJECT
      end select
      call endDEFAULT(myOBJ,tag,error)
      end subroutine endOPTIONS

!--------------------------------------------------------
      subroutine startINTERFACE(tag,attribs,error)
! handle START tags in the INTERFACE section
!--------------------------------------------------------
      character(len=*)                  :: tag
      character(len=*), dimension(:,:)  :: attribs
      logical                           :: error
      character(LEN_ID) :: CLSNAME
      character(LEN_ID) :: IDNAME
      character(LEN_NAME) :: CNAME
      select case(trim(tag))
      case('SPECTROMETER')
        call readClassAttrib(attribs,CLSNAME,IDNAME,CNAME,error)
        if (error) return
        if (REDEFINE) then
          myOBJ=CREATE_CLASS(CLSNAME,IDNAME,CNAME)
        else
          myOBJ=PCLASS(INSOBJ)
        endif
        XHND_ICLS=myOBJ%ICLS
      !  write(*,*) ' startINTERFACE ',trim(myOBJ%ID),' ',XHND_ICLS
        if (myOBJ%ICLS.gt.0) then
          call CLASSES_NAMESPACE(myOBJ%ICLS,XHND_NAMES) ! define namespace for SPECTROMETER
        else
          call MSG_ERROR('XMLSIMRES','Unknown interface '//trim(CLSNAME)//' name='//trim(IDNAME),0,1)
          error=.true.
        endif
      case default
        call startDEFAULT(tag,attribs,error)
      end select
      end subroutine startINTERFACE

!--------------------------------------------------------
      subroutine endINTERFACE_(tag,error)
! handle END tags in the FRAME section
!--------------------------------------------------------
      character(len=*)                  :: tag
      logical  :: error
      integer :: i
      select case(trim(tag))
      case('SPECTROMETER')
        if (REDEFINE) i=AddInterface(myOBJ%P_INS)
      !  write(*,*) ' endINTERFACE ',trim(myOBJ%ID),' ',XHND_ICLS
        call CLEAR_TMP_OBJECT
      end select
      call endDEFAULT(myOBJ,tag,error)
      end subroutine endINTERFACE_

!--------------------------------------------------------
      subroutine startCOMPONENTS(tag,attribs,error)
! handle START tags in the PRIMARY,SECONDARY and SPECIMEN sections
!--------------------------------------------------------
      character(len=*)                  :: tag
      character(len=*), dimension(:,:)  :: attribs
      logical                           :: error
      character(LEN_ID) :: CLSNAME
      character(LEN_ID) :: IDNAME
      character(LEN_NAME) :: CNAME
      select case(trim(tag))
      case('FRAME')
        call readClassAttrib(attribs,CLSNAME,IDNAME,CNAME,error)
        if (error) return
        if (REDEFINE) then
          myOBJ=CREATE_CLASS(CLSNAME,IDNAME,CNAME)
        else
          myOBJ=getClassByID(trim(IDNAME))
        endif
        XHND_ICLS=myOBJ%ICLS
        if (myOBJ%ICLS.gt.0) then
          call XHNDPass('FRAME')
          call CLASSES_NAMESPACE(myOBJ%ICLS,XHND_NAMES)
          ! write(*,*) ' startCOMPONENTS ',trim(myOBJ%ID),' ',myOBJ%ICLS
          ! write(*,*) '  names=',trim(XHND_NAMES)
        else
          call MSG_ERROR('XMLSIMRES','Unknown component class: '//trim(CLSNAME),0,1)
          error=.true.
        endif
      case default
        call startDEFAULT(tag,attribs,error)
      end select
      end subroutine startCOMPONENTS

!--------------------------------------------------------
      subroutine endCOMPONENTS(tag,error)
! handle END tags in the PRIMARY,SECONDARY and SPECIMEN sections
!--------------------------------------------------------
      character(len=*)                  :: tag
      logical                           :: error
      select case(trim(tag))
      case('PRIMARY','SECONDARY','SPECIMEN')
      !  write(*,*) ' endCOMPONENTS ',trim(tag),' ',XHND_ICLS
        call XHNDReturn
      end select
      end subroutine endCOMPONENTS

!--------------------------------------------------------
      subroutine startINSTRUMENT(tag,attribs,error)
! handle START tags in the INSTRUMENT section
!--------------------------------------------------------
      character(len=*)                  :: tag
      character(len=*), dimension(:,:)  :: attribs
      logical                           :: error
      select case(trim(tag))
      case('PRIMARY','SECONDARY','SPECIMEN')
        select case(trim(tag))
        case('PRIMARY')
          ibeam=1
        case('SECONDARY')
          ibeam=2
        case('SPECIMEN')
          ibeam=0
        end select
        call XHNDPass('COMPONENTS')
  ! FRAME components can be outside PRIMARY or SECONDARY sections when REDEFINE=false
      case('FRAME')
        if (.not.REDEFINE) then
          call XHNDPass('COMPONENTS')
          call startCOMPONENTS(tag,attribs,error)
        endif
      case('INTERFACE')
        call XHNDPass('INTERFACE')
      case('CFGTITLE','MONOCHROMATORS','ANALYZERS')
       ! allow them, handled by dataMAIN
      case default
       ! nothing else to be handled
        ! call startDEFAULT(tag,attribs,error)
      end select
      end subroutine startINSTRUMENT

!--------------------------------------------------------
      subroutine endINSTRUMENT(tag,error)
! handle END tags in the INSTRUMENT section
!--------------------------------------------------------
      character(len=*)                  :: tag
      logical                           :: error
      integer :: ierr
      select case(trim(tag))
      case('INSTRUMENT')
        !write(*,*) 'endINSTRUMENT '
        call INST_ADJUST(ierr,'endINSTRUMENT')
        !write(*,*) 'endINSTRUMENT adjusted, ierr=',ierr
        if (ierr.ne.0) then
          error=.true.
          call MSG_ERROR('XMLSIMRES','Can''t adjust instrument to the imported parameters',0,1)
        endif
      end select
      call endDEFAULT(myOBJ,tag,error)
      end subroutine endINSTRUMENT

!--------------------------------------------------------
      subroutine startSIMRES(tag,attribs,error)
! handle START tags in the SIMRES section
!--------------------------------------------------------
      character(len=*)                  :: tag
      character(len=*), dimension(:,:)  :: attribs
      logical                           :: error
      select case(trim(tag))
      case('INSTRUMENT')
        IBEAM=0
        REDEFINE=readAttribBoolean(attribs,'redefine',.true.)
        if (REDEFINE) call DeleteComponents
        !write(*,*) 'startSIMRES REDEFINE=',REDEFINE
        call XHNDPass('INSTRUMENT')
      case('OPTIONS')
        call XHNDPass('OPTIONS')
      case('COMMANDS')
        call XHNDPass('COMMANDS')
      case('HTML')
          ! do nothing, used by GUI for comment output
      end select
      end subroutine startSIMRES

!--------------------------------------------------------
      subroutine endSIMRES(tag,error)
! handle END tags in the SIMRES section
!--------------------------------------------------------
      character(len=*)                  :: tag
      logical                           :: error
      select case(trim(tag))
      case('SIMRES')
        call CLEAR_XMLHANDLER
      case('HTML')
          ! do nothing, used by GUI for comment output
      end select
      if (error) DOPARSING=.false.
      end subroutine endSIMRES

!--------------------------------------------------------
      subroutine startMAIN(tag,attribs,error)
! dispatch to handlers for particular section
!--------------------------------------------------------
      character(len=*)                  :: tag
      character(len=*), dimension(:,:)  :: attribs
      logical                           :: error
      if (DOPARSING) then
  ! dispatch to child's handlers
    ! ignored tags:
        select case(trim(tag))
        case('GROUP')
          return
        end select
        select case(XHND)
        case('SIMRES')
          call startSIMRES(tag,attribs,error)
        case('OPTIONS')
          call startOPTIONS(tag,attribs,error)
        case('INSTRUMENT')
          call startINSTRUMENT(tag,attribs,error)
        case('INTERFACE')
          call startINTERFACE(tag,attribs,error)
        case('COMPONENTS')
          call startCOMPONENTS(tag,attribs,error)
        case('FRAME')
          call startFRAME(tag,attribs,error)
        case('COMMANDS')
          call startCOMMANDS(tag,attribs,error)
        end select
      else if (trim(tag).eq.'SIMRES') then
  ! SIMRES section starts
        call CLEAR_XMLHANDLER
        DOPARSING=.true.
        XHND='SIMRES'
      endif
      if (error) DOPARSING=.false.
      end subroutine startMAIN

!--------------------------------------------------------
      subroutine endMAIN(tag,error)
! dispatch to handlers for particular section
!--------------------------------------------------------
      character(len=*)                  :: tag
      logical                           :: error
      if (DOPARSING) then
      !write(*,*) 'endMAIN ',trim(XHND),' ',trim(tag)
    ! ignored tags:
        select case(trim(tag))
        case('GROUP')
          return
        end select
        select case(XHND)
        case('SIMRES')
          call endSIMRES(tag,error)
        case('OPTIONS')
          call endOPTIONS(tag,error)
        case('INSTRUMENT')
          call endINSTRUMENT(tag,error)
        case('INTERFACE')
          call endINTERFACE_(tag,error)
        case('COMPONENTS')
          call endCOMPONENTS(tag,error)
        case('FRAME')
          call endFRAME(tag,error)
        case('COMMANDS')
          call endCOMMANDS(tag,error)
        end select
      endif
      if (error) DOPARSING=.false.
      end subroutine endMAIN

!--------------------------------------------------------
      subroutine dataMAIN(tag,xmldata,error)
! main data handler
!--------------------------------------------------------
      character(len=*)                  :: tag
      character(len=*)                  :: xmldata
      logical                           :: error
    ! process only non-empty data
      call PARSE_INPUT(xmldata)
      if ((.not.DOPARSING).or.(len_trim(xmldata).le.0)) RETURN

      if (trim(xmldata).eq.'undefined') then
        write(*,*) 'dataMAIN undefined ',trim(XHND)
      endif
      select case(XHND)
      case('INSTRUMENT')
          select case(trim(tag))
          case('CFGTITLE')
            CFGTITLE=trim(xmldata)
          case('MONOCHROMATORS')
            MONOCHROMATORSID=trim(xmldata)
          case('ANALYZERS')
            ANALYZERSID=trim(xmldata)
          case default
           ! INSTRUMENT has no other data
          end select
      case('INTERFACE')
        ! MONOCHROMATORS and ANALYZERS ae handled also here for backward compatibility
          select case(trim(tag))
          case('MONOCHROMATORS')
            MONOCHROMATORSID=trim(xmldata)
          case('ANALYZERS')
            ANALYZERSID=trim(xmldata)
          case default
             call dataDEFAULT(tag,xmldata,error)
          end select
      case('FRAME','OPTIONS','COMMANDS')
          call dataDEFAULT(tag,xmldata,error)
      case('HTML')
          ! do nothing
      end select
      if (error) DOPARSING=.false.
      end subroutine dataMAIN

      end module XMLSIMRES
