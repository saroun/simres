!//////////////////////////////////////////////////////////////////////
!////  $Id: xmlwriter.f,v 1.47 2019/08/15 15:02:09 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2012, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.47 $
!////     $Date: 2019/08/15 15:02:09 $
!//////////////////////////////////////////////////////////////////////
!////
!////  XML writer for SIMRES
!////  Procedures needed to write SIMRES configuration in XML file
!////  Includes option for non-XML format (console)
!//////////////////////////////////////////////////////////////////////
      MODULE XMLWRITER
    !  use FIELDS
      use xmlparse
      USE CLASSES
      use COMMANDS
      USE COMPONENTS
      use VCTABLE
      use VCTABLE_COMMANDS
      USE FILETOOLS
      use INSTCONTROL
      use XMLINFO
      use XMLUTILS
      implicit none

      private dbgreport
      contains

      subroutine dbgreport
        TYPE(TCLASS) :: OBJ
        character(LEN_ID) :: CLSID,OBJID
        character(LEN_NAME) :: OBJNAME
        OBJ=PCLASS(SAMOBJ)
        call GetClassAttrib(OBJ,CLSID,OBJID,OBJNAME)
        write(*,*) ' XMLWRITER sample=',OBJ%ICLS,OBJ%TCLS,trim(CLSID),' ',trim(OBJID),' ',trim(OBJNAME)
      end subroutine dbgreport


C---------------------------------------------------------
      SUBROUTINE save_config(fname)
C write the whole instrument to XML
C---------------------------------------------------------
      character(*)     :: fname
      character(MAX_FNAME_LENGTH) :: FRES
      integer :: IRES
      type(XML_PARSE) :: info
      logical :: std=.false.
      std=(trim(fname).eq.'STDIN')
      if (.not.std) write(*,*) 'save_config to '//trim(fname)

      if (std) then
        call DLG_FILEOPEN(FNAME,' ','xml',0,1,IRES,FRES)
      else
        call DLG_FILEOPEN(FNAME,OUTPATH,'xml',0,1,IRES,FRES)
      endif
      if (IRES.gt.0) then
      !  write(*,*) 'save_config OK 1'
        call INST_ADJUST(IRES,'save_config')
        call xml_open(info,trim(FRES), .false.)
        if (.not.info%error) call xml_options(info,report_lun=XML_STDOUT,report_details=.false.)
        if (.not.info%error) call EXPORT2XML(info,(.not.std))
        call xml_close(info)
        if (std) write(*,*) 'cfg updated'
        if ((.not.info%error).and.(.not.std)) call MSG_INFO('Instrument setup saved to '//trim(FRES),1)
      endif
      end subroutine save_config

C---------------------------------------------------------
      SUBROUTINE save_component(fname,id)
C write a component to XML file
C---------------------------------------------------------
      INCLUDE 'config.inc'
      character(*),intent(in)     :: fname,id
      character(MAX_FNAME_LENGTH) :: FRES
      integer :: IRES
      type(XML_PARSE) :: info
      character(len=128), dimension(2,64)    :: attribs
      character(32) :: GROUPTAG
      TYPE(TCLASS) :: OBJ
      TYPE(TCOMPID) :: COMPID
      logical :: std=.false.
      std=(trim(fname).eq.'STDIN')
      call PCLASS_CLEAR(OBJ)
      OBJ=getClassByID(trim(id))
      ! write(*,*) 'save_component ',trim(id)
      SELECT CASE(OBJ%TCLS)
      case(TCLS_COM,TCLS_SAM,TCLS_INS)
        if (std) then
          call DLG_FILEOPEN(FNAME,' ','xml',0,1,IRES,FRES)
        else
          call DLG_FILEOPEN(FNAME,OUTPATH,'xml',0,1,IRES,FRES)
        endif
        !write(*,*) 'save_component to ',trim(FRES), ires,std

        if (IRES.gt.0) then
          call INST_ADJUST(IRES,'save_component')
          call xml_open(info,trim(FRES), .false.)
          if (.not.info%error) call xml_options(info,report_lun=XML_STDOUT,report_details=.false.)
          info%error=.false.
          info%tag='SIMRES'
          attribs(1,1)='version'
          attribs(2,1)=trim(PACKAGE_VERSION)
          call xml_put(info,attribs,1,' ','open')
    ! group start tag
          info%tag='INSTRUMENT'
          attribs(1,1)='redefine'
          attribs(2,1)='no'
          call xml_put(info,attribs,1,' ','open')
          GROUPTAG=' '
          select case(OBJ%TCLS)
          case(TCLS_COM)
            call GetComponentByID(OBJ%P_COM%ID,COMPID)
            select case (COMPID%IBEAM)
              case(1)
                GROUPTAG='PRIMARY'
              case(2)
                GROUPTAG='SECONDARY'
            end select
          case(TCLS_SAM)
            GROUPTAG='SPECIMEN'
          case(TCLS_INS)
             GROUPTAG='INTERFACE'
          END select
          info%tag=trim(GROUPTAG)
          call xml_put(info,attribs,0,' ','open')
   ! print component
          if (.not.info%error) then
          ! exclude DIST when saving a single component
            if ((OBJ%TCLS==TCLS_COM).or.(OBJ%TCLS==TCLS_SAM)) then
              call OBJ2XML(info,OBJ,.true.,'|DIST|')
            else
              call OBJ2XML(info,OBJ,.true.,'')
            endif
          endif
   ! group end tag
          info%tag=trim(GROUPTAG)
          call xml_put(info,attribs,0,' ','close')
          info%tag='INSTRUMENT'
          call xml_put(info,attribs,0,' ','close')
          info%tag='SIMRES'
          call xml_put(info,attribs,0,' ','close')
          call xml_close(info)
          if ((.not.info%error).and.(.not.std)) call MSG_INFO(trim(id)//' saved to '//trim(FRES),1)
        endif
      end select
      end subroutine save_component

C---------------------------------------------------------
      SUBROUTINE save_commands(fname)
C write commands data to XML
C---------------------------------------------------------
      character(*)     :: fname
      character(MAX_FNAME_LENGTH) :: FRES
      integer :: IRES
      type(XML_PARSE) :: info
      call DLG_FILEOPEN(FNAME,OUTPATH,'xml',0,1,IRES,FRES)
      if (IRES.gt.0) then
        call xml_open(info,trim(FRES), .false.)
        if (.not.info%error) call xml_options(info,report_lun=XML_STDOUT,report_details=.false.)
        if (.not.info%error) call EXPORTCMDS(info)
        call xml_close(info)
        if (.not.info%error) call MSG_INFO('Commands saved to '//trim(FRES),1)
      endif
      end subroutine save_commands

C---------------------------------------------------------
      SUBROUTINE CLASS_PRN(IU,OBJ)
C write class to XML
C---------------------------------------------------------
      integer,intent(in) :: IU
      type(TCLASS),intent(in) :: OBJ
      type(XML_PARSE) :: info
      info%lun=IU
      call OBJ2XML(info,OBJ,.false.,'')
      end subroutine CLASS_PRN

!-----------------------------------------------------------------------
      subroutine simxml_put(isXML,info, attribs, no_attribs,xmldata,task)
! print component data
! isXML=true   .. as XML (call xml_put)
! isXML=false  .. as a text on console
!-----------------------------------------------------------------------
      logical,intent(in) :: isXML
      type(XML_PARSE),  intent(inout)               :: info
      character(len=*), intent(in), dimension(:,:)  :: attribs
      integer,          intent(in)                  :: no_attribs
      character(len=*), intent(in)                  :: xmldata
      character(len=*), intent(in)                  :: task
      integer :: i
      TYPE(PCLSDEF) :: CLS
      character(32) :: LVAL,NVAL
      logical :: error

1     format(50('-'),/,'parameters for ',a,' [',a,']',/,50('-'))
2     format('  parameters for ',a,':')
5     format(50('-'))

      if (isXML) then
        call xml_put(info, attribs, no_attribs,xmldata,task)
      else
        select case(task)
        case('open')
      ! update tagpath and parenttag
          call update_tagpath(info,'open')
          select case(trim(info%tag))
      ! base class of a family => print ID string and class name
          case('FRAME','SPECTROMETER','TRACING')
            do i=1,no_attribs
              select case(trim(attribs(1,i)))
              case('id')
                write(info%lun,1) trim(adjustl(attribs(2,i))),trim(adjustl(info%tag))
              end select
            enddo
      ! descendant classes
          case default
            call readAttrib(attribs,'length',LVAL,error)
            if (.not.error) then
            ! list field
              call readAttrib(attribs,'name',NVAL,error)
              write(info%lun,"(a10,' = list of ',a,' items',20x,a)") adjustl(info%tag),trim(LVAL),trim(NVAL)
            else
              write(info%lun,2) trim(adjustl(info%tag))
            endif
          end select
        case('elem')
          call element_put(isXML,info, attribs, no_attribs,xmldata)
        case('close')
          select case(trim(info%tag))
          case('FRAME','SPECTROMETER','TRACING')
            write(*,5)
          end select
      ! update tagpath and parenttag
          call update_tagpath(info,'close')
        end select
      endif
      end subroutine simxml_put



!---------------------------------------------------------
      SUBROUTINE FIELDS2XML(info,OBJ,CLS,isXML,exclude)
! Print to XML: group of fields from class definition CLS and data from OBJ
! Print the whole inheritance path, the base ancestor first
! exclude ... string in format |str1|str2|... lists field ID's to be excluded from export
!---------------------------------------------------------
      type(XML_PARSE) :: info
      TYPE(TCLASS) :: OBJ
      type(PCLSDEF) :: CLS
      character(*)     :: exclude
      logical,intent(in) :: isXML
      character(len=128), dimension(2,2)    :: attribs
      type(PCLSDEF) :: CLS1
      TYPE(TFIELD):: FIELD
      character(LEN_ID+1) :: prefix
      character(LEN_ID) :: clstag,id
      integer :: i,j,NA
1     format(' <!-- FIELDS2XML  ',a,'  -->')
2     format(' <!-- FIELDS2XML  ',a,':  ',a,'  -->')
3     format(' <!-- FIELDS2XML  ',a,':  ',I5,'  -->')
      if (.not.associated(CLS%C)) return
      !write(*,1) trim(CLS%C%ID)
      do i=1,CLS%C%NF
        if (associated(CLS%C%FIELDS(i)%F)) then
          FIELD%DEF=>CLS%C%FIELDS(i)%F
          id=trim(FIELD%DEF%ID)
          !write(*,2) trim(id),trim(FIELD%DEF%TID)
          if (index(exclude,'|'//trim(id)//'|')<=0) then
          !write(*,3) 'dim',FIELD%DEF%DIM
          !if (ALLOC_FIELD(FIELD,FIELD%DEF%DIM)) then
            !write(*,2) 'allocated',trim(id)
            select case (trim(FIELD%DEF%TID))
            case(ID_FIELD_CLASS)
        ! class field
              write(*,1) 'CLASSOBJ '//trim(FIELD%DEF%TID)
              call GetClassDef(trim(FIELD%DEF%CID),CLS1)
              if (associated(CLS1%C).and.(CLS1%C%NF.gt.0)) then
                attribs(1,1)='id'
                attribs(2,1)=trim(FIELD%DEF%ID)
                attribs(1,2)='name'
                attribs(2,2)=trim(FIELD%DEF%NAME)
                prefix=trim(FIELD%DEF%ID)//'.'
                clstag=trim(FIELD%DEF%ID)
                info%tag=clstag
                call simxml_put(isXML,info,attribs,2,' ','open')
                do j=1,CLS1%C%NF
                  if (associated(CLS1%C%FIELDS(j)%F)) then
                    FIELD%DEF=>CLS1%C%FIELDS(j)%F
                    if (ALLOC_FIELD(FIELD,FIELD%DEF%DIM)) then
                      !if (INDEX(clstag,'VFOC').gt.0) then
                      !  write(*,*) 'FIELDS2XML 1 [',trim(prefix)//trim(FIELD%DEF%ID),']'
                      !endif
                      call PCLASS_OUT(OBJ,trim(prefix)//trim(FIELD%DEF%ID),FIELD,NA)
                      if (NA.gt.0) call writeFieldXML(FIELD,info,isXML)
                    endif
                  endif
                enddo
                info%tag=clstag
                call simxml_put(isXML,info,attribs,0,' ','close')
              endif
        ! standard field
            case default
              !if (INDEX(FIELD%DEF%ID,'VFOC').gt.0) then
              !    write(*,*) 'FIELDS2XML 2 [',trim(FIELD%DEF%ID),']'
             ! endif
              if (ALLOC_FIELD(FIELD,FIELD%DEF%DIM)) then
                call PCLASS_OUT(OBJ,FIELD%DEF%ID,FIELD,NA)
                if (NA.gt.0) call writeFieldXML(FIELD,info,isXML)
              endif
            end select
            call DISPOSE_FIELD(FIELD)
          !endif
          endif
        endif
      enddo
      end subroutine FIELDS2XML

C---------------------------------------------------------
      SUBROUTINE OBJ2XML(info,OBJ,isXML,exclude)
C Print object to info%lun in XML format
! Print the whole inheritance path, the base ancestor first
C---------------------------------------------------------
      type(XML_PARSE) :: info
      TYPE(TCLASS) :: OBJ
      logical,intent(in) :: isXML
      character(*)     :: exclude
      character(len=128), dimension(2,3)    :: attribs
      integer :: noattribs
      character(LEN_ID) :: CLSID,OBJID,TOPID
      character(LEN_NAME) :: OBJNAME
      integer :: PCLS(1:10)
      INTEGER :: IP,NCLS
      type(PCLSDEF) :: CLS
! collect all ancestors IDs in PCLS
      call GetClassParents(OBJ%ICLS,PCLS,NCLS)
      if (NCLS.le.0) then
        info%error=.true.
        return
      endif
      call GetClassAttrib(OBJ,CLSID,OBJID,OBJNAME)
  ! get attribs for the base class
    ! class
      attribs(1,1)='class'
      attribs(2,1)=trim(CLSID)
    ! ID
      attribs(1,2)='id'
      attribs(2,2)=trim(OBJID)
    ! name
      attribs(1,3)='name'
      attribs(2,3)=trim(OBJNAME)
      noattribs=3
  ! scan the inheritance path, starting with the base class
      IP=NCLS
      TOPID=' '
      DO WHILE(IP.GT.0)
    ! print start tag for the class
        call GetClassDef(PCLS(IP),CLS)
        if (associated(CLS%C)) then
          if (len_trim(TOPID).eq.0) TOPID=trim(CLS%C%ID)
      ! COMMAND and OPTION are abstract types => use them as the top tag names
          select case(trim(TOPID))
          case('COMMAND','OPTION')
            if (IP.eq.NCLS-1) then
              info%tag=trim(TOPID)
            else
              info%tag=trim(CLS%C%ID)
            endif
          case default
            info%tag=trim(CLS%C%ID)
          end select
          if (CLS%C%NF.gt.0) then
    ! write class header
            call simxml_put(isXML,info,attribs,noattribs,' ','open')
    ! scan fields
            call FIELDS2XML(info,OBJ,CLS,isXML,exclude)
            noattribs=0 ! no attributes for subclasses
          endif
        endif
        IP=IP-1
      enddo
    ! close tags for all subclasses
      do IP=1,NCLS
        call GetClassDef(PCLS(IP),CLS)
        if (associated(CLS%C)) then
          select case(trim(TOPID))
          case('COMMAND','OPTION')
            if (IP.eq.NCLS-1) then
              info%tag=trim(TOPID)
            else
              info%tag=trim(CLS%C%ID)
            endif
          case default
            info%tag=trim(CLS%C%ID)
          end select
          if (CLS%C%NF.gt.0) then
            call simxml_put(isXML,info,attribs,0,' ','close')
            noattribs=0
          endif
        endif
      enddo
      END subroutine OBJ2XML

C---------------------------------------------------------
      SUBROUTINE EXPORT2XML(info,redefine)
C write simres configuration to XML
C---------------------------------------------------------
      INCLUDE 'config.inc'
      type(XML_PARSE) :: info
      logical :: redefine
      TYPE(TCLASS) :: OBJ
      integer :: i
      character(len=128), dimension(2,64)    :: attribs
      info%error=.false.
      info%tag='SIMRES'
      attribs(1,1)='version'
      attribs(2,1)=trim(PACKAGE_VERSION)
      call xml_put(info,attribs,1,' ','open')
! COMMANDS
      info%tag='COMMANDS'
      call xml_put(info,attribs,0,' ','open')
      do i=1,COMMANDS_NUM
        OBJ=PCLASS(getCmdObj(i))
        call OBJ2XML(info,OBJ,.true.,'')
      enddo
      info%tag='COMMANDS'
      call xml_put(info,attribs,0,' ','close')
! OPTIONS
      info%tag='OPTIONS'
      call xml_put(info,attribs,0,' ','open')
      OBJ=PCLASS(POPTOBJ(TROPT))
      call OBJ2XML(info,OBJ,.true.,'')
      OBJ=PCLASS(POPTOBJ(REPOPT))
      call OBJ2XML(info,OBJ,.true.,'')
      info%tag='OPTIONS'
      call xml_put(info,attribs,0,' ','close')
! INSTRUMENT
      info%tag='INSTRUMENT'
      if (redefine) then
        call xml_put(info,attribs,0,' ','open')
      else
        attribs(1,1)='redefine'
        attribs(2,1)='no'
        call xml_put(info,attribs,1,' ','open')
      endif
      info%tag='CFGTITLE'
      call xml_put(info,attribs,0,trim(CFGTITLE),'elem')
      if (MONO_N.gt.0) then
        info%tag='MONOCHROMATORS'
        call xml_put(info,attribs,0,trim(MONOCHROMATORSID),'elem')
      endif
      if (ANAL_N.gt.0) then
        info%tag='ANALYZERS'
        call xml_put(info,attribs,0,trim(ANALYZERSID),'elem')
      endif
! INTERFACE
      info%tag='INTERFACE'
      call xml_put(info,attribs,0,' ','open')
      if (INSOBJ%ICLS.gt.0) then
        OBJ=PCLASS(INSOBJ)
        call OBJ2XML(info,OBJ,.true.,'')
      endif
      info%tag='INTERFACE'
      call xml_put(info,attribs,0,' ','close')
! SPECIMEN
      if (SAMOBJ%ICLS.gt.0) then
        info%tag='SPECIMEN'
        call xml_put(info,attribs,0,' ','open')
        OBJ=PCLASS(SAMOBJ)
        call OBJ2XML(info,OBJ,.true.,'')
        info%tag='SPECIMEN'
        call xml_put(info,attribs,0,' ','close')
      endif
! BEAMLINE
      if (BEAMLINE_NC(1).gt.0) then
        info%tag='PRIMARY'
        call xml_put(info,attribs,0,' ','open')
        do i=1,BEAMLINE_NC(1)
          OBJ=PCLASS(getComByIndex(i,1))
          call OBJ2XML(info,OBJ,.true.,'')
        enddo
        info%tag='PRIMARY'
        call xml_put(info,attribs,0,' ','close')
      endif
      if (BEAMLINE_NC(2).gt.0) then
        info%tag='SECONDARY'
        call xml_put(info,attribs,0,' ','open')
        do i=1,BEAMLINE_NC(2)
          OBJ=PCLASS(getComByIndex(i,2))
          call OBJ2XML(info,OBJ,.true.,'')
        enddo
        info%tag='SECONDARY'
        call xml_put(info,attribs,0,' ','close')
      endif
      info%tag='INSTRUMENT'
      call xml_put(info,attribs,0,' ','close')
      info%tag='SIMRES'
      call xml_put(info,attribs,0,' ','close')
      END subroutine EXPORT2XML


C---------------------------------------------------------
      SUBROUTINE EXPORTCMDS(info)
C write commands input data to XML
C---------------------------------------------------------
      INCLUDE 'config.inc'
      type(XML_PARSE) :: info
      TYPE(TCLASS) :: OBJ
      integer :: i
      character(len=128), dimension(2,64)    :: attribs
      info%error=.false.
      info%tag='SIMRES'
      attribs(1,1)='version'
      attribs(2,1)=trim(PACKAGE_VERSION)
      call xml_put(info,attribs,1,' ','open')
! COMMANDS
      info%tag='COMMANDS'
      call xml_put(info,attribs,0,' ','open')
      do i=1,COMMANDS_NUM
        OBJ=PCLASS(getCmdObj(i))
        call OBJ2XML(info,OBJ,.true.,'')
      enddo
      info%tag='COMMANDS'
      call xml_put(info,attribs,0,' ','close')
      info%tag='SIMRES'
      call xml_put(info,attribs,0,' ','close')
      END subroutine EXPORTCMDS


!-----------------------------------------------------------------------
      subroutine element_put(isXML,info, attribs, no_attribs,xmldata)
! print a parameter value as single element
! isXML=true   .. as XML (call xml_put)
! isXML=false  .. as a text on console
!-----------------------------------------------------------------------
      logical,intent(in) :: isXML
      type(XML_PARSE),  intent(inout)               :: info
      character(len=*), intent(in), dimension(:,:)  :: attribs
      integer,          intent(in)                  :: no_attribs
      character(len=*), intent(in)                  :: xmldata
      character(32) :: USTR,NSTR,VSTR,ESTR,LSTR
      integer :: i,IVAL,j,IL
      logical :: IsInteger
! list items
2     format(5x,a,'(',a,')',10x,a)
! enum
1     format(a10,' = ',a33,'  ',a)
! other
3     format(a10,' = ',a20,' ',a12,'  ',a)
! units
4     format('[',a,']')
      if (isXML) then
        call xml_put(info, attribs, no_attribs,xmldata,'elem')
      else
        USTR=' '
        NSTR=' '
        ESTR=' '
        VSTR=' '
        IVAL=-1
        IL=-1
        do i=1,no_attribs
          select case(trim(attribs(1,i)))
          case('units')
            write(USTR,4) trim(attribs(2,i))
          case('name')
            NSTR=trim(attribs(2,i))
          case('enum')
            ESTR=trim(attribs(2,i))
          case('value')
            VSTR=trim(attribs(2,i))
            if (IsInteger(VSTR,j)) IVAL=j
          case('length')
            LSTR=trim(attribs(2,i))
            if (IsInteger(LSTR,j)) IL=j
          end select
        enddo
!        write(*,*) 'element_put ',trim(NSTR),' [',trim(xmldata),']'
  !      write(info%lun,3) trim(info%tag),trim(VSTR),trim(USTR),trim(NSTR)
      ! list elements
        if (trim(info%tag).eq.'ITEM') then
          write(info%lun,2) trim(NSTR),trim(VSTR),trim(adjustl(xmldata))
      ! enumerator
        else if ((ESTR.ne.' ').and.(VSTR.ne.' ')) then
          call GetEnumString(ESTR,IVAL,USTR)
          write(info%lun,1) info%tag,adjustl(USTR),trim(NSTR)
        else
       ! other
          write(info%lun,3) info%tag,trim(adjustl(xmldata)),USTR,trim(NSTR)
        endif
       endif
      end subroutine element_put


!--------------------------------------------------------------------------
      subroutine writeFieldXML(FIELD,info,isXML)
! collect data for XML output of a parameter of OBJ and write it to info
!--------------------------------------------------------------------------
      type(TFIELD), intent(in) :: FIELD
      type(XML_PARSE),intent(out) :: info
      logical,intent(in)  :: isXML
      integer :: noattribs
      character(LEN_LINE)                  :: xmldata
      character(len=128), dimension(2,64)    :: attribs
      integer :: NA
      character(32) :: CNUM
      TYPE(TENUMDEF) :: EDEF
!      write(*,*) 'writeFieldXM ',trim(FIELD%ID),' tid=',trim(FIELD%DEF%TID)
  ! handle RANGE variables separately
      if (trim(FIELD%DEF%TID).eq.ID_FIELD_RANGE) then
        call writeRangeXML(FIELD,info,isXML)
        return
      else if (trim(FIELD%DEF%TID).eq.ID_FIELD_SELECT) then
        call writeSelectXML(FIELD,info,isXML)
        return
      else if (trim(FIELD%DEF%TID).eq.ID_FIELD_TABLE) then
        call writeTableXML(FIELD,info,isXML)
        return
      endif
  ! create xml attributes and data
      xmldata=' '
      noattribs=0
      info%error=.true.
  ! name attribute for all types
      attribs(1,1)="name";attribs(2,1)=trim(FIELD%DEF%NAME)
      noattribs=1
      select case(trim(FIELD%DEF%TID))
      case(ID_FIELD_FLOAT)
      ! units attribute for FLOAT variables
        attribs(1,2)="units"
        attribs(2,2)=trim(FIELD%DEF%UNITS)
        noattribs=2
      case(ID_FIELD_ENUM)
! NOTE: enum attribute is not used by parsers, but it is required by
! element_put in non-XML regime.
      ! enum attribute = enum type
        call GetEnumDef(FIELD%DEF%ENUM,EDEF)
        attribs(1,2)="enum";attribs(2,2)=trim(EDEF%ID)
      ! value attribute = numerical value
        call INT2STR(FIELD%EORD,CNUM)
        attribs(1,3)="value";attribs(2,3)=trim(CNUM)
        noattribs=3
      end select
      call FIELD2STR(FIELD,0,xmldata,NA)
  ! set tag = parameter ID string
      info%tag=trim(FIELD%ID)
  ! all is O.K.
      info%error=.false.
      call element_put(isXML,info, attribs, noattribs,xmldata)
      end subroutine writeFieldXML

!--------------------------------------------------------------------------
      subroutine writeRangeXML(FIELD,info,isXML)
! collect data for XML output of a parameter of OBJ and write it to info
!--------------------------------------------------------------------------
      type(TFIELD) :: FIELD
      type(XML_PARSE),intent(out) :: info
      logical,intent(in)  :: isXML
      integer :: noattribs
      character(LEN_LINE)                  :: xmldata
      character(len=128), dimension(2,64)    :: attribs
      integer :: NA,i
      character(32) :: CNUM
!      write(*,*) 'writeRangeXML ',trim(FIELD%ID),' assoc=',associated(FIELD%DEF),' NP=',FIELD%NP
      xmldata=' '
      info%error=.false.
      attribs(1,1)="name";attribs(2,1)=trim(FIELD%DEF%NAME)
      call INT2STR(FIELD%NP,CNUM)
      attribs(1,2)="length";attribs(2,2)=trim(CNUM)
      noattribs=2
      info%tag=trim(FIELD%ID)
!      write(*,*) 'writeRangeXML 2'
      call simxml_put(isXML,info, attribs, noattribs,xmldata,'open')
!      write(*,*) 'writeRangeXML 3'
      info%tag='ITEM'
!      write(*,*) 'writeRangeXML ',trim(FIELD%ID),' NP=',FIELD%NP
      do i=1,FIELD%NP
        attribs(1,1)="name";attribs(2,1)=trim(FIELD%ID)
        call INT2STR(i,CNUM)
        attribs(1,2)="value";attribs(2,2)=trim(CNUM)
        call FIELD2STR(FIELD,i,xmldata,NA)
  !      write(*,*) 'writeRangeXML ',i,' ',trim(xmldata),' NA=',NA
        call element_put(isXML,info, attribs,noattribs,xmldata)
      enddo
      info%tag=trim(FIELD%ID)
      call simxml_put(isXML,info, attribs, noattribs,xmldata,'close')
      end subroutine writeRangeXML

!--------------------------------------------------------------------------
      subroutine writeSelectXML(FIELD,info,isXML)
! write SELECT type data to XML
!--------------------------------------------------------------------------
      type(TFIELD) :: FIELD
      type(XML_PARSE),intent(out) :: info
      logical,intent(in)  :: isXML
      integer :: noattribs
      character(LEN_LINE)                  :: xmldata
      character(len=128), dimension(2,64)    :: attribs
      integer :: i
      character(32) :: CNUM
      xmldata=' '
      info%error=.false.
      attribs(1,1)="name";attribs(2,1)=trim(FIELD%DEF%NAME)
      call INT2STR(FIELD%NP,CNUM)
      attribs(1,2)="length";attribs(2,2)=trim(CNUM)
      call INT2STR(FIELD%EORD,CNUM)
      attribs(1,3)="selected";attribs(2,3)=trim(CNUM)
      noattribs=3
      info%tag=trim(FIELD%ID)
      call simxml_put(isXML,info, attribs, noattribs,xmldata,'open')
      info%tag='ITEM'
      do i=1,FIELD%NP
      !  write(*,*) 'writeSelectXML ITEM ',i,trim(FIELD%SFIELD(i))
        call element_put(isXML,info, attribs,0,FIELD%SFIELD(i))
      enddo
      info%tag=trim(FIELD%ID)
      call simxml_put(isXML,info, attribs, noattribs,xmldata,'close')
      end subroutine writeSelectXML


!--------------------------------------------------------------------------
      subroutine writeTableXML(FIELD,info,isXML)
! write SELECT type data to XML
!--------------------------------------------------------------------------
      type(TFIELD) :: FIELD
      type(XML_PARSE),intent(out) :: info
      logical,intent(in)  :: isXML
      integer :: noattribs
      character(LEN_LINE)                  :: xmldata
      character(len=128), dimension(2,64)    :: attribs
      integer :: i,na
      character(32) :: CNUM
      character(LEN_LINE) :: LINE
      xmldata=' '
      !write(*,*) 'writeTableXML ',FIELD%NP

      info%error=.false.
      attribs(1,1)="name";attribs(2,1)=trim(FIELD%DEF%NAME)
      call INT2STR(FIELD%NP,CNUM)
      attribs(1,2)="rows";attribs(2,2)=trim(CNUM)
      noattribs=2
      info%tag=trim(FIELD%ID)
      call simxml_put(isXML,info, attribs, noattribs,xmldata,'open')
      info%tag='ITEM'
      do i=1,FIELD%NP
        call ROW2STR(FIELD,i,0,LINE,na)
        !write(*,*) 'writeTableXML ',trim(LINE)
        call element_put(isXML,info, attribs,0,trim(LINE))
      enddo
      info%tag=trim(FIELD%ID)
      call simxml_put(isXML,info, attribs, noattribs,xmldata,'close')
      end subroutine writeTableXML


      end module XMLWRITER