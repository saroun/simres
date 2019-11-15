!//////////////////////////////////////////////////////////////////////
!////  $Id: vctable_components.f,v 1.35 2019/08/15 15:02:06 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2011, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.35 $
!////     $Date: 2019/08/15 15:02:06 $
!//////////////////////////////////////////////////////////////////////
!////
!////   Virtual classes table for COMPONENTS
!////  - serves to emulation of polymorphism in Fortran
!//////////////////////////////////////////////////////////////////////
      MODULE VCTABLE_COMPONENTS
      use CLASSES
      use CLASSDEF,ONLY:GetICLS
      use FRAMES
      use FRAMES_TRACE
      use SOURCES
      use SOURCES_TRACE
      use DETECTORS
      use DETECTORS_TRACE
      use GUIDES
      use GUIDES_TRACE
      use SGUIDES
      use SGUIDES_TRACE
      use CRYSTALS
      use CRYSTALS_TRACE
      use XTALS
      use XTALS_TRACE
      use DCHOPPERS
      use DCHOPPERS_TRACE
      use MONITORS
      use MONITORS_TRACE
      ! use GENERATOR
      implicit none

! TCOMOBJ - abstract type with pointer to a component
      type TCOMOBJ
        integer :: ICLS   ! index to the list of all class definitions from CLSDEF (see classes.xml)
        integer :: CCLS   ! class ID from CCLS_NAMES (see classes.f)
        integer :: INST   ! instance index
        character(LEN_ID) :: ID ! id name for simple access
      end type TCOMOBJ

      character(1024) :: CCLS_NAMES


      contains

!--------------------------------------------------------
      subroutine CLASSES_INIT
! register component names, order must agree with the definition of the CCLS_XXX
! in classes.f
!--------------------------------------------------------
        CCLS_NAMES='FRAME'
        CCLS_NAMES=trim(CCLS_NAMES)//':'//'SOURCE'
        CCLS_NAMES=trim(CCLS_NAMES)//':'//'DETECTOR'
        CCLS_NAMES=trim(CCLS_NAMES)//':'//'GUIDE'
        CCLS_NAMES=trim(CCLS_NAMES)//':'//'CRYSTAL'
        CCLS_NAMES=trim(CCLS_NAMES)//':'//'XTAL'
        CCLS_NAMES=trim(CCLS_NAMES)//':'//'DCHOPPER'
        CCLS_NAMES=trim(CCLS_NAMES)//':'//'MONITOR'
        CCLS_NAMES=trim(CCLS_NAMES)//':'//'SGUIDE'
        CCLS_NAMES=trim(CCLS_NAMES)//':'//'IGROUP'
      end subroutine CLASSES_INIT


!--------------------------------------------------------
      integer function PCOM_ADD(CCLS,INST)
! add a new component
! CCLS ... class ID number
! INST is the instance index of previously created component
! New component of the same tye is created and data copied from INST
!--------------------------------------------------------
      integer,intent(in) :: CCLS,INST
      integer :: ICOMP
      ICOMP=0
      select case(CCLS)
      case(CCLS_FRAME)
        ICOMP=AddFRAME(INST)
      case(CCLS_DETECTOR)
        ICOMP=AddDETECTOR(INST)
      case(CCLS_SOURCE)
        ICOMP=AddSOURCE(INST)
      case(CCLS_GUIDE)
        ICOMP=AddGUIDE(INST)
      case(CCLS_SGUIDE)
        ICOMP=AddSGUIDE(INST)
      case(CCLS_CRYSTAL)
        ICOMP=AddCRYSTAL(INST)
      case(CCLS_XTAL)
        ICOMP=AddXTAL(INST)
      case(CCLS_DCHOPPER)
        ICOMP=AddDCHOPPER(INST)
      case(CCLS_MONITOR)
        ICOMP=AddMONITOR(INST)
      end select
      PCOM_ADD=ICOMP
      end function PCOM_ADD


!--------------------------------------------------------
! nullifies all pointers in TCOMOBJ
!--------------------------------------------------------
      subroutine pcom_clear(OBJ)
      type(TCOMOBJ) :: OBJ
        OBJ%ICLS=0
        OBJ%CCLS=0
        OBJ%INST=0
        OBJ%ID='undefined'
      end subroutine pcom_clear


!-------------------------------------------------------------------------
! True if CID is class ID of a valid component
!-------------------------------------------------------------------------
      LOGICAL function PCOM_IS_COMPONENT(CID)
      character(*) :: CID
      logical :: LOG1
      LOG1=.false.
      LOG1=(LOG1.or.(CID.eq.trim('FRAME')))
      LOG1=(LOG1.or.(CID.eq.trim('SOURCE')))
      LOG1=(LOG1.or.(CID.eq.trim('DETECTOR')))
      LOG1=(LOG1.or.(CID.eq.trim('GUIDE')))
      LOG1=(LOG1.or.(CID.eq.trim('SGUIDE')))
      LOG1=(LOG1.or.(CID.eq.trim('CRYSTAL')))
      LOG1=(LOG1.or.(CID.eq.trim('XTAL')))
      LOG1=(LOG1.or.(CID.eq.trim('DCHOPPER')))
      LOG1=(LOG1.or.(CID.eq.trim('MONITOR')))
      PCOM_IS_COMPONENT=LOG1
      end function PCOM_IS_COMPONENT

!-------------------------------------------------------------------------
! Creator for COMOBJ, return instance number
! Instance can be obtained by XTAL_INSTANCE(IDX)
!-------------------------------------------------------------------------
      type(TCOMOBJ) function PCOM_CREATE(CID,ID,NAMESTR)
      character(*) :: CID,ID,NAMESTR
      type(TCOMOBJ) :: RES
    !  FRAME:SOURCE:DETECTOR:GUIDE:CRYSTAL:XTAL
      call pcom_clear(RES)
      select case(trim(CID))
      case('FRAME')
        RES%INST=FRAME_CREATE(ID,NAMESTR)
        RES%CCLS=CCLS_FRAME
      case('SOURCE')
        RES%INST=SOURCE_CREATE(ID,NAMESTR)
        RES%CCLS=CCLS_SOURCE
      case('DETECTOR')
        RES%INST=DETECTOR_CREATE(ID,NAMESTR)
        RES%CCLS=CCLS_DETECTOR
      case('GUIDE')
        RES%INST=GUIDE_CREATE(ID,NAMESTR)
        RES%CCLS=CCLS_GUIDE
      case('SGUIDE')
        RES%INST=SGUIDE_CREATE(ID,NAMESTR)
        RES%CCLS=CCLS_SGUIDE
      case('CRYSTAL')
        RES%INST=CRYSTAL_CREATE(ID,NAMESTR)
        RES%CCLS=CCLS_CRYSTAL
      case('XTAL')
        RES%INST=XTAL_CREATE(ID,NAMESTR)
        RES%CCLS=CCLS_XTAL
      case('DCHOPPER')
        RES%INST=DCHOPPER_CREATE(ID,NAMESTR)
        RES%CCLS=CCLS_DCHOPPER
      case('MONITOR')
        RES%INST=MONITOR_CREATE(ID,NAMESTR)
        RES%CCLS=CCLS_MONITOR
      case('IGROUP')
        RES%INST=XTAL_CREATE(ID,NAMESTR)
        RES%CCLS=CCLS_IGROUP
      end select
      if (RES%CCLS.gt.0) then
        RES%ICLS=getICLS(RES%CCLS,CCLS_NAMES)
        RES%ID=trim(ID)
      endif
    !  write(*,*) 'PCOM_CREATE '//trim(RES%ID),' INST=',RES%INST
      PCOM_CREATE=RES
      end function PCOM_CREATE

!-------------------------------------------------------------------------
      SUBROUTINE PCOM_DISPOSE(OBJ)
!-------------------------------------------------------------------------
      type(TCOMOBJ) :: OBJ
    !  write(*,*) 'PCOM_DISPOSE '//trim(OBJ%ID),' INST=',OBJ%INST
      select case(OBJ%CCLS)
        case(CCLS_FRAME)
          call FRAME_DISPOSE(OBJ%INST)
        case(CCLS_SOURCE)
          call SOURCE_DISPOSE(OBJ%INST)
        case(CCLS_DETECTOR)
          call DETECTOR_DISPOSE(OBJ%INST)
        case(CCLS_GUIDE)
          call GUIDE_DISPOSE(OBJ%INST)
        case(CCLS_SGUIDE)
          call SGUIDE_DISPOSE(OBJ%INST)
        case(CCLS_CRYSTAL)
          call CRYSTAL_DISPOSE(OBJ%INST)
        case(CCLS_XTAL)
          call XTAL_DISPOSE(OBJ%INST)
        case(CCLS_DCHOPPER)
          call DCHOPPER_DISPOSE(OBJ%INST)
        case(CCLS_MONITOR)
          call MONITOR_DISPOSE(OBJ%INST)
      end select
      ! write(*,*) 'PCOM_DISPOSE OK 1'
      call pcom_clear(OBJ)
      ! write(*,*) 'PCOM_DISPOSE OK 2'
      end SUBROUTINE PCOM_DISPOSE

!-------------------------------------------------------------------------
      subroutine PCOM_DISPOSE_ALL
!-------------------------------------------------------------------------
        call FRAME_DISPOSE_ALL
        call SOURCE_DISPOSE_ALL
        call DETECTOR_DISPOSE_ALL
        call GUIDE_DISPOSE_ALL
        call SGUIDE_DISPOSE_ALL
        call CRYSTAL_DISPOSE_ALL
        call XTAL_DISPOSE_ALL
        call DCHOPPER_DISPOSE_ALL
        call MONITOR_DISPOSE_ALL
      end subroutine PCOM_DISPOSE_ALL

!-------------------------------------------------------------------------
! Call GO on given component
!-------------------------------------------------------------------------
      logical function PCOM_GO(OBJ)
      type(TCOMOBJ) :: OBJ
      logical :: RES
      ! TYPE(PFRAME) :: F
      RES=.false.
      select case(OBJ%CCLS)
        case(CCLS_FRAME)
          RES=FRAME_GO(OBJ%INST)
        case(CCLS_SOURCE)
          RES=SOURCE_GO(OBJ%INST)
        case(CCLS_DETECTOR)
          RES=DETECTOR_GO(OBJ%INST)
        case(CCLS_GUIDE)
          RES=GUIDE_GO(OBJ%INST)
        case(CCLS_SGUIDE)
          RES=SGUIDE_GO(OBJ%INST)
        case(CCLS_CRYSTAL)
          RES=CRYSTAL_GO(OBJ%INST)
        case(CCLS_XTAL)
          RES=XTAL_GO(OBJ%INST)
        case(CCLS_DCHOPPER)
          RES=DCHOPPER_GO(OBJ%INST)
          !if (RES) then
            !call PCOMOBJ_GET(OBJ,F)
            !ID=ADCHOPPERS(OBJ%INST)%X%FRAME%ID
            !write(*,*) 'PCOMOBJ_GO  ',trim(ID),' CNT=',ADCHOPPERS(OBJ%INST)%X%FRAME%COUNT,' INST=',OBJ%INST
         ! endif
        case(CCLS_MONITOR)
          RES=MONITOR_GO(OBJ%INST)
        end select
      PCOM_GO=RES
      end function PCOM_GO

!---------------------------------------------------------
       SUBROUTINE PCOMOBJ_GET(CLS,F)
! get PFRAME data from a component
!---------------------------------------------------------
      TYPE(TCOMOBJ),intent(in) :: CLS
      TYPE(PFRAME),intent(out) :: F
        F%IDX=0
        NULLIFY(F%X)
        select case(CLS%CCLS)
        case(CCLS_FRAME)
          call FRAME_GET(CLS%INST,F)
        case(CCLS_SOURCE)
          call SOURCE_GET_FRAME(CLS%INST,F)
        case(CCLS_DETECTOR)
          call DETECTOR_GET_FRAME(CLS%INST,F)
        case(CCLS_GUIDE)
          call GUIDE_GET_FRAME(CLS%INST,F)
        case(CCLS_SGUIDE)
          call SGUIDE_GET_FRAME(CLS%INST,F)
        case(CCLS_CRYSTAL)
          call CRYSTAL_GET_FRAME(CLS%INST,F)
        case(CCLS_XTAL)
          call XTAL_GET_FRAME(CLS%INST,F)
        case(CCLS_DCHOPPER)
          call DCHOPPER_GET_FRAME(CLS%INST,F)
        case(CCLS_MONITOR)
          call MONITOR_GET_FRAME(CLS%INST,F)
        end select
      end SUBROUTINE PCOMOBJ_GET

!---------------------------------------------------------
       SUBROUTINE PCOMOBJ_DEFAULT(CLS)
! dispatch calls to *_DEFAULT procedures for each object
!---------------------------------------------------------
      TYPE(TCOMOBJ) :: CLS
        select case(CLS%CCLS)
        case(CCLS_FRAME)
          call FRAME_DEFAULT(CLS%INST)
        case(CCLS_SOURCE)
          call SOURCE_DEFAULT(CLS%INST)
        case(CCLS_DETECTOR)
          call DETECTOR_DEFAULT(CLS%INST)
        case(CCLS_GUIDE)
          call GUIDE_DEFAULT(CLS%INST)
        case(CCLS_SGUIDE)
          call SGUIDE_DEFAULT(CLS%INST)
        case(CCLS_CRYSTAL)
          call CRYSTAL_DEFAULT(CLS%INST)
        case(CCLS_XTAL)
          call XTAL_DEFAULT(CLS%INST)
        case(CCLS_DCHOPPER)
          call DCHOPPER_DEFAULT(CLS%INST)
        case(CCLS_MONITOR)
          call MONITOR_DEFAULT(CLS%INST)
        end select
      end SUBROUTINE PCOMOBJ_DEFAULT

!---------------------------------------------------------
       SUBROUTINE PCOMOBJ_INIT(CLS)
! dispatch calls to *_INIT procedures for each object
!---------------------------------------------------------
      TYPE(TCOMOBJ) :: CLS
        select case(CLS%CCLS)
        case(CCLS_FRAME)
          call FRAME_INIT(CLS%INST)
        case(CCLS_SOURCE)
          call SOURCE_INIT(CLS%INST)
        case(CCLS_DETECTOR)
          call DETECTOR_INIT(CLS%INST)
        case(CCLS_GUIDE)
          call GUIDE_INIT(CLS%INST)
        case(CCLS_SGUIDE)
          call SGUIDE_INIT(CLS%INST)
        case(CCLS_CRYSTAL)
          call CRYSTAL_INIT(CLS%INST)
        case(CCLS_XTAL)
          call XTAL_INIT(CLS%INST)
        case(CCLS_DCHOPPER)
          call DCHOPPER_INIT(CLS%INST)
        case(CCLS_MONITOR)
          call MONITOR_INIT(CLS%INST)
        end select
      end SUBROUTINE PCOMOBJ_INIT

!---------------------------------------------------------
       SUBROUTINE PCOMOBJ_PREPARE(CLS,IERR)
! dispatch calls to *_PREPARE procedures for each object
!---------------------------------------------------------
      TYPE(TCOMOBJ) :: CLS
      integer,intent(out) :: IERR
        IERR=0
        select case(CLS%CCLS)
        case(CCLS_FRAME)
          call FRAME_PREPARE(CLS%INST,IERR)
        case(CCLS_SOURCE)
          call SOURCE_PREPARE(CLS%INST,IERR)
        case(CCLS_DETECTOR)
          call DETECTOR_PREPARE(CLS%INST,IERR)
        case(CCLS_GUIDE)
          call GUIDE_PREPARE(CLS%INST,IERR)
        case(CCLS_SGUIDE)
          call SGUIDE_PREPARE(CLS%INST,IERR)
        case(CCLS_CRYSTAL)
          call CRYSTAL_PREPARE(CLS%INST,IERR)
        case(CCLS_XTAL)
          call XTAL_PREPARE(CLS%INST,IERR)
        case(CCLS_DCHOPPER)
          call DCHOPPER_PREPARE(CLS%INST,IERR)
        case(CCLS_MONITOR)
          call MONITOR_PREPARE(CLS%INST,IERR)
        end select
      end SUBROUTINE PCOMOBJ_PREPARE


!---------------------------------------------------------
       SUBROUTINE PCOMOBJ_RESET(CLS)
! clear event counter of a component
!---------------------------------------------------------
      TYPE(TCOMOBJ),intent(in) :: CLS
      TYPE(PFRAME) :: F
      call PCOMOBJ_GET(CLS,F)
      if (associated(F%X)) then
        !write(*,*) 'PCOMOBJ_RESET ',trim(F%X%ID),' CNT=',F%X%COUNT,' INST=',CLS%INST
        !write(*,*) '               REGISTRY=',F%X%REGISTRY
        F%X%COUNT=0
      endif
      ! clear monitors
      select case(CLS%CCLS)
        case(CCLS_MONITOR)
          call MONITOR_CLR(CLS%INST)
        case(CCLS_DETECTOR)
          call DETECTOR_CLR(CLS%INST)
        case(CCLS_SGUIDE)
          call SGUIDE_CLR(CLS%INST)
        case(CCLS_GUIDE)
          call GUIDE_CLR(CLS%INST)
      end select
      end SUBROUTINE PCOMOBJ_RESET

!---------------------------------------------------------
       SUBROUTINE PCOMOBJ_FINALIZE(CLS)
! dispatch calls to *_FINALIZE procedures for each object
!---------------------------------------------------------
      TYPE(TCOMOBJ) :: CLS
        select case(CLS%CCLS)
        case(CCLS_XTAL)
          call XTAL_FINALIZE(CLS%INST)
        end select
      end SUBROUTINE PCOMOBJ_FINALIZE

!---------------------------------------------------------
       SUBROUTINE PCOMOBJ_ENDTRACE(CLS,NEU)
! dispatch calls to *_ENDTRACE procedures for each object
! specify what components should do after succsessful tracing
! of a single event
!---------------------------------------------------------
      TYPE(TCOMOBJ) :: CLS
      TYPE(NEUTRON) :: NEU

        select case(CLS%CCLS)
        case(CCLS_MONITOR)
          call MONITOR_ENDTRACE(CLS%INST,NEU%CNT,NEU%P)
        case(CCLS_DETECTOR)
          call DETECTOR_ENDTRACE(CLS%INST,NEU%CNT,NEU%P)
        case(CCLS_SGUIDE)
          call SGUIDE_ENDTRACE(CLS%INST,NEU%CNT,NEU%P)
        case(CCLS_GUIDE)
          call GUIDE_ENDTRACE(CLS%INST,NEU%CNT,NEU%P)
        case(CCLS_XTAL)
          call XTAL_ENDTRACE(CLS%INST,NEU%CNT,NEU%P)
        end select
      end SUBROUTINE PCOMOBJ_ENDTRACE

!---------------------------------------------------------
       SUBROUTINE PCOMOBJ_ADJUST(CLS,IERR)
! dispatch calls to *_ADJUST procedures for each object
! IERR<>0 if there is a problem with positioning
!         (e.g. unmatched Bragg condition)
!---------------------------------------------------------
      TYPE(TCOMOBJ) :: CLS
      integer,intent(out) :: IERR
        IERR=0
        select case(CLS%CCLS)
        case(CCLS_FRAME)
          call FRAME_ADJUST(CLS%INST)
        case(CCLS_SOURCE)
          call SOURCE_ADJUST(CLS%INST,IERR)
        case(CCLS_DETECTOR)
          call DETECTOR_ADJUST(CLS%INST,IERR)
        case(CCLS_GUIDE)
          call GUIDE_ADJUST(CLS%INST,IERR)
        case(CCLS_SGUIDE)
          call SGUIDE_ADJUST(CLS%INST,IERR)
        case(CCLS_CRYSTAL)
          call CRYSTAL_ADJUST(CLS%INST,IERR)
        case(CCLS_XTAL)
          call XTAL_ADJUST(CLS%INST,IERR)
        case(CCLS_DCHOPPER)
          call DCHOPPER_ADJUST(CLS%INST,IERR)
        case(CCLS_MONITOR)
          call MONITOR_ADJUST(CLS%INST,IERR)
        end select
      end SUBROUTINE PCOMOBJ_ADJUST


!---------------------------------------------------------
       SUBROUTINE PCOMOBJ_UPDATE_LAB(CLS,RLAB,TLAB)
! Update transformation matrix for laboratory system for given object.
! Must be called in sequence so that RLAB and TLAB contain values corresponding
! to the preceding component.
! Updates also matrices owned by the components:
! OBJ%TLAB: translation to incident coord. in laboratory system
! OBJ%RLAB: rotation from incident coord. to previous component exit coord.
!---------------------------------------------------------
      TYPE(TCOMOBJ) :: CLS
      REAL(KIND(1.D0)),intent(inout) :: RLAB(3,3),TLAB(3)
      TYPE(PFRAME) :: F
        select case(CLS%CCLS)
        case(CCLS_GUIDE)
          call GUIDE_UPDATE_LAB(CLS%INST,RLAB,TLAB)
        case(CCLS_SGUIDE)
          call SGUIDE_UPDATE_LAB(CLS%INST,RLAB,TLAB)
        case DEFAULT
          call PCOMOBJ_GET(CLS,F)
          call FRAME_UPDATE_LAB(F%X,RLAB,TLAB)
        end select
      end SUBROUTINE PCOMOBJ_UPDATE_LAB


!---------------------------------------------------------
       SUBROUTINE PCOMOBJ_INP(CLS,PNAME,ARG,NARG)
! dispatch calls to *_INP procedures for each object
!---------------------------------------------------------
      TYPE(TCOMOBJ) :: CLS
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
        select case(CLS%CCLS)
        case(CCLS_FRAME)
          call FRAME_INP(CLS%INST,PNAME,ARG,NARG)
        case(CCLS_SOURCE)
          call SOURCE_INP(CLS%INST,PNAME,ARG,NARG)
        case(CCLS_DETECTOR)
          call DETECTOR_INP(CLS%INST,PNAME,ARG,NARG)
        case(CCLS_GUIDE)
          call GUIDE_INP(CLS%INST,PNAME,ARG,NARG)
        case(CCLS_SGUIDE)
          call SGUIDE_INP(CLS%INST,PNAME,ARG,NARG)
        case(CCLS_CRYSTAL)
          call CRYSTAL_INP(CLS%INST,PNAME,ARG,NARG)
        case(CCLS_XTAL)
          call XTAL_INP(CLS%INST,PNAME,ARG,NARG)
        case(CCLS_DCHOPPER)
          call DCHOPPER_INP(CLS%INST,PNAME,ARG,NARG)
        case(CCLS_MONITOR)
          call MONITOR_INP(CLS%INST,PNAME,ARG,NARG)
        end select
      end SUBROUTINE PCOMOBJ_INP

!---------------------------------------------------------
       SUBROUTINE PCOMOBJ_OUT(CLS,PNAME,ARG,NARG)
! dispatch calls to *_OUT procedures for each object
!---------------------------------------------------------
      TYPE(TCOMOBJ) :: CLS
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
        select case(CLS%CCLS)
        case(CCLS_FRAME)
          call FRAME_OUT(CLS%INST,PNAME,ARG,NARG)
        case(CCLS_SOURCE)
          call SOURCE_OUT(CLS%INST,PNAME,ARG,NARG)
        case(CCLS_DETECTOR)
          call DETECTOR_OUT(CLS%INST,PNAME,ARG,NARG)
        case(CCLS_GUIDE)
          call GUIDE_OUT(CLS%INST,PNAME,ARG,NARG)
        case(CCLS_SGUIDE)
          call SGUIDE_OUT(CLS%INST,PNAME,ARG,NARG)
        case(CCLS_CRYSTAL)
          call CRYSTAL_OUT(CLS%INST,PNAME,ARG,NARG)
        case(CCLS_XTAL)
          call XTAL_OUT(CLS%INST,PNAME,ARG,NARG)
        case(CCLS_DCHOPPER)
          call DCHOPPER_OUT(CLS%INST,PNAME,ARG,NARG)
        case(CCLS_MONITOR)
          call MONITOR_OUT(CLS%INST,PNAME,ARG,NARG)
        end select
      end SUBROUTINE PCOMOBJ_OUT

      end MODULE VCTABLE_COMPONENTS
