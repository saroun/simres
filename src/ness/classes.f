!//////////////////////////////////////////////////////////////////////
!////  $Id: classes.f,v 1.39 2012/04/09 22:35:04 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.39 $
!////     $Date: 2012/04/09 22:35:04 $
!//////////////////////////////////////////////////////////////////////
!////
!////   Definitions of all classes and parameters parsed from classes.xml
!////  - class IDs, indexes, parameter IDs, names, units, hints, enum. definitions etc.
!////
!//////////////////////////////////////////////////////////////////////
      MODULE CLASSES
      implicit none
      save

! ELEMENTS DEFINED BY SIMRES:
! only those with named indices

! class IDs used internally in dispatch procedures
      integer,parameter :: TCLS_COM=1  ! component
      integer,parameter :: TCLS_SAM=2  ! sample
      integer,parameter :: TCLS_INS=3  ! instrument interface
      integer,parameter :: TCLS_OPT=4  ! option group
      integer,parameter :: TCLS_CMD=5  ! command
      integer,parameter :: TCLS_IGP=6  ! interacting group of components

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
      integer, PARAMETER :: CCLS_XTAL = 6
      integer, PARAMETER :: CCLS_DCHOPPER = 7
      integer, PARAMETER :: CCLS_MONITOR = 8
      integer, PARAMETER :: CCLS_SGUIDE = 9
      integer, PARAMETER :: CCLS_IGROUP = 10
      integer, parameter :: CCLS_NUMBER = CCLS_IGROUP
  ! class identificator (keep the same order as class ID !)
  !    character(56),parameter :: CCLS_NAMES='FRAME:SOURCE:DETECTOR:GUIDE:CRYSTAL:XTAL:DCHOPPER:IGROUP'

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


! BEAMLINE (keep it here, it is neded to define dimensions by many other modules)
      integer, parameter :: BEAMLINE_DIM=128  ! max. number of components

! CLASSES


  ! partition table for parameter arrays

! local variables useful during parsing
!      integer,private :: CLS_INDEX
!      integer,private :: ENUM_INDEX
!      CHARACTER(LEN_ID),private  :: CLASS_ID
!      CHARACTER(LEN_ID),private  :: PARAM_ID
!      CHARACTER(LEN_NAME),private  :: CLASS_NAME



  !    interface getEnumTypeOrd
  !      module procedure getEnumTypeOrd_1
  !      module procedure getEnumTypeOrd_2
  !    end interface

!      public OCLS_TRACING,OCLS_REPORTS,OCLS_NAMES
!      public SCLS_SAMPLE,SCLS_SINGLE,SCLS_POLY,SCLS_NAMES
!      public CCLS_FRAME,CCLS_GUIDE,CCLS_DETECTOR,CCLS_SOURCE,CCLS_CRYSTAL,CCLS_NAMES,CCLS_XTAL
!      public DCLS_TAS,DCLS_PWDIF,DCLS_SPEC,DCLS_NAMES
!      public TCLS_COM,TCLS_SAM,TCLS_OPT,TCLS_INS,TCLS_CMD,TCLS_NAMES
!      public KFMODE_FLAT,KFMODE_TOF
!      public FRAME_SHAPE_BOX,FRAME_SHAPE_CYLINDER,FRAME_SHAPE_DISC, FRAME_SHAPE_ELLIPSOID
!      public BEAMLINE_DIM





      end MODULE CLASSES