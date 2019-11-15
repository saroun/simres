!//////////////////////////////////////////////////////////////////////
!////  $Id: reportsopt.f,v 1.9 2019/08/15 17:24:08 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.9 $
!////     $Date: 2019/08/15 17:24:08 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Ray-tracing reporting options
!////
!//////////////////////////////////////////////////////////////////////
      MODULE REPORTSOPT
      use FIELDDEF
      use FIELDDATA
      use CLASSES
      implicit none

! structure with reporting option data
      type TREPOPT;sequence
        INTEGER :: CLASS                  ! class ID
        character(LEN_ID)   :: ID         ! component ID string
        CHARACTER(LEN_NAME) :: NAME       ! descriptive name
        logical :: ELOG       ! switch event logs
        integer :: NRAYS      ! number of rays to be logged
        logical :: STATS      ! reporting of tracing statistics
        logical :: VOL        ! reporting of sampling volume
        logical :: VARI       ! reporting of variance reduction progress
        logical :: PROG       ! reporting of tracing progress
        logical :: RES        ! reporting basic results (intensity, E-spread)
        logical :: MAXV       ! report max. value optimization
        LOGICAL :: MUTE       ! switch all reports on/off
        REAL(KIND(1.D0)) :: TIME       !  progress bar update time in [s]
      end type TREPOPT

      type(TREPOPT) :: REPOPT

      contains

C--------------------------------------------------------
      SUBROUTINE REPOPT_DEFAULT(OBJ)
C Set default settings for tracing module
C These can be ovewritten by user before the call of TRACING_INI
C--------------------------------------------------------
      TYPE(TREPOPT) :: OBJ
      OBJ%CLASS=OCLS_REPORTS
      OBJ%ELOG=.FALSE.
      OBJ%NRAYS=40
      OBJ%STATS=.TRUE.
      OBJ%VOL=.TRUE.
      OBJ%VARI=.TRUE.
      OBJ%PROG=.TRUE.
      OBJ%RES=.TRUE.
      OBJ%MAXV=.TRUE.
      OBJ%TIME=1  ! time in [s]
      OBJ%MUTE=.FALSE.
      END SUBROUTINE REPOPT_DEFAULT

C--------------------------------------------------------
      SUBROUTINE REPOPT_MUTE(OBJ)
C Set settings for mute behaviour
C--------------------------------------------------------
      TYPE(TREPOPT) :: OBJ
      OBJ%CLASS=OCLS_REPORTS
      OBJ%ELOG=.FALSE.
      OBJ%NRAYS=40
      OBJ%STATS=.FALSE.
      OBJ%VOL=.FALSE.
      OBJ%VARI=.FALSE.
      ! OBJ%PROG=.TRUE.
      OBJ%RES=.FALSE.
      OBJ%MAXV=.FALSE.
      OBJ%TIME=1  ! time in [s]
      OBJ%MUTE=.TRUE.
      END SUBROUTINE REPOPT_MUTE

!---------------------------------------------------------
      SUBROUTINE REPOPT_INP(OBJ,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      type(TREPOPT) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      LR=0
      NARG=0
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call FIELD2ARRAY(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
          call REPOPT_INP_R(OBJ,trim(ARG%ID),NUM,LR)
        case default
          write(*,*) 'REPOPT_INP: undefined field type '//trim(ARG%DEF%TID)
      end select
      NARG=LR
      END SUBROUTINE REPOPT_INP


!---------------------------------------------------------
      SUBROUTINE REPOPT_OUT(OBJ,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      type(TREPOPT) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      LR=0
      NARG=0
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call REPOPT_OUT_R(OBJ,trim(ARG%ID),NUM,LR)
          call ARRAY2FIELD(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
        case default
          write(*,*) 'REPOPT_OUT: undefined field type '//trim(ARG%DEF%TID)
      end select
      NARG=LR
      END SUBROUTINE REPOPT_OUT

C---------------------------------------------------------
      SUBROUTINE REPOPT_INP_R(OBJ,PNAME,ARG,NARG)
C input data from ARG to OBJ for TREPOPT class
C INPUT
C    PNAME   .... parameter name
C    ARG     .... REAL*8 parameter values
C INPUT
C    OBJ     .... TREPOPT structure
C    NARG    .... number of items read from ARG
C---------------------------------------------------------
      TYPE(TREPOPT),intent(out) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(in) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR
      LR=1
      SELECT CASE (trim(PNAME))
      CASE('ELOG'); OBJ%ELOG=(NINT(ARG(1)).eq.1)
      CASE('NRAYS'); OBJ%NRAYS=NINT(ARG(1))
      CASE('STATS'); OBJ%STATS=(NINT(ARG(1)).eq.1)
      CASE('VOL'); OBJ%VOL=(NINT(ARG(1)).eq.1)
      CASE('VARI'); OBJ%VARI=(NINT(ARG(1)).eq.1)
      CASE('PROG');   OBJ%PROG=(NINT(ARG(1)).eq.1)
      CASE('RES'); OBJ%RES=(NINT(ARG(1)).eq.1)
      CASE('MAXV'); OBJ%MAXV=(NINT(ARG(1)).eq.1)
      CASE('TIME'); OBJ%TIME=ARG(1)
      END SELECT
      NARG=LR
      END SUBROUTINE REPOPT_INP_R

C---------------------------------------------------------
      SUBROUTINE REPOPT_OUT_R(OBJ,PNAME,ARG,NARG)
C output data from OBJ to ARG for TTROPT class
C INPUT
C    PNAME   .... parameter name
C    ARG     .... REAL*8 parameter values
C INPUT
C    OBJ     .... TREPOPT structure
C    NARG    .... number of items read from ARG
C---------------------------------------------------------
      TYPE(TREPOPT),intent(in) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(out) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR
      LR=1
      SELECT CASE (trim(PNAME))
      CASE('ELOG'); ARG(1)=0; if (OBJ%ELOG) ARG(1)=1.D0
      CASE('NRAYS'); ARG(1)=OBJ%NRAYS
      CASE('STATS'); ARG(1)=0; if (OBJ%STATS) ARG(1)=1.D0
      CASE('VOL'); ARG(1)=0; if (OBJ%VOL) ARG(1)=1.D0
      CASE('VARI'); ARG(1)=0; if (OBJ%VARI) ARG(1)=1.D0
      CASE('PROG');   ARG(1)=0; if (OBJ%PROG) ARG(1)=1.D0
      CASE('RES'); ARG(1)=0; if (OBJ%RES) ARG(1)=1.D0
      CASE('MAXV'); ARG(1)=0; if (OBJ%MAXV) ARG(1)=1.D0
      CASE('TIME'); ARG(1)=OBJ%TIME
      END SELECT
      NARG=LR
      END SUBROUTINE REPOPT_OUT_R



      end MODULE REPORTSOPT
