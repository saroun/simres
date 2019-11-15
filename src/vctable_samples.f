!//////////////////////////////////////////////////////////////////////
!////  $Id: vctable_samples.f,v 1.43 2019/08/15 15:02:06 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.43 $
!////     $Date: 2019/08/15 15:02:06 $
!//////////////////////////////////////////////////////////////////////
!////
!////   Virtual classes table for SAMPLES
!////  - serves to emulation of polymorphism in Fortran
!//////////////////////////////////////////////////////////////////////
      MODULE VCTABLE_SAMPLES
      use CLASSDEF,ONLY:GetICLS
      USE SAMPLE_SINGLE
      USE SAMPLE_POLY
      use SAMPLES_TRACE
      use XMLINFO
      use FRAMES_TRACE
      use NSTORE
      implicit none

! TSAMOBJ - abstract type with pointer to a component
      type TSAMOBJ
        integer :: ICLS  ! class ID from CLS_IDS (all defined in classes.xml)
        integer :: SCLS  ! class ID from SCLS_NAMES
        integer :: ISTORE(2) ! reference to the event storage for initial and final neutron states
        type(TFRAME),pointer :: P_FRAME ! parent class frame
        type(TSAMPLE),pointer :: P_SAMPLE ! parent class
        type(TPCRYST),pointer :: P_POLY
        type(TSCRYST),pointer :: P_SINGLE
      end type TSAMOBJ

! overloaded POBJ
      interface PSAMOBJ
        module procedure psamobj_sample_
        module procedure psamobj_poly_
        module procedure psamobj_single_
      end interface

! temporary storage for neutron states
      type(NEUTRON) :: TMPN1,TMPN2

! global instance of TSAMOBJ
      type(TSAMOBJ) :: SAMOBJ

      private psamobj_sample_
      private psamobj_poly_
      private psamobj_single_
      private TMPN1,TMPN2
      private SAMPLE_GO

      contains

!--------------------------------------------------------
! nullifies all pointers in TSAMOBJ
!--------------------------------------------------------
      subroutine PSAMOBJ_CLEAR(OBJ)
      type(TSAMOBJ) :: OBJ
        OBJ%ICLS=0
        OBJ%SCLS=0
        NULLIFY(OBJ%P_FRAME)
        NULLIFY(OBJ%P_SAMPLE)
        NULLIFY(OBJ%P_POLY)
        NULLIFY(OBJ%P_SINGLE)
      end subroutine PSAMOBJ_CLEAR

!-------------------------------------------
! overloaded PSAMOBJ: define for each object type
!-------------------------------------------
      type(TSAMOBJ) function psamobj_sample_(OBJ)
        type(TSAMPLE),target :: OBJ
        type(TSAMOBJ) :: RES
        call PSAMOBJ_CLEAR(RES)
        RES%SCLS=SCLS_SAMPLE
        RES%P_SAMPLE => OBJ
        RES%ICLS=getICLS(RES%SCLS,SCLS_NAMES)
        RES%P_FRAME => OBJ%FRAME
        psamobj_sample_=RES
      end function psamobj_sample_
      type(TSAMOBJ) function psamobj_poly_(OBJ)
        type(TPCRYST),target :: OBJ
        type(TSAMOBJ) :: RES
        call PSAMOBJ_CLEAR(RES)
        RES%SCLS=SCLS_POLY
        RES%P_POLY => OBJ
        RES%ICLS=getICLS(RES%SCLS,SCLS_NAMES)
        RES%P_FRAME => OBJ%SAM%FRAME
        RES%P_SAMPLE => OBJ%SAM
        psamobj_poly_=RES
      end function psamobj_poly_
      type(TSAMOBJ) function psamobj_single_(OBJ)
        type(TSCRYST),target :: OBJ
        type(TSAMOBJ) :: RES
        call PSAMOBJ_CLEAR(RES)
        RES%SCLS=SCLS_SINGLE
        RES%P_SINGLE => OBJ
        RES%ICLS=getICLS(RES%SCLS,SCLS_NAMES)
        RES%P_FRAME => OBJ%SAM%FRAME
        RES%P_SAMPLE => OBJ%SAM
        psamobj_single_=RES
      end function psamobj_single_


!  ******             DISPATCH PROCEDURES FOR VIRTUAL METHODS           ******

!---------------------------------------------------------
       SUBROUTINE PSAMOBJ_ASSIGN(OBJ,SRC)
! copy object data from SRC to OBJ
!---------------------------------------------------------
      TYPE(TSAMOBJ),intent(out) :: OBJ
      TYPE(TSAMOBJ),intent(in) :: SRC
      if (associated(SRC%P_SAMPLE)) then
        call psamobj_clear(OBJ)
        select case(SRC%SCLS)
        case(SCLS_SAMPLE)
          OBJ%P_SAMPLE=SRC%P_SAMPLE
        case(SCLS_POLY)
          OBJ%P_POLY=SRC%P_POLY
          OBJ%P_SAMPLE => OBJ%P_POLY%SAM
        case(SCLS_SINGLE)
          OBJ%P_SINGLE=SRC%P_SINGLE
          OBJ%P_SAMPLE => OBJ%P_SINGLE%SAM
        case default
          return
        end select
        OBJ%ICLS=SRC%ICLS
        OBJ%SCLS=SRC%SCLS
        OBJ%P_FRAME => OBJ%P_SAMPLE%FRAME
      endif
      end SUBROUTINE PSAMOBJ_ASSIGN

!---------------------------------------------------------
       SUBROUTINE PSAMOBJ_INP(CLS,PNAME,ARG,NARG)
! dispatch calls to *_INP procedures for each object
!---------------------------------------------------------
      TYPE(TSAMOBJ) :: CLS
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
        select case(CLS%SCLS)
        case(SCLS_SAMPLE)
          call SAMPLE_INP(CLS%P_SAMPLE,PNAME,ARG,NARG)
        case(SCLS_POLY)
          call PCRYST_INP(CLS%P_POLY,PNAME,ARG,NARG)
        case(SCLS_SINGLE)
          call SCRYST_INP(CLS%P_SINGLE,PNAME,ARG,NARG)
        end select
      end SUBROUTINE PSAMOBJ_INP

!---------------------------------------------------------
       SUBROUTINE PSAMOBJ_OUT(CLS,PNAME,ARG,NARG)
! dispatch calls to *_OUT procedures for each object
!---------------------------------------------------------
      TYPE(TSAMOBJ) :: CLS
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
        select case(CLS%SCLS)
        case(SCLS_SAMPLE)
          call SAMPLE_OUT(CLS%P_SAMPLE,PNAME,ARG,NARG)
        case(SCLS_POLY)
          call PCRYST_OUT(CLS%P_POLY,PNAME,ARG,NARG)
        case(SCLS_SINGLE)
          call SCRYST_OUT(CLS%P_SINGLE,PNAME,ARG,NARG)
        end select
      end SUBROUTINE PSAMOBJ_OUT

!---------------------------------------------------------
       SUBROUTINE PSAMOBJ_DEFAULT(CLS)
! dispatch calls to *_DEFAULT procedures for each object
!---------------------------------------------------------
      TYPE(TSAMOBJ) :: CLS
        select case(CLS%SCLS)
        case(SCLS_SAMPLE)
          call SAMPLE_DEFAULT(CLS%P_SAMPLE)
        case(SCLS_POLY)
          call PCRYST_DEFAULT(CLS%P_POLY)
        case(SCLS_SINGLE)
          call SCRYST_DEFAULT(CLS%P_SINGLE)
        end select
      end SUBROUTINE PSAMOBJ_DEFAULT


!---------------------------------------------------------
       SUBROUTINE PSAMOBJ_PREPARE(CLS,IERR)
! dispatch calls to *_PREPARE procedures for each object
!---------------------------------------------------------
      TYPE(TSAMOBJ) :: CLS
      integer,intent(out) :: IERR
        select case(CLS%SCLS)
        case(SCLS_SAMPLE)
          call SAMPLE_PREPARE(CLS%P_SAMPLE,IERR)
        case(SCLS_POLY)
          call PCRYST_PREPARE(CLS%P_POLY,IERR)
        case(SCLS_SINGLE)
          call SCRYST_PREPARE(CLS%P_SINGLE,IERR)
        end select
      end SUBROUTINE PSAMOBJ_PREPARE


!---------------------------------------------------------
       SUBROUTINE PSAMOBJ_GETQE(CLS,Q,E)
! return |Q| for given sample
!---------------------------------------------------------
      TYPE(TSAMOBJ) :: CLS
      real(kind(1.D0)),intent(out) :: Q,E
        select case(CLS%SCLS)
        case(SCLS_SAMPLE)
          call SAMPLE_GETQE(CLS%P_SAMPLE,Q,E)
        case(SCLS_POLY)
          call PCRYST_GETQE(CLS%P_POLY,Q,E)
        case(SCLS_SINGLE)
          call SCRYST_GETQE(CLS%P_SINGLE,Q,E)
        end select
      end SUBROUTINE PSAMOBJ_GETQE

!---------------------------------------------------------
       SUBROUTINE PSAMOBJ_DEPENDENCES(CLS,IERR)
! update sample data for given wavelength
!---------------------------------------------------------
      TYPE(TSAMOBJ) :: CLS
      integer,intent(out) :: IERR
        select case(CLS%SCLS)
        case(SCLS_POLY)
          if (CLS%P_SAMPLE%TYP.eq.isam_POWDER) then
            call PCRYST_DEPENDENCES(CLS%P_POLY,IERR)
          ENDIF
        case(SCLS_SINGLE)
          call SCRYST_DEPENDENCES(CLS%P_SINGLE,IERR)
        end select
      end SUBROUTINE PSAMOBJ_DEPENDENCES

!-----------------------------------------------------------------
      SUBROUTINE PSAMOBJ_INIT(CLS)
! Adjust sample - dispatch procedure
! Arguments are stored in ARG(:) array
! ensure also correct sample type, compatible with CLS type
!-----------------------------------------------------------------------
      TYPE(TSAMOBJ) :: CLS
      if (associated(CLS%P_SAMPLE)) then
      ! check sample type
        select case(CLS%SCLS)
        case(SCLS_SAMPLE)
          call SAMPLE_INIT(CLS%P_SAMPLE)
          select case(CLS%P_SAMPLE%TYP)
          case(isam_INELASTIC,isam_ELASTIC,isam_VANAD)
            ! all is OK
          case default
            write(*,*) 'PSAMOBJ_INIT Incorrect type: ',CLS%P_SAMPLE%TYP
            write(*,*) '             allowed = ',isam_INELASTIC,isam_ELASTIC,isam_VANAD
            CLS%P_SAMPLE%TYP=isam_INELASTIC
            call MSG_WARN('Incorrect type for generic sample: set to inelastic',1)
          end select
        case(SCLS_POLY)
          call PCRYST_INIT(CLS%P_POLY)
          select case(CLS%P_SAMPLE%TYP)
          case(isam_ELASTIC,isam_POWDER)
            ! all is OK
          case default
            write(*,*) 'PSAMOBJ_INIT Incorrect type: ',CLS%P_SAMPLE%TYP
            write(*,*) '             allowed = ',isam_ELASTIC,isam_POWDER
            CLS%P_SAMPLE%TYP=isam_ELASTIC
            call MSG_WARN('Incorrect type for polycrystalline sample: set to elastic',1)
          end select
        case(SCLS_SINGLE)
          call SCRYST_INIT(CLS%P_SINGLE)
          select case(CLS%P_SAMPLE%TYP)
          case(isam_ELASTIC,isam_INELASTIC,isam_PHONON,isam_VANAD)
            ! all is OK
          case default
            write(*,*) 'PSAMOBJ_INIT Incorrect type: ',CLS%P_SAMPLE%TYP
            write(*,*) '             allowed = ',isam_ELASTIC,isam_INELASTIC,isam_PHONON,isam_VANAD
            CLS%P_SAMPLE%TYP=isam_INELASTIC
            call MSG_WARN('Incorrect type for single crystal sample: set to inelastic',1)
          end select
        end select
      ! set coordinate system type
        select case(CLS%P_SAMPLE%TYP)
        case(isam_INELASTIC,isam_ELASTIC,isam_VANAD,isam_PHONON)
          CLS%P_SAMPLE%COORD=csam_SANS
        case(isam_POWDER)
          CLS%P_SAMPLE%COORD=csam_PWD
        end select
      endif
      end SUBROUTINE PSAMOBJ_INIT

!-----------------------------------------------------------------
      SUBROUTINE PSAMOBJ_FINALIZE(CLS)
! Finalization subroutine for samples
! called after instrument tracing
!-----------------------------------------------------------------------
      TYPE(TSAMOBJ) :: CLS
    ! if (associated(CLS%P_SAMPLE)) then
    !  endif
      end SUBROUTINE PSAMOBJ_FINALIZE

!---------------------------------------------------------
       SUBROUTINE PSAMOBJ_ENDTRACE(ibeam,NEU)
! Save initial and final neutron states
!---------------------------------------------------------
      integer,intent(in) :: ibeam
      TYPE(NEUTRON) :: NEU
      if (associated(SAMOBJ%P_FRAME)) then
        select case(ibeam)
        case (1)
          if (SAMOBJ%ISTORE(1)>0) then
            ! transform to incident axis coordinates
            call FRAME_LOCAL(-1,SAMOBJ%P_FRAME,TMPN1)
            TMPN1%P=NEU%P     ! set the final weight to the event
            TMPN1%CNT=NEU%CNT ! must set counter to neutron before calling NSTORE_SETEVENT !!
            TMPN1%REF=NEU%REF ! set reference to the initial neutron !!
            call NSTORE_SETEVENT(SAMOBJ%ISTORE(1),TMPN1)
          endif
        case(2)
          if (SAMOBJ%ISTORE(2)>0) then
            ! transform to incident axis coordinates
            call FRAME_LOCAL(-1,SAMOBJ%P_FRAME,TMPN2)
            TMPN2%P=NEU%P     ! set the final weight to the event
            TMPN2%CNT=NEU%CNT ! must set counter to neutron before calling NSTORE_SETEVENT !!
            TMPN2%REF=NEU%REF ! set reference to the initial neutron !!
            call NSTORE_SETEVENT(SAMOBJ%ISTORE(2),TMPN2)
          endif
        end select
      endif
      end SUBROUTINE PSAMOBJ_ENDTRACE

!--------------------------------------------------------------------
      LOGICAL FUNCTION PSAMOBJ_GO()
!--------------------------------------------------------------------
! handle non-associated pointers
      if ((.not.associated(SAMOBJ%P_SAMPLE)).or.(.not.associated(SAMOBJ%P_FRAME))) then
        PSAMOBJ_GO=.false.
      else
        PSAMOBJ_GO=SAMPLE_GO(SAMOBJ)
      endif
      END FUNCTION PSAMOBJ_GO


!--------------------------------------------------------------------
      LOGICAL FUNCTION PSAMOBJ_CONTAINER_GO()
!--------------------------------------------------------------------
! handle non-associated pointers
      logical :: LOG1
      LOG1=.false.
      if ((.not.associated(SAMOBJ%P_SAMPLE)).or.(.not.associated(SAMOBJ%P_FRAME))) then
        LOG1=.false.
      else
        LOG1=CONTAINER_GO(SAMOBJ%P_FRAME)
        if (LOG1) then
          ! save initial state
          TMPN1=NEUT
        endif
      endif
      PSAMOBJ_CONTAINER_GO=LOG1
      END FUNCTION PSAMOBJ_CONTAINER_GO

!--------------------------------------------------------------------
      LOGICAL FUNCTION SAMPLE_GO(CLS)
! Top-level procedure for tracing through sample
! Works in LOCAL coordinates at the begining.
! It must follow previous call to CONTAINER_GO, which moves neutron to
!    a random scattering point inside the sample container and leaves
!    NEUT on the local coordinates
! On exit, convert to the exit coordinates.
!--------------------------------------------------------------------
  ! TOF to the scattering point, sampled from exp(-x/L)
  !      DT=T1-LOG(1.D0-P1*RAN1())*LPATH
      TYPE(TSAMOBJ) :: CLS
      LOGICAL :: LOG1
      REAL(kind(1.D0)) :: T1,T2,PATH0,PATH1,PATH2,P0,PP,SIGMA
      REAL(kind(1.D0)) :: R(3),K(3),t(4),dt(2)
      REAL(kind(1.D0)), parameter :: EPS=1.D-8
      integer :: it
      logical :: dbg=.false.
1     format(a,': ',6(G10.4,1x),I6)
! handle non-associated pointers
      if ((.not.associated(CLS%P_SAMPLE)).or.(.not.associated(CLS%P_FRAME))) then
        SAMPLE_GO=.false.
        return
      endif
! save initial state
      TMPN1=NEUT
      PP=1.D0
! get cross-section with container borders (also for hollow shapes)


      !dbg=((NEUT%R(3)>0).and.(CLS%P_FRAME%COUNT<30))

      it = CONTAINER_BORDER(CLS%P_FRAME,NEUT%R,NEUT%K,t,dt)
      LOG1 = (it>0)
! get path in the material before and after the current position along incident beam
      T1 = dt(1)
      T2 = dt(2)
      ! LOG1=(FRAME_BORDER(CLS%P_FRAME,NEUT%R,NEUT%K,T1,T2))
      PATH0=NEUT%K0*(T2+T1) ! path along incident beam through the whole sample
      LOG1=(LOG1.and.(PATH0>1.D-4)) ! neglect thickness < 0.1 um
  ! reject condition
  ! allow continuation for transmission samples
      if (.not.LOG1) then
        LOG1=(CLS%P_SAMPLE%TRANS.and.(.not.CLS%P_SAMPLE%SCATT))
        goto 100
      endif
      PATH1=NEUT%K0*T1
      PATH2=NEUT%K0*T2
      P0=NEUT%P
      SIGMA=0.D0
  ! convert to Q-coordinates
      R=NEUT%R
      K=NEUT%K
      CALL M3XV3(1,CLS%P_SAMPLE%QMAT,R,NEUT%R)
      CALL M3XV3(1,CLS%P_SAMPLE%QMAT,K,NEUT%K)
      LOG1=.false.
      SELECT case (CLS%P_SAMPLE%TYP)
        case(isam_POWDER)
          if (associated(CLS%P_POLY)) LOG1=PCRYST_GO_INST(CLS%P_POLY,PATH0,SIGMA)
!       case(isam_PHONON)
!            if (associated(CLS%P_SINGLE)) LOG1=SCRYST_GO(CLS%P_SINGLE,NEUT)
!       case(isam_VANAD)
!            LOG1=SAM_VANAD(CLS%P_SAMPLE,NEUT)
        case default
          LOG1=SAMPLE_GO_INST(CLS%P_SAMPLE,PATH0,SIGMA)
      end SELECT
! don't loose time with absorbed events
      if (.not.LOG1) goto 100
! handle transmitted neutrons if allowed
      !if (CLS%P_SAMPLE%TRANS) then
      !  if (SIGMA.le.0.D-10) then
      !    NEUT%P=NEUT%P/PATH0
      !    goto 100
      !  endif
      !  Z=SIGMA*PATH0
      !  NEUT%P=NEUT%P/(1.D0-exp(-Z))
      !endif

  ! convert from Q-coordinates back to local
      R=NEUT%R
      K=NEUT%K
      CALL M3XV3(-1,CLS%P_SAMPLE%QMAT,R,NEUT%R)
      CALL M3XV3(-1,CLS%P_SAMPLE%QMAT,K,NEUT%K)
    ! weight by beam attenuation
      if (SIGMA>0.D0) then
        ! get cross-section with container borders for scattered ray
        it = CONTAINER_BORDER(CLS%P_FRAME,NEUT%R,NEUT%K,t,dt)
        if (it>0) then
          PATH2=NEUT%K0*dt(2)
        endif
        !IF (FRAME_BORDER(CLS%P_FRAME,NEUT%R,NEUT%K,T1,T2)) THEN
        !  PATH2=NEUT%K0*T2
        !endif
        PP=PP*exp(-SIGMA*(PATH1+PATH2))
        LOG1=(PP.GT.1.D-9)
      endif
      call AddEventLog(CLS%P_FRAME,NEUT)
! save final state
      TMPN2=NEUT

! convert to exit coordinates
      call SLIT_POST(CLS%P_FRAME)

100   if (LOG1) then
        NEUT%P=NEUT%P*PP
        call incCounter(CLS%P_FRAME%COUNT)
      else
        NEUT%P=0.D0
      endif
      SAMPLE_GO=LOG1
      END FUNCTION SAMPLE_GO



      end MODULE VCTABLE_SAMPLES

