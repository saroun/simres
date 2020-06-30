!//////////////////////////////////////////////////////////////////////
!////  $Id: components.f,v 1.108 2019/08/15 15:02:08 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2011, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.108 $
!////     $Date: 2019/08/15 15:02:08 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Definition of TAS components - former structures
!////  Provides abstraction layer for components
!////
!//////////////////////////////////////////////////////////////////////
      MODULE COMPONENTS
      use CLASSES
      use FRAMES
      USE DETECTORS
      USE CRYSTALS
      USE SAMPLES
      USE GUIDES
      USE SGUIDES
      use XTALS_TRACE
      use VCTABLE_COMPONENTS
      USE VCTABLE_SAMPLES
      USE VCTABLE_OPTIONS
      USE VCTABLE_INSTRUMENTS
      use VCTABLE
      use IGROUPS
      use NSTORE
      use MESSAGES
      use MIRRLOG
      use GENERATOR
      implicit none

! TCOMPID - component reference in BEAMLINE arrays
      TYPE  TCOMPID
        TYPE(TCOMOBJ) :: OBJ  ! object reference
        integer :: IBEAM      ! beamline number: (1) primary (2) secondary
        integer :: IORD       ! component order index in the BEAMLINE array
        integer :: ISTORE     ! reference to events storage (= index on the list of registered storages in NSTORE)
        integer :: IREG       ! registered as monitor IREG (= index on the list of monitors, BMONITORS
      end type TCOMPID

! This is the array describing actual beamline
! 2nd index stands for PRIMARY or SECONDARY part
      integer :: BEAMLINE_NC(2) ! actual number of components in the beamline
      TYPE(TCOMPID) :: BEAMLINE(BEAMLINE_DIM,2)

! References to the components acting as monochromator and analyzers
      integer,parameter :: MAX_MONO=16
      integer :: MONO_N,ANAL_N
      TYPE(TCOMOBJ) :: MONOCHROMATORS(MAX_MONO)
      TYPE(TCOMOBJ) :: ANALYZERS(MAX_MONO)

      TYPE(TCOMOBJ) :: BEAMSOURCE
      character(4*(LEN_ID+1)) :: MONOCHROMATORSID,ANALYZERSID
! transformation arrays for laboratory coord. system
      real(KIND(1.D0)), private  :: RLAB(3,3),TLAB(3)

      contains


!--------------------------------------------------------
      subroutine REGISTER_MONO
! scan IDs of monochromators or analyzers and register them
!--------------------------------------------------------
      integer :: i,IC
      MONO_N=0
      ANAL_N=0
    ! monochromators
      do i=1,BEAMLINE_NC(1)
        call GETPARINDX(MONOCHROMATORSID,':',trim(BEAMLINE(i,1)%OBJ%ID),IC)
        if ((IC.gt.0).and.(MONO_N.lt.MAX_MONO)) then
          MONO_N=MONO_N+1
          MONOCHROMATORS(MONO_N)=BEAMLINE(i,1)%OBJ
        endif
      enddo
    ! analyzers
      do i=1,BEAMLINE_NC(2)
        call GETPARINDX(ANALYZERSID,':',trim(BEAMLINE(i,2)%OBJ%ID),IC)
        if ((IC.gt.0).and.(ANAL_N.lt.MAX_MONO)) then
          ANAL_N=ANAL_N+1
          ANALYZERS(ANAL_N)=BEAMLINE(i,2)%OBJ
        endif
      enddo
      end subroutine REGISTER_MONO


!--------------------------------------------------------
      subroutine REGISTER_LOGGERS
! register bounce loggers
!--------------------------------------------------------
      integer :: i
      TYPE(TCOMOBJ) :: OBJ  ! object reference
      call MLOGS_FREEALL
      do i=1,BEAMLINE_NC(1)
        OBJ=BEAMLINE(i,1)%OBJ
        select case(OBJ%CCLS)
          case(CCLS_GUIDE)
            call GUIDE_REGMLOG(OBJ%INST)
          case(CCLS_SGUIDE)
            call SGUIDE_REGMLOG(OBJ%INST)
        end select
      enddo
      end subroutine REGISTER_LOGGERS

!-------------------------------------------------------------------------
! Get component localization on the beamline
! INPUT:
!   IDSTR  .. ID string
! OUTPUT:
!   COMPID   .. record with component identification
!-------------------------------------------------------------------------
      subroutine GetComponentByID(IDSTR,COMPID)
      character(*) ,intent(in) :: IDSTR
      type(TCOMPID), intent(out) :: COMPID
      logical :: found
      integer :: i,j
      COMPID%IBEAM=0
      COMPID%IORD=0
      call PCOM_CLEAR(COMPID%OBJ)
      if (len_trim(IDSTR).le.0) return
      found=.false.
      i=1
      do j=1,2
        i=1
        do while ((.not.found).and.(i.le.BEAMLINE_NC(j)))
          if (BEAMLINE(i,j)%OBJ%ID.eq.trim(IDSTR)) then
            found=.true.
            COMPID%IBEAM=j
            COMPID%IORD=i
            COMPID%OBJ=BEAMLINE(i,j)%OBJ
          endif
          i=i+1
        enddo
      enddo
      end subroutine GetComponentByID

!-------------------------------------------------------------
! get component data reference from beamline index
!-------------------------------------------------------------
      type(TCOMOBJ) function getComByIndex(INDX,IBEAM)
        integer,intent(in) :: INDX,IBEAM
        type(TCOMOBJ) :: RES
        call PCOM_CLEAR(RES)
        if (INDX.gt.0) then
          select case (IBEAM)
          case(1,2)
            if (INDX.le.BEAMLINE_NC(IBEAM)) RES=BEAMLINE(INDX,IBEAM)%OBJ
          end select
        endif
        GetComByIndex=RES
      end function GetComByIndex

C-------------------------------------------------------------
! get component pointer by ID string
C-------------------------------------------------------------
      type(TCOMOBJ) function getComByID(IDSTR)
        character(*) :: IDSTR
        type(TCOMPID) :: COMPID
        call GetComponentByID(IDSTR,COMPID)
        getComByID=getComByIndex(COMPID%IORD,COMPID%IBEAM)
      end function getComByID

!-------------------------------------------------------------
! get TCLASS by ID string
!-------------------------------------------------------------
      type(TCLASS) function getClassByID(IDSTR)
        character(*) :: IDSTR
        type(TCLASS) :: RES
        character(LEN_ID) :: CLSID,OBJID
        character(LEN_NAME) :: OBJNAME
        call PCLASS_CLEAR(RES)
! try options
        if (RES%ICLS.le.0) then
          RES=PCLASS(getOptByID(IDSTR))
        endif
! try interface
        if ((RES%ICLS.le.0).and.(INSOBJ%ICLS.gt.0)) then
          RES=PCLASS(INSOBJ)
          call GetClassAttrib(RES,CLSID,OBJID,OBJNAME)
          if (trim(IDSTR).ne.trim(OBJID)) call PCLASS_CLEAR(RES)
        endif
! try sample
        if ((RES%ICLS.le.0).and.(SAMOBJ%ICLS.gt.0)) then
          RES=PCLASS(SAMOBJ)
          call GetClassAttrib(RES,CLSID,OBJID,OBJNAME)
          if (trim(IDSTR).ne.trim(OBJID)) call PCLASS_CLEAR(RES)
        endif
! try a beamline object
        if (RES%ICLS.le.0) then
          RES=PCLASS(getComByID(IDSTR))
        endif
        getClassByID=RES
      end function getClassByID

!-------------------------------------------
! delete all compoents from beamlines and deallocate associated menmory
!-------------------------------------------
      subroutine DeleteComponents
      integer :: i,j
        do i=1,MAX_MONO
          call PCOM_CLEAR(MONOCHROMATORS(i))
          call PCOM_CLEAR(ANALYZERS(i))
        enddo
        MONO_N=0
        ANAL_N=0
        MONOCHROMATORSID=' '
        ANALYZERSID=' '
        call PSAMOBJ_CLEAR(SAMOBJ)
        call PINSOBJ_CLEAR(INSOBJ)
        RLAB=0.D0
        do i=1,3
          RLAB(i,i)=1.D0
          TLAB(i)=0.D0
        enddo
        call XTALS_REF_DISPOSE
        call PCOM_DISPOSE_ALL
        call IGROUP_DISPOSE_ALL
        DO j=1,2
          BEAMLINE_NC(j)=0
          do i=1,BEAMLINE_DIM
            call PCOM_CLEAR(BEAMLINE(i,j)%OBJ)
            BEAMLINE(i,j)%IORD=0
            BEAMLINE(i,j)%IBEAM=0
          enddo
        enddo
      end subroutine DeleteComponents

!-------------------------------------------------------------
! Call _PREPARE procedure on all components
! It should be called after each change of components parameters
! so that the components can update their internal fields
! like transformation matrices etc.
!--------------------------------------------------------------
      subroutine PrepareBeamline
      integer :: IERR
      integer :: i,j,im
      TYPE(TCOMOBJ) :: OBJ
      character(128) :: S
      IERR=0
      do j=1,2
        do i=1,BEAMLINE_NC(j)
          OBJ=BEAMLINE(i,j)%OBJ
          call PCOMOBJ_PREPARE(OBJ,IERR)
          if (IERR.ne.0) then
            S='Can''t prepare internal fields for the component '//trim(OBJ%ID)
            im= ADD_MESSAGE(m_error,'PCOMOBJ_PREPARE',trim(S),m_low)
          endif
        enddo
      enddo
      if (SAMOBJ%ICLS.gt.0) then
        call PSAMOBJ_PREPARE(SAMOBJ,IERR)
        if (IERR.ne.0) then
          S='Can''t prepare internal fields for the sample'
          im= ADD_MESSAGE(m_error,'PSAMOBJ_PREPARE',trim(S),m_low)
        endif
      endif
      !write(*,*) 'PrepareBeamline done.'

      end subroutine PrepareBeamline


!-------------------------------------------------------------
! Get coordinates of a neutron in global coordinates for logging
! Assumes NEUT is in the state after adjust procedure, i.e.
! in the exit axis coordinates of a component.
! Input:
! OBJ: component reference as TCOMOBJ
!-------------------------------------------------------------
      subroutine GetNCoord(F, NGLOB, DIST)
      TYPE(TFRAME), intent(in) :: F
      TYPE(NEUTRON), intent(out) :: NGLOB
      real(KIND(1.D0)) :: DIST
      TYPE(NEUTRON) :: NEU
      NEU = NEUT
      call FRAME_AXIS(1,F,NEU)
      call Ax2Glob(F,NEU,NGLOB)
      DIST= F%DIST
      end subroutine GetNCoord


!-------------------------------------------------------------
! Adjust position and orientation of all components except samples
! Adjust transformation matrices for laboratory system for the whole beamline
! Lab. coordinates = incident axis coordinates of the 1st component
! if IO>0, write component table to this unit:
!   ID, x, y, z, distance, tof
! all in global coordinates, mm, us
!--------------------------------------------------------------
      subroutine AdjustBeamline(IERR, IO)
      integer,intent(out) :: IERR
      integer :: i,IO
      logical :: tmp_g
      real(KIND(1.D0)) :: ddet,dsam,dt0, aux, d
      TYPE(TCOMOBJ) :: OBJ
      TYPE(PSOURCE) :: PSRC
      TYPE(NEUTRON) :: NGLOB
      TYPE(PFRAME) :: F
      character(32) :: CNUM,CNUM1
      real(KIND(1.D0)) :: RL(3,3), TL(3)

      !real(kind(1.D0)) :: K0
1     format(a,': ',7(G13.6,1x))
11    format(a,1x,2(1x,G13.6),3(1x,G15.7))
! initialize RLAB, TLAB arrays for conversion to LAB coord
      !call getRotYXY(PSRC%X%FRAME%REXI,RL)
      call UNIMAT(RLAB,3,3)
      TLAB=0.D0
      ierr=0
      if (IO>0) then
        write(IO,11) '# ID, x, y, z, distance, tof'
      endif
      TRACING_UP=.false.
! don't use gravity for adjustment
      tmp_g=use_gravity
      use_gravity=.false.
      call SOURCE_GET(BEAMSOURCE%INST,PSRC)
      PULSED_SRC=(PSRC%X%TYP.eq.src_pulsed)
      SRC_OVERLAP=PSRC%X%OVERLAP
      SRC_PERIOD=HOVM/PSRC%X%FRQ
! set reference times to zero before tracing along the beam axis
! TOF_ZERO is adjusted during following scan by the T0 chopper, if it exists
      TOF_ZERO=0.D0
      TOF_ZERO2=0.D0
! generate event at the origin of the beam axis
      call GENERATE_NOM(INSOBJ%P_SPEC%KI,NEUT)
      !write(*,1) 'AdjustBeamline R,K',NEUT%R, NEUT%K
      do i=1,BEAMLINE_NC(1)
        OBJ=BEAMLINE(i,1)%OBJ
        aux = NEUT%K0*NEUT%T
        call PCOMOBJ_ADJUST(OBJ,IERR)
        if (ierr.ne.0) goto 99
        call PCOMOBJ_GET(OBJ,F)
        !write(*,1) 'adjusted '//trim(F%X%ID),NEUT%R, NEUT%K
        !write(*,1) 'adjust '//trim(F%X%ID),F%X%TEXI, F%X%REXI(1,3), NEUT%T
        !write(*,1) 'move by :  ',aux, NEUT%K0*NEUT%T,NEUT%K0*NEUT%T-aux
        call PCOMOBJ_UPDATE_LAB(OBJ,RLAB,TLAB)

        if (IO>0) then
          call PCOMOBJ_GET(OBJ,F)
          call GetNCoord(F%X, NGLOB, d)
          !write(IO,11) OBJ%ID, NGLOB%R, d, NGLOB%T/HOVM
          CALL M3XV3(1,PSRC%X%FRAME%REXI,F%X%TLAB,TL)
          write(IO,11) OBJ%ID, TL, d, NGLOB%T*NEUT%K0
        endif

      enddo
      ! combinded with McStas: add extra path before
      IF (TROPT%MCSTAS.and.(TROPT%PATHIN.gt.0)) then
        !aux = NEUT%K0*NEUT%T
        NEUT%T = NEUT%T + TROPT%PATHIN/NEUT%K0
        !write(*,1) 'AdjustBeamline, added before [mm]:  ',TROPT%PATHIN
        !write(*,1) 'move by :  ',aux, NEUT%K0*NEUT%T
      endif
 ! set TOF_SAMPLE at the end of primary beamline if there is no sample
      TOF_SAMPLE=NEUT%T
      NEUI=NEUT
      if (SAMOBJ%ICLS.gt.0) then
      ! adjust sample, adjust NEUI relative to the secondary axis using sample data
        call SAMPLE_ADJUST(SAMOBJ%P_SAMPLE,IERR)
        TOF_SAMPLE=NEUT%T
        call FRAME_AXIS(1,SAMOBJ%P_FRAME,NEUI)
        call FRAME_UPDATE_LAB(SAMOBJ%P_FRAME,RLAB,TLAB)
        if (IO>0) then
          call GetNCoord(SAMOBJ%P_FRAME, NGLOB, d)
          !write(IO,11) SAMOBJ%P_FRAME%ID, NGLOB%R, d, NGLOB%T/HOVM
          CALL M3XV3(1,PSRC%X%FRAME%REXI,SAMOBJ%P_FRAME%TLAB,TL)
          write(IO,11) SAMOBJ%P_FRAME%ID, TL, d, NGLOB%T*NEUT%K0
        endif
      else
      ! no sample, only adjust NEUI relative to the secondary axis using instrument data
        !write(*,1) 'SPEC_ADJUST from :  ',NEUT%K0*NEUT%T, INSOBJ%P_SPEC%FRAME%DIST
        call SPEC_ADJUST(INSOBJ%P_SPEC,IERR)
        !write(*,1) 'SPEC_ADJUST to 1:  ',NEUT%K0*NEUT%T
        TOF_SAMPLE=NEUT%T
        call FRAME_AXIS(1,INSOBJ%P_SPEC%FRAME,NEUI)
        !write(*,1) 'SPEC_ADJUST to 2:  ',NEUT%K0*NEUT%T
        call FRAME_UPDATE_LAB(INSOBJ%P_SPEC%FRAME,RLAB,TLAB)
        !write(*,1) 'SPEC_ADJUST to 3:  ',NEUT%K0*NEUT%T
      endif
      NEUF=NEUT
      ! combinded with McStas: add extra path after
      IF (TROPT%MCSTAS.and.(TROPT%PATHOUT.gt.0)) then
        !aux = NEUT%K0*NEUT%T
         NEUT%T = NEUT%T + TROPT%PATHOUT/NEUT%K0
        !write(*,1) 'AdjustBeamline, added after [mm]:  ',TROPT%PATHIN
        !write(*,1) 'move by :  ',aux, NEUT%K0*NEUT%T
      endif
      !write(*,1) 'starting secondary at :  ',NEUT%K0*NEUT%T
      do i=1,BEAMLINE_NC(2)
        OBJ=BEAMLINE(i,2)%OBJ
        !aux = NEUT%K0*NEUT%T
        call PCOMOBJ_ADJUST(OBJ,IERR)
        if (ierr.ne.0) goto 99
        call PCOMOBJ_UPDATE_LAB(OBJ,RLAB,TLAB)
        if (IO>0) then
          call PCOMOBJ_GET(OBJ,F)
          call GetNCoord(F%X, NGLOB, d)
          !write(IO,11) OBJ%ID, NGLOB%R, d, NGLOB%T/HOVM
          CALL M3XV3(1,PSRC%X%FRAME%REXI,F%X%TLAB,TL)
          write(IO,11) OBJ%ID, TL, d, NGLOB%T*NEUT%K0
        endif
        !write(*,1) 'move by :  ',aux, NEUT%K0*NEUT%T
      enddo
      TOF_DETECTOR=NEUT%T
      call AdjustChoppers
      dt0=NEUI%K0*TOF_ZERO
      dsam=NEUI%K0*TOF_SAMPLE
      ddet=dsam+NEUF%K0*(TOF_DETECTOR-TOF_SAMPLE)
      if (IO>0) then
        if (EVENTGEN%KPROF.ne.0) then
          call FLOAT2STR(TWOPI/EVENTGEN%KMAX,CNUM)
          call FLOAT2STR(TWOPI/EVENTGEN%KMIN,CNUM1)
          write(IO,11) '# Wavelength range set to '//trim(CNUM)//' .. '//trim(CNUM1)
        else
          call FLOAT2STR(TWOPI/EVENTGEN%K0,CNUM)
          write(IO,11) '# Wavelength set to:  '//trim(CNUM)
        endif
        write(IO,11) '# TOF  T0, sample, detector [us] ',TOF_ZERO/HOVM,TOF_SAMPLE/HOVM,TOF_DETECTOR/HOVM
        write(IO,11) '# DIST T0, sample, detector [mm] ',dt0,dsam,ddet
      endif
      !write(*,*) 'AdjustBeamline done.'

99    if (ierr.ne.0) then
        call FLOAT2STR(NEUT%K0,CNUM)
        call MSG_ERROR('AdjustBeamline','can''t adjust '//trim(OBJ%ID)//' for |K|='//trim(CNUM),0,1)
      endif
      use_gravity=tmp_g
      end subroutine AdjustBeamline

!-----------------------------------------------------------------------
! Adjust chopper phases
! For choppers with AUTOADJ=true,
! set phase=0 for the nominal wavelength and moderator time = 0
! lock phase to TOF_ZERO if required by lockT0 flag
!-----------------------------------------------------------------------
      subroutine AdjustChoppers
      integer :: i,j
      TYPE(TCOMOBJ) :: OBJ
      do j=1,2
        do i=1,BEAMLINE_NC(j)
          OBJ=BEAMLINE(i,j)%OBJ
          if (OBJ%CCLS.eq.CCLS_DCHOPPER) then
            call DCHOPPER_SETPHASE(OBJ%INST)
          endif
        enddo
      enddo
      end subroutine AdjustChoppers


C/// Following subroutines add individual components to the beamline IBEAM
C/// There is no polymorphism in Fotran, so we have to work around
C/// The objects in arguments have already be prepared, i.e. with actual
C/// position, orientation, size etc.


!-------------------------------------------------------------------------
! Add a component: TSAMOBJ
! SAMPLE is not part of a beamline
!-------------------------------------------------------------------------
      INTEGER function AddSpecimen(OBJ)
      type(TSAMOBJ),intent(in) :: OBJ
      integer :: i
        i=0
        if (associated(OBJ%P_SAMPLE)) then
          i=1
          select case(OBJ%SCLS)
          case(SCLS_SAMPLE)
            SAMPLE=OBJ%P_SAMPLE
            SAMOBJ=PSAMOBJ(SAMPLE)
          case(SCLS_POLY)
            PCRYST=OBJ%P_POLY
            SAMOBJ=PSAMOBJ(PCRYST)
          case(SCLS_SINGLE)
            SCRYST=OBJ%P_SINGLE
            SAMOBJ=PSAMOBJ(SCRYST)
          case default
            i=0
          end select
        endif
        if (i.gt.0) call PSAMOBJ_INIT(SAMOBJ)
        AddSpecimen=i
      end function AddSpecimen

!-------------------------------------------------------------------------
! Add a component: TINSOBJ
!-------------------------------------------------------------------------
      INTEGER function AddInterface(OBJ)
      type(TINSOBJ),intent(in) :: OBJ
      integer :: i
        i=0
        if (associated(OBJ%P_SPEC)) then
          i=1
          select case(OBJ%DCLS)
          case(DCLS_SPEC)
            if (associated(OBJ%P_SPEC)) then
              SPEC=OBJ%P_SPEC
              INSOBJ=PINSOBJ(SPEC)
              i=1
            endif
      !  case(DCLS_PWDIF)
      !  case(DCLS_TAS)
          case default
            i=0
          end select
        endif
        if (i.gt.0) call PINSOBJ_INIT(INSOBJ)
        AddInterface=i
      end function AddInterface


!-------------------------------------------------------------------------
! Add a component: TOPTOBJ
!-------------------------------------------------------------------------
      INTEGER function AddOptions(OBJ)
      type(TOPTOBJ),intent(in) :: OBJ
      integer :: i
        i=0
        select case(OBJ%OCLS)
        case(OCLS_TRACING)
          TROPT=OBJ%P_TRA
          TRAOBJ=POPTOBJ(TROPT)
        case(OCLS_REPORTS)
          REPOPT=OBJ%P_REP
          REPOBJ=POPTOBJ(REPOPT)
        case default
          i=0
        end select
        AddOptions=i
  !      write(*,*) 'AddOptions ',associated(OBJ%P_TRA),OPTOBJ%ICLS,OPTOBJ%OCLS,associated(OPTOBJ%P_TRA)
      end function AddOptions

!-------------------------------------------------------------------------
! Add TCOMOBJ to the BEAMLINE array
! INPUT:
!   OBJ    .. TCOMOBJ data
!   IBEAM  .. beam index
!-------------------------------------------------------------------------
      INTEGER function AddComponent(OBJ,IBEAM)
      type(TCOMOBJ) :: OBJ
      integer,intent(in) :: IBEAM
      integer :: NC,ICOMP
      !TYPE(PFRAME) :: F
      ICOMP=0
      AddComponent=0
      if ((IBEAM.LE.0).or.(IBEAM.GT.2)) return
      if (BEAMLINE_NC(IBEAM).ge.BEAMLINE_DIM) return
      ICOMP=PCOM_ADD(OBJ%CCLS,OBJ%INST)
  ! add OBJ to a beamline
      if (ICOMP.gt.0) then
        NC=BEAMLINE_NC(IBEAM)+1
        BEAMLINE_NC(IBEAM)=NC
        BEAMLINE(NC,IBEAM)%OBJ=OBJ
        BEAMLINE(NC,IBEAM)%IORD=NC
        BEAMLINE(NC,IBEAM)%IBEAM=IBEAM
        BEAMLINE(NC,IBEAM)%OBJ%INST=ICOMP
        BEAMLINE(NC,IBEAM)%ISTORE=0
        BEAMLINE(NC,IBEAM)%IREG=0
        ! assign INST to the objects frame data
        !call PCOMOBJ_GET(OBJ,F)
        !if (associated(F%X)) then
        !  F%X%INST=ICOMP
        !endif
        !write(*,*) 'AddComponent ',NC,IBEAM,trim(OBJ%ID)
        ! register source component
        if (OBJ%CCLS.eq.CCLS_SOURCE) BEAMSOURCE=BEAMLINE(NC,IBEAM)%OBJ
      endif
      AddComponent=ICOMP
      end function AddComponent

!--------------------------------------------------------
      logical function GetDetector(DET)
! get PDETECTOR object, return true if successful:
! Detector must be the last component of the secondary bemline !!!
!--------------------------------------------------------
      type(PDETECTOR) :: DET
      type(TCOMOBJ) :: OBJ
      logical :: LOG1
      OBJ=BEAMLINE(BEAMLINE_NC(2),2)%OBJ
      if (OBJ%CCLS.eq.CCLS_DETECTOR) then
        call DETECTOR_GET(OBJ%INST,DET)
        LOG1=.true.
      else
        LOG1=.false.
      endif
      GetDetector=LOG1
      end function GetDetector

!--------------------------------------------------------
      SUBROUTINE AddIGROUP(INST,IBEAM)
! add group of interacting components to the beamline
!--------------------------------------------------------
      INTEGER :: INST
      integer,intent(in) :: IBEAM
      integer :: NC
      type(PIGROUP) :: PIG
      if ((IBEAM.LE.0).or.(IBEAM.GT.2)) return
      if (BEAMLINE_NC(IBEAM).ge.BEAMLINE_DIM) return
      if (IGROUP_isValid(INST)) then
        call IGROUP_GET(INST,PIG)
        NC=BEAMLINE_NC(IBEAM)+1
        BEAMLINE_NC(IBEAM)=NC
        BEAMLINE(NC,IBEAM)%OBJ%ID=PIG%X%FRAME%ID
        BEAMLINE(NC,IBEAM)%OBJ%ICLS=0
        BEAMLINE(NC,IBEAM)%OBJ%CCLS=CCLS_IGROUP
        BEAMLINE(NC,IBEAM)%OBJ%INST=PIG%IDX
        BEAMLINE(NC,IBEAM)%IORD=NC
        BEAMLINE(NC,IBEAM)%ISTORE=0
        BEAMLINE(NC,IBEAM)%IREG=0
        BEAMLINE(NC,IBEAM)%IBEAM=IBEAM
      endif
      end SUBROUTINE AddIGroup


!--------------------------------------------------------
      SUBROUTINE COMPONENTS_CLR(IBEAM)
! clear event counters of each component on the beamline IBEAM
!--------------------------------------------------------
      integer, intent(in) :: IBEAM
      integer:: i
      select case(IBEAM)
      case(1,2)
        do i=1,BEAMLINE_NC(IBEAM)
          select case(BEAMLINE(i,IBEAM)%OBJ%CCLS)
          case(CCLS_IGROUP)
            call IGROUP_RESET(BEAMLINE(i,IBEAM)%OBJ%INST)
          case default
            call PCOMOBJ_RESET(BEAMLINE(i,IBEAM)%OBJ)
          end select
        enddo
      end select
      end SUBROUTINE COMPONENTS_CLR




!-------------------------------------------------------------------------
! Call GO on given component
!-------------------------------------------------------------------------
      logical function COMPONENT_GO(OBJ)
      type(TCOMOBJ) :: OBJ
      logical :: RES
      RES=.false.
      select case(OBJ%CCLS)
        case(CCLS_IGROUP)
          RES=IGROUP_GO(OBJ%INST)
        case default
          RES=PCOM_GO(OBJ)
        end select
      COMPONENT_GO=RES
      end function COMPONENT_GO

!-------------------------------------------------------------------------
! Call INIT procedue on each component on the beamline IBEAM
! This is called just before tracing starts, after adjustment
!-------------------------------------------------------------------------
      subroutine COMPONENTS_INIT(IBEAM)
      integer, intent(in) :: IBEAM
      integer:: i
      select case(IBEAM)
      case(1,2)
        do i=1,BEAMLINE_NC(IBEAM)
          select case(BEAMLINE(i,IBEAM)%OBJ%CCLS)
          case(CCLS_IGROUP)
            call IGROUP_INIT(BEAMLINE(i,IBEAM)%OBJ%INST)
          casE default
            call PCOMOBJ_INIT(BEAMLINE(i,IBEAM)%OBJ)
          end select
        enddo
      end select
      end SUBROUTINE COMPONENTS_INIT


!-------------------------------------------------------------------------
! Call FINALIZE procedue on each component on the beamline IBEAM
!-------------------------------------------------------------------------
      subroutine COMPONENTS_FINALIZE(IBEAM)
      integer, intent(in) :: IBEAM
      integer:: i
      select case(IBEAM)
      case(1,2)
        do i=1,BEAMLINE_NC(IBEAM)
          select case(BEAMLINE(i,IBEAM)%OBJ%CCLS)
         case(CCLS_IGROUP)
            call IGROUP_FINALIZE(BEAMLINE(i,IBEAM)%OBJ%INST)
          case DEFAULT
            call PCOMOBJ_FINALIZE(BEAMLINE(i,IBEAM)%OBJ)
          end select
        enddo
      end select
      end SUBROUTINE COMPONENTS_FINALIZE


!-------------------------------------------------------------------------
! Call FINALIZE procedue on each component on the beamline IBEAM
!-------------------------------------------------------------------------
      subroutine COMPONENTS_ENDTRACE(IBEAM,NEU)
      type(NEUTRON) :: NEU
      integer, intent(in) :: IBEAM
      integer:: i
      ! write(*,*) 'COMPONENTS_ENDTRACE ',IBEAM,NEU%CNT
      select case(IBEAM)
      case(1,2)
        do i=1,BEAMLINE_NC(IBEAM)
          select case(BEAMLINE(i,IBEAM)%OBJ%CCLS)
          case(CCLS_IGROUP)
            call IGROUP_ENDTRACE(BEAMLINE(i,IBEAM)%OBJ%INST,NEU)
          case DEFAULT
            call PCOMOBJ_ENDTRACE(BEAMLINE(i,IBEAM)%OBJ,NEU)
          end select
        enddo
      end select
      end SUBROUTINE COMPONENTS_ENDTRACE



      end MODULE COMPONENTS


