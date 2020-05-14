!//////////////////////////////////////////////////////////////////////
!////  $Id: instcontrol.f,v 1.65 2019/01/11 17:35:58 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.65 $
!////     $Date: 2019/01/11 17:35:58 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Basic interface for instrument control
!////  High-level control of beamline components
!////
!//////////////////////////////////////////////////////////////////////
      MODULE INSTCONTROL
      use DATASETS
      use COMPONENTS
      use COMPONENTS_IO
      use EVENTMONITOR
      use VCTABLE
      USE GENERATOR
      use SOURCES_TABLE
      use MESSAGES
      implicit none

      contains

!--------------------------------------------------------
      subroutine SELECT_DSET(IX)
! select given dataset
! copy parameters from DATASETS module to BPAR_ARRAY
! update beamline by calling BEAM_INP
!--------------------------------------------------------
      integer,intent(in) :: IX
      integer :: IERR
      logical :: needadj
      if (DSET_WRITE(IX,BPAR_ARRAY(1),BPAR_NP)) then
        needadj=(DSET_ISEL.ne.IX)
        DSET_ISEL=IX
        !write(*,*) 'SELECT_DSET before BEAM_INP'
        call BEAM_INP
        !write(*,*) 'SELECT_DSET after BEAM_INP'
        if (needadj) call INST_ADJUST(IERR,'SELECT_DSET')
      endif
      end subroutine SELECT_DSET

!--------------------------------------------------------
      subroutine DEFINE_DSET(IX)
! define given dataset using current setup
! adjust spectrometer
! update BPAR_ARRAY by calling BEAM_OUT
!--------------------------------------------------------
      integer,intent(in) :: IX
      character(32) :: CNUM
      integer :: IERR,IS
      logical :: needadj
      TYPE(PDATASET) :: DTS
      IS=IX
      if (IX.le.0) IS=1
      needadj=(DSET_ISEL.ne.IS)
      !write(*,*) 'DEFINE_DSET DSET_ISEL,IX,IS ',DSET_ISEL,IX,IS
      if (IS.gt.DSET_N) then
      ! destroy previously registered stores
        if (DSET_N.eq.0) call NSTORE_UNREGISTER_ALL
      ! new instrument is defined ...
        call REGISTER_BMONITORS
        !write(*,*) 'DEFINE_DSET REGISTER_MONITORS'
        call INST_ADJUST(IERR,'DEFINE_DSET 1')
        !write(*,*) 'DEFINE_DSET INST_ADJUST'
        call BEAM_OUT
        if (DSET_N.eq.0) call DSET_MKREF(BPAR_ARRAY(1),BPAR_NP)
        call DSET_ADD
        !write(*,*) 'DEFINE_DSET DSET_ADD ',DSET_N
        DSET_ISEL=IS
        DTS=DSET_DATA(DSET_ISEL)
        DTS%P%REG1=BMONITORS(BMONITOR_PRI)%ISTORE
        DTS%P%REG2=BMONITORS(BMONITOR_SEC)%ISTORE
        call INT2STR(DSET_N,CNUM)
        if (XMLOUT.le.0) call MSG_INFO('Created dataset '//trim(CNUM),1)
        needadj=.false.
      ENDIF
      call DSET_READ(IS,BPAR_ARRAY(1),BPAR_NP)
      call DSET_SETPOS(IS,INSOBJ%P_SPEC%NOM)
      !write(*,*) 'DEFINE_DSET DSET_SETPOS'
      if (needadj) then
        call INST_ADJUST(IERR,'DEFINE_DSET 2')
        call BEAM_OUT
      endif
      !write(*,*) 'DEFINE_DSET OK'
      end subroutine DEFINE_DSET

!--------------------------------------------------------
! Apply settings to beamline components
! set LOGMSG and DBGM=true to track the calling subroutine
!--------------------------------------------------------
      subroutine INST_ADJUST(IERR,LOGMSG)
      integer,intent(out) :: IERR
      logical, PARAMETER :: DBGM=.false.
      character*(*) :: LOGMSG
      logical :: hasSample
      REAL(KIND(1.D0)) :: Q0,EN,QAB(3)
      character(128) :: S
      integer :: im
1     format('INST_ADJUST ',a,': ',6(G12.5,1x))
      ierr=1
      if (DBGM) then
        write(*,*) 'INST_ADJUST '//trim(LOGMSG)
      endif
! scan IDs of monochromators or analyzers and register them
      call REGISTER_MONO
! nothing to adjust without a defined instrument interface
      if (.not.associated(INSOBJ%P_SPEC)) return


! consolidate input to make KI,KF,EN,Q,TTHETA consistent if possible
      call SPEC_VALIDATE_POS(INSOBJ%P_SPEC,IERR)

! copy input Q,E values
      Q0=INSOBJ%P_SPEC%Q
      EN=INSOBJ%P_SPEC%EN
      !write(*,1) 'INST_ADJUST Q,E',Q0,EN

! exchange parameters with SAMPLE settings if required
      hasSample=associated(SAMOBJ%P_SAMPLE)
      if (hasSample) then
        TROPT%MCSTAS=.false.
        TROPT%PATHIN=0.D0
        TROPT%PATHOUT=0.D0
      ! first initialize sample QE by SPEC data
        SAMOBJ%P_SAMPLE%Q0=Q0
        SAMOBJ%P_SAMPLE%SGN=INSOBJ%P_SPEC%SS
        SAMOBJ%P_SAMPLE%ESC=EN
      ! prepare sample, e.g. UB matrix
        call PSAMOBJ_PREPARE(SAMOBJ,IERR)
      ! if auto-adjust, get Q,E from the sample
        if (INSOBJ%P_SPEC%ADJSPEC) then
          call PSAMOBJ_GETQE(SAMOBJ,Q0,EN)
        ! keep SPEC.Q=const when adjusting instrument
          !write(*,1) 'INST_ADJUST Q,E',Q0,EN
          INSOBJ%P_SPEC%QCONST=.true.
        endif
      else
        INSOBJ%P_SPEC%ADJSPEC=.false.
      endif

! sample may have changed Q,E => adjust SPEC fields again
      call SPEC_SET_EN(INSOBJ%P_SPEC,EN,IERR)
      if (ierr>0) then
        S='Can''t adjust Ki or Kf for given energy transfer: '//trim(LOGMSG)
        im= ADD_MESSAGE(m_error,'INST_ADJUST',trim(S),m_low)
        return
      endif
      call SPEC_SETQ(INSOBJ%P_SPEC,Q0,ierr)
      if (ierr>0) then
        S='Can''t adjust Ki,Kf for given Q: '//trim(LOGMSG)
        im= ADD_MESSAGE(m_error,'INST_ADJUST',trim(S),m_low)
        return
      endif

! Prepare for instrument adjustment
      call INST_ADJUST_PREPARE
      call PrepareBeamline

! at this point, SPEC should have consistent KI,KF,Q,EN,LAMBDA,TTHETA
! update sample data if required (scattering power, UB matrix, Bragg angle etc ...)
      if (hasSample) then
        SAMOBJ%P_SAMPLE%KI0=INSOBJ%P_SPEC%KI
        SAMOBJ%P_SAMPLE%KF0=INSOBJ%P_SPEC%KF
        SAMOBJ%P_SAMPLE%Q0=INSOBJ%P_SPEC%Q
        SAMOBJ%P_SAMPLE%SGN=INSOBJ%P_SPEC%SS
        SAMOBJ%P_SAMPLE%ESC=INSOBJ%P_SPEC%EN
        call PSAMOBJ_DEPENDENCES(SAMOBJ,ierr)
        if (ierr>0) then
          S='Can''t calculate sample parameters (e.g. scattering angle) for given instrument setting.'
          im= ADD_MESSAGE(m_error,'INST_ADJUST',trim(S),m_low)
          return
        endif
! from now on, SAMPLE should have all its internal fields updated
! for example the UB marrix should be up to date for single crystals
      ! single crystals only: adjust instrument off-plane Q component when required
        if (INSOBJ%P_SPEC%ADJSPEC) then
          if (INSOBJ%P_SPEC%OffPlane.and.(SAMOBJ%SCLS==SCLS_SINGLE)) then
            call SCRYST_GETQAB(SAMOBJ%P_SINGLE,QAB)
            INSOBJ%P_SPEC%PSI=-asin(QAB(3)/INSOBJ%P_SPEC%KF)
            !write(*,1)  'Calculated PSI [deg]:',INSOBJ%P_SPEC%PSI/deg
            INSOBJ%P_SPEC%KFMODE=KFMODE_FLAT
          else
            INSOBJ%P_SPEC%PSI=0.D0
          endif
        endif
      endif
! calculate orientation of the KF axis
! SAMPLE.FRAME.AX or INSOBJ%P_SPEC%FRAME%AX, and corresponding rotation matrix
      call INST_ADJUST_AXIS(ierr)
      if (ierr>0) return
! adjust sample goniometer
      if (hasSample) then
        call INST_ADJUST_SAMPLE(ierr)
        if (ierr>0) return
      endif

! update dataset structure, SPEC%NOM
      INSOBJ%P_SPEC%NOM%KI=INSOBJ%P_SPEC%KI
      if (hasSample) then
      ! has sample: get coordinates from the sample
        call M3xV3(-SAMOBJ%P_FRAME%MEXI,SAMOBJ%P_FRAME%REXI,(/0.D0,0.D0,INSOBJ%P_SPEC%KF/),INSOBJ%P_SPEC%NOM%KF)
        INSOBJ%P_SPEC%NOM%MCL=SAMOBJ%P_SAMPLE%MCL
        INSOBJ%P_SPEC%NOM%MFL=SAMOBJ%P_FRAME%REXI
      else
      ! no sample: get coordinates from the instrument setting
        call M3xV3(-1,INSOBJ%P_SPEC%FRAME%REXI,(/0.D0,0.D0,INSOBJ%P_SPEC%KF/),INSOBJ%P_SPEC%NOM%KF)
        call UNIMAT(INSOBJ%P_SPEC%NOM%MCL,3,3) ! MCL makes no sense without sample
        INSOBJ%P_SPEC%NOM%MFL=INSOBJ%P_SPEC%FRAME%REXI
      endif

! adjust position of other components
      call AdjustBeamline(IERR, 0)
      end subroutine INST_ADJUST

!--------------------------------------------------------
      subroutine INST_ADJUST_SAMPLE(IERR)
! adjust sample goniometer
!--------------------------------------------------------
      integer,intent(out) :: IERR
      ! orient single crystals when required
      if ((SAMOBJ%SCLS==SAMOBJ%SCLS).and.INSOBJ%P_SPEC%ORISAM) then
        call SAMPLE_ORIENT(SAMOBJ%P_SAMPLE,INSOBJ%P_SPEC%OffPlane)
      ! otherwise keep GON values and update transform. matrices only
      else
        call SetGonioAngles(SAMOBJ%P_FRAME,SAMOBJ%P_FRAME%GON)
      endif
      ! now kf-axis and sample gonio are in correct position
      ! update various transformation matrices related to the sample Q-vector
      call SAMPLE_SETQSC(SAMOBJ%P_SAMPLE)
      end subroutine INST_ADJUST_SAMPLE

!--------------------------------------------------------
! Preparatory tasks for instrument adjustment.
! No positioning is done here, only preparatory tasks, namely:
! 1) PRIMARY beamline: set generator constraints
! 2) SECONDARY beamline: set FCONE=true for the 1st analyzer crystal on the beamline in FLATCONE mode
! Requires KI - call only after KI is finally defined !!
!--------------------------------------------------------
      subroutine INST_ADJUST_PREPARE
      TYPE(PSOURCE) :: PSRC
      TYPE(PCRYSTAL) :: CR
      TYPE(PXTAL) :: XT
      REAL(KIND(1.D0)) :: K0,kmin,kmax
1     format(a,': ',6(G12.5,1x))
      !REAL(KIND(1.D0)) :: dtdef,t0,tmin,tmax
      character(32) :: CNUM,CNUM1
! PRIMARY beam
! adjust Generator k-range if required (no monochromator)
      call SOURCE_GET(BEAMSOURCE%INST,PSRC)
      K0=INSOBJ%P_SPEC%KI ! always use instrument setting to define nominal KI
      if (MONO_N.le.0) then
        if (associated(PSRC%X)) then
          call FLUX_SET_RANGE(K0,PSRC%X%LAMW)
          call FLUX_RANGE(PSRC%X%TEMP,kmin,kmax)
          ! write(*,1) 'INST_ADJUST_PREPARE ',K0,PSRC%X%LAMW,kmin,kmax
        else
          call FLUX_SET_RANGE(K0,0.D0)
          call FLUX_RANGE(300.D0,kmin,kmax)
        endif
        call GENERATOR_SETK(K0,kmin,kmax)
        call FLOAT2STR(TWOPI/kmax,CNUM)
        call FLOAT2STR(TWOPI/kmin,CNUM1)
        write(*,*) 'No monochromator: wavelength range set to ',trim(CNUM),' .. ',trim(CNUM1)
      else
    ! no clip range if there is a monochromator on the beam
        call FLUX_SET_RANGE(K0,0.D0)
        call GENERATOR_SETK(K0,0.D0,0.D0)
      endif
! SECONDARY beam
! set FCONE flag for the 1st analyzer according to the INSOBJ%P_SPEC%FLATCONE value
      if (ANAL_N.gt.0) then
        select case(ANALYZERS(1)%CCLS)
        case(CCLS_CRYSTAL)
    ! CRYSTAL analyzers
          call CRYSTAL_GET(ANALYZERS(1)%INST,CR)
          CR%X%FCONE=(INSOBJ%P_SPEC%FLATCONE.ne.0)
        case(CCLS_XTAL)
    ! XTAL analyzers
          call XTAL_GET(ANALYZERS(1)%INST,XT)
          XT%X%FCONE=(INSOBJ%P_SPEC%FLATCONE.ne.0)
        end select
      endif
      end subroutine INST_ADJUST_PREPARE

!--------------------------------------------------------
! Adjust sample exit axis angles for given |Q| in [A^-1]
! The sample gonio angles remain unchanged.
!--------------------------------------------------------
      subroutine INST_ADJUST_AXIS(IERR)
      integer,intent(out) :: IERR
      real(kind(1.D0)) :: KI,KF,Q,PSI
      integer :: KFMODE,im
      character(128) :: S
      ierr=1
      KI=INSOBJ%P_SPEC%KI
      KF=INSOBJ%P_SPEC%KF
      Q=INSOBJ%P_SPEC%Q
      PSI=INSOBJ%P_SPEC%PSI
      KFMODE=INSOBJ%P_SPEC%KFMODE
      if (SAMOBJ%ICLS.gt.0) then
        call SAMPLE_ADJAX(SAMOBJ%P_SAMPLE,KI,KF,Q,PSI,KFMODE,IERR)
        INSOBJ%P_SPEC%FRAME%AX = SAMOBJ%P_SAMPLE%FRAME%AX
        INSOBJ%P_SPEC%FRAME%REXI = SAMOBJ%P_SAMPLE%FRAME%REXI
        INSOBJ%P_SPEC%FRAME%MEXI = SAMOBJ%P_SAMPLE%FRAME%MEXI
      else
        call SPEC_ADJAX(INSOBJ%P_SPEC,KI,KF,Q,PSI,KFMODE,IERR)
      endif
      if (ierr>0) then
        S='Can''t adjust Kf axis for given KI,KF,Q,PSI'
        im= ADD_MESSAGE(m_error,'INST_ADJUST_AXIS',trim(S),m_low)
      endif

      end subroutine INST_ADJUST_AXIS

!--------------------------------------------------------
! set scattering angle, adjust dhkl for a polycryst. sample
!--------------------------------------------------------
      subroutine INST_SET_THETAS(THETAS,IERR)
      real(KIND(1.D0)) :: THETAS
      integer,intent(out) :: IERR
      ierr=0
  ! nothing to do without sample, this is not an error
      if (.not.associated(SAMOBJ%P_SAMPLE)) return
      ierr=1
  ! adjust dhkl
      select case(SAMOBJ%SCLS)
      case(SCLS_POLY)
        SAMOBJ%P_POLY%DHKL=PI/INSOBJ%P_SPEC%KI/sin(THETAS/2.D0)
      end select
      call INST_ADJUST_SAMPLE(IERR)
      end subroutine INST_SET_THETAS

      end module INSTCONTROL

