!//////////////////////////////////////////////////////////////////////
!////  $Id $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2011, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.30 $
!////     $Date: 2019/08/15 15:02:07 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Ray-tracing subroutines for SOURCE
!////
!////////////////////////////////////////////////////////////////////////
      MODULE SOURCES_TRACE
      use TRACELIB
      USE TRACINGDATA
      USE TRACINGOPT
      USE SOURCES
      use SOURCES_TABLE
      use EVENTLOGS
      USE FRAMES_TRACE
      use RNDGEN
      implicit none
      private

      logical dbg
      integer :: idbg
      public SOURCE_GO,SOURCE_ADJUST,SOURCE_INIT,SOURCE_PULSE_RANGE

      contains

!-------------------------------------------------------------
      logical function SOURCE_GO(INST)
! call GO on given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
        SOURCE_GO=SOURCE_GO_INST(ASOURCES(INST)%X)
      end function SOURCE_GO


!---------------------------------------------------------------------
      SUBROUTINE SOURCE_INIT(INST)
!----------------------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(SOURCE),POINTER :: OBJ
      if (.not.SOURCE_isValid(INST)) return
      OBJ => ASOURCES(INST)%X
      CALL FRAME_INIT(OBJ%FRAME)
      dbg=.false.
      idbg=0
      END SUBROUTINE SOURCE_INIT

!-------------------------------------------------
      SUBROUTINE SOURCE_ADJUST(INST,IERR)
! Update position of the component using nominal trajectory represented by NEUT
! Call FRAME_INIT_MAT first.
! Then use SLIT_PRE & SLIT_POST procedures to propagate NEUT
!-------------------------------------------------
      integer,intent(in) :: INST
      integer,intent(out) :: IERR
      TYPE(SOURCE),POINTER :: OBJ
      if (.not.SOURCE_isValid(INST)) return
      OBJ => ASOURCES(INST)%X
      IERR=0
      call SLIT_PRE(OBJ%FRAME)
      call SLIT_TURN(OBJ%FRAME)
      call SLIT_POST(OBJ%FRAME)
      end SUBROUTINE SOURCE_ADJUST

!-------------------------------------------------
      LOGICAL FUNCTION SOURCE_GO_INST(OBJ)
!-------------------------------------------------
      TYPE(SOURCE) :: OBJ
      LOGICAL :: LOG1
      REAL(KIND(1.D0)) :: PP,tin,tout,DT,w
      integer :: SS,it,NTAB
1     format(a,': ',8(1x,G12.5))

      dbg=(dbg.and.(idbg<5))
      if (dbg) idbg=idbg+1

! move to the entry plane and transform to local coordinates
      CALL SLIT_PRE(OBJ%FRAME)
      PP=0.D0
      SS=1
      if (TRACING_UP) SS=-1


! get the number of source tables
      NTAB=SOURCE_TABLE_NUM()
      if (NTAB.le.1) then
        it=NTAB
        w=1.D0
      else
! for multiple tables, select one in random with appropriate weight factor
        it=GENINT(1,NTAB)
        w=1.D0*NTAB
      endif

      if (dbg) write(*,1) 'SOURCE_GO_INST',NEUT%CNT, OBJ%FRAME%COUNT,TRACING_UP, NEUT%P, w

! get intersection times for the source
      LOG1=FRAME_BORDER(OBJ%FRAME,NEUT%R,NEUT%K,tin,tout)
      if (dbg) write(*,1) '     start',NEUT%R(1:3),NEUT%K(1:3),NEUT%T,LOG1
      if (.not.LOG1) GOTO 10

! ensure we are at z=0 (there might be an additional rotation or shift of the component)
!      call SLIT_MIDDLE(.false.)


! Move neutron to the moderator surface
! Neutron is now at the z=0 surface
! The true surface is defined by the FRAME and by the
! surfaces in the source table.
! NOTE: the FRAME must envelope any additional surface.
! Ony positive flight time is allowed for downstream tracing
! and negative for up-stream

! emission surfaces are defined?
      if (LOG1.and.SURFACE_TABLE_DEF(it)) then
        ! move to the entry point
        call TransportG(tin)
        ! get time to the nearest surface (only DT>0 is allowed)
        LOG1=FLUX_SURFACE_TIME(it,NEUT%R,NEUT%K,DT,1)
        ! check that intersectin withthe surface is inside the source frame
        LOG1=(LOG1.and.(DT<(tout-tin)))
        ! reject failed events
        if (dbg) write(*,1) '     entry',NEUT%R(1:3),NEUT%K(1:3),NEUT%T,LOG1
        if (.not.LOG1) GOTO 10
      else
        ! emission surface is defined by the front face of the frame
        if (SS<0) then
          DT=tin
        else
          DT=tout
        endif
      endif
      ! move to the emission surface
      call TransportG(DT)
      if (dbg) write(*,1) '   surface',NEUT%R(1:3),NEUT%K(1:3),NEUT%T,LOG1

! adjust time so that T0 corresponds to the emission surface
! WARNING: this works correctly only if the source is aligned with the origin (no shift, no rotation)
      DT = NEUT%R(3)/NEUT%K(3)
      NEUT%T=NEUT%T-DT

  ! weight by flux distribution
      PP=SOURCE_FLUX(it,OBJ)
      LOG1=(PP.GT.0.D0)
      if (dbg) write(*,1) '   flux',PP,LOG1

! log event at the exit of a source
      IF (LOG1) THEN
        !dbg=((idbg<20).and.(NEUT%R(3)<80))
        !if (dbg) idbg=idbg+1
        !if (dbg) write(*,1) '   transported ',NEUT%R(1:3),NEUT%T,SS*tout

        !if (trace_cnt<10) then
        !  write(*,1) 'SOURCE_GO_INST', it,w,TWOPI/NEUT%K0,NEUT%P,PP
        !endif
        NEUT%P=PP*w*NEUT%P
        call AddEventLog(OBJ%FRAME,NEUT)
        CALL SLIT_POST(OBJ%FRAME)
        call incCounter(OBJ%FRAME%COUNT)
      ENDIF
10    SOURCE_GO_INST=LOG1
      END FUNCTION SOURCE_GO_INST

!----------------------------------------------------
      REAL(KIND(1.D0)) FUNCTION SOURCE_FLUX(it,OBJ)
! it = index to a source table, can be zero
!----------------------------------------------------
      integer :: it
      TYPE(SOURCE) :: OBJ
      REAL(KIND(1.D0)) :: RES,K(3),PW,T,PHASE,kmin,kmax,RH,TN
1     format(a,': ',6(1x,G12.5))
      !logical :: DBG=.true.
      RES=0.D0
! handle frame overlaps
      if (OBJ%TYP.eq.src_pulsed) then
        TN=NEUT%T/HOVM
        PHASE=OBJ%FRQ*TN
        if (OBJ%OVERLAP) then
          T=TN-NINT(PHASE)/OBJ%FRQ-OBJ%DELAY
        else if (abs(PHASE).gt.0.5D0) then
          SOURCE_FLUX=0.D0
          return
        else
          T=TN-OBJ%DELAY
        endif
      endif
      PW=1.D0
! apply wavelength limits
      call FLUX_RANGE(OBJ%TEMP,kmin,kmax)
      if (kmax>kmin) then
        if ((NEUT%K0<kmin).or.(NEUT%K0>kmax)) then
          SOURCE_FLUX=0.D0
          return
        endif
      endif
! spallation source table
      if (SPS_TABLE_DEF(it)) then
        if (OBJ%TYP.eq.src_pulsed) then
          RES=SPS_TAB(it,OBJ%FLUX,T,NEUT%K0)
        else
          RES=OBJ%FRQ*SPS_TAB_AVE(it,OBJ%FLUX,NEUT%K0)
          if (.NOT.TOF_MODE) RES=1.D6*RES
        endif
        ! RES is per microsec. if TOF_MODE
      else
! other flux table
        if (FLUX_TABLE_DEF(it)) then
          RES=FLUX_TAB(it,OBJ%FLUX,NEUT%K0)
        else
! Maxwell source
          RES=FLUX_MAXW(OBJ%FLUX,OBJ%TEMP,NEUT%K0)
        endif
        ! RES is per sec
    ! pulse is defined by:
        if (OBJ%is_pulsed) then
      ! pulse table
          if (PULSE_TABLE_DEF(it)) then
            if (OBJ%TYP.eq.src_pulsed) then
              PW=PULSE_TAB(it,T)
            else
              PW=OBJ%FRQ*1.D6*GET_FLXPULSE_AVE(it) ! flux/pulse ... convert to pulse/sec
            endif
      ! pulse width
          else if (OBJ%PULSW>0.D0) then
            if (OBJ%TYP.eq.src_steady) PW=OBJ%FRQ*OBJ%PULSW
          endif
        endif
        ! TOF mode: integration is done over time in [us]
        ! but RES is per sec
        if (TOF_MODE) PW=PW*1.D-6
      endif

    ! RES returns differential flux dPHI/dK^3 in [1e14/s/cm^2*Ang^3]
    ! tr_lspace requires dPHI/dLambda/dOmega in [1e14/s/cm^2*Ang/sr]
    !  if (TROPT%PSMODE==tr_lspace) then
    !    cosa=abs(NEUT%K(3))/SQRT(NEUT%K(1)**2+NEUT%K(3)**2)
    !    RES=RES*NEUT%K0**4/TWOPI/cosa
    !  endif
      ! clip time if defined
      if (OBJ%PULSW>0.D0) then
        if (abs(T) > OBJ%PULSW/2.D0) PW=0.D0
      endif
      RES=RES*PW
! correct for source inhomogeneity
      if (HOMOG_TABLE_DEF(it)) then
        K=NEUT%K
        if (TRACING_UP)  K=-NEUT%K
        RH=HOMOG_TAB(it,NEUT%R,K,NEUT%K0)
        if (dbg) write(*,1) 'SOURCE_FLUX LAM,FLX,X,RH', TWOPI/NEUT%K0,RES,NEUT%R(1),RH
        RES=RES*RH
        !if (dbg) then
        !  write(*,1) 'SOURCE_FLUX ',NEUT%K0,RH,RES,FLXNR
        !  dbg=.false.
        !endif
      endif
        !if (dbg) then
        !  write(*,1) 'SOURCE_FLUX ',NEUT%K0,RH,RES,FLXNR
        !  dbg=.false.
        !endif
      SOURCE_FLUX=RES
      end FUNCTION SOURCE_FLUX

!---------------------------------------------------------
      SUBROUTINE SOURCE_PULSE_RANGE(OBJ,tmin,tmax)
! Get time range [m/h units] in pulse mode
!---------------------------------------------------------
      TYPE(SOURCE) :: OBJ
      REAL(KIND(1.D0)),intent(out) :: tmin,tmax
      tmin=0.D0
      tmax=0.D0
  ! lookup table
      if (OBJ%TYP.eq.src_pulsed) then
        call PULSE_RANGE(tmin,tmax)
        if (OBJ%PULSW>0.D0) then
          if (tmax>tmin) then
            tmin=max(tmin,-0.5D0*OBJ%PULSW*HOVM)
            tmax=min(tmax,+0.5D0*OBJ%PULSW*HOVM)
          else
            tmin=-0.5D0*OBJ%PULSW*HOVM
            tmax=+0.5D0*OBJ%PULSW*HOVM
          endif
        endif
        ! add (-1;+2) of source period if overlap is allowed
        !if (OBJ%OVERLAP) then
        !  TSRC = HOVM/OBJ%FRQ
        !  tmin = tmin - 1*TSRC
        !  tmax = tmax + 2*TSRC
        !endif
      endif
      end SUBROUTINE SOURCE_PULSE_RANGE

!---------------------------------------------------------
      REAL(KIND(1.D0)) FUNCTION FLUX_MAXW(FLX,TEMP,K)
! Return maxwellian flux for given temperature and |K|-value
! FLX is the integral neutron flux [1e14/s/cm^2]
! Returns dPhi/dK/dOmega in [1e14/s/cm^2/ster*Ang]
! CHANGED 9/10/2007, ver. 1.21:
!     Returns dPhi/d^3K in [1e14/s/cm^2*Ang^3]
!---------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: FLX,TEMP,K
      REAL(KIND(1.D0)) :: K2,VK2,RES
      ! 2*mn*kB/h_bar^2 [Angstrom^-2*K^-1]
      REAL(KIND(1.D0)),parameter :: CMAXW=4.15871716D-2
      K2=K**2
      VK2=CMAXW*TEMP
      if ((K2.lt.25.D0*VK2).and.(K.gt.25.D-4*VK2)) then
        RES=FLX*K*EXP(-K2/VK2)/(2.D0*PI*VK2**2)
      else
        RES=0.D0
      endif
      FLUX_MAXW=RES
      END FUNCTION FLUX_MAXW

      end MODULE SOURCES_TRACE