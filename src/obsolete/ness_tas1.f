
C--------------------------------------------------------
      SUBROUTINE DOUBLE_INI
C Define a beamline with double monochromator
C--------------------------------------------------------
c      USE COMPONENTS
      USE TRACING
      INCLUDE 'const.inc'
      INCLUDE 'inout.inc'
      INCLUDE 'ness_common.inc'
      real(KIND(1.D0)) :: k0,dkx,dky,dkz

! set defaults first
      call TRACING_DEFAULTS
! uncomment the following to override default settings:
  ! down-stream tracing (set 1 for up-stream)
c      TRACING_MODE=1
  ! switch event logs
      EventLogs=.TRUE.
  ! switch variance reduction
c      TRACING_VARI=.TRUE.
  ! swicth max. value optimization
c      TRACING_MAXV=.FALSE.
  ! reporting of tracing statistics
c      TRACING_REPORT_STAT=.FALSE.
  ! reporting of sampling volume
c      TRACING_REPORT_VOL=.FALSE.
  ! reporting of variance reduction progress
c      TRACING_REPORT_VARI=.FALSE.
  ! reporting of tracing progress
c      TRACING_REPORT_PROG=.FALSE.
  ! reporting basic results (intensity, E-spread)
c      TRACING_REPORT_RES=.FALSE.
  ! report max. value optimization
c      TRACING_REPORT_MAXV=.FALSE.

      EFFICIENCY_LIMIT=1.D5

 ! prepare generator
      call GENERATOR_SETDEFAULT(STP.KI)

  ! clean beamline
      call DeleteComponents

  ! add components
      call AddSOURCE(SOU)
      call AddGUIDE(GDEa)
      call AddGUIDE(GUIDE)
      call AddGUIDE(SOL1)
      call AddCRYSTAL(MON)
      call AddGUIDE(SOL3)
      call AddCRYSTAL(ANA)
      call AddGUIDE(SOL4)
      call AddSLIT(DET%FRAME)

! if there is no crystal in the beam:
      if (CRYSTALS_NC.le.0) then
        call SOURCE_RANGE(k0,dkx,dky,dkz)
        call GENERATOR_SETDEFAULT(k0)
      endif

! progress bar update
      PROGRES_UPDATE_INT=3.0  ! time in [s]

      END SUBROUTINE DOUBLE_INI


C--------------------------------------------------------
      SUBROUTINE TAS1_INI
C Define beamline primary TAS
C only a generator and crystal
C--------------------------------------------------------
c      USE COMPONENTS
      USE TRACING
      INCLUDE 'const.inc'
      INCLUDE 'inout.inc'
      INCLUDE 'ness_common.inc'
      real(KIND(1.D0)) :: k0,dkx,dky,dkz

! set defaults first
      call TRACING_DEFAULTS
! uncomment the following to override default settings:
  ! down-stream tracing (set 1 for up-stream)
c      TRACING_MODE=1
  ! switch event logs
c      EventLogs=.TRUE.
  ! switch variance reduction
c      TRACING_VARI=.TRUE.
  ! swicth max. value optimization
c      TRACING_MAXV=.FALSE.
  ! reporting of tracing statistics
c      TRACING_REPORT_STAT=.FALSE.
  ! reporting of sampling volume
c      TRACING_REPORT_VOL=.FALSE.
  ! reporting of variance reduction progress
c      TRACING_REPORT_VARI=.FALSE.
  ! reporting of tracing progress
c      TRACING_REPORT_PROG=.FALSE.
  ! reporting basic results (intensity, E-spread)
c      TRACING_REPORT_RES=.FALSE.
  ! report max. value optimization
c      TRACING_REPORT_MAXV=.FALSE.

 ! prepare generator
      call GENERATOR_SETDEFAULT(STP.KI)

  ! clean beamline
      call DeleteComponents

  ! add TAS1 components
      call AddSOURCE(SOU)
      call AddGUIDE(GDEa)
      call AddGUIDE(GUIDE)
      call AddGUIDE(SOL1)
      call AddCRYSTAL(MON)
      call AddGUIDE(SOL2a)
      call AddGUIDE(SOL2)
      SAM%FRAME%AXY=0.D0 ! ignore scattering angle
      call AddSLIT(SAM%FRAME)

! if there is no crystal in the beam:
      if (CRYSTALS_NC.le.0) then
        call SOURCE_RANGE(k0,dkx,dky,dkz)
        call GENERATOR_SETDEFAULT(k0)
      endif

      END SUBROUTINE TAS1_INI


