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
!////     $Revision: 1.16 $
!////     $Date: 2019/08/16 17:16:25 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Ray-tracing subroutines for CRYSTAL
!////
!////////////////////////////////////////////////////////////////////////
      MODULE CRYSTALS_TRACE
      use TRACELIB
      use RNDGEN
      USE TRACINGDATA
      USE CRYSTALS
      use EVENTLOGS
      use FRAMES_TRACE
      implicit none
      private
      SAVE

      INTEGER,PARAMETER :: MSEG=128
      INTEGER :: NSEG
      REAL(KIND(1.D0)) :: TSEG(0:MSEG),TSEG0(0:MSEG),ALPHA(0:MSEG),GRAD(0:MSEG)
      REAL(KIND(1.D0)) :: PSEG(0:MSEG),PSEG0(0:MSEG),PATH,QKIN
      integer :: ISEG(3,MSEG)


      public CRYSTAL_GO,CRYSTAL_ADJUST,CRYSTAL_INIT

      logical :: dbg=.false.
      integer :: ndbg

      contains

!-------------------------------------------------------------
      logical function CRYSTAL_GO(INST)
! call GO on given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
        CRYSTAL_GO=CRYST_GO(ACRYSTALS(INST)%X)
      end function CRYSTAL_GO

!-------------------------------------------------------------------------
      subroutine CRYSTAL_ADJUST(INST,IERR)
! Update position of the component using nominal trajectory represented by NEUT
! Call FRAME_INIT_MAT first.
! Then use SLIT_PRE & SLIT_POST procedures to propagate NEUT
! Set wavelength and update Bragg angle
! return IERR<>0 if lambda>2*dhkl
!--------------------------------------------------------------------------
      integer,intent(in) :: INST
      integer,intent(out) :: IERR
      TYPE(CRYSTAL),POINTER :: OBJ
      if (.not.CRYSTAL_isValid(INST)) return
      OBJ => ACRYSTALS(INST)%X
      ierr=1
      call FRAME_INIT_MAT(OBJ%FRAME)
      OBJ%LAMBDA=twopi/NEUT%K0
! G-vector
      OBJ%GTOT=2*PI/OBJ%DHKL
      OBJ%G(1)=OBJ%GTOT*SIN(OBJ%CHI)*COS(OBJ%PSI)
      OBJ%G(2)=-SIN(OBJ%PSI)
      OBJ%G(3)=OBJ%GTOT*COS(OBJ%CHI)*COS(OBJ%PSI)
! adjust crystal position
      if (OBJ%AUTOADJ) then
        call CRYSTAL_ORIENT(OBJ,0.D0,ierr)
      else
        ierr=0
      endif
! adjust focal distances
      if (ierr.eq.0) call CRYSTAL_AUTOFOCUS(OBJ)
! propagate NEUT
      call SLIT_PRE(OBJ%FRAME)
      ! define scattering vector
      NEUT%K=NEUT%K+OBJ%G
      !call SLIT_TURN(OBJ%FRAME)
      call SLIT_POST(OBJ%FRAME)
!1     format(a,': ',6(G12.5,2x))

      end subroutine CRYSTAL_ADJUST


! ---------------------------------------
      SUBROUTINE CRYSTAL_INIT(INST)
! initialization of crystal class
! REQUIRES previous call to CRYST_ADJUST
! ---------------------------------------
      integer,intent(in) :: INST
      TYPE(CRYSTAL),POINTER :: OBJ
      real(kind(1.D0)) ::  Z,GETQKIN,GETREFDYN,GETMI,tout,tin
      INTEGER :: I,J
      if (.not.CRYSTAL_isValid(INST)) return
      OBJ => ACRYSTALS(INST)%X

1     format(a,': ',7(G11.4,1x))

  ! define scattering vector
      OBJ%GTOT=TWOPI/OBJ%DHKL
      OBJ%G(1)=OBJ%GTOT*SIN(OBJ%CHI)*COS(OBJ%PSI)
      OBJ%G(2)=-SIN(OBJ%PSI)
      OBJ%G(3)=OBJ%GTOT*COS(OBJ%CHI)*COS(OBJ%PSI)

      CALL FRAME_INIT(OBJ%FRAME)

! G-gradient
      DO I=1,3
        OBJ%MAPG(I)=.FALSE.
        DO J=1,3
          OBJ%DG_DR(I,J)=0.D0
        ENDDO
      ENDDO
      if (OBJ%RH.NE.0.D0) then
        if ((OBJ%TYP.eq.ctyp_gradient).or.(OBJ%TYP.eq.ctyp_bent)) then
          OBJ%DG_DR(1,1)=-COS(OBJ%CHI)*OBJ%GTOT*OBJ%RH
          OBJ%DG_DR(1,3)=SIN(OBJ%CHI)*OBJ%GTOT*OBJ%RH
          OBJ%DG_DR(3,1)=SIN(OBJ%CHI)*OBJ%GTOT*OBJ%RH
          OBJ%DG_DR(2,2)=0.D0    ! no vertical bending
          OBJ%DG_DR(3,3)=-OBJ%POI*COS(OBJ%CHI)*OBJ%GTOT*OBJ%RH
          OBJ%MAPG(1)=.TRUE.
          OBJ%MAPG(3)=.TRUE.
        endif
      endif
      if ((OBJ%DGR.NE.0.D0).and.(OBJ%TYP.eq.ctyp_gradient)) then
          Z=1.D-4*OBJ%GTOT*OBJ%DGR
          OBJ%DG_DR(1,1)=OBJ%DG_DR(1,1)+Z*cos(OBJ%DGA+OBJ%CHI)
          OBJ%DG_DR(1,3)=OBJ%DG_DR(1,3)-Z*sin(OBJ%DGA+OBJ%CHI)
          OBJ%DG_DR(3,1)=OBJ%DG_DR(3,1)-Z*sin(OBJ%DGA+OBJ%CHI)
          OBJ%DG_DR(3,3)=OBJ%DG_DR(3,3)-Z*cos(OBJ%DGA+OBJ%CHI)
          OBJ%MAPG(1)=.TRUE.
          OBJ%MAPG(3)=.TRUE.
      ENDIF
! unit vector |- to G
      OBJ%gama(1)=COS(OBJ%CHI)*COS(OBJ%PSI)
      OBJ%gama(3)=-SIN(OBJ%CHI)*COS(OBJ%PSI)
      OBJ%gama(2)=0
! extinction length
      if (OBJ%LAMBDA.LE.0.D0) OBJ%LAMBDA=OBJ%DHKL
! kin. reflectivity
      OBJ%QHKL=GETQKIN(OBJ,OBJ%lambda)
      ! write(*,*) ' CRYSTAL_INIT DEXT' ,OBJ%DHKL,OBJ%LAMBDA,OBJ%QML
      OBJ%dext=OBJ%DHKL/OBJ%LAMBDA*SQRT(4*PI/OBJ%QML)
! absorption coefficient
      OBJ%MI=GETMI(OBJ,OBJ%lambda,300.D0)
! peak reflectivity (no mosaic)
      OBJ%REF=GETREFDYN(OBJ,OBJ%LAMBDA)
! Darwin box width
      !  write(*,*) ' CRYSTAL_INIT DELTA' ,OBJ%QHKL,OBJ%Dext,PI
      OBJ%DELTA=OBJ%QHKL*OBJ%Dext*1D-4/PI

! special
      OBJ%Ext1=1.D0
      select case(OBJ%TYP)
      case (ctyp_mosaic,ctyp_gradient)
  ! mosaicity is always larger than Darwin box
!    write(*,*) ' CRYSTAL_INIT ' ,OBJ%HMOS,OBJ%DELTA
        OBJ%HMOS=MAX(OBJ%HMOS,OBJ%DELTA)
  ! primary extinction
        Z=OBJ%dlam/OBJ%Dext
        IF (Z.GT.1D-5) OBJ%Ext1=TANH(Z)/Z
      END select

! set temporary TRACING_UP flag
! TRACE_UP is used by SLIT_PRE, but it should always be false when *_INIT is called
!      TRACING_UP=(TRACING_DIR.lt.0)

! get OBJ%DEPTH, mean transmission depth
      NEUT%K0=2.D0*PI/OBJ%LAMBDA
      NEUT%R=(/0.D0,0.D0,0.D0/)
      NEUT%K=(/0.D0,0.D0,NEUT%K0/)
      call SLIT_PRE(OBJ%FRAME)
      call BORDER_BOX(NEUT%R,NEUT%K,OBJ%FRAME%SIZE,tin,tout)
      OBJ%DEPTH=(tout-tin)*NEUT%K0

! get OBJ%P0, reflection probability for the middle ray, without absorption
      NEUT%R=(/0.D0,0.D0,0.D0/)
      NEUT%K=(/0.D0,0.D0,NEUT%K0/)
      QKIN=OBJ%QHKL/10.*OBJ%DW*OBJ%Ext1
      call SLIT_PRE(OBJ%FRAME)    ! convert to local coordinates
    ! for random walk model: move to the entry point in the crystal
      if (ARRAY3D_ENTRY(OBJ%A3D,NEUT%R,NEUT%K,tin)) then
  ! move to the crystal entry
        call Transport(tin)
        call CR_SCAN_SEGMENTS(OBJ,1,0.D0,Z)

      !call SLIT_ENTRY(OBJ%FRAME)  ! go to entry surface
      !call CR_BORDERS(OBJ)
      !call CR_SCATT(OBJ,1,0.D0)
        OBJ%P0=PSEG(NSEG)
      else
        OBJ%P0=1.D0
      endif
      ndbg=0
      dbg=.false.


      !write(*,1) 'CRYSTAL_INIT'//trim(OBJ%REFNAME), OBJ%LAMBDA
      !write(*,1) '    QHKL,dext,delta, ext1',OBJ%QHKL,OBJ%dext,OBJ%DELTA,OBJ%Ext1
      !write(*,1) '    MU, DW ',OBJ%MI,OBJ%DW
      !write(*,1) '    QKIN, P0',QKIN,OBJ%P0

      END SUBROUTINE CRYSTAL_INIT


!-------------------------------------------
      LOGICAL FUNCTION CRYST_GO(CR)
! tracing procedure for all crystal types
!-------------------------------------------
      TYPE (CRYSTAL) :: CR
  !    TYPE(NEUTRON) :: NAUX
1     format(a,': ',6(G10.4,1x))

      !ndbg=ndbg+1
      !dbg=(ndbg<5)
      call ARRAY3D_SETDBG(dbg)
! handle transparent objects
      if (CR%FRAME%TRANSPARENT.EQ.1) then
        CRYST_GO=NONE_GO(CR%FRAME)
        return
      endif
      if (dbg)  write(*,1) trim(CR%FRAME%ID)//' entry axis  ',NEUT%R(1),NEUT%R(3),NEUT%K(1),NEUT%K(3),NEUT%P

      CALL SLIT_PRE(CR%FRAME)

! polarization condition
! quite simplified, should be improved
      ! write(*,*) 'CRYST_GO POL=',NEUT%S(1),' MAG=',CR%MAG
      if (CR%MAG*NEUT%S(1).LT.0.D0) then
        GOTO 99
      endif
! trace through the  array of crystals
      if (dbg)  write(*,*)
      if (dbg)  write(*,1) 'CRYST_GO ',NEUT%CNT
      if (dbg)  write(*,1) trim(CR%FRAME%ID)//' entry local  ',NEUT%R(1),NEUT%R(3),NEUT%K(1),NEUT%K(3),NEUT%P
      select case(CR%TYP)
      case(ctyp_filter)
        call CRYST_FILTER_GO(CR)
      case default
        call CRYSTAL_ARRAY_GO(CR)
      end select
      if (dbg)  write(*,1) trim(CR%FRAME%ID)//' exit local  ',NEUT%R(1),NEUT%R(3),NEUT%K(1),NEUT%K(3),NEUT%P
      if (dbg)  write(*,*)

  ! reject invalid events
      if (NEUT%P.LE.0.D0) goto 99
  ! increment event counter
      call incCounter(CR%FRAME%COUNT)

! log event at the exit
      call AddEventLog(CR%FRAME,NEUT)
! convert to exit coordinates
      CALL SLIT_POST(CR%FRAME)

  !    write(*,1) trim(CR%FRAME%ID)//' post  ',NEUT%R(1),NEUT%R(3),NEUT%K(1),NEUT%K(3),NEUT%P

10    CRYST_GO=.TRUE.
      RETURN

99    NEUT%P=0
      CRYST_GO=.FALSE.
      END FUNCTION CRYST_GO

!-----------------------------------------------------------
      subroutine CRYST_FILTER_GO(CR)
! Treat crystal as a single-crystal filter
! no Bragg reflections, only absorption.
! works in local coordinates
!-----------------------------------------------------------
      TYPE(CRYSTAL) :: CR
      integer :: i
      REAL(KIND(1.D0)) :: MI,GETMI,tin

      ! call SLIT_ENTRY(CR%FRAME)
      if (.not.ARRAY3D_ENTRY(CR%A3D,NEUT%R,NEUT%K,tin)) goto 99
  ! move to the crystal entry
      call Transport(tin)
      ! call ARRAY3D_BORDERS(CR%A3D,NEUT%R,NEUT%K,TSEG0,TSEG,ISEG,NSEG)
      call ARRAY3D_BORDERS_NEW(CR%A3D,NEUT%R,NEUT%K,TSEG0,TSEG,ISEG,NSEG)
  ! count total flight path in the material
      PATH=0.D0
      do i=1,NSEG
        PATH=PATH+TSEG(i)-TSEG0(i)
      enddo
      PATH=PATH*NEUT%K0
  ! calculate absorption for actual wavelength
      MI=GETMI(CR,2.D0*PI/NEUT%K0,300.D0)/10.0
  ! weight event
  !    write(*,*) 'CRYST_FILTER_GO ',NEUT%P,MI,PATH,EXP(-MI*PATH)
  !    read(*,*)
      NEUT%P=NEUT%P*EXP(-MI*PATH)
      return
99    NEUT%P=0.D0
      end subroutine CRYST_FILTER_GO

!-------------------------------------------------------------
      subroutine CRYSTAL_ARRAY_GO(CR)
! Main procedure for tracing in CRYSTAL.
! Starts and ends in local coordinates !!
! INPUT:
!   CR    .. crystal data
! implements only random walk mode
!-------------------------------------------------------------
      TYPE (CRYSTAL),intent(inout) :: CR
      REAL(KIND(1.D0)) :: PIN,MI,TIN,TOUT,G(3)
      REAL(KIND(1.D0)), PARAMETER :: EPS=1D-3
1     format(a,' :',6(1x,G12.5))

  ! remember initial weight
      PIN=NEUT%P
  ! calculate values, which are constant in the course of the iteration
      MI=CR%MI/10.
  ! accumulate path length through the material in PATH
      PATH=0.D0
! make the 1st step
  ! get time to crystal entry, go to end if entry is missed
  ! permit back step if R is inside the array, to ensure the start at the array entry
      if (.not.ARRAY3D_ENTRY(CR%A3D,NEUT%R,NEUT%K,tin)) goto 99
  ! move to the crystal entry
      call Transport(tin)
      if (dbg) write(*,1) 'CRYSTAL_ARRAY_GO entry   ',NEUT%R,NEUT%K


! ****  Begin of multiple scattering iteration cycle  ****

!// generate random walk step in K direction
!---------------------------------------------------
50    CALL CR_WALKSTEP(CR,1,G,TOUT)
!      write(*,*) 'CR_WALKSTEP K  ',TOUT,NEUT%R(3),G
      if (NEUT%P.LE.EPS*PIN) goto 99   ! don't follow trajectories with low probability
  ! log event
      call AddEventLog(CR%FRAME,NEUT)
      call CR_REFLECT(G)
      if (dbg) write(*,1)   '           reflected',NEUT%R,NEUT%K

!// generate 2nd random walk step in K+G direction
!-------------------------------------------------
51    CALL CR_WALKSTEP(CR,-1,G,TOUT)
!      write(*,*) 'CR_WALKSTEP K+G',TOUT,NEUT%R(3),G
      if (NEUT%P.LE.EPS*PIN) GOTO 99
  ! log event
      call AddEventLog(CR%FRAME,NEUT)
    !  if (CR%model.eq.0) TOUT=0.D0
      if (TOUT.LE.0) THEN
  ! reached exit
        NEUT%P=NEUT%P*EXP(-MI*PATH)
      else
  ! continue with random walk
        call CR_REFLECT(G)
        if (dbg) write(*,1) '      back-reflected',NEUT%R,NEUT%K
        GOTO 50
      endif
! end of multiple scattering loop
!-----------------------------------
! ****  End of multiple scattering cycle  ****
      if (dbg) write(*,1) 'CRYSTAL_ARRAY_GO exit',PATH,MI,NEUT%P/PIN
      RETURN
99    NEUT%P=0.D0
      if (dbg) write(*,1) 'CRYSTAL_ARRAY_GO stop',NEUT%P
      end subroutine CRYSTAL_ARRAY_GO


!--------------------------------------------------------------
      SUBROUTINE CR_SCAN_SEGMENTS(CR,DIR,PHI,TOUT)
! Scan through segments along K and calculate cross times, scattering
! probabilities and other fields required for random walk tracing.
! TOUT is the time to the crystal exit point
! Return TOUT=0 if there is no intersection with any segment
!--------------------------------------------------------------
      TYPE (CRYSTAL) :: CR
      integer,intent(in) :: DIR
      REAL(KIND(1.D0)), intent(in) :: PHI
      REAL(KIND(1.D0)), intent(out) :: TOUT
      REAL(KIND(1.D0)) :: TR(3),RLOC(3,3)
      integer :: i,MLOC
      TYPE(NEUTRON) :: NEU0
1     format(a,': ',6(G12.5,2x))
      NEU0=NEUT
! get crossing times for segments
      ! call ARRAY3D_BORDERS(CR%A3D,NEUT%R,NEUT%K,TSEG0,TSEG,ISEG,NSEG)
      call ARRAY3D_BORDERS_NEW(CR%A3D,NEUT%R,NEUT%K,TSEG0,TSEG,ISEG,NSEG)
      IF (NSEG.LE.0) then
        TOUT=0.D0
        return
      endif
      if (dbg) write(*,1) 'CR_SCAN_SEGMENTS ',NEUT%K(3)*TSEG0(1:NSEG)
! scan crossed segments and accumulate probabilities for all segments
      PSEG(0)=0.D0
      do i=1,NSEG
  ! transform NEUT to segment local coordinates
        call ARRAY3D_SEG_TRANS(CR%A3D,ISEG(1,i),TR,RLOC,MLOC)
        call NEUT_TO_SEGMENT(TR,RLOC,MLOC,NEU0,NEUT)
  ! move to the segment entry
        call Transport(TSEG0(i))
    !    if (dbg) then
    !      call R_TO_ARRAY(TR,RLOC,MLOC,NEUT%R,R)
    !      write(*,1) '    iseg, z',i,R(3)
    !    endif
    ! Calculate PSEG(i) ... cummulative scattering probability
    ! and associated fields: PSEG0(i),alpha(i),grad(i)
        CALL CR_SEGSCATT(CR,i,DIR,NEUT%R,NEUT%K,PHI,TSEG(i)-TSEG0(i))
    !    if (dbg) write(*,1) '    PSEG',PSEG0(i),PSEG(i),NEUT%K(3)*TSEG0(i),NEUT%K(3)*TSEG(i)

      enddo
      TOUT=TSEG(NSEG) ! time to the last segment exit
      NEUT=NEU0
      end SUBROUTINE CR_SCAN_SEGMENTS

!--------------------------------------------------------------
      SUBROUTINE CR_WALKSTEP(CR,DIR,G,TOUT)
! Move to the reflection point in the crystal - one step of random walk
! INPUT:
!   DIR  .. 1,-1 for directions K,K+G
! OUTPUT:
!   G(3) .. Local diffraction vector at R
!   TOUT .. Time to reach assembly exit
! NOTE: DIR=-1 AND TOUT=0 stops random walk
!--------------------------------------------------------------
      TYPE (CRYSTAL) :: CR
      integer,intent(in) :: DIR
      REAL(KIND(1.D0)),intent(out) :: G(3),TOUT
      REAL(KIND(1.D0)), PARAMETER :: EPS=1D-3
      INTEGER :: I,M,MLOC
      REAL(KIND(1.D0)) :: Z,PHI,P,DT
      REAL(KIND(1.D0)) :: TR(3),RLOC(3,3),R(3),G1(3)
1     format(a,': ',6(G12.5,2x))

! add random vertical mosaic angle
      PHI=CR%VMOS*GASDEV1(0.D0,3.D0)
      call CR_SCAN_SEGMENTS(CR,DIR,PHI,TOUT)
      IF (TOUT.LE.0) GOTO 99
      P=PSEG(NSEG)   ! scattering probability up to the array exit
      IF (DIR.LT.0) P=1. ! for K+G direction, consider the possibility of leaving the crystal
      IF (P.LT.EPS*CR%P0) GOTO 99 ! don't follow trajectories with low probability
! Select the next segment according to scattering probabilities
      Z=P*RAN1()
      M=1
      DO WHILE (M.LE.NSEG.AND.Z.GE.PSEG(M))
          M=M+1
      END DO
      if (dbg) write(*,1) '     DIR,M,NSEG',DIR,M,NSEG
      if (dbg) write(*,1) '           ISEG',ISEG(1:3,M)
      if (dbg) write(*,1) '           PSEG',PSEG(1:NSEG)
      if (dbg) write(*,1) '           TSEG',TSEG0(1:NSEG)
! No reflection in any segment
      IF(M.GT.NSEG) THEN
        IF(DIR.LT.0) THEN
          GOTO 90 ! go to the last segment exit and return
        ELSE
          GOTO 99 ! event stopped (no reflection)
        ENDIF
      ENDIF
      ! get NEUT in the M-th segment coordinates
      call ARRAY3D_SEG_TRANS(CR%A3D,ISEG(1,M),TR,RLOC,MLOC)
      !if (dbg) write(*,1) 'CR_WALKSTEP M, TR,RLOC',M, TR,RLOC(1:3,1)

! Find point of reflection inside the M-th segment
      CALL CR_SEGTRACE(CR,NEUT%K0,M,DT)
      !if (dbg) write(*,1) 'CR_WALKSTEP DT, TSEG, TSEG0',DT,TSEG(M),TSEG0(M)
      IF (DT.LE.0) GOTO 99

! accumulate flight-path through the material
      DO I=1,M-1
        PATH=PATH+(TSEG(I)-TSEG0(I))*NEUT%K0
      ENDDO
      PATH=PATH+DT*NEUT%K0
      if (dbg) write(*,1) '              TOF',TSEG0(M),DT
! Move to the point of reflection
      call Transport(DT+TSEG0(M))
      if (dbg) write(*,1) '      transported',NEUT%R(1),NEUT%R(3),NEUT%K(1),NEUT%K(3)
! Calculate local G-vector
      call R_TO_SEGMENT(TR,RLOC,MLOC,NEUT%R,R)
      CALL CR_LOCALG(CR,DIR,R,ALPHA(M)+GRAD(M)*DT,PHI,G1)
      CALL M3XV3(-MLOC,RLOC,G1,G)
      NEUT%P=NEUT%P*P
      RETURN

! exit when going along K+G (DIR<0)
90    call Transport(TSEG(NSEG))
      ! add path
      DO I=1,NSEG
          PATH=PATH+(TSEG(I)-TSEG0(I))*NEUT%K0
      ENDDO
      if (dbg) write(*,1) '  transported out',NEUT%R(1),NEUT%R(3),NEUT%K(1),NEUT%K(3)
      TOUT=0.D0
      RETURN

! Left without reflection
99    NEUT%P=0.D0
      TOUT=0.D0

      END SUBROUTINE CR_WALKSTEP


!-----------------------------------------------------------------------------
      SUBROUTINE CR_SEGSCATT(CR,J,DIR,R,K,PHI,DT)
! Calculate scattering cross-section probability for J-th segment on the path DT along K(3)
! INPUT:
!   J    .. segment index
!   DIR  .. 1,-1 for directions K,K+G
!   R,K  .. starting neutron coordinates
!   PHI  .. vertical mosaic angle
!   DT   .. time-of-flight along K(3) through the segment
! OUTPUT:
!   stored in /CRBORDERS/
!-----------------------------------------------------------------------------
      TYPE (CRYSTAL) :: CR
      integer,intent(in) :: J,DIR
      real(KIND(1.D0)), intent(in) :: R(3),K(3),PHI,DT
      real(KIND(1.D0)) :: TAU,sigma,K0
1     format(a,': ',7(G11.4,1x))
      IF (DT.LE.0) THEN
          PSEG0(J)=0.D0
          PSEG(J)=PSEG(J-1)
          RETURN
      ENDIF
      K0=SQRT(K(1)**2+K(2)**2+K(3)**2)

      select case(CR%TYP)
! MOSAIC CRYSTAL and MOSAIC CRYSTAL WITH G-GRADIENT
      case(ctyp_mosaic,ctyp_gradient)
        call CR_BRAGGDEV(CR,DIR,R,K,PHI,alpha(J),grad(J))
! ONLY G-GRADIENT
      case(ctyp_bent)
      ! no double reflection in the same segment
        IF ((TSEG0(J).NE.0.D0).or.(DIR.GT.0)) then
          call CR_FINDREF(CR,DIR,R,K,TAU,grad(J))
          alpha(j)=-grad(j)*TAU
        endif
      END select
      sigma=K0*QKIN*CR_SEGPOWER(CR,J,DIR,DT)
    !  if (dbg) write(*,1) 'CR_SEGSCATT  K0, QKIN',K0,QKIN
    !  if (dbg) write(*,1) 'CR_SEGSCATT  J ... sigma',J,alpha(J),grad(J),K0*DT,sigma
!// get scattering probability PSEG0 (from this segment) and PSEG (cumulative)
      IF(sigma.gt.14) then
        PSEG0(J)=1.D0
        PSEG(J)=1.D0
      else
        PSEG0(J)=1.D0-exp(-sigma)
        PSEG(J)=1.D0-(1.D0-PSEG(J-1))*(1.D0-PSEG0(J))
      endif

      END SUBROUTINE CR_SEGSCATT


!-----------------------------------------------------------------------------
      real(KIND(1.D0)) FUNCTION CR_SEGPOWER(CR,J,DIR,DT)
! Calculate SIGMA=scattering power for J-th segment on the path DT along K(3)
! for unit |K| and Qkin, i.e. scattering probability = 1-exp(-K*Qkin*SIGMA)
! Requires ALPHA(J),GRAD(J) arrays to be ready for use
! INPUT:
!   J    .. segment index
!   DIR  .. 1,-1 for directions K,K+G
!   DT   .. time-of-flight along K(3) through the segment
!-----------------------------------------------------------------------------
      TYPE (CRYSTAL) :: CR
      integer,intent(in) :: J,DIR
      real(KIND(1.D0)), intent(in) :: DT
      real(KIND(1.D0)) :: Z1,Z2,TAU,sigma
      real(KIND(1.D0)) :: HMOS_DIST,ERF_CR
      sigma=0.D0
      select case(CR%TYP)
! MOSAIC CRYSTAL
      case(ctyp_mosaic)
        sigma=HMOS_DIST(alpha(J)/CR%HMOS)/CR%HMOS*DT
! MOSAIC CRYSTAL WITH G-GRADIENT
      case(ctyp_gradient)
        Z1=ERF_CR((alpha(J)+grad(J)*DT)/CR%HMOS,0)
        Z2=ERF_CR(alpha(J)/CR%HMOS,0)
        sigma=(Z1-Z2)/grad(J)
! ONLY G-GRADIENT
      case(ctyp_bent)
        sigma=0.D0
      ! no double reflection in the same segment
        IF ((TSEG0(J).NE.0.D0).or.(DIR.GT.0)) then
          TAU=-alpha(j)/grad(j)
          IF ((TAU.GT.0).AND.(TAU.LT.DT)) sigma=1.D0/abs(grad(J))
        endif
      END select
      CR_SEGPOWER=sigma

      END FUNCTION CR_SEGPOWER


!-------------------------------------------------------------------------
      SUBROUTINE CR_SEGTRACE(CR,K0,M,DT)
! return random step size to a point of reflection inside the M-th segment
!-------------------------------------------------------------------------
      TYPE (CRYSTAL) :: CR
      INTEGER*4 M
      REAL*8 ksi,Z,Z1,sigma,P0,DT,AA
      REAL*8 K0
      REAL*8 ERF_CR,HMOS_DIST
1     format(a,': ',6(G12.5,2x))
      DT=0.D0

! MOSAIC CRYSTAL
      !if (dbg) write(*,1) 'CR_SEGTRACE ',M,GRAD(M),ALPHA(M),PSEG0(M)
      !if (dbg) write(*,1) 'CR_SEGTRACE ',K0,QKIN
      IF(ABS(GRAD(M)).LT.1.D-6*K0*QKIN) THEN
          ksi=RAN1()
          P0=PSEG0(M)
          AA=ALPHA(M)/CR%HMOS
          sigma=K0*QKIN*HMOS_DIST(AA)/CR%HMOS
          !if (dbg) write(*,1) 'CR_SEGTRACE ',P0,AA,CR%HMOS,QKIN,sigma

          DT=-Log(1-ksi*P0)/sigma
      ELSE
! MOSAIC CRYSTAL WITH D-GRADIENT
        IF(CR%HMOS.GT.SEC) THEN
          ksi=RAN1()
          P0=PSEG0(M)
          Z=GRAD(M)*Log(1-ksi*P0)/(K0*QKIN)
          AA=ALPHA(M)/CR%HMOS
!  ERF_CR shoud be used only on (-inf.;0) because of num. precision
          IF(AA.GT.0) THEN
            Z1=ERF_CR(-AA,0)
            IF(1.-Z1-Z.GT.0.5D0) THEN
              DT=-ERF_CR(Z1+Z,1)-AA
            ELSE
              DT=ERF_CR(1.-Z1-Z,1)-AA
            ENDIF
          ELSE
            Z1=ERF_CR(AA,0)
            IF(Z1-Z.GT.0.5D0) THEN
              DT=-ERF_CR(1.-Z1+Z,1)-AA
            ELSE
              DT=ERF_CR(Z1-Z,1)-AA
            ENDIF
          ENDIF
           DT=DT*CR%HMOS/GRAD(M)
! ONLY D-GRADIENT
        ELSE
          DT=-ALPHA(M)/GRAD(M)
        ENDIF
      ENDIF

      END SUBROUTINE CR_SEGTRACE



!-----------------------------------------------------------------------------
      SUBROUTINE CR_BRAGGDEV(CR,DIR,R,K,PHI,THETA,THETA1)
! Calculate local angular deviation from Bragg condition and its gradient
! INPUT:
!   DIR  .. 1,-1 for directions K,K+G
!   PHI  .. vertical mosaic angle
!   R,K  .. neutron coordinates
! OUTPUT:
!   THETA  .. angular deviation
!   THETA1 .. dTheta/dT, T is the flight path in units of [h_bar/m]
!-----------------------------------------------------------------------------
      TYPE (CRYSTAL) :: CR
      integer,intent(in) :: DIR
      real(KIND(1.D0)), intent(in) :: R(3),K(3),PHI
      real(KIND(1.D0)), intent(out) :: THETA,THETA1

      INTEGER*4 I,I1
      real(KIND(1.D0)) :: a,b,Z,KK,GABS
      real(KIND(1.D0)) :: KaG(3),G(3)

      KK=0.
      DO I=1,3
         KK=KK+K(I)**2
      END DO

! get local G vector, including vertical mosaic angle
      CALL CR_LOCALG(CR,DIR,R,0.D0,PHI,G)

! get angular deviation from the Bragg condition
      GABS=0.
      DO I=1,3
         KaG(I)=K(I)+G(I)
         GABS=GABS+G(I)**2
      END DO
      GABS=SQRT(GABS)
      a=KaG(1)**2+ KaG(2)**2+KaG(3)**2-KK    ! (K+G)^2-K^2
      b=DIR*GABS*(KaG(1)*CR%GAMA(1)+KaG(2)*CR%GAMA(2)+KaG(3)*CR%GAMA(3))

! get angular deviation from the Bragg condition (=alpha)
      THETA=-a/(2.*b)
! get gradient of the angular deviation
      THETA1=0.D0
      DO I=1,3
        IF (CR%MAPG(I)) THEN
          Z=0.D0
          DO I1=1,3
            Z=Z+DIR*CR%DG_DR(I,I1)*K(I1)
          ENDDO
          THETA1=THETA1+KaG(I)*Z
        ENDIF
      END DO
      THETA1=-THETA1/b

!// add random angle within Darwin box width
!        THETA=THETA+(RAN1(1)-0.5)*Q*CR%DEXT*1D-3/PI

      END SUBROUTINE CR_BRAGGDEV


! -------------------------------------------------------------
      SUBROUTINE CR_LOCALG(CR,DIR,R,ETA,PHI,G)
! Calculate local G-vector at R
! R is in segment coordinates
! ETA .. horizontal tilt angle of mosaic domain
! PHI .. vertical tilt angle of mosaic domain
! -------------------------------------------------------------
      TYPE (CRYSTAL) :: CR
      integer, intent(in) :: DIR
      real(KIND(1.D0)), intent(in) :: R(3),ETA,PHI
      real(KIND(1.D0)), intent(out) :: G(3)
      INTEGER :: I,J
      real(KIND(1.D0)) :: W(3),G0(3)
      real(KIND(1.D0)) :: Z,GABS,GII,GT

! Calculate local G-vector (only G-gradient)
      DO I=1,3
        G0(I)=CR%G(I)
        IF (CR%MAPG(I)) THEN
          W(I)=0.D0
          DO J=1,3
            W(I)=W(I)+CR%DG_DR(I,J)*R(J)
          ENDDO
          G0(I)=G0(I)+W(I)
        ENDIF
      END DO
      GII=SQRT(CR%G(1)**2+CR%G(3)**2)
      GT=CR%G(2)
      GABS=SQRT(G0(1)**2+G0(2)**2+G0(3)**2)

! Add the angle of the mosaic block
      G(1)=G0(1)+G0(3)*ETA - GT*CR%G(1)/GII*PHI
      G(2)=G0(2)+GII*PHI
      G(3)=G0(3)-G0(1)*ETA - GT*CR%G(3)/GII*PHI
! Renormalize
      Z=GABS/SQRT(G(1)**2+G(2)**2+G(3)**2)
      DO i=1,3
         G(i)=DIR*G(i)*Z
      ENDDO
      END SUBROUTINE CR_LOCALG


!-----------------------------------------------------------------------------
      SUBROUTINE CR_FINDREF(CR,DIR,R,K,DT,GR)
! Calculate TOF to the point of reflection in BENT CRYSTAL
! INPUT:
!   DIR  .. 1,-1 for directions K,K+G
!   R,K  .. neutron coordinates
! OUTPUT:
!   DT  .. time of flight in [h_bar/m]
!   GR  .. gradient of Bragg angle = (K+G).DG_DR.K /[G(K+G).gamma]
!-----------------------------------------------------------------------------
      TYPE (CRYSTAL) :: CR
      integer,intent(in) :: DIR
      real(KIND(1.D0)), intent(in) :: R(3),K(3)
      real(KIND(1.D0)), intent(out) :: DT,GR
      real(KIND(1.D0)),parameter :: EPS=1.D-20
      INTEGER*4 I,I1
      real(KIND(1.D0)) :: a,b,c,Z,KK,GABS,signb,det
      real(KIND(1.D0)) :: KaG(3),G(3)

1     format(a,' :',6(G10.4,2x))
      KK=0.
      DO I=1,3
         KK=KK+K(I)**2
      END DO

! get local G vector
      CALL CR_LOCALG(CR,DIR,R,0.D0,0.D0,G)

! get local deviation from the Bragg condition
      GABS=0.
      DO I=1,3
         KaG(I)=K(I)+G(I)
         GABS=GABS+G(I)**2
      END DO
      GABS=SQRT(GABS)
      c=KaG(1)**2+ KaG(2)**2+KaG(3)**2-KK    ! (K+G)^2-K^2
! get linear term, b=(K+G).DG_DR.K
! and quadratic term, a=|DG_DR.K|^2
      b=0.D0
      a=0.D0
      DO I=1,3
        IF (CR%MAPG(I)) THEN
          Z=0.D0
          DO I1=1,3
            Z=Z+CR%DG_DR(I,I1)*K(I1)
          ENDDO
          a=a+Z**2
          b=b+KaG(I)*Z
        ENDIF
      ENDDO
      GR=-b/GABS/(KaG(1)*CR%GAMA(1)+KaG(2)*CR%GAMA(2)+KaG(3)*CR%GAMA(3))
      b=DIR*b
  ! solve quadratic equation: get the smaller solution
      signb=SIGN(1.D0,b)
      det=b*b-a*c
      if ((a.gt.EPS).and.(det.gt.EPS)) then
        DT=signb*(-abs(b)+SQRT(det))/a
      else
        DT=1.D30
      endif

      END SUBROUTINE CR_FINDREF

!-----------------------------------------------------------
      subroutine CR_REFLECT(G)
! Do reflection K -> K+G, ensure energy conservation
!-----------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: G(3)
      REAL(KIND(1.D0)) :: Z
      integer :: i
  ! set K=K+G
      DO I=1,3
        NEUT%K(I)=NEUT%K(I)+G(I)
      ENDDO
  ! renormalize
      Z=NEUT%K0/SQRT(NEUT%K(1)**2+NEUT%K(2)**2+NEUT%K(3)**2)  ! set |K+G| = |K|
      DO I=1,3
        NEUT%K(I)=NEUT%K(I)*Z
      ENDDO
      end subroutine CR_REFLECT

      end MODULE CRYSTALS_TRACE