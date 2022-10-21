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
!////     $Revision: 1.20 $
!////     $Date: 2019/08/15 15:02:07 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Ray-tracing subroutines for choppers
!////
!////////////////////////////////////////////////////////////////////////
      MODULE DCHOPPERS_TRACE
      !use CONSTANTS
      !use CLASSES
    !  use FIELDS
      !use FIELDDATA
      !use FRAMES
      use TRACELIB
      USE TRACINGDATA
      USE DCHOPPERS
      use EVENTLOGS
      USE FRAMES_TRACE
      implicit none
      private
      logical :: dbg=.false.

      public DCHOPPER_GO,DCHOPPER_ADJUST,DCHOPPER_INIT,DCHOPPER_SETPHASE

      contains


!----------------------------------------------------
      LOGICAL FUNCTION DCHOPPER_GO(INST)
! Assume thin rotating disc in the middle between entrance and exit appertures
! The shape, distance and size appertures is defined by the FRAME%SIZE and SHAPE
!----------------------------------------------------
      integer,intent(in) :: INST
      TYPE(DCHOPPER),POINTER :: OBJ
      logical :: LOG1,passed
      integer :: IC,i,NMAX,IFRAME
      real(kind(1.D0)) :: R1(3),K1(3),PHASE,PHASE0
      integer,parameter :: OR(2) = (/2,1/)

      OBJ => ADCHOPPERS(INST)%X
! Convert to local coordinate system
      CALL SLIT_PRE(OBJ%FRAME)
! move to the middle plane (chopper disc)
      call SLIT_MIDDLE(.false.)
! log event at the entry of each object
      call AddEventLog(OBJ%FRAME,NEUT)
      PHASE0=OBJ%FRQ*NEUT%T/HOVM+OBJ%PHI
      IFRAME=NINT(PHASE0)
    ! check frame overlap
      if (.not.OBJ%OVERLAP) then
        !PHASE=OBJ%FRQ*(NEUT%T-TOF_ZERO)/HOVM+OBJ%PHI
        !if (OBJ%FRAME%COUNT>9000) write(*,*) 'DCHOPPER_GO PHASE=',PHASE,(NEUT%T-TOF_ZERO)/HOVM
        if (IFRAME.ne.0) then
          DCHOPPER_GO=.false.
          return
        endif
      endif

! check passage through the entry and exit windows
      LOG1=TLIB_PASS_FRAME(NEUT%R,NEUT%K,OBJ%FRAME%SIZE,OBJ%FRAME%SHAPE)
! treat stopped choppers always in transparent position, regardless of actual phase
      if (.not.OBJ%ACTIVE) goto 10
! check passage through the chopper disc
      ! chopper orientation axis (2 for vertical, 1 for horizontal)
      IC=OR(OBJ%ORI+1)
      !NMAX=1
      !if (OBJ%OVERLAP) NMAX=OBJ%NW
      NMAX=OBJ%NW
      passed=.false.
      i=0
      do while (LOG1.and.(i<NMAX).and.(.not.passed))
        i=i+1
          call TLIB_ROT_POS(NEUT%R,NEUT%K,OBJ%FRQ,OBJ%PHI+OBJ%PHASES(i),3,OBJ%RAD,IC,NEUT%T,R1,K1)
          passed=DCHOPPER_PASS(OBJ,i,R1(1:2))
!         !dbg=(NEUT%T/HOVM>250000)
!            if (dbg.and.passed.and.(OBJ%FRQ<200*1.D-6)) then
!            if (dbg) then
!1     format(a,7(1x,G13.5))
!              !PHASE=OBJ%FRQ*(NEUT%T-TOF_ZERO)/HOVM+OBJ%PHI
!              PHASE=PHASE0+OBJ%PHASES(i)
!              if (i==1) write(*,*)
!              write(*,1) '    frame        i,R,K =',i,R1,K1
!              write(*,1) '    OM*T, PHASE =',PHASE0,PHASE
!              write(*,1) '    T,T0,TSAM =',NEUT%T/HOVM,TOF_ZERO/HOVM,TOF_SAMPLE/HOVM
!              dbg=.false.
!            endif
      enddo
      LOG1=passed
! post-transformation
10    if (LOG1) then
        CALL SLIT_POST(OBJ%FRAME)
        call incCounter(OBJ%FRAME%COUNT)
      ! store T0 when passing through a T0 chopper
        if (OBJ%isT0) then
          NEUT%T0=TOF_ZERO  ! don't do this: +IFRAME*HOVM/OBJ%FRQ
        endif
      else
        NEUT%P=0.D0
      endif

      DCHOPPER_GO=LOG1
      END FUNCTION DCHOPPER_GO

!-------------------------------------------------
      SUBROUTINE DCHOPPER_ADJUST(INST,IERR)
! Update position of the component using nominal trajectory represented by NEUT
! Call FRAME_INIT_MAT first.
! Then use SLIT_PRE & SLIT_POST procedures to propagate NEUT
!-------------------------------------------------
      integer,intent(in) :: INST
      integer,intent(out) :: IERR
      TYPE(DCHOPPER),POINTER :: OBJ
      ! character(32) :: CNUM
1     format(a,': ',6(1x,G12.6))
      IERR=0
      if (.not.DCHOPPER_isValid(INST)) return
      OBJ => ADCHOPPERS(INST)%X
      IERR=0
      call SLIT_PRE(OBJ%FRAME)
      OBJ%T0=NEUT%T
  ! set T0 chopper
      if (OBJ%isT0.and.(OBJ%FRQ.ne.0.D0)) then
      ! first running T0 chopper:
        if (TOF_ZERO.eq.0.D0) then
          TOF_ZERO=NEUT%T  ! don't subtract window phase: -OBJ%PHASES(1)*HOVM/OBJ%FRQ
  ! for two T0 choppers: use the average
        else if (TOF_ZERO2.eq.0.D0) then
          TOF_ZERO2=NEUT%T  ! don't subtract window phase: -OBJ%PHASES(1)*HOVM/OBJ%FRQ
          ! call FLOAT2STR(TOF_ZERO,CNUM)
          !write(*,1) 'DCHOPPER_ADJUST '//trim(OBJ%FRAME%ID)//' T02',TOF_ZERO,TOF_ZERO2,(TOF_ZERO+TOF_ZERO2)/2.D0
          TOF_ZERO=(TOF_ZERO+TOF_ZERO2)/2.D0
        else
          call MSG_WARN(trim(OBJ%FRAME%ID)//': can''t define more than two T0 choppers',1)
        endif
      endif
      call DCHOPPER_READ_FRAMES(OBJ)
      call SLIT_POST(OBJ%FRAME)
      !write(*,1) 'Distance '//trim(OBJ%FRAME%ID),NEUT%K0*NEUT%T
      end SUBROUTINE DCHOPPER_ADJUST


!------------------------------------
      SUBROUTINE DCHOPPER_INIT(INST)
! Called just before tracing
!------------------------------------
      integer,intent(in) :: INST
      TYPE(DCHOPPER),POINTER :: OBJ
      integer :: i
      real(kind(1.D0)) :: SGN, wx, wy, w
	  integer,parameter :: IX(2)=(/1,2/)
      integer,parameter :: IY(2)=(/2,1/)
      if (.not.DCHOPPER_isValid(INST)) return
      OBJ => ADCHOPPERS(INST)%X
      do i=1,OBJ%NW
        SGN=1.D0
        if (OBJ%WIN*OBJ%WIDTHS(i)>0.5D0) SGN=-1.D0
        OBJ%TANWS(i)=SGN*abs(tan(0.5D0*OBJ%WIN*OBJ%WIDTHS(i)*TWOPI))
      enddo
	  ! determine if the chopper is active
	  OBJ%ACTIVE=.true.
	  if ((OBJ%FRQ==0.D0).and.(OBJ%PHI+OBJ%PHASES(1)==0.D0)) then
	    wx = OBJ%FRAME%SIZE(IX(OBJ%ORI+1)) ! chopper frame width
		wy = OBJ%FRAME%SIZE(IY(OBJ%ORI+1))
	    w = (2*OBJ%RAD-wy)*OBJ%TANWS(1) ! 1st window width at the bottom edge
		OBJ%ACTIVE = (w<=wx) ! not active if the window width is bigger than the frame width
		if (OBJ%ACTIVE) then
			write(*,*) 'active chopper: '//trim(OBJ%FRAME%ID)
		endif
	  endif
	  
      END SUBROUTINE DCHOPPER_INIT


!====================================================
! component specific code
!====================================================

!-------------------------------------------------
      SUBROUTINE DCHOPPER_SETPHASE(INST)
! Set chopper phase automatically if required.
! Must be called after all choppers have executed DCHOPPER_ADJUST !!
!-------------------------------------------------
      integer,intent(in) :: INST
      TYPE(DCHOPPER),POINTER :: OBJ
      if (.not.DCHOPPER_isValid(INST)) return
      OBJ => ADCHOPPERS(INST)%X
      if (OBJ%AUTOADJ) then
        if (OBJ%lockT0) then
          OBJ%PHI=-OBJ%FRQ*TOF_ZERO/HOVM
        else
          OBJ%PHI=-OBJ%FRQ*OBJ%T0/HOVM
        endif
      endif
      end SUBROUTINE DCHOPPER_SETPHASE


!----------------------------------------------------
      LOGICAL FUNCTION DCHOPPER_PASS(OBJ,IW,R)
! Check passage through the chopper window
! R is the position in rotating frame
!----------------------------------------------------
      TYPE(DCHOPPER),intent(in) :: OBJ
      integer,intent(in) :: IW
      real(kind(1.D0)),intent(in) :: R(2)
      integer,parameter :: IX(2)=(/1,2/)
      integer,parameter :: IY(2)=(/2,1/)
      real(kind(1.D0)) :: W,X,Y,SGN
      logical :: LOG1
      Y=R(IY(OBJ%ORI+1))

      SGN=SIGN(1.D0,OBJ%TANWS(IW))
    ! SGN>0 if window width < 0.5 period, else SGN<0
    ! for SGN>0, allow only upper position (avoid PI-period)
    ! for SGN<0, allow all events in upper position
      if (SGN*(Y+OBJ%RAD)<0.D0) then
        if (dbg) write(*,*) 'DCHOPPER_PASS 1 ',SGN,Y+OBJ%RAD
        DCHOPPER_PASS=(SGN<0)
        return
      endif
      X=R(IX(OBJ%ORI+1))
      select case(OBJ%PROF)
    ! U-shape
      case(1)
        W=abs((OBJ%RAD)*OBJ%TANWS(IW))
    ! V-shape
      case default
        Y=R(IY(OBJ%ORI+1))
        W=abs((OBJ%RAD+Y)*OBJ%TANWS(IW))
      end select
      LOG1=(ABS(X)*SGN < W*SGN)
      if (dbg) write(*,*) 'DCHOPPER_PASS  ',LOG1,SGN,X,W

      DCHOPPER_PASS=LOG1
      end FUNCTION DCHOPPER_PASS


      end MODULE DCHOPPERS_TRACE
