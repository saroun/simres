!//////////////////////////////////////////////////////////////////////
!////  $Id: guides_trace.f,v 1.36 2019/08/15 17:24:08 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2016, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.36 $
!////     $Date: 2019/08/15 17:24:08 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Ray-tracing subroutines for CRYSTAL
!////
!////////////////////////////////////////////////////////////////////////
      MODULE GUIDES_TRACE
      use TRACELIB
      USE TRACINGDATA
      USE GUIDES
      use EVENTLOGS
      use FRAMES_TRACE
      use RNDGEN
      use NSTORE
      use TABLE_CRYSTALS
      use MIRROR_TABLE
      implicit none
      private

      type(NEUTRON) :: TMPN(GUIDES_DIM),DBGN

      public GUIDE_GO,GUIDE_ADJUST,GUIDE_INIT
      public GUIDE_CLR,GUIDE_ENDTRACE

      contains


!-------------------------------------------------------------
      logical function GUIDE_GO(INST)
! call GO on given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
        GUIDE_GO=GUIDE_GO_INST(AGUIDES(INST)%X,INST)
      end function GUIDE_GO


!------------------------------------
      SUBROUTINE GUIDE_INIT(INST)
!------------------------------------
      integer,intent(in) :: INST
      TYPE(GUIDE),POINTER :: OBJ
      ! write(*,*) 'GUIDE_INIT '
      if (.not.GUIDE_isValid(INST)) return
      OBJ => AGUIDES(INST)%X
      select case (OBJ%TYP)
      case(-1)
        CALL FRAME_INIT(OBJ%FRAME)
      case (2,3)
        CALL PGUIDE_INIT(OBJ)
      case (4)
        CALL EGUIDE_INIT(OBJ)
      case DEFAULT
        CALL GDE_INIT(OBJ)
      END select
      ! deselect logger = it will be reallocated later, before simulation starts (TRACING)
      OBJ%MLOG=0
      GUIDES_idbg=0
      GUIDES_dbg=.false.
      END SUBROUTINE GUIDE_INIT

!-------------------------------------------------------------
      SUBROUTINE GUIDE_CLR(INST)
! return reference to given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(GUIDE),pointer :: OBJ
      ! write(*,*) 'GUIDE_CLR '
        if (.not.GUIDE_isValid(INST)) return
        OBJ => AGUIDES(INST)%X
        OBJ%FRAME%COUNT=0
        ! clear nstore registry
         call NSTORE_CLEAR(OBJ%FRAME%REGISTRY)
        ! clear bounce logger)
        call MLOGS_REC_CLR(OBJ%MLOG)
      end SUBROUTINE GUIDE_CLR

!--------------------------------------------------------------------------
      subroutine GUIDE_ENDTRACE(INST,CNT,P)
! Add previously registered event to the monitor storage, using given weight
! update internal oger of wall bounces
!--------------------------------------------------------------------------
      integer,intent(in) :: INST,CNT
      REAL(kind(1.D0)),intent(in)  :: P
      TYPE(PGUIDE) :: OBJ
1     format(a,6(1x,G12.5))
        OBJ%X => AGUIDES(INST)%X
        if (OBJ%X%FRAME%REGISTRY>0) then
          ! transform to incident axis coordinates
          call FRAME_LOCAL(-1,OBJ%X%FRAME,TMPN(INST))
          TMPN(INST)%P=P     ! set the final weight to the event
          TMPN(INST)%CNT=CNT ! must set counter to neutron before calling NSTORE_SETEVENT !!
          call NSTORE_SETEVENT(OBJ%X%FRAME%REGISTRY,TMPN(INST))
        endif
        if (OBJ%X%LOGBNC>0) then
          ! write(*,1) 'GUIDE_ENDTRACE ',INST,P
          call  MLOGS_REC_UPD(OBJ%X%MLOG,P)
        endif
      end subroutine GUIDE_ENDTRACE

!-------------------------------------------------
      SUBROUTINE GUIDE_ADJUST(INST,IERR)
! Update position of the component using nominal trajectory represented by NEUT
! Call FRAME_INIT_MAT first.
! Then use SLIT_PRE & SLIT_POST procedures to propagate NEUT
!-------------------------------------------------
      integer,intent(in) :: INST
      integer,intent(out) :: IERR
      integer :: SS
      REAL(KIND(1.D0)) :: ROPT
      TYPE(GUIDE),POINTER :: OBJ
1     format(a,': ',7(1x,G12.6))
      if (.not.GUIDE_isValid(INST)) return
      ! write(*,*) 'GUIDE_ADJUST '
      OBJ => AGUIDES(INST)%X
      IERR=0
      ! adjust curvature of a curved guide if requiered
      if ((OBJ%TYP == 1).and.(OBJ%NODIR).and.(OBJ%FRAME%SIZE(3)>0.D0)) then
        if (abs(OBJ%ROH)>0.D0) then
          SS=NINT(sign(1.D0,OBJ%ROH))
          ROPT=4.D0*(OBJ%FRAME%SIZE(1)+OBJ%W2)/OBJ%FRAME%SIZE(3)**2
          if (abs(OBJ%ROH)<ROPT) OBJ%ROH=SS*ROPT
        else if (abs(OBJ%ROV)>0.D0) then
          SS=NINT(sign(1.D0,OBJ%ROV))
          ROPT=4.D0*(OBJ%FRAME%SIZE(2)+OBJ%H2)/OBJ%FRAME%SIZE(3)**2
          if (abs(OBJ%ROV)<ROPT) OBJ%ROV=SS*ROPT
        endif
      endif
      if (OBJ%TYP<=0) OBJ%LOGBNC=0
      call FRAME_INIT_MAT(OBJ%FRAME)
      call SLIT_PRE(OBJ%FRAME)
      call SLIT_POST(OBJ%FRAME)
      ! translucent material is allowed only for straight guides
      if (OBJ%TYP.ne.1) OBJ%MATER=0
      ! set material to absorbing if 1/MU>1 um*A
      select case (OBJ%MATER)
      case(0)
        OBJ%MU=1.D4
      case(2)
        OBJ%MU=GET_CR_ABS(MU_Si,NEUT%K0)/TWOPI
      case(3)
        OBJ%MU=GET_CR_ABS(MU_Al2O3,NEUT%K0)/TWOPI
      end select

      if (OBJ%MONITOR) then
        write(*,1) trim(OBJ%FRAME%ID)//' is monitor'
      endif
      !write(*,1) 'Distance '//trim(OBJ%FRAME%ID),NEUT%K0*NEUT%T
      end SUBROUTINE GUIDE_ADJUST


!----------------------------------------------------
      LOGICAL FUNCTION GUIDE_GO_INST(OBJ,INST)
!----------------------------------------------------
      TYPE(GUIDE) :: OBJ
      integer,intent(in) :: INST
      INTEGER :: I
      REAL(KIND(1.D0)) :: dj,OSCD,OSCA,DUM
      COMMON /OSCBEND/ OSCD,OSCA
1     format(a,6(1x,G12.5))
! handle transparent objects
      if (OBJ%FRAME%TRANSPARENT.EQ.1) then
        GUIDE_GO_INST=NONE_GO(OBJ%FRAME)
        return
      endif
! generate misalignent using the value of OBJ%MISALIGN
      call GUIDE_MISALIGN(OBJ)
! Convert to local coordinate system
      CALL FRAME_PRE(OBJ%FRAME)
! transparent collimator:
      IF ((OBJ%TYP.LT.0).OR.(OBJ%FRAME%SIZE(3).LE.0)) GOTO 700
! move to the bender entry window
      call FRAME_ENTRY(OBJ%FRAME)
! log event at the entry of each object
      call AddEventLog(OBJ%FRAME,NEUT)
! Make a random shift to simulate oscillating collimator
      OSCD=0
      OSCA=0
      IF (OBJ%OSCILLATE) THEN
        dj=(RAN1()-0.5)/OBJ%NLH
        OSCA=dj*(OBJ%W2-OBJ%FRAME%SIZE(1))/OBJ%FRAME%SIZE(3)
        OSCD=dj*OBJ%FRAME%SIZE(1)
      ENDIF
! reset bounce recorder
      call MLOGS_RESET(OBJ%MLOG)
      !GUIDES_dbg=.false.
      if (OBJ%FRAME%SHAPE.eq.FRAME_SHAPE_DISC) then
        CALL TUBE_GO(OBJ)
      else
        select case(OBJ%TYP)
        case(0)
          CALL SOLLER_GO(OBJ)
        case(1)
          CALL BENDER_GO(OBJ)
        case(2,3,4)
          CALL FGUIDE_GO(OBJ)
        case DEFAULT
          NEUT%STATE=1 ! absorb
          NEUT%P=0.D0
          GOTO 300
        END select
      endif
! event accepted
700   continue

!--------------------------
! finalize
! absorbed
300   IF (NEUT%STATE>0) then
  ! stop unless monitoring of absorbed events is required
        if ((.not.OBJ%MONITOR).or.(NEUT%P<=0.D0)) then
          NEUT%P=0.D0
          GUIDE_GO_INST=.FALSE.
          if (GUIDES_dbg)  write(*,1) 'GUIDE_GO_INST STOPPED ',NEUT%CNT
          return
        endif
      ! label neutron with its state number for monitoring
        NEUT%LABEL=NEUT%STATE
      else
        if (OBJ%CLOSED) then
          NEUT%P=0.D0
          GUIDE_GO_INST=.FALSE.
          if (GUIDES_dbg)  write(*,1) 'GUIDE_GO_INST CLOSED'
          return
        endif
      endif
! passed
      if (GUIDES_dbg)  write(*,1) 'GUIDE_GO_INST PASSED ',NEUT%CNT
! renormalize K due to numerical errors
      DUM=SQRT(NEUT%K(1)**2+NEUT%K(2)**2+NEUT%K(3)**2)
      DO I=1,3
        NEUT%K(I)=NEUT%K(I)*NEUT%K0/DUM
      ENDDO
  ! log event at the exit of a guide
      call AddEventLog(OBJ%FRAME,NEUT)
  ! register event
      TMPN(INST)=NEUT
      if (TRACING_UP) TMPN(INST)%K=-TMPN(INST)%K
  ! post-transformation to exit coordinates, correct for deflection
      ! DBGN=NEUT
      CALL FRAME_POST(OBJ%FRAME)
      call incCounter(OBJ%FRAME%COUNT)
      GUIDE_GO_INST=.TRUE.


! DEBUG
      !return
      !if (abs(DBGN%K(1)/DBGN%K(3))>0.005) then
      !if (abs(NEUT%R(1))<1.5) then
      !  GUIDES_idbg=GUIDES_idbg+1
      !  GUIDES_dbg=(GUIDES_idbg<30)
      !endif
      !endif
      IF (GUIDES_dbg) write(sout,1)  'GUIDE_GO_INST EXIT: ',OBJ%FRAME%COUNT
      IF (GUIDES_dbg.and.(OBJ%MLOG>0))  call MLOGS_PRINT_TMP(OBJ%MLOG,1,6)
      IF (GUIDES_dbg) write(sout,1)  '           R,K: ',DBGN%R,DBGN%K
      !IF (GUIDES_dbg) write(sout,1)  '      POST R,K: ',NEUT%R,NEUT%K
      GUIDES_dbg=.false.
!-------------------------------


      END FUNCTION GUIDE_GO_INST


!--------------------------------------------------------------
      SUBROUTINE SOLLER_GO(OBJ)
! GO procedure for a simple collimator (non reflecting, TYP=0)
! INPUT:   assume R,K at the entry in local coordinates !
! RETURN:  R,K at the exit in local coordinates, P=P*transmission, T=T+passage time
!--------------------------------------------------------------
      TYPE(GUIDE) :: OBJ
      LOGICAL*4 LOG1, BENDER_PASS
      INTEGER*4 IH,IV,IH1,IV1
      real(kind(1.D0)) :: DT
1     format(a,': ',7(G10.4,1x))
!  check passage through the entry
      LOG1=BENDER_PASS(OBJ,NEUT%R,.true.,IH,IV)
      !if (NEUT_DBG) then
      !  write(*,1) 'SOLLER_GO start: ', LOG1, IH, IV
      !endif
      IF (.NOT.LOG1) then
        NEUT%P=0.D0
        GOTO 100
      endif
!  move to the exit
      if (TRACING_UP) then
        DT=-NEUT%R(3)/NEUT%K(3)
      else
        DT=(OBJ%FRAME%SIZE(3)-NEUT%R(3))/NEUT%K(3)
      endif
      call TransportG(DT)
!  check passage through the same slit at the exit

      !if (NEUT_DBG) then
      !  I = BENDER_INLAM(OBJ,NEUT%R,.true.,IH1,IV1)
      !  write(*,1) 'SOLLER_GO I, IH1, IH2: ', I, IH1, IV1
      !endif

      LOG1=(LOG1.AND.BENDER_PASS(OBJ,NEUT%R,.true.,IH1,IV1))
      IF ((.NOT.LOG1).OR.(IH1.NE.IH).OR.(IV1.NE.IV)) GOTO 100
      RETURN

! no passage
100   NEUT%STATE=1
      END SUBROUTINE SOLLER_GO


!--------------------------------------------------------------
      SUBROUTINE TUBE_GO(OBJ)
! GO procedure for a coarse circular collimator (SHAPE=0)
!--------------------------------------------------------------
      TYPE(GUIDE) :: OBJ
      LOGICAL :: LOG1
      real(kind(1.D0)) :: DT
! up-stream
      if (TRACING_UP) then
    !  check passage through the exit window
        LOG1=(((NEUT%R(1)/OBJ%W2)**2+(NEUT%R(2)/OBJ%H2)**2).lt.2.5D-1)
        IF (LOG1) then
    !  move to the entry
          DT=-NEUT%R(3)/NEUT%K(3)
          call TransportG(DT)
    !  check passage through the entry
          LOG1=(((NEUT%R(1)/OBJ%FRAME%SIZE(1))**2+(NEUT%R(2)/OBJ%FRAME%SIZE(2))**2).lt.2.5D-1)
        endif
! down-stream
      else
    !  check passage through the entry
        LOG1=((NEUT%R(1)/OBJ%FRAME%SIZE(1))**2+(NEUT%R(2)/OBJ%FRAME%SIZE(2))**2).lt.2.5D-1
        IF (LOG1) then
    !  move to the exit
          DT=(OBJ%FRAME%SIZE(3)-NEUT%R(3))/NEUT%K(3)
          call TransportG(DT)
    !  check passage through the exit
          LOG1=(((NEUT%R(1)/OBJ%W2)**2+(NEUT%R(2)/OBJ%H2)**2).lt.2.5D-1)
        endif
      endif
      if (.not.LOG1) then
        NEUT%P=0.D0
        NEUT%STATE=1
      endif
      END SUBROUTINE TUBE_GO

!----------------------------------------------------
      SUBROUTINE FGUIDE_GO(OBJ)
! GO procedure for elliptic or parabolic guide guide (TYP=2,3,4)
! INPUT:   NEUT at the entry in local coordinates
! RETURN:  NEUT at the exit in local coordinates
! On exit, NEUT%STATE must have one of these values:
! STATE=0 .. passed
! STATE=1 .. escaped by missing the guide, transition to substrate, etc.
! STATE=2 .. absorbed in the coating
! NOTE: absorption in coating is considered only if MONITOR option is on.
! Set NEUT%P=0 if you wish to stop the event (=> no escape monitoring).
!----------------------------------------------------
      TYPE(GUIDE) :: OBJ
      REAL(kind(1.D0)),parameter :: Q2Ni=0.5D0*QTNi
      LOGICAL :: BENDER_PASS
      INTEGER :: IH,IV,IC,NL,PROC
      REAL(kind(1.D0)) :: DT,Q,N(3),KG(3),P

! Parabolic guide with varying lamellae lengths, traced from the wider side
      call PGUIDE_MOVE_ENTRY(OBJ)
!  check passage through the entry
      IF (.NOT.BENDER_PASS(OBJ,NEUT%R,.false.,IH,IV)) GOTO 200
      NL=0
! iterate through reflections
      do while (NEUT%STATE.eq.0)
        NL=NL+1
        if (OBJ%TYP==4) then
          CALL EGUIDE_CON(OBJ,NEUT%R,NEUT%K,IH,IV,IC,DT,Q,N,KG,P)
        else
          CALL PGUIDE_CON(OBJ,NEUT%R,NEUT%K,IH,IV,IC,DT,Q,N,KG,P)
        endif
        IF (Q.LT.0.D0) GOTO 200 ! signal from EGUIDE_CON that event is stopped
  ! move to the next node
        call TransportG(DT)
  ! check for free passage
        IF (IC.LE.0) goto 100
  ! handle interaction with mirror
        call GUIDE_MIRROR_ACTION(OBJ,IC,0.5D0*Q,NEUT%K0,NEUT%S(1),PROC)
        select case(PROC)
      ! reflection
        case(0)
          NEUT%K=KG
          NEUT%P=NEUT%P*P
          if (OBJ%MLOG>0)  call MLOGS_REC_TMP(OBJ%MLOG,IC,NEUT%R,NEUT%K,Q*Q2Ni)
          call AddEventLog(OBJ%FRAME,NEUT)
      ! escape
        case(1)
          NEUT%STATE=1
      ! absorption in coating
        case(2)
          NEUT%STATE=2
      ! scattering in coating
        case(4)
          NEUT%STATE=4
        end select
      enddo

100   CONTINUE
      return
! stop event
200   NEUT%P=0.D0
      NEUT%STATE=1
      END SUBROUTINE FGUIDE_GO

!----------------------------------------------------
      SUBROUTINE BENDER_GO(OBJ)
! GO procedure for a guide with flat walls (TYP=1)
! INPUT:   NEUT at the entry in local coordinates
! RETURN:  NEUT at the exit in local coordinates
! On exit, NEUT%STATE must have one of these values:
! STATE=0 .. passed
! STATE=1 .. escaped by missing the guide, transition to substrate, etc.
! STATE=2 .. absorbed in the coating
! NOTE: absorption in coating is considered only if lamellas are translucent
! or if the MONITOR option is on.
! Set NEUT%P=0 if you wish to stop the event (=> no escape monitoring).
!----------------------------------------------------
      TYPE(GUIDE) :: OBJ
      REAL(kind(1.D0)),parameter :: Q2Ni=0.5D0/TWOPI/GammaNi
      INTEGER ::  BENDER_INLAM
      INTEGER :: IH,IV,IX,I,IC,IN_LAM,ILOOP,PROC
      REAL(KIND(1.D0)) :: DT,DDT,Q,N(3),KG(3),pd,MU2PI,P,QT
      logical :: TRANS,CAN_TRANSIT
      REAL(KIND(1.D0)) :: RN
      integer :: ILAM(2),NLAM(2),IE

      integer :: ISTATE=1 ! neutron position: (1) in the gap (2) inside lamella
      character*7, parameter :: CSTATE(2)=(/'channel','lamella'/)

! Set numerical precission:
! NOTE: longitudinal position may shift at the reflection point by < EPSD/EPSA !
! However, trajectory should still be sufficiently accurate
    ! distance from lamella surface [mm]
      REAL*8,parameter :: EPSD=1.D-5
    ! angle w.r.t. mirror surface (lower angles may be rejected) [rad]
      REAL*8,parameter :: EPSA=1.D-6

! depending on IC as an index, the following flags determine:
  ! the side of a lamella, IL: left and top (1) , right and bottom (0)
      integer,parameter :: IL(4)=(/1,0,1,0/)
  ! and its normal direction, ID: horizontal(1), vertical(2)
      integer,parameter :: ID(4)=(/1,1,2,2/)
! table indexes and m=-values for the 4 surfaces
      integer :: ITAB(4)
      REAL(KIND(1.D0)) :: MTAB(4)
1     format(a,': ',9(1x,G13.6))

      !GUIDES_dbg = (trim(OBJ%FRAME%ID).eq.'GCA1'.and.(NINT(trace_tot)==5))
      if (GUIDES_dbg) GUIDES_idbg=GUIDES_idbg+1
      !GUIDES_dbg=(NINT(trace_tot)==5)
      ITAB=(/OBJ%NHLU,OBJ%NHRU,OBJ%NVT,OBJ%NVB/)
      MTAB=(/OBJ%GHLU,OBJ%GHRU,OBJ%GVT,OBJ%GVB/)/GammaNi
      NLAM=(/OBJ%NLH,OBJ%NLV/)

      if (GUIDES_dbg) write(*,1) 'GDE_GO '//trim(OBJ%FRAME%ID),NEUT%R,NEUT%K,NEUT%T
! check passage through the entry
! if inner lamella is hit, IN_LAM will be > 0
! if we are outside the guide profile, IN_LAM < 0
      IN_LAM=BENDER_INLAM(OBJ,NEUT%R,.false.,IH,IV)
! start inside a lamella?
      TRANS=(IN_LAM.GT.0)

! filter neutrons which may enter the guide
    ! entry outside the front window
      if (IN_LAM<0) then
        NEUT%STATE=1
        NEUT%P=0.D0 ! no continuation, no event monitoring
    ! entry through lamellas crossing (assumed to be absorbing)
      else if (IN_LAM>2) then
        NEUT%STATE=1
    ! entry through a lamella front face
      else if (TRANS) then
        if ((OBJ%MATER==0).or.(OBJ%TRFRONT==0)) NEUT%STATE=1
      endif
      if (NEUT%STATE>0) goto 100

! IC determines wall: left (1) right (2) top (3) or bottom (4)
      IC=1
      IF (TRANS) IC=2*IN_LAM ! 2 or 4, lamella can be either RIGHT or BOTTOM only

! absorption coefficient in cm^-1.A^-1)
! transmission probability = EXP(-MU2PI*DT)
      select case (OBJ%MATER)
      case(1)
        MU2PI=OBJ%MU*TWOPI
      case(2)
        MU2PI=GET_CR_ABS(MU_Si,NEUT%K0)
      case(3)
        MU2PI=GET_CR_ABS(MU_Al2O3,NEUT%K0)
      case default
        MU2PI=1.D10
      end select

      DDT=0.D0
      ILOOP=0
      IX=ID(IC)
      ILAM=(/IH,IV/)
      ISTATE=1
      IF (TRANS) ISTATE=2
!  iterate through the guide
! --------------------------------------------------------------
      do while (IC.GT.0)
    ! NOTE the difference:
    ! GUIDE_TRANS takes IC for INPUT and OUTPUT
    ! GUIDE_CON takes IC only for OUTPUT

       if (GUIDES_dbg)  write(*,1) 'GDE_GO TRANS ',TRANS,IC,IH,IV,CSTATE(istate)

! calculate next touch point and related info
        if (TRANS) then
          if (GUIDES_dbg)  write(*,1) 'GDE_GO STATE',NEUT%STATE
          call GUIDE_TRANS(OBJ,NEUT%R,NEUT%K,IH,IV,IC,DT,Q,N,KG,P)
          ! absorption in the lamella (MU is in mm^-1.A^-1)
          ! transmission probability = EXP(-OBJ%MU*2*PI*(DT+DDT))
          ! sample penetration depth
          RN=1.D0*RAN1()
          pd=-LOG(RN)/MU2PI
          if (GUIDES_dbg)  write(*,1) '    IC,STATE,P',IC,NEUT%STATE,P
          if (GUIDES_dbg)  write(*,1) '    N,KG',N,KG
          if (GUIDES_dbg)  write(*,1) '    RN,PD,DT,DDT',RN,PD,DT,DDT

          if (pd<DT+DDT) then
            DT=pd-DDT
            NEUT%STATE=1
          endif
        else
          CALL GUIDE_CON(OBJ,NEUT%R,NEUT%K,IH,IV,IC,DT,Q,N,KG,P)
        endif
        QT=0.5D0*Q
        if (GUIDES_dbg)  write(*,1) '    STATE,IC,DT,Q',NEUT%STATE,IC,DT,Q*Q2Ni
! move neutron to the next node
        call TransportG(DT)

        if (GUIDES_dbg)  write(*,1) '    NEXT R, KG, N',NEUT%R, KG, N

! check for free passage
        IF (IC.LE.0) goto 100
! reject, too low reflection angle
        IF (Q.LE.EPSA*NEUT%K0) NEUT%STATE=1
! if absorbed => finish
        if (NEUT%STATE>0) goto 100
! record event
        call AddEventLog(OBJ%FRAME,NEUT)
        IX=ID(IC)
        DDT=0.D0
! handle interaction at the touch point
        call GUIDE_MIRROR_ACTION(OBJ,IC,QT,NEUT%K0,NEUT%S(1),PROC)

        if (GUIDES_dbg)  write(*,1) '    PROC ',PROC

! reflection:
        IF (PROC==0) THEN
          ! Set K=K+Q
          NEUT%K=KG
          NEUT%P=NEUT%P*P
          if (OBJ%MLOG>0) call MLOGS_REC_TMP(OBJ%MLOG,IC,NEUT%R,NEUT%K,Q*Q2Ni)
          goto 10 ! continue
        endif
! if there is no reflection, decide about transition into next sector (lamella or channel)
    ! can transit only to free space or material with defined absorption
        CAN_TRANSIT=(TRANS.or.(OBJ%MATER>0))
    ! avoid transition into outer walls, only reflection is then permitted
      ! IE = index of the interacting lamella
        IE=ILAM(IX)+IL(IC)
        CAN_TRANSIT=(CAN_TRANSIT.and.(IE.GT.0).AND.(IE.LT.NLAM(IX)))
! absorption or scattering in coating
        if ((PROC==2).or.(PROC==4)) then
          CAN_TRANSIT=.false.
          NEUT%STATE=PROC
          goto 100 ! no reflection, absorb
! if transmission is not allowed, escape
        else IF (.NOT.CAN_TRANSIT) THEN
          NEUT%STATE=1
          goto 100 ! no reflection, absorb
        ELSE
! else transmitted to the next channel or lamella
          ! in the lamella: for right/bottom wall, decrease IH/IV by 1
          ! in the channel: for left/top wall, increase IH/IV by 1
          I=IL(IC)
          IF (TRANS) I=IL(IC)-1
          select case (IX)
          case(1)
            IH=IH+I
          case(2)
            IV=IV+I
          end select
          ILAM=(/IH,IV/)
          TRANS=(.NOT.TRANS)
          ISTATE=1
          IF (TRANS) ISTATE=2
          if (GUIDES_dbg)  write(*,1) '    TRANSMITTED ',CAN_TRANSIT, TRANS, IC, IH

        ENDIF
  ! move a little away from the surface (avoid problems with num. precission)
10      if (Q>0.D0) then
          DDT=EPSD/Q
          call TransportG(DDT)
          if (GUIDES_dbg)  write(*,1) '    DDT',DDT,EPSD,Q
        endif
        ILOOP=ILOOP+1
        IF (ILOOP.GT.100000) goto 999  ! debug: prevent infinite loops
      ENDDO
!  end of iteration
! --------------------------------------------------------------

100   CONTINUE
      return

! too many iterations -> BUG?
999   write(*,*) OBJ%FRAME%NAME,OBJ%FRAME%COUNT
      write(*,*) 'too many iterations through guide, z=',NEUT%R(3)
      write(*,*) 'program exits'
      END SUBROUTINE BENDER_GO

!---------------------------------------------------
      subroutine GUIDE_MIRROR_ACTION(OBJ,IC,QT,K0,S,PROC)
! handle neutron interaction at the mirror surface
! IC = wall index, (1-4)=(left,right,top,bottom)
! QT = 0.5 of the scattering vector = 2PI/lambda*theta
! K0 = 2PI/lambda
! S  = spin projection
! RETURN:
! PROC=0 ... reflected
! PROC=1 ... transmitted under coating
! PROC=2 ... absorbed in the coating
! PROC=4 ... scattered in the coating
!---------------------------------------------------
      TYPE(GUIDE) :: OBJ
      integer,intent(in) :: IC
      real(kind(1.D0)),intent(in) :: QT,K0,S
      integer,intent(out) :: PROC
      real(kind(1.D0)) :: MC
      integer :: NR
      call GUIDE_GET_MTAB(IC,OBJ,NR,MC)
      call MIRROR_REF_EX(NR,OBJ%MONITOR.or.(OBJ%MATER>0),QT,S,MC,K0,PROC)
      end subroutine GUIDE_MIRROR_ACTION

      end MODULE GUIDES_TRACE
