!//////////////////////////////////////////////////////////////////////
!////  $Id $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2012, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.38 $
!////     $Date: 2019/08/16 17:16:25 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Ray-tracing subroutines for SGUIDE
!////
!////////////////////////////////////////////////////////////////////////
      MODULE SGUIDES_TRACE
      use SIMATH
      use TRACELIB
      USE TRACINGDATA
      USE SGUIDES
      use EVENTLOGS
      use FRAMES_TRACE
      use MIRROR_TABLE
      use RNDGEN
      use mtmod
      use NSTORE
      implicit none
      private

      logical :: dbg=.false.
      integer :: dbgcnt=0
      REAL(kind(1.D0)),parameter :: EPS_DX=1.D-7 ! tolerance for wall crossing [mm]
      REAL(kind(1.D0)),parameter :: PMIN=1.D-5   ! threshold for reflection probability
      logical :: TRCONT



      type(NEUTRON) :: TMPN(SGUIDES_DIM)

      public SGUIDE_GO,SGUIDE_ADJUST,SGUIDE_INIT,SGUIDE_ENDTRACE,SGUIDE_CLR

      contains


!-------------------------------------------------------------
      logical function SGUIDE_GO(INST)
! call GO on given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
        SGUIDE_GO=SGUIDE_GO_INST(ASGUIDES(INST)%X,INST)
      end function SGUIDE_GO


!------------------------------------
      SUBROUTINE SGUIDE_INIT(INST)
!------------------------------------
      integer,intent(in) :: INST
      TYPE(SGUIDE),POINTER :: OBJ
      if (.not.SGUIDE_isValid(INST)) return
      OBJ => ASGUIDES(INST)%X
      CALL FRAME_INIT(OBJ%FRAME)
      ! deselect logger = it will be reallocated later, before simulation starts (TRACING)
      OBJ%MLOG=0
      OBJ%DX(1:2,0:OBJ%NSEG)=0.D0
      OBJ%DA(1:2,1:OBJ%NSEG)=0.D0
      END SUBROUTINE SGUIDE_INIT


!-------------------------------------------------------------
      SUBROUTINE SGUIDE_CLR(INST)
! return reference to given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(SGUIDE),pointer :: OBJ
        if (.not.SGUIDE_isValid(INST)) return
        OBJ => ASGUIDES(INST)%X
        OBJ%FRAME%COUNT=0
        call NSTORE_CLEAR(OBJ%FRAME%REGISTRY)
        ! clear bounce logger)
        call MLOGS_REC_CLR(OBJ%MLOG)
      end SUBROUTINE SGUIDE_CLR

!--------------------------------------------------------------------------
      subroutine SGUIDE_ENDTRACE(INST,CNT,P)
! Add previously registered event to the monitor storage, using given weight
!--------------------------------------------------------------------------
      integer,intent(in) :: INST,CNT
      REAL(kind(1.D0)),intent(in)  :: P
      TYPE(PSGUIDE) :: OBJ
        OBJ%X => ASGUIDES(INST)%X
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
      end subroutine SGUIDE_ENDTRACE


!-------------------------------------------------
      SUBROUTINE SGUIDE_ADJUST(INST,IERR)
! Update position of the component using nominal trajectory represented by NEUT
! Call FRAME_INIT_MAT first.
! Then use SLIT_PRE & SLIT_POST procedures to propagate NEUT
!-------------------------------------------------
      integer,intent(in) :: INST
      integer,intent(out) :: IERR
      integer :: i,is
      TYPE(SGUIDE),POINTER :: OBJ
1     format(a,': ',7(1x,G12.6))
      if (.not.SGUIDE_isValid(INST)) return
      OBJ => ASGUIDES(INST)%X
      ! call SGUIDE_ADJUST_SEGMENTS(INST)
      call SGUIDE_SET_PROFILE(INST,1)
      call SGUIDE_SET_PROFILE(INST,2)
      call SGUIDE_CALC_PARAM(INST)
      call FRAME_INIT_MAT(OBJ%FRAME)
      call SLIT_PRE(OBJ%FRAME)
      call SLIT_POST(OBJ%FRAME)
      dbgcnt=0
      is=0
      do i=1, 4
        if (OBJ%MONITOR(i)) is = is+1
      enddo
      if (is>0) then
        write(*,1) trim(OBJ%FRAME%ID)//' has monitors: ',is
      endif

      !if (OBJ%dbg_NSUM(3)>0) then
      !  write(*,1) 'Events sum, '//trim(OBJ%FRAME%ID)
      !  WRITE(*,1) 'reflected: ', OBJ%dbg_SUM(1), OBJ%dbg_NSUM(1)
      !  WRITE(*,1) 'transmitted: ', OBJ%dbg_SUM(2), OBJ%dbg_NSUM(2)
      !  WRITE(*,1) 'captured: ', OBJ%dbg_SUM(3), OBJ%dbg_NSUM(3)
      !  WRITE(*,1) 'scattered: ', OBJ%dbg_SUM(4), OBJ%dbg_NSUM(4)
      !endIF
      OBJ%dbg_SUM=0.D0
      OBJ%dbg_NSUM=0
      end SUBROUTINE SGUIDE_ADJUST

!---------------------------------------------------------------------------
      SUBROUTINE SEG_ENTRY(OBJ,ISEG)
! move to the segment entry point
!-----------------------------------------------------------------------------
      TYPE(SGUIDE) :: OBJ
      integer,intent(in) :: ISEG
      REAL(KIND(1.D0)) :: Z0,DT
! trace up-stream: move to z=LENGTH
      if (TRACING_UP) then
        Z0=OBJ%WIN(3,ISEG)
! trace down-stream: move to z=0
      else
        Z0=OBJ%WIN(3,ISEG-1)
      endif
      if (abs(NEUT%K(3)).gt.1.D-20) then
        DT=(Z0-NEUT%R(3))/NEUT%K(3)
      else
        DT=1.D30
        NEUT%P=0.D0
      endif
      call TransportG(DT)
      END subroutine SEG_ENTRY

!----------------------------------------------------
      LOGICAL*4 FUNCTION SGUIDE_SEG_INSIDE(OBJ,ISEG)
! check that neutron is within entry window of the segment
!----------------------------------------------------
      TYPE(SGUIDE) :: OBJ
      integer,intent(in) :: ISEG
      REAL(kind(1.D0)) :: SZ(2),R0(2)
      logical :: INS
      integer :: ISG
1     format(a,7(1x,G12.5))
      if (TRACING_UP) then
      ! start at the end of the segment
        ISG = ISEG
      else
      ! start at the entry of the segment
        ISG = ISEG-1
      endif
      ! position wrt entry window center
      if (ISG==0) then
        SZ=obj%frame%size(1:2)
      ELSE
        SZ=(/OBJ%W(ISG),OBJ%H(ISG)/)
      ENDIF
      R0(1)=NEUT%R(1)-OBJ%WIN(1,ISG)-OBJ%DX(1,ISG)
      R0(2)=NEUT%R(2)-OBJ%WIN(2,ISG)-OBJ%DX(2,ISG)
      if (dbg) write(*,1) 'SGUIDE_SEG_INSIDE ISG, SZ, R0: ',ISG, SZ, R0
      INS=TLIB_INSIDE(R0(1:2),SZ(1:2),OBJ%FRAME%SHAPE)
      SGUIDE_SEG_INSIDE = INS
      end FUNCTION SGUIDE_SEG_INSIDE


!----------------------------------------------------
      LOGICAL*4 FUNCTION SGUIDE_GO_INST(OBJ,INST)
!----------------------------------------------------
      TYPE(SGUIDE) :: OBJ
      integer,intent(in) :: INST
      integer :: DI,NL,ISEG,JC,i
      logical :: LOG1,INS

1     format(a,8(1x,G14.7))
      JC=0
      !dbgcnt=dbgcnt+1
      !dbg = ((dbgcnt<1050).and.(dbgcnt>1000))
      ! dbg = (trim(OBJ%FRAME%ID).eq.'GEX1'.and.(NINT(trace_tot)==5))
! Convert to local coordinate system
      CALL FRAME_PRE(OBJ%FRAME)
! move to the guide entry window
      call FRAME_ENTRY(OBJ%FRAME)
      ! call SGUIDE_PRE(OBJ)
! log event at the entry of each object
      call AddEventLog(OBJ%FRAME,NEUT)
      if (dbg) write(*,1) 'SGUIDE_GO_INST '//trim(OBJ%FRAME%ID),NEUT%R,NEUT%K

! generate segments misalignment if required
      call SGUIDE_MISALIGN(OBJ)

!............................
      if (TRACING_UP) then
    ! start at the end of the last segment
        NL=OBJ%NSEG
        DI=-1
        !SZ=(/OBJ%W(NL),OBJ%H(NL)/)
        ! position wrt entry window center
        !R0(1)=NEUT%R(1)-OBJ%WIN(1,NL)-OBJ%DX(1,NL)
        !R0(2)=NEUT%R(2)-OBJ%WIN(2,NL)-OBJ%DX(2,NL)
      else
    ! start at the entry of the 1st segment
        NL=1
        DI=1
        !SZ=OBJ%FRAME%SIZE(1:2)
        ! position wrt entry window center
        !R0(1)=NEUT%R(1)-OBJ%WIN(1,0)-OBJ%DX(1,0)
        !R0(2)=NEUT%R(2)-OBJ%WIN(2,0)-OBJ%DX(2,0)
      endif


      if (dbg) then
      !if (dbgcnt==1001) then
        write(*,1) 'SGUIDE_GO_INST list, NL=',NL
        do i=1,OBJ%NSEG
          write(*,1) '   ',i, OBJ%WIN(3,i), OBJ%W(i), OBJ%H(i)
        enddo
      !endif
      endif

      INS = SGUIDE_SEG_INSIDE(OBJ,NL)

    ! check passage through the 1st window
      if (OBJ%FRONTMASK) then
        !INS=TLIB_INSIDE(R0(1:2),SZ(1:2),OBJ%FRAME%SHAPE)
        if (dbg) write(*,1) 'SGUIDE_GO_INST INSIDE',INS
! missed front area ?
        if (.not.INS) then
          NEUT%STATE=1
          goto 300
        endif
      endif
! trace through all guide segments
      ISEG=NL
      TRCONT=.true.
      if (dbg) write(*,1) 'SGUIDE_GO_INST ISEG ',ISEG,NEUT%STATE
! reset bounce recorder
      call MLOGS_RESET(OBJ%MLOG)
      do while (TRCONT)
          ! LOG1=.true.
          call SEG_ENTRY(OBJ,ISEG)
          LOG1 = SGUIDE_SEG_INSIDE(OBJ,ISEG)
          if (.not.LOG1) then
            ! stoped: missed segment entry
            NEUT%STATE=1
            if (dbg) write(*,1) 'SGUIDE_GO_INST MISSED SEG: ',ISEG,NEUT%R(1:3), ' z0=',OBJ%WIN(3,ISEG-1)
          else
            ! trace reflection within one segment
            do while (LOG1)
              LOG1=SGUIDE_CROSS(OBJ,ISEG,JC)
              if (dbg) write(*,1) 'SGUIDE_GO_INST JC ',JC,NEUT%R(3),LOG1
              if (LOG1) call AddEventLog(OBJ%FRAME,NEUT)
            enddo
            if (dbg) write(*,1) 'SGUIDE_GO_INST ISEG ',ISEG,NEUT%R(3),NEUT%P
            ISEG=ISEG+DI
          endif
          INS=((ISEG>0).and.(ISEG<=OBJ%NSEG))
          TRCONT=(TRCONT.and.INS.and.(NEUT%STATE==0))
      enddo

!--------------------------
! Check whether the neutron should be stopped
300   IF (NEUT%STATE>0) then
! Absorbed: stop unless monitoring of absorbed events is required
            LOG1 = (JC<=0) ! missed front window area ?
            if (JC>0) LOG1 = (.not.OBJ%MONITOR(JC)) ! monitoring not requested
            LOG1 = (LOG1.or.(NEUT%P<=0.D0)) ! absorbed neutron - zero probability
            if (LOG1) then
                  if (dbg) write(*,1) 'SGUIDE_GO_INST STOPPED',JC,OBJ%MONITOR(JC),NEUT%P
            endif
      else
! Transmitted: stop if guide is closed
            LOG1 = OBJ%CLOSED
      endif
      if (LOG1) then
            NEUT%P=0.D0
            SGUIDE_GO_INST=.FALSE.
            return
      endif
      ! label neutron with its state number for monitoring
      NEUT%LABEL=NEUT%STATE
! passed
      if (dbg) write(*,1) 'SGUIDE_GO_INST PASSED '
  ! log event at the exit of a guide
      call AddEventLog(OBJ%FRAME,NEUT)
  ! register event
      TMPN(INST)=NEUT
      if (TRACING_UP) TMPN(INST)%K=-TMPN(INST)%K
  ! post-transformation to exit coordinates, correct for deflection
      CALL FRAME_POST(OBJ%FRAME)
      call incCounter(OBJ%FRAME%COUNT)
      SGUIDE_GO_INST=.TRUE.
      END FUNCTION SGUIDE_GO_INST


!----------------------------------------------------
      subroutine SGUIDE_GOTOEND(OBJ)
! transport to the guide end, free flight, no walls
!----------------------------------------------------
      TYPE(SGUIDE) :: OBJ
      REAL(kind(1.D0)) :: z0,kz,zout,tout
      z0=NEUT%R(3)
      kz=NEUT%K(3)
    ! reverse tracing: exchange exit and entry
      if (TRACING_UP) then
        zout=0.D0   ! position of the guide entry
      else
        zout=OBJ%WIN(3,OBJ%NSEG)     ! position of the guide exit
      endif
      tout=(zout-z0)/kz   ! time to the segment exit
      call TransportG(tout)
      end subroutine SGUIDE_GOTOEND


!----------------------------------------------------
      LOGICAL FUNCTION SGUIDE_CROSS(OBJ,ISEG,JC)
! trace neutron to the cross-section with any wall
! or to the exit of the segment ISEG
! return : true if there is a contact with a wall.
! If there is a contact, perform corresponding reflection
! with reflection probability
! return reflectin wall index
!----------------------------------------------------
      TYPE(SGUIDE) :: OBJ
      integer,intent(in) :: ISEG
      integer,intent(out) :: JC
      integer,parameter :: SIDE(4)=(/1,1,2,2/)
      !integer,parameter :: SIDEINV(4)=(/2,2,1,1/)
      ! WSIDE contains indices of perpendicular walls corresponding
      ! to give wall index
      integer,parameter :: WSIDE(8)=(/3,4,3,4,1,2,1,2/)
      integer,parameter :: SG(4)=(/1,-1,1,-1/)
      integer,parameter :: SPOL(4)=(/2,2,1,1/)
      integer,parameter :: IXY(4)=(/1,2,1,2/)
      integer,parameter :: II(0:4) = (/1, 2, 3, -1, 4/)
      integer :: j,k,k1,k2,NR,IC1,IC2,PROC
      REAL(kind(1.D0)) :: gx,ang,Q,QT,N(3),d(3),NW(3),KG(3)
      REAL(kind(1.D0)) :: z,z0,zin,zout,kz,z0seg
      REAL(kind(1.D0)) :: t1,t2,tc,tin,tout
      REAL(kind(1.D0)) :: x,kx,H,P,dx,dkx
      REAL(kind(1.D0)),parameter :: Q2Ni=0.5D0/TWOPI/GammaNi
      logical :: touch
1     format(a,7(1x,G14.7))

      z0=NEUT%R(3)
      kz=NEUT%K(3)

    ! reverse tracing: exchange exit and entry
      if (TRACING_UP) then
        zin=OBJ%WIN(3,ISEG)      ! position of the segment entry
        zout=OBJ%WIN(3,ISEG-1)   ! position of the segment exit
      else
        zin=OBJ%WIN(3,ISEG-1)    ! position of the segment entry
        zout=OBJ%WIN(3,ISEG)     ! position of the segment exit
      endif
      z0seg=OBJ%WIN(3,ISEG-1)    ! position of the segment entry in forward direction
      tin=(z0-zin)/kz     ! time from the segment entry
      tout=(zout-z0)/kz   ! time to the segment exit

      TC=TLIB_INF
      JC=0
      touch=.false.
      ! find cross-sections with all active walls and select the nearest one
      do j=1,4
      if (OBJ%ACTIVE(j)) then
        d(1:3)=OBJ%D(2:4,j,ISEG)
        k=SIDE(j)
        dx=OBJ%DX(k,ISEG-1)+OBJ%DA(k,ISEG)*(z0-z0seg)
        dkx=kz*OBJ%DA(k,ISEG)
        x=NEUT%R(k)-dx
        kx=NEUT%K(k)-dkx
        gx=NEUT_GR(k)
      ! for smoothed model or gravity, solve the quadratic equation
        if (OBJ%SMOOTH(j).or.abs(gx)>0.D0) then
          call TLIB_CROSS(kx,kz,x,z0-z0seg,d,gx,t1,t2)
          if (t1<=TLIB_EPS) t1=t2
        else
      ! flat surface = one solution
          dkx=kx-d(2)*kz
          if (abs(dkx)>TLIB_EPS) then
            t1=(d(3)-x+d(2)*(z0-z0seg))/dkx
          else
            t1=TLIB_INF
          endif
        endif
      ! take the minimum positive solution
        if ((t1>TLIB_EPS).and.(t1<TC)) then
          TC=t1
          JC=j
        endif
      endif
      enddo
    ! is there a cross-section ?
    ! touch condition - contact must be within segment end points
      if (JC>0) touch=(TC<tout)
    ! free passage = transport to the end of the segment

      if (dbg) write(*,1) 'SGUIDE_CROSS touch ',touch,JC,TC
      if (.not.touch) then
        call TransportG(tout-TLIB_EPS)
        goto 10
      endif
    ! transport to the next touch point
      call TransportG(TC-TLIB_EPS)
      if (dbg) write(*,1) 'SGUIDE_CROSS touched at ',NEUT%R, NEUT%K
    ! check the top and bottom segment boundaries.
    ! allow free passage if the top/bottom walls are not active
        k=SIDE(jc)
        !iy=SIDEINV(jc)
        z0=NEUT%R(3)-z0seg
        IC1=WSIDE(2*(jc-1)+1)
        IC2=WSIDE(2*(jc-1)+2)
        k1=SIDE(IC1)
        if (.not.OBJ%ACTIVE(IC1)) then
          H=SGUIDE_POS(0,OBJ%D(1:4,IC1,ISEG),z0)
          H=H+OBJ%DX(k1,iseg-1)+OBJ%DA(k1,iseg)*z0
          touch=(touch.and.(NEUT%R(k1)<H))
        endif
        k2=SIDE(IC2)
        if (.not.OBJ%ACTIVE(IC2)) then
          H=SGUIDE_POS(0,OBJ%D(1:4,IC2,ISEG),z0)
          H=H+OBJ%DX(k2,iseg-1)+OBJ%DA(k2,iseg)*z0
          touch=(touch.and.(NEUT%R(k2)>H))
        endif
        if (.not.touch) goto 10
    ! check for gaps
        z=NEUT%R(3)
        if ((abs(z-zin)<0.5D0*OBJ%GAP).or.(abs(z-zout)<0.5D0*OBJ%GAP)) then
           if (OBJ%GAPABS) then
          ! absorbing gap
            goto 30
          else
          ! allow leakage
            goto 10
          endif
        endif
    ! hit active wall - try reflection
      ! perform reflection
      ! get 1st derivative of the surface
        ANG=SGUIDE_POS(1,OBJ%D(1:4,jc,ISEG),z0)
        ANG=ANG+OBJ%DA(k,iseg)
        ! get surface normal unit vector
        N=0.D0
        N(3)=SG(jc)*sin(ANG)
        N(k)=-SG(jc)*SQRT(1.D0-N(3)**2)

        !mirror_dbg = dbg
        !rndgen_dbg = dbg
        if (dbg) write(*,1) 'SGUIDE_CROSS ngen1 ',mtmod_ngen
        call MIRROR_REFLECT(OBJ%WAV,NEUT%K,N,jc,KG,NW,Q,P)
        if (dbg) write(*,1) 'SGUIDE_CROSS ngen2 ',mtmod_ngen
        if (dbg) write(*,1) 'SGUIDE_CROSS K,  N ',NEUT%K, N
        if (dbg) write(*,1) 'SGUIDE_CROSS KG, NW',KG, NW
        if (dbg) write(*,1) 'SGUIDE_CROSS Q, P',Q, P
        if (P<=0.D0) then
          goto 40
        endif
        !write(*,1) 'SGUIDE_CROSS P=',Q , P, NEUT%P*P
        NR=OBJ%NR(jc,ISEG)
        QT=0.5D0*Q
      ! handle interaction with mirror
        call MIRROR_REF_EX(NR,OBJ%MONITOR(jc),QT,NEUT%S(SPOL(jc)),OBJ%MC(jc,ISEG),NEUT%K0,PROC)
        select case(PROC)
      ! reflection
        case(0)
          NEUT%K=KG
          NEUT%P=NEUT%P*P
          if (OBJ%MLOG>0)  call MLOGS_REC_TMP(OBJ%MLOG,jc,NEUT%R,NEUT%K,Q*Q2Ni)
          goto 20
      ! escape=transmitted
        case(1)
          goto 30
      ! absorption in coating
        case(2)
          goto 40
        case(4)
          goto 50
        end select

! free passage
10    SGUIDE_CROSS=.false.
      if (dbg) write(*,1) 'SGUIDE_CROSS free '
      return
! reflect and continue
20    SGUIDE_CROSS=.true.
      if (dbg) write(*,1) 'SGUIDE_CROSS reflect ',NEUT%K
      return
! escape=transmitted
30    NEUT%STATE=1 ! escape
      if (dbg) write(*,1) 'SGUIDE_CROSS escape '
      SGUIDE_CROSS=.false.
      return
! absorb in coating
40    NEUT%STATE=2 ! absorb
      if (dbg) write(*,1) 'SGUIDE_CROSS absorb '
      SGUIDE_CROSS=.false.
      return
! scatter in coating
50    NEUT%STATE=4 ! scatter
      if (dbg) write(*,1) 'SGUIDE_CROSS scatter '
      SGUIDE_CROSS=.false.
      end FUNCTION SGUIDE_CROSS

      end MODULE SGUIDES_TRACE
