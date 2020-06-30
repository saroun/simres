!//////////////////////////////////////////////////////////////////////
!////  $Id: xtals_trace.f,v 1.30 2019/08/16 17:16:26 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2009, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.30 $
!////     $Date: 2019/08/16 17:16:26 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Describes component class: XTAL
!////  implements ray-tracing algorithm
!////////////////////////////////////////////////////////////////////////
      MODULE XTALS_TRACE
      USE CONSTANTS
      use IO
      use XTALS
      use XTALS_REF
      use TRACINGDATA
      use EVENTLOGS
      use XMLINFO
      use ARRAY3D
      use FRAMES_TRACE
      use RNDGEN
      use NSTORE
      implicit none
      private
      save

      type(NEUTRON) :: TMPN(XTALS_DIM),DBGN

      logical :: DBG=.false.
      logical :: myDBG=.false.
      INTEGER,parameter :: DBGREC=0
      INTEGER :: DBGIT=0
      integer,parameter :: MAXITER=100
      integer :: ITER
      real(kind(1.D0)) :: XTALS_TRACE_REC(6,MAXITER),TMP_TRACE_REC(6,MAXITER)
      public XTAL_GO,PRINT_REC,XTAL_FINALIZE,XTAL_ADJUST,XTAL_INIT,XTAL_REPORT_ABS,XTAL_ENDTRACE

      contains

!-------------------------------------------------------------
      logical function XTAL_GO(INST)
! call GO on given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
      integer :: i
      i=AXTALS(INST)%X%FRAME%COUNT
       ! XTAL_GO=.false.
      !  XTAL_GO=XTAL_GO_INST(AXTALS(INST)%X)
        !myDBG=(i<10)
        !if (i>498) then
        !  NEUT_DBG=.true.
        !  myDBG=NEUT_DBG
        !endif
        XTAL_GO=XTAL_ARRAY_GO(AXTALS(INST)%X,INST)
        if (myDBG) write(*,*) '   XTAL_GO ',XTAL_GO,i

      end function XTAL_GO

!--------------------------------------------------------------------------
      subroutine XTAL_ENDTRACE(INST,CNT,P)
! Add previously registered event to the monitor storage, using given weight
! update internal oger of wall bounces
!--------------------------------------------------------------------------
      integer,intent(in) :: INST,CNT
      REAL(kind(1.D0)),intent(in)  :: P
      TYPE(PXTAL) :: OBJ
        OBJ%X => AXTALS(INST)%X
        if (OBJ%X%FRAME%REGISTRY>0) then
        if (TMPN(INST)%P>0.D0) then
          ! transform to incident axis coordinates
          call FRAME_LOCAL(-1,OBJ%X%FRAME,TMPN(INST))
          TMPN(INST)%P=P     ! set the final weight to the event
          TMPN(INST)%CNT=CNT ! must set counter to neutron before calling NSTORE_SETEVENT !!
          call NSTORE_SETEVENT(OBJ%X%FRAME%REGISTRY,TMPN(INST))
        endif
        endif
      end subroutine XTAL_ENDTRACE

! ---------------------------------------
      SUBROUTINE XTAL_FINALIZE(INST)
! finalization of Xtal class
! called after instrument tracing
! ---------------------------------------
      integer,intent(in) :: INST
      TYPE(XTAL),POINTER :: OBJ
      character(32) :: CNUM
      if (.not.XTAL_isValid(INST)) return
      OBJ => AXTALS(INST)%X
    !  write(*,*) 'XTAL_FINALIZE DBGREC,ICR=',DBGREC,OBJ%ICR,OBJ%MAX_ITER_COUNT
      if (OBJ%MAX_ITER_COUNT.gt.0) then
        call INT2STR(OBJ%MAX_ITER_COUNT,CNUM)
        call MSG_WARN('Iteration limit reached in '//trim(OBJ%FRAME%ID)//', '//trim(CNUM)//' times',1)
        call PRINT_REC(OBJ)
      endif
      END SUBROUTINE XTAL_FINALIZE

!-------------------------------------------------------------
      SUBROUTINE PRINT_REC(CR)
!-------------------------------------------------------------
      TYPE (XTAL),intent(in) :: CR
      integer :: j
      character(128) CNUM1,CNUM2
      if (DBGREC.eq.CR%ICR) then
        write(*,*) 'Iteration limit reached in '//trim(CR%FRAME%ID)
        do j=1,MAXITER
          call ARRAY2STR(XTALS_TRACE_REC(5,j),3,CNUM1)
          call ARRAY2STR(XTALS_TRACE_REC(6,j),3,CNUM2)
          write(*,"(I3,2x,a,2x,a)") j,trim(CNUM1),trim(CNUM2)
        enddo
      endif
      end SUBROUTINE PRINT_REC


!-------------------------------------------------------------------------
      subroutine XTAL_ADJUST(INST,IERR)
! Update position of the component using nominal trajectory represented by NEUT
! Call FRAME_INIT_MAT first.
! Then use SLIT_PRE & SLIT_POST procedures to propagate NEUT
! return IERR<>0 if lambda>2*dhkl
! calculates also following fields in CR:
! BMAT,DHKL,G
!--------------------------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(XTAL),POINTER :: OBJ
      ! real(kind(1.D0)),intent(in) :: LAMBDA
      integer,intent(out) :: IERR
      TYPE(TCR_ITEM) :: CDATA
      real(kind(1.D0)) :: KI(3),G(3),A(3)
      integer :: i
      if (.not.XTAL_isValid(INST)) return
      OBJ => AXTALS(INST)%X
      OBJ%LAMBDA=twopi/NEUT%K0
      ierr=2
  ! get crystal data
      if (OBJ%ITAB.le.0) return
      ierr=0
      call GET_CR_ITEM(OBJ%ITAB,CDATA)
  ! calculate B-matrix, DHKL, G
      call CalBMatrix(CDATA%LAT,OBJ%A,OBJ%B,OBJ%CHI,OBJ%PSI,OBJ%T,OBJ%BMAT)
  ! cal. coordinate transformations for actual gonio and axis angles
  ! TODO: this call is probably duplicit, already done in _PREPARE:
      call FRAME_INIT_MAT(OBJ%FRAME)
  ! adjust axis and gonio
      if (OBJ%AUTOADJ) call XTAL_ORIENT(OBJ,ierr)
      call XTAL_AUTOFOCUS(OBJ)

  ! find available primary reflections
      !NEUT%K=(/0.D0,0.D0,2.D0*PI/LAMBDA/)
      !NEUT%K0=2.D0*PI/LAMBDA
      call SLIT_PRE(OBJ%FRAME)
      KI=NEUT%K
      call CAL_PREF(OBJ,NEUT%K)
  ! report reflections (debug)
    !  call XREF_REPORT(OBJ%ICR)
      call RLNODE_REPORT(OBJ%ICR)
      ! get nominal diffracton vector in local coordinates
      do i=1,3
        A(i)=1.D0*OBJ%HKL(i)
      enddo
      call M3xV3(1,OBJ%BMAT,A,G)
      NEUT%K=KI+G
      ! convert to exit coordinates
      call SLIT_POST(OBJ%FRAME)
    ! update LAMBDA of the crystal object
    ! it will correspond to the selected reflecction
      ! LOG1=XTAL_ADJUST_LAMBDA(INST)

    !  write(*,*) 'XTAL_ADJUST LAMBDA,MI',LAMBDA,XTAL_GETMI(OBJ,NEUT%K0)
    !  write(*,*) '  SIGABS,SIGINC,SIGTDS ',OBJ%SIGABS,OBJ%SIGINC,OBJ%SIGTDS(1:3)
      ! write(*,*) 'XTAL_ADJUST ICR=',OBJ%ICR
      end subroutine XTAL_ADJUST

! ---------------------------------------
      SUBROUTINE XTAL_INIT(INST)
! initialization of Xtal class
! called before instrument tracing
! ---------------------------------------
      integer,intent(in) :: INST
      TYPE(XTAL),POINTER :: OBJ
      TYPE(TCR_ITEM) :: CDATA
      if (.not.XTAL_isValid(INST)) return
      DBGIT=0
      OBJ => AXTALS(INST)%X
    ! call inherited method
        call FRAME_INIT(OBJ%FRAME)
    ! clear counters
        OBJ%MAX_ITER_COUNT=0
    ! get crystal data
        call GET_CR_ITEM(OBJ%ITAB,CDATA)
    ! calculate beam attenuation cross-sections
    ! assumption: LAMBDA has been defined in XTAL_ADJUST
        call GET_CR_SIGMA(CDATA,OBJ%LAMBDA,OBJ%T,OBJ%SIGABS,OBJ%SIGINC,OBJ%SIGTDS)
      END SUBROUTINE XTAL_INIT


!---------------------------------------------------------
      subroutine XTAL_REPORT_ABS(INST)
!---------------------------------------------------------
      integer,intent(in) :: INST
      character(128) :: CNUM1
      TYPE(XTAL),POINTER :: OBJ
      TYPE(TCR_ITEM) :: CDATA
      REAL(KIND(1.D0)) :: Z,lam2

1     format('Attenuation parameters:')
21    format(a10,' = ',a12,'   [',a,']')

      if (.not.XTAL_isValid(INST)) return
      call XTAL_INIT(INST)
      OBJ => AXTALS(INST)%X
      call GET_CR_ITEM(OBJ%ITAB,CDATA)
      write(*,1)
      call FLOAT2STR(OBJ%LAMBDA,CNUM1)
      write(*,21) 'lambda',trim(CNUM1),'A'
      call FLOAT2STR(OBJ%T,CNUM1)
      write(*,21) 'T',trim(CNUM1),'K'
      call FLOAT2STR(XTAL_GETMI(OBJ,TWOPI/OBJ%LAMBDA)*10.D0,CNUM1)
      write(*,21) 'mu',trim(CNUM1),'1/cm'
      call FLOAT2STR(CDATA%LAT%VOL,CNUM1)
      write(*,21) 'cell vol.',trim(CNUM1),'A^3'
      call FLOAT2STR(OBJ%SIGABS*10.D0,CNUM1)
      write(*,21) 'sigma_abs',trim(CNUM1),'1/cm/A'
      call FLOAT2STR(OBJ%SIGINC*10.D0,CNUM1)
      write(*,21) 'sigma_inc',trim(CNUM1),'1/cm'
      call FLOAT2STR(OBJ%SIGTDS(1)*OBJ%LAMBDA*10.D0,CNUM1)
      write(*,21) 'sigma_sph',trim(CNUM1),'1/cm/A'
      lam2=OBJ%LAMBDA**2
      Z=OBJ%SIGTDS(3)/lam2 + OBJ%SIGTDS(4)/lam2**2
      call FLOAT2STR(OBJ%SIGTDS(2)*(1.D0-exp(-Z))*10.D0,CNUM1)
      write(*,21) 'sigma_mph',trim(CNUM1),'1/cm'
      call FLOAT2STR(OBJ%SIGTDS(2)*10.D0,CNUM1)
      write(*,21) 'sigma_fa',trim(CNUM1),'1/cm'
      call FLOAT2STR(OBJ%SIGTDS(3),CNUM1)
      write(*,21) 'B',trim(CNUM1),'A^2'
      call FLOAT2STR(OBJ%SIGTDS(4),CNUM1)
      write(*,21) 'D',trim(CNUM1),'A^4'
      end subroutine XTAL_REPORT_ABS


!-------------------------------------------------------------
      logical function XTAL_ARRAY_GO(CR,INST)
! Main procedure for tracing in XTAL
! INPUT:
!   CR    .. crystal data
!-------------------------------------------------------------
      TYPE (XTAL),intent(inout) :: CR
      integer,intent(in) :: INST
  ! segment data
      integer,parameter :: MSEG=128
      REAL(KIND(1.D0)) :: TSEG(0:MSEG),TSEG0(0:MSEG)
      integer :: ISEG(3,MSEG),NSEG
  ! tracing nodes
      integer,parameter :: MSEL=MSEG*MAX_NODES ! max. MAX_NODES primary reflections per segment
      integer :: NSEL,SSEL(MSEL),SEL_MLOC(MSEL),RSEL(MSEL)
      REAL(KIND(1.D0)) :: PSEL(0:MSEL),TSEL(MSEL)
      REAL(KIND(1.D0)) :: SEL_TR(3,MSEL),SEL_RLOC(3,3,MSEL)
  ! reflection data
      REAL(KIND(1.D0)) :: TREF(MAX_NODES+1),PREF(MAX_NODES)
      integer :: IREF(MAX_NODES)
  ! other local variables
      REAL(KIND(1.D0)) :: tin,TOF,MI,TR(3),RLOC(3,3),Z
      integer :: INODE,SELNODE,i,j,k,MLOC,NP
      logical :: transmitted
      TYPE(NEUTRON) :: NEU0
      TYPE (PRLNODE) :: NODE
1     format(a,': ',7(1x,G12.5))

      myDBG=.false.
      TMPN(INST)%P=0.D0

      ! if (myDBG) write(*,1) 'XTAL_ARRAY_GO TIN ',NEUT%T, CR%FRAME%COUNT
      CALL SLIT_PRE(CR%FRAME)
    !  write(*,*) 'XTAL_ARRAY_GO ',CR%FRAME%COUNT
      XTAL_ARRAY_GO=.false.

      !myDBG=((DBGIT>1000).and.(DBGIT<1020))
      !myDBG=(NEUT%R(2)>55.D0)

      !call ARRAY3D_SETDBG(myDBG)
      !DBGIT=DBGIT+1
      if (myDBG) write(*,*)
      if (myDBG) write(*,1) '--------------------------------------------------------'
      if (myDBG) write(*,1) 'XTAL_ARRAY_GO R,K',NEUT%R, NEUT%K
      if (myDBG) write(*,1) 'bounding box',CR%A3D%SA

! get time to crystal entry, go to end if entry is missed
! permit back step if R is inside the array, to ensure the start at the array entry
      if (.not.ARRAY3D_ENTRY(CR%A3D,NEUT%R,NEUT%K,tin)) goto 98
! move to the crystal entry
      call Transport(tin)

      ITER=0    ! iteration counter
      TOF=0.D0  ! accumulated time of flight
      INODE=1   ! starting node in rec. lattice (1=direct beam)
      PSEL(0)=0.D0  ! this is always zero
      MI=XTAL_GETMI(CR,NEUT%K0)  ! absorption coeficient
      SELNODE=RLNODE_GET(CR%ICR,CR%HKL) ! target node index
      NEU0=NEUT  ! NEU0 is neutron in array local coordinates
      transmitted=.false.
  !    if (myDBG) write(*,1) 'tin',tin
  !    if (myDBG) write(*,1) 'array entry R,K',NEUT%R, NEUT%K
! BEGIN LOOP
      do while(.not.transmitted)
!----------------------------------------------------
      ! get crossing times for segments
      call ARRAY3D_BORDERS_NEW(CR%A3D,NEU0%R,NEU0%K,TSEG0,TSEG,ISEG,NSEG)
  !    if (myDBG) write(*,*) 'NSEG=',NSEG,' tout=',TSEG(NSEG)
    !  exit
! scan crossed segments and accumulate probabilities for all segments and reflections
      NSEL=0
      ! get current node data
      if (NSEG.gt.0) call RLNODE_GET_ITEM(CR%ICR,INODE,NODE)
      do i=1,NSEG
      ! transform NEUT to segment local coordinates
        call ARRAY3D_SEG_TRANS(CR%A3D,ISEG(1,i),TR,RLOC,MLOC)
        call NEUT_TO_SEGMENT(TR,RLOC,MLOC,NEU0,NEUT)
      ! move to the segment entry
        call Transport(TSEG0(i))
        if (myDBG) write(*,1) 'local iseg, R,K',i,NEUT%R, NEUT%K
      ! collect acessible paths from the current node
        call XTAL_GET_PATHS(CR%ICR,TSEG(i)-TSEG0(i),NODE,TREF,PREF,IREF,NP)
      ! accumulate total reflection probability and store other info for each path
        do j=1,NP
          !if (myDBG) write(*,*) 'nseg=',i,'  iref=',j,' tof =',TREF(j),' p=',PREF(j)
          NSEL=NSEL+1
          PSEL(NSEL)=PSEL(NSEL-1)+(1.D0-PSEL(NSEL-1))*PREF(j)
          TSEL(NSEL)=TREF(j) ! record reflection times relative to segment entry
          RSEL(NSEL)=IREF(j) ! record reflection indices for i-th segment
          SSEL(NSEL)=i ! record segment indices
          ! store transformation data
          SEL_MLOC(NSEL)=MLOC
          do k=1,3
            SEL_TR(k,NSEL)=TR(k)
            SEL_RLOC(1:3,k,NSEL)=RLOC(1:3,k)
          enddo
        enddo
      enddo

    ! play Russian rulette to select segment and reflection
      if (NSEL.gt.0) then
        Z=RAN1()
        I=1
        do while ((I.LE.NSEL).and.(Z.gt.PSEL(I)))
          I=I+1
        enddo
        if (I.gt.NSEL) transmitted=.true.
    ! no reflection on the way => exit crystal array
      else
        transmitted=.true.
      endif
    ! transmitted => move to array exit point
    ! this should be the only way out of the cycle
      if (transmitted) then
      ! get TOF through material
        do j=1,NSEG
          TOF=TOF+TSEG(j)-TSEG0(j)
        enddo
      ! return to the local coordinates and transport to the end of the last segment
        NEUT=NEU0
        call Transport(TSEG(NSEG))
        if (myDBG) write(*,1) '   transmitted r,k',NEUT%R,NEUT%K
        if (myDBG) write(*,1) '   path',TOF*NEUT%K0,TOF*NEUT%K(3)

      else
        !myDBG=((DBGIT>=0).and.(DBGIT<50))
        !DBGIT=DBGIT+1
        ! if (myDBG) write(*,1) '   reflected in seg',ISEG(1:3,i)
        if (myDBG) write(*,1) '   reflected from hkl',NODE%N%HKL
        ! transform NEUT to selected segment local coordinates
        call NEUT_TO_SEGMENT(SEL_TR(1,i),SEL_RLOC(1,1,i),SEL_MLOC(i),NEU0,NEUT)
        ! if (myDBG) write(*,*) '    iter=',ITER,' sel=',i,' seg=',SSEL(i),' R=',NEUT%R
      ! get next node index
        INODE=NODE%N%NEXT(RSEL(I))
        if (myDBG) write(*,1) '   node selected, next',SELNODE,INODE
      ! get TOF through material
        do j=1,SSEL(i)-1
          TOF=TOF+TSEG(j)-TSEG0(j)
        enddo
        TOF=TOF+TSEL(I)
        if (myDBG) write(*,1) '   dist. to reflection',(TSEG0(SSEL(i))+TSEL(i))*NEUT%K0
        ! transport to ref. point and reflect (must be done in segment local coord.)
        call XTAL_REFLECT(CR%ICR,NODE,TSEG0(SSEL(i))+TSEL(i),RSEL(i))
        ! if (myDBG) write(*,1) '   reflected at local r,k',NEUT%R,NEUT%K
        !if (myDBG) write(*,*) '    reflect ',' Z=',NEUT%K(3)
      ! transform NEUT back to the array local coordinates
        call NEUT_TO_ARRAY(SEL_TR(1,i),SEL_RLOC(1,1,i),SEL_MLOC(i),NEUT,NEU0)
        ! if (myDBG) write(*,*) '    trans   ',' Z=',NEUT%K(3)
      ! register event
        TMPN(INST)=NEU0
        if (TRACING_UP) TMPN(INST)%K=-TMPN(INST)%K
        call AddEventLog(CR%FRAME,NEU0)
        if (DBGREC.eq.CR%ICR) then
          do j=1,3
            TMP_TRACE_REC(j,ITER+1)=NODE%N%HKL(j)
            TMP_TRACE_REC(j+3,ITER+1)=NEUT%R(j)
          enddo
        endif
        if (myDBG) write(*,1) '   reflected at array r,k',NEU0%R,NEU0%K
      endif

      ITER=ITER+1
      if (ITER.GT.MAXITER) then
        CR%MAX_ITER_COUNT=CR%MAX_ITER_COUNT+1
        if (DBGREC.eq.CR%ICR) then
          do j=1,MAXITER
            XTALS_TRACE_REC(1:6,j)=TMP_TRACE_REC(1:6,j)
          enddo
        endif
        goto 99
      endif
! END LOOP
      enddo

!----------------------------------------------------
      call AddEventLog(CR%FRAME,NEUT)
    ! end node is not the selected one ?
      if (INODE.ne.SELNODE) goto 99
  ! weight by transmission probability
      NEUT%P=NEUT%P*EXP(-MI*TOF*NEUT%K0)
      if (myDBG) write(*,1) 'XTAL_ARRAY_GO finished at r,k,p',NEUT%R,NEUT%K,EXP(-MI*TOF*NEUT%K0)
  ! increment event counter
      call incCounter(CR%FRAME%COUNT)
      CALL SLIT_POST(CR%FRAME)
      ! if (myDBG) write(*,1) ' R,KF=',NEUT%R,NEUT%K
      XTAL_ARRAY_GO=.TRUE.
    !  if (CR%FRAME%COUNT.gt.10) stop
      RETURN

! missed array
98    if (SELNODE.eq.1) then
    ! selected node is (0,0,0) = transmission
  ! increment event counter
        call incCounter(CR%FRAME%COUNT)
        CALL SLIT_POST(CR%FRAME)
        XTAL_ARRAY_GO=.TRUE.
!        myDBG=((DBGIT>0).and.(DBGIT<20).and.(NEUT%R(1)<-4))
!        DBGIT=DBGIT+1
!        if (myDBG) TMPN(INST)=NEU0
        if (myDBG) write(*,1) 'XTAL_ARRAY_GO missed, R, tin',NEUT%R,tin
        return
      endif
99    NEUT%P=0
      if (myDBG) write(*,1) 'XTAL_ARRAY_GO no entry, stopped, tin',tin
      XTAL_ARRAY_GO=.FALSE.
      end function XTAL_ARRAY_GO

!-------------------------------------------------------------
      SUBROUTINE XTAL_REFLECT(ICR,NODE,TR,IR)
! transfer to the point of reflection and turn direction to a new node
! ICR ... crystal index
! NODE ... current reflection node
! TR ... time to reflection point
! IR ... index to reflection data (index to NODE%PATH)
!-------------------------------------------------------------
      integer,intent(in) :: ICR,IR
      TYPE (PRLNODE),intent(in) :: NODE
      real(kind(1.D0)),intent(in) :: TR
      TYPE (PXREF) :: PXR
      real(kind(1.D0)) :: Z,AUX(3)
      integer :: j
1     format(a,': ',7(1x,G12.5))
      ! transfer to the point of reflection
        call Transport(TR)
      ! turn direction to a new node
      ! get reflection data
        call XREF_GET_ITEM(ICR,NODE%N%PATH(IR),PXR)
      ! Calculate local G-vector deviation
        call XTAL_SEG_GLOC(PXR%X,NEUT%R,AUX)
        if (myDBG) write(*,1) 'XTAL_REFLECT hkl, G',PXR%X%HKL,AUX
        ! call M3xV3(1,PXR%X%DG_DR,NEUT%R,AUX)
      ! turn direction to a new node
        NEUT%K=NEUT%K+AUX
      ! renormalize
        Z=NEUT%K0/SQRT(NEUT%K(1)**2+NEUT%K(2)**2+NEUT%K(3)**2)  ! set |K+G| = |K|
        DO j=1,3
          NEUT%K(j)=NEUT%K(j)*Z
        ENDDO
      ! add a small shift along new direction to avoid locking at the reflection point
        call Transport(1.D-5/NEUT%K0)
      end SUBROUTINE XTAL_REFLECT

!-----------------------------------------------------------------------------
      SUBROUTINE XTAL_GET_PATHS(ICR,TOUT,NODE,TREF,PREF,IREF,NP)
! Calculate reflection probabilities for all accessible reflections
! for given crystal and given rec. lattice node
! INPUT:
!   CR  .. crystal data
!   NODE .. reciprocal lattice node (item from RLNODE_TAB)
! OUTPUT:
!   TREF   .. flight times [h_bar/m] of reflections
!   PREF   .. corresponding cumulative probabilities
!   IREF   .. pointers to selected reflections (ref. to RLNODE_TAB)
!   NP     .. number of elements written to TREF,PREF arrays
! NOTE: dimension of TREF must be at least NP+1
! TREF(NP+1) stores time to the frame exit
!-----------------------------------------------------------------------------
    !  TYPE (XTAL),intent(in) :: CR
      TYPE (PRLNODE),intent(in) :: NODE
      real(kind(1.D0)),intent(in) :: TOUT
      integer,intent(in) :: ICR
      real(kind(1.D0)),intent(out) :: TREF(:),PREF(:)
      integer,intent(out) :: IREF(:)
      integer,intent(out) :: NP
      integer :: I,J,K
      real(kind(1.D0)) :: DT,P,TMIN
      TYPE (PXREF) :: PXR
    ! logical :: FRAME_BORDER
        NP=0
    !  if (FRAME_BORDER(CR%FRAME,NEUT%R,NEUT%K,tin,tout)) then
      ! lower limit for TOF to the reflection point
      ! required to avoid inf. cycles due to num. precission
        TMIN=1.D-3/NEUT%K0
        J=0
        do I=1,NODE%N%NP
          call XREF_GET_ITEM(ICR,NODE%N%PATH(I),PXR)
          call XTAL_FINDREF(PXR%X,DT,P)
          if ((DT.gt.TMIN).and.(DT.le.TOUT).and.(P.gt.1.D-10)) then
        ! keep ascending order in TREF
          ! append
            if (NP.eq.0) then
                  J=NP+1
            else if (TREF(NP).le.DT) then
                  J=NP+1
          ! insert at appropriate position
            else
              J=1
              do while ((J.le.NP).and.(TREF(J).le.DT))
                J=J+1
              enddo
            ! make space, push the stack
              K=NP
              do while (K.ge.J)
                TREF(K+1)=TREF(K)
                PREF(K+1)=PREF(K)
                IREF(K+1)=IREF(K)
                K=K-1
              enddo
            endif
            TREF(J)=DT
            PREF(J)=P
            IREF(J)=I
            NP=NP+1
          endif
        enddo
      ! calculate correct cumulative probabilities
        do K=2,NP
          PREF(K)=PREF(K-1)+(1.D0-PREF(K-1))*PREF(K)
        enddo
        ! store time to frame exit in TREF(NP+1)
        TREF(NP+1)=TOUT
      ! if there is no cross-section with FRAME, set tout=0
      !  if (NP.eq.0) TREF(1)=0.D0
      end SUBROUTINE XTAL_GET_PATHS


!-----------------------------------------------------------------------------
      SUBROUTINE XTAL_FINDREF(XR,DT,PREF)
! Calculate TOF to the point of reflection in BENT CRYSTAL
! Works on NEUT, assumed to be in local coordinates
! INPUT:
!   XR  .. reflection data
! OUTPUT:
!   DT   .. time of flight in [h_bar/m]
!   PREF .. peak reflectivity
!-----------------------------------------------------------------------------
      TYPE (XREF),intent(in) :: XR
      real(KIND(1.D0)), intent(out) :: DT,PREF
      real(KIND(1.D0)),parameter :: EPS=1.D-20
      real(KIND(1.D0)),parameter :: DKKMIN=1.D-6
      INTEGER*4 I,J
      real(KIND(1.D0)) :: a,b,c,Z,det,KG2,GABS,DT1,DT2,RTD
      real(KIND(1.D0)) :: KaG(3),G(3),R(3),dz,dkk,DT0
      logical :: db=.false.
1     format(a,': ',7(1x,G12.5))

      ! db=(DBGIT<20)
      ! db=myDBG
      ! db=((DBGIT>100000).and.(DBGIT<100020))
      if (db) write(*,1) 'XTAL_FINDREF R',DBGIT,NEUT%R(1:3),XR%HKL
      dkk=1.D30
      ITER=0
      DT=0.D0
      DT0=0.D0
      R=NEUT%R
      PREF=0.D0
      GABS=0.D0
      b=0.D0

      do while ((DT.ge.-1.D1).and.(abs(dkk)>DKKMIN).and.(ITER<3))
        ITER=ITER+1
        call XTAL_SEG_GLOC(XR,R,G)

        if (db) write(*,1) '   local G',G
! Get local deviation from the Bragg condition, c=(K+G)^2-K^2
        GABS=0.D0
        DO I=1,3
          KaG(I)=NEUT%K(I)+G(I)
          GABS=GABS+G(I)**2
        ENDDO
        GABS=SQRT(GABS)
        KG2=KaG(1)**2+ KaG(2)**2+KaG(3)**2
        c=KG2-NEUT%K0**2
      ! delta_k/k from the Bragg point
        dkk=c/GABS/NEUT%K0
        if (db) write(*,1) '   iter, dkk, dz, DT',iter, dkk, dz, DT
        if (abs(dkk)>DKKMIN) then
! get linear term, b=(K+G).DG_DR.K
! and quadratic term, a=|DG_DR.K|^2
          b=0.D0
          a=0.D0
          DO I=1,3
            Z=0.D0
            DO J=1,3
              Z=Z+XR%DG_DR(I,J)*NEUT%K(J)
            ENDDO
            a=a+Z**2
            b=b+KaG(I)*Z
          ENDDO
  ! solve quadratic equation: get the smaller positive solution, if any
          det=b*b-a*c
          if (db) write(*,1) '        a,b,det:',a,b,det
          if ((a.gt.EPS).and.(det.gt.EPS)) then
            RTD=SQRT(det)
            if (b>0) then
              DT=(-b+RTD)/a+DT0
            else
              DT=(-b-RTD)/a+DT0
            endif

          !  DT1=(-b+RTD)/a+DT0
          !  DT2=(-b-RTD)/a+DT0
          !  if ((DT1.gt.0).and.(DT2.gt.0)) then
          !    DT=min(DT1,DT2)
          !  else
          !    DT=max(DT1,DT2)
          !  endif
            R=R+NEUT%K*(DT-DT0)
            DT0=DT
          else
            DT=-1.D30
          endif
          if (db) write(*,1) '        dt1, dt2, dt:',DT1,DT2,DT
        endif
      enddo
  ! get peak reflectivity
  ! r = 1 - exp(-Qkin/|beta|), where beta=b/(G.k^2.cos(thetaB))
  ! Qkin=PI.QML/(k^2.d.cos(thetaB))
  ! hence Qkin/|beta| = PI.G.QML/(DHKL.b)
      if (DT.ge.0.D0) then
        if (abs(b).lt.1.D-10) then
          PREF=1.D0
        else
          PREF = 1.D0 - exp(-1.D-1*XR%QML*PI*GABS/(XR%DHKL*abs(b)))
        endif
      ELSE
        DT=-1.D30
      endif
      END SUBROUTINE XTAL_FINDREF


!=================================
! OBSOLETE
!==================================


!-------------------------------------------------------------
      logical function XTAL_GO_INST(CR)
! obsolete, call XTAL_ARRAY_GO instead
! Main procedure for tracing in XTAL, single slab
! INPUT:
!   CR    .. crystal data
!-------------------------------------------------------------
      TYPE (XTAL),intent(inout) :: CR
      real(kind(1.D0)) :: TOF,MI,tin,tout
      integer :: INODE,SELNODE,i
    !  DBG=((CR%FRAME%COUNT.gt.90).and.(CR%FRAME%COUNT.lt.100))
    !  DBG=(trace_ctot.eq.100)

      CALL SLIT_PRE(CR%FRAME)
      if (DBG) call PrintEvent(smes,CR%FRAME,' ')
    ! if there is no cross-section with FRAME, set tout=0
      if (.not.FRAME_BORDER(CR%FRAME,NEUT%R,NEUT%K,tin,tout)) goto 99
    ! move to the crystal entry
      call Transport(tin)
! trace through the crystal
  ! accumulate time-of-flight in TOF
      TOF=0.D0  ! accumulated time of flight
      INODE=1   ! starting node in rec. lattice
      ITER=0    ! iteration counter
      MI=XTAL_GETMI(CR,NEUT%K0)  ! absorption coeficient
      SELNODE=RLNODE_GET(CR%ICR,CR%HKL) ! target node index
      do while ((ITER.LT.MAXITER).and.(.not.XTAL_STEP(CR,INODE,TOF)))
        ITER=ITER+1
      enddo
      if (DBG) write(*,*) 'XTAL_GO TOF=',TOF,' INODE=',INODE
    ! too many iterations - warning on console, result=false
      if (ITER.ge.MAXITER) then
      ! report iteration limit
      !  write(*,*) 'XTAL_GO: too many iterations ',CR%MAX_ITER_COUNT
        CR%MAX_ITER_COUNT=CR%MAX_ITER_COUNT+1
        if (DBGREC.eq.CR%ICR) then
          do i=1,MAXITER
            XTALS_TRACE_REC(1:6,i)=TMP_TRACE_REC(1:6,i)
          enddo
        endif
        goto 99
      endif
    ! no intersection with the crystal
      if (TOF.le.0) goto 99
    ! end node is not the selected one
      if (DBG) write(*,*) 'XTAL_GO selnode=',SELNODE,' hkl=',CR%HKL
      if (INODE.ne.SELNODE) goto 99
  ! weight by transmission probability
      NEUT%P=NEUT%P*EXP(-MI*TOF*NEUT%K0)
    !  DBG=.true.
      if (DBG) call PrintEvent(smes,CR%FRAME,'passed')
  ! reject invalid events
    !  if (NEUT%P.LE.0.D0) goto 99
  ! increment event counter
      call incCounter(CR%FRAME%COUNT)
! convert to exit coordinates
      CALL SLIT_POST(CR%FRAME)
      XTAL_GO_INST=.TRUE.
      RETURN

99    NEUT%P=0
      XTAL_GO_INST=.FALSE.
      END FUNCTION XTAL_GO_INST


!-------------------------------------------------------------
      logical function XTAL_STEP(CR,INODE,DT)
! obsolete, using XTAL_ARRAY_GO instead
! Make one step of random walk in the crystal with NEUT
! INPUT:
!   CR    .. crystal data
!   INODE .. reciprocal lattice node index (ref. to RLNODE_TAB)
!   DT    .. cumulative tie of flight
! RETURNS:
!   true, if neutron exits the crystal without further reflection
!-------------------------------------------------------------
      TYPE (XTAL),intent(in) :: CR
      REAL(KIND(1.D0)),intent(inout) :: DT
      integer,intent(inout) :: INODE
      REAL(KIND(1.D0)) :: Z,TREF(MAX_NODES+1),PREF(MAX_NODES),tin,tout
      integer :: i,j,IREF(MAX_NODES),NP
      TYPE (PRLNODE) :: NODE
      logical :: transmitted
      call RLNODE_GET_ITEM(CR%ICR,INODE,NODE)

      if (FRAME_BORDER(CR%FRAME,NEUT%R,NEUT%K,tin,tout)) then
        call XTAL_GET_PATHS(CR%ICR,tout,NODE,TREF,PREF,IREF,NP)
      else
        NP=0
        TREF(1)=0.D0
      endif
      if (NP.gt.0) then
    ! play Russian rulette to select reflection or transmission
        ! Z=RAN1()*PREF(NP) ! this was a bug !!
        Z=RAN1()
        I=1
        do while ((I.LE.NP).and.(Z.gt.PREF(I)))
          I=I+1
        enddo
      endif
    ! reflection
      if ((NP.gt.0).and.(I.LE.NP)) then
        INODE=NODE%N%NEXT(IREF(I))
        transmitted=.false.
      ! transfer to the point of reflection
        DT=DT+TREF(I)
      ! turn direction to a new node
        call XTAL_REFLECT(CR%ICR,NODE,TREF(I),IREF(I))
        if (DBGREC.eq.CR%ICR) then
          do j=1,3
            TMP_TRACE_REC(j,ITER+1)=NODE%N%HKL(j)
            TMP_TRACE_REC(j+3,ITER+1)=NEUT%R(j)
          enddo
        endif
    ! transmission
      else
        transmitted=.true.
        DT=DT+tout
        call Transport(tout)
      endif
      XTAL_STEP=transmitted
      ! log event
      call AddEventLog(CR%FRAME,NEUT)
      end function XTAL_STEP


      end module XTALS_TRACE