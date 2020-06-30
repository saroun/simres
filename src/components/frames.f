!//////////////////////////////////////////////////////////////////////
!////  $Id: frames.f,v 1.46 2019/06/22 20:48:04 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.46 $
!////     $Date: 2019/06/22 20:48:04 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Describes FRAME, the base object (ancestor) of all beamline objects
!////
!////////////////////////////////////////////////////////////////////////
      MODULE FRAMES
      use SIMATH
      use CLASSES
      use FIELDDEF
      use FIELDDATA
      use XMLINFO
      USE MESSAGES
      implicit none

      integer, parameter :: top_gonio=0
      integer, parameter :: top_stage=1


C SLIT:
      TYPE TFRAME; SEQUENCE
        INTEGER :: CLASS ! class ID
        character(LEN_ID)   :: ID                 ! component ID string (e.g. MON)
        CHARACTER(LEN_NAME) :: NAME               ! descriptive name of the component (e.g. Monochromator)
        REAL(KIND(1.D0)) :: DIST                      ! distance
        REAL(KIND(1.D0)) :: AX(3)                     ! exit axis angles
        REAL(KIND(1.D0)) :: SIZE(3),STA(3),GON(3)     ! dimensions, linear stage, gonio
        REAL(KIND(1.D0)) :: VEL(3)                    ! velocity (assumed constant)
        ! REAL(KIND(1.D0)) :: CTR(3)                    ! position of body center in local coordinates
        REAL(KIND(1.D0)) :: RLOC(3,3)           ! rotation from incident axis to local coordinates
        REAL(KIND(1.D0)) :: REXI(3,3)           ! rotation from incident to exit axis coordinates
        REAL(KIND(1.D0)) :: RLAB(3,3)           ! rotation matrix: incident axis to lab. coordinates
        real(kind(1.D0)) :: RDEFL(3,3)          ! deflection rotation matrix
        REAL(KIND(1.D0)) :: TEXI(3)             ! translation to exit frame origin (rel. to incident axis)
        REAL(KIND(1.D0)) :: TLAB(3)             ! transition vector: incident axis to lab. coordinates
        real(kind(1.D0)) :: TDEFL(3)            ! deflection translation vector
        REAL(KIND(1.D0)) :: DM(3)               ! misalignment [mm]
        REAL(KIND(1.D0)) :: INSIZE(3)           ! inner size for hollow shapes
        logical :: ISHOLLOW                     ! true if hollow shape
        logical :: ISDEFL                       ! true if deflecting
        !integer :: INST    ! instance number of the component
        integer :: REGISTRY                     ! storage index, for components that store neutrons
        integer :: ORDER                  ! top_gonio or top_stage (irder of stage and gonio)
        INTEGER :: MLOC,MEXI              ! mask (=0 if RLOC or REXI are not required, else 1)
        INTEGER :: COUNT                  ! event counter
        INTEGER :: SHAPE            ! shape ID
        INTEGER :: TRANSPARENT      ! flag indicating a transparent component
        INTEGER :: IRNDP            ! indices to random numbers array (penetration depth)
        INTEGER :: ISMOVED          ! 1 if FRAME_LOCAL is required (STA,GON or VEL <>0)
      END TYPE TFRAME


      private FRAME_INP_INST,FRAME_INP_OBJ

! overloaded INP
      interface FRAME_INP
        module procedure FRAME_INP_INST
        module procedure FRAME_INP_OBJ
      end interface
! overloaded OUT
      interface FRAME_OUT
        module procedure FRAME_OUT_INST
        module procedure FRAME_OUT_OBJ
      end interface

      integer, parameter :: FRAMES_DIM=128
! pointer type
      type PFRAME
        integer :: IDX  ! instance index
        TYPE(TFRAME),pointer :: X
      end type PFRAME
! instances of XTAL. AXTALS(0) is always unallocated
      integer :: FRAMES_NC
      type(PFRAME) :: AFRAMES(1:FRAMES_DIM)


      private FRAME_INIT_I,FRAME_INIT_O
! overloaded INIT
      interface FRAME_INIT
        module procedure FRAME_INIT_I
        module procedure FRAME_INIT_O
      end interface

      contains


!-------------------------------------------------------------------------
! Creator for FRAME, return instance number
! Memory is allocated on the first free element of AFRAMES array
!-------------------------------------------------------------------------
      integer function FRAME_CREATE(ID,NAMESTR)
      character(*) :: ID,NAMESTR
      integer :: ierr,i
      character(32) :: S
        ierr=1
        i=1
        do while ((i<FRAMES_DIM).and.(AFRAMES(i)%IDX>0))
          i=i+1
        enddo
        if (AFRAMES(i)%IDX.le.0) then
          allocate(AFRAMES(i)%X,STAT=ierr)
          if(ierr.eq.0) then
            FRAMES_NC=FRAMES_NC+1
            call FRAME_DEFAULT(i)
            AFRAMES(i)%IDX=i
            AFRAMES(i)%X%ID=trim(ID)
            AFRAMES(i)%X%NAME=trim(NAMESTR)
          endif
        endif
        if (ierr.eq.0) then
          FRAME_CREATE=i
        else
          call INT2STR(FRAMES_DIM,S)
          call MSG_ERROR('FRAME_CREATE','Reached maximum number or instances ('//trim(S)//')',1,1)
          FRAME_CREATE=0
        endif
      end function FRAME_CREATE

!---------------------------------------------------------
      subroutine FRAME_DISPOSE(INST)
!---------------------------------------------------------
        integer,intent(in) :: INST
        if ((inst.gt.0).and.(inst.le.FRAMES_DIM)) then
          if (associated(AFRAMES(inst)%X)) then
            DEallocate(AFRAMES(inst)%X)
            nullify(AFRAMES(inst)%X)
            FRAMES_NC=FRAMES_NC-1
          endif
          AFRAMES(inst)%IDX=0
        endif
      end subroutine FRAME_DISPOSE



!----------------------------------------------------
      SUBROUTINE FRAME_INIT_I(INST)
! initialize object FRAME:
! - clear event counter
!----------------------------------------------------
      integer,intent(in) :: INST
      TYPE(TFRAME),POINTER :: OBJ
      if (.not.FRAME_isValid(INST)) return
      OBJ => AFRAMES(INST)%X
      call FRAME_INIT_O(OBJ)
      END SUBROUTINE FRAME_INIT_I

!----------------------------------------------------
      SUBROUTINE FRAME_INIT_O(OBJ)
! initialize object FRAME:
! - clear event counter
!----------------------------------------------------
      TYPE(TFRAME),intent(out) :: OBJ
      REAL(KIND(1.D0)) :: AUX(3)
      logical :: rep
      rep=.false.
      !write(*,*) 'FRAME_INIT ' ,trim(OBJ%ID),' ',OBJ%COUNT
! clear counter
      OBJ%COUNT=0
      OBJ%IRNDP=0
! clear misalignment
      OBJ%DM=0.0
! LINEAR STAGE
  ! The stage is ON TOP of the goniometer, with axes // to the local coordinates
! debug report
      if (rep) then
        call XML_RSXDUMP(6,' ',1)
        call XML_ARRAY(6,3,'position of '//trim(OBJ%ID),'x|y|z',OBJ%STA)
        AUX(1:3)=OBJ%AX(1:3)/deg
        call XML_ARRAY(6,3,'exit axis of '//trim(OBJ%ID),'y|x|z',AUX(1))
        AUX(1:3)=OBJ%GON(1:3)/deg
        call XML_ARRAY(6,3,'goniometer of '//trim(OBJ%ID),'y|x|y1',AUX(1))
        call XML_MATRIX(6,3,3,'gonio matrix for '//trim(OBJ%ID),'x|y|z','x|y|z',OBJ%RLOC(1,1),3)
        call XML_MATRIX(6,3,3,'axis matrix for '//trim(OBJ%ID),'x|y|z','x|y|z',OBJ%REXI(1,1),3)
        call XML_RSXDUMP(6,' ',0)
      endif
      END SUBROUTINE FRAME_INIT_O

!---------------------------------------------------------
      subroutine FRAME_DISPOSE_ALL
!---------------------------------------------------------
        integer :: i
        do i=1,FRAMES_DIM
          call FRAME_DISPOSE(i)
        enddo
      end subroutine FRAME_DISPOSE_ALL

!-------------------------------------------------------------------------
! Add a component
!-------------------------------------------------------------------------
      INTEGER function AddFRAME(INST)
      integer,intent(in) :: INST
      integer :: i
      i=0
      if (FRAME_isValid(INST)) then
        i=FRAME_CREATE(AFRAMES(INST)%X%ID,AFRAMES(INST)%X%NAME)
        if (i.gt.0) then
          AFRAMES(i)%X=AFRAMES(INST)%X
          AFRAMES(i)%IDX=i
        !  write(*,*) 'AddFRAME ',trim(AFRAMES(i)%X%ID),i
        endif
      endif
      AddFRAME=i
      end function AddFRAME

!-------------------------------------------------------------
      logical function FRAME_isValid(INST)
! check index and association
!-------------------------------------------------------------
      integer,intent(in) :: INST
      FRAME_isValid= ((INST.gt.0).and.(INST.le.FRAMES_DIM).and.associated(AFRAMES(INST)%X))
      end function FRAME_isValid

!-------------------------------------------------------------
      SUBROUTINE FRAME_GET(INST,OBJ)
! return reference to given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(PFRAME),intent(out) :: OBJ
        if (FRAME_isValid(INST)) then
          OBJ%IDX=AFRAMES(INST)%IDX
          OBJ%X => AFRAMES(INST)%X
        else
          NULLIFY(OBJ%X)
          OBJ%IDX=0
        endif
      end SUBROUTINE FRAME_GET

!------------------------------------
      SUBROUTINE FRAME_DEFAULT(INST)
! set default parameters
!------------------------------------
      integer,intent(in) :: INST
      TYPE(TFRAME),POINTER :: OBJ
        if (.not.FRAME_isValid(INST)) return
        OBJ => AFRAMES(INST)%X
        call FRAME_CLEAR(OBJ)
        OBJ%ORDER=top_gonio
      END SUBROUTINE FRAME_DEFAULT

!------------------------------------
      SUBROUTINE FRAME_PREPARE(INST,IERR)
! set default parameters
!------------------------------------
      integer,intent(in) :: INST
      integer,intent(out) :: IERR
      TYPE(TFRAME),POINTER :: OBJ
        IERR=0
        if (.not.FRAME_isValid(INST)) return
        OBJ => AFRAMES(INST)%X
        call FRAME_INIT_MAT(OBJ)
      END SUBROUTINE FRAME_PREPARE


!------------------------------------
      SUBROUTINE FRAME_CLEAR(OBJ)
! set default parameters
!------------------------------------
      TYPE(TFRAME),intent(out)  :: OBJ
        OBJ%SHAPE=FRAME_SHAPE_BOX
        OBJ%CLASS=CCLS_FRAME
        OBJ%DIST=0.D0
        OBJ%AX=0.D0
        OBJ%SIZE=0.D0
        OBJ%STA=0.D0
        OBJ%GON=0.D0
        OBJ%VEL=0.D0
        OBJ%INSIZE=0.D0
        OBJ%ISMOVED=0
        OBJ%ISDEFL=.false.
        OBJ%ISHOLLOW=.false.
        OBJ%REGISTRY=0
        OBJ%TRANSPARENT=0
        OBJ%DM=0.0
        OBJ%TEXI=0.0
        OBJ%TDEFL=0.0
        OBJ%TLAB=0.0
        call UNIMAT(OBJ%RLOC,3,3)
        call UNIMAT(OBJ%REXI,3,3)
        call UNIMAT(OBJ%RDEFL,3,3)
        call UNIMAT(OBJ%RLAB,3,3)
      END SUBROUTINE FRAME_CLEAR

!-----------------------------------------------------------
! Set exit axis angles and update transformation matrices.
! Incident axis is (0,0,1)
! INPUT:
!   AX   .. (y,x,y) Euler angles (!! changed from y,x,z system !!)
!-----------------------------------------------------------
      subroutine SetAxisAngles(OBJ,AX)
      type(TFRAME) :: OBJ
      real(kind(1.D0)),intent(in) :: AX(3)
      real(kind(1.D0)) :: R1(3,3),R2(3,3),R3(3,3),R(3,3)
      integer :: i
      OBJ%AX(1:3)=AX(1:3)
  ! rot. matrix for lower y-axis (horizontal, usually scattering angle)
      !call MK_ROT3(2,OBJ%AX(1),R1)
  ! rot. matrix for sagital angle
      !call MK_ROT3(1,OBJ%AX(2),R2)
  ! rot. matrix for the upper y-axis (attached to the next component)
      !call MK_ROT3(2,OBJ%AX(3),R3)
  ! total rot. matrix from incident to exit axis coordinates
      !CALL M3XM3(1,R3,R2,R)
      !CALL M3XM3(1,R,R1,OBJ%REXI)

! use SIMATH library:
      call getRotYXY(AX,OBJ%REXI)
! flag for avoiding unnecessary calls of transformation subroutines
      OBJ%MEXI=0
      do i=1,3
        if (OBJ%REXI(i,i).NE.1.D0) OBJ%MEXI=1
      enddo
      end subroutine SetAxisAngles


!-----------------------------------------------------------
! Set goniometer angles defined as (y,x,y) Euler angles
! Update also transformation matrix.
! Incident axis is (0,0,1)
! INPUT:
!   GON(3) .. the (y,x,y) Euler angles of the goniometer
!-----------------------------------------------------------
      subroutine SetGonioAngles(OBJ,GON)
      type(TFRAME),intent(out) :: OBJ
      real(kind(1.D0)),intent(in) :: GON(3)
      OBJ%GON(1:3)=GON(1:3)
      call CalRotMatrix(GON,OBJ%RLOC,OBJ%MLOC)
      end subroutine SetGonioAngles

!-----------------------------------------------------------
      subroutine CalRotMatrix(GON,RLOC,MLOC)
! Get rotation matrix for given segment angles
! RLOC=R3.R2.R1, where Ri are partical rotation matrices
! return MLOC=0 if there is no rotation, else MLOC=1
!-----------------------------------------------------------
      real(kind(1.D0)),intent(in) :: GON(3)
      real(kind(1.D0)),intent(out) :: RLOC(3,3)
      integer,intent(out) :: MLOC
      real(kind(1.D0)) :: R1(3,3),R2(3,3),R3(3,3),R(3,3)
      call UNIMAT(R1,3,3)
      MLOC=0

      !! this is equivalent to
      !call getRotYXY(GON,RLOC)
      !do i=1,3
      !  if (RLOC(i,i).NE.1.D0) MLOC=1
      !enddo
      !! but we keep following as more efficient

  ! RLOC=R3.R2.R1
      if (GON(1).ne.0.D0) then
        call MK_ROT3(2,GON(1),R1)
        MLOC=1
      endif
      if (GON(2).ne.0.D0) then
        call MK_ROT3(1,GON(2),R2)
        CALL M3XM3(1,R2,R1,R)
        MLOC=1
      else
        R=R1
      endif
      if (GON(3).ne.0.D0) then
        call MK_ROT3(2,GON(3),R3)
        CALL M3XM3(1,R3,R,RLOC)
        MLOC=1
      else
        RLOC=R
      endif
      end subroutine CalRotMatrix


!-----------------------------------------------------------
! Calculate exit axis angles of a component for given scattering angle.
! The direction of the exit axis is defined by THETA and PHI
! Incident axis is (0,0,1)
! Exit axis is defined by YXY Euler angles
! INPUT:
!   THETA scattering angle
!   PHI   rotation of Kf around Ki or sagital angle (depends on KFMODE)
!   KFMODE  .. defines the meaning PHI
!-----------------------------------------------------------
      subroutine CalAxisAngles(THETA,PHI,KFMODE,AX)
      real(kind(1.D0)),intent(in) :: THETA,PHI
      integer,intent(in) :: KFMODE
      real(kind(1.D0)),intent(inout) :: AX(3)
      real(kind(1.D0)) :: RPHI(3,3),RTHETA(3,3),R(3,3),CPHI,COMEGA
      select case(KFMODE)
! PHI is sagital angle (X-craddle)
      case(KFMODE_FLAT)
        CPHI=cos(PHI)
        if (CPHI.gt.0) then
          COMEGA=cos(THETA)/CPHI
          if (abs(COMEGA).le.1.D0) then
            AX=(/SIGN(1.D0,THETA)*acos(COMEGA),PHI,0.D0/)
          else
            AX=(/0.D0,THETA,0.D0/)
          endif
        else
          AX=(/THETA,0.D0,0.D0/)
        endif
! PHI is azimuthal angle (rotation of Kf around Ki)
      case(KFMODE_TOF)
        call MK_ROT3(3,PHI,RPHI)
        call MK_ROT3(2,THETA,RTHETA)
        call M3xM3(1,RTHETA,RPHI,R)
        call getEulerYXY(R,AX)
      end select
      end subroutine CalAxisAngles

!--------------------------------------------------------
! Calculate (y,x,y) Euler angles for a goniometer
! Serves e.g. for orienting crystals in diffraction position.
! We want to adjust single crystal in reflection position,
! i.e. Q=G. It can be done in multiple ways and therefore we
! have to define sutiable constraints.
! This transformation is performed in 2 steps:
!   1) orient G along z, using minimum tilt strategy
!   2) Adjust z along Q.
! To do (2), we define rotation from incident coordinates to a system with
! z' || Q. There are two suitable modes, which differ only if Q has a component
! out of the scattering plane.:
!   a) x-axis stays in the horizontal plane
!   b) x-axis stays in the scattering plane
! Incident axis is (0,0,1)
! INPUT:
!   Q(3)      .. scattering vector (in incident coord., only direction is important)
!   G(3)      .. a vector in component's local coord. to be aligned with Q
!   XStayHorizontal .. use mode a) if true
! OUTPUT:
!   GON(3) .. the (y,x,y) Euler angles of the goniometer
!--------------------------------------------------------
      subroutine CalGonioAngles(Q,G,XStayHorizontal,GON)
      real(kind(1.D0)),intent(in) :: Q(3),G(3)
      logical,intent(in) :: XStayHorizontal
      real(kind(1.D0)),intent(out) :: GON(3)
      real(kind(1.D0)) :: QQ,GG,AA,EE,A(3),E1(3,3),E(3,3),R(3,3)
      integer :: i,j,k,im
      REAL(KIND(1.D0)),parameter :: EPS=1.D-10
      character(128) :: S

! set default (no rotation)
      GON=0.D0
! Q and G must have non-zero length
      QQ=SQRT(Q(1)**2+Q(2)**2+Q(3)**2)
      GG=SQRT(G(1)**2+G(2)**2+G(3)**2)
      if ((QQ.LE.0.D0).or.(GG.LE.0.D0)) return
! E are the base vectors of a system with z' || G, x' || (y x G),
! in local coordinates. (Reminder: y is vertical in local coordinates)

      call V3cV3((/0.D0,1.D0,0.D0/),G,A)
      AA=SQRT(A(1)**2+A(2)**2+A(3)**2)
      if (AA.GT.EPS) then
        DO i=1,3
          E(i,1)=A(i)/AA
          E(i,3)=G(i)/GG
        enddo
        call V3cV3(E(1,3),E(1,1),E(1,2))
      else
      ! G is parallel to y => set x'=x, y'=-z, z'=y
        E(1:3,1)=(/1.D0,0.D0,0.D0/)
        E(1:3,2)=(/0.D0,0.D0,-1.D0/)
        E(1:3,3)=(/0.D0,1.D0,0.D0/)
      endif

! Adjust z along Q. To do so, define rotation from incident coordinates
! to a system with z' || Q. There are two suitable modes, which differ
! only if Q has a non-zero vertical component:
      if (XStayHorizontal) then
! a) Keep x horizontal. This mode is useful e.g.
! for a sample in flat-cone mode with lifted multianalyzer, which then
! scans in a plane at nearly constant height above the reference plane
        ! set E1 columns so that z"//Q, x"//(y x Q)
        E1(1:3,2)=(/0.D0,1.D0,0.D0/)
        E1(1:3,3)=Q(1:3)/QQ
        call V3cV3(E1(1,2),E1(1,3),E1(1,1))
        EE=SQRT(E1(1,1)**2+E1(2,1)**2+E1(3,1)**2)
        if (EE.lt.EPS) then
          S='Scattering vector can''t be vertical'
          im= ADD_MESSAGE(m_error,'CalGonioAngles',trim(S),m_low)
          return
        endif
        do i=1,3
          E1(i,1)=E1(i,1)/EE
        enddo
        call V3cV3(E1(1,3),E1(1,1),E1(1,2))
      else
! b) Keep x-axis in the scattering plane (Ki x Q). We set y' || (Ki x Q).SS.
! SS is the sign of (Ki x Q)(2), so that y' always points to the upper hemisphere.
! This mode should be used e.g. when orienting multianalyzers in horizontal plane.
        ! A is the normal to the scattering plane (Ki,Q)
        call V3cV3((/0.D0,0.D0,1.D0/),Q,A)
        AA=SIGN(1.D0,A(2))*SQRT(A(1)**2+A(2)**2+A(3)**2)
        ! if Q//z, then no rotation is needed
        if (abs(AA).lt.EPS) then
          call UNIMAT(E1,3,3)
        else
          DO i=1,3
            E1(i,2)=A(i)/AA
            E1(i,3)=Q(i)/QQ
          enddo
          call V3cV3(E1(1,2),E1(1,3),E1(1,1))
        endif
      endif
! Now construct the rotation matrix from Incident to Local coordinates
! of the crystal in reflection position. R=E.E1^T
      do i=1,3
        do j=1,3
          R(i,j)=0.D0
          do k=1,3
            R(i,j)=R(i,j)+E(i,k)*E1(j,k)
          enddo
        enddo
      enddo
! Get Euler angles from trans. matrix
      call getEulerYXY(R,GON)
      end subroutine CalGonioAngles

!-----------------------------------------------------------------
      SUBROUTINE GetKfAngles(KI0,KF0,Q0,SGN,PHI,KFMODE,AXIS,IERR)
! Get sample axis angles for given KI0,KF0,Q0,PHI,KFMODE
! 1) Rotates AX to match the scattering angle
! 2) Does not move gonio, only updates transformation matrices
! Following fields of OBJ are updated here:
! FRAME%AX,KI0,KF0,KI0,ESC + dependences in SETQSC
!------------------------------------------- ----------------------
      real(kind(1.D0)),intent(in) :: KI0,KF0,Q0,PHI
      integer,intent(in) :: SGN,KFMODE
      integer,intent(out) :: ierr
      real(kind(1.D0)),intent(out) :: AXIS(3)
      real(kind(1.D0)) :: CTHETA,STHETA,THETA
      ierr=1
      if ((KI0.gt.0).and.(KF0.gt.0)) then
        ierr=2
    ! check triangle
        if (Q0<=0.D0) then
          call CalAxisAngles(0.D0,PHI,KFMODE,AXIS)
          ierr=0
        else
          CTHETA=(KI0**2+KF0**2-Q0**2)/2.D0/KI0/KF0
          if (abs(CTHETA).lt.1.D0) then
            STHETA=SGN*SQRT(1.D0-CTHETA**2)
            THETA=atan2(STHETA,CTHETA)
            call CalAxisAngles(THETA,PHI,KFMODE,AXIS)
            ierr=0
          endif
        endif
      endif
      select case(ierr)
      case(1)
        call MSG_ERROR('GetKfAngles','Ki,Kf must be > 0',0,1)
      case(2)
        call MSG_ERROR('GetKfAngles','Can''t close scattering triangle',0,1)
      end select
      end SUBROUTINE GetKfAngles

!-------------------------------------------------
      SUBROUTINE FRAME_INIT_MAT(OBJ)
! initialize transformation matrices
! for beam deflecting components, call FRAME_DEFLECTION afterwards
!-------------------------------------------------
      TYPE(TFRAME),intent(out) :: OBJ
      REAL(KIND(1.D0)) :: Z
      ! set EXIT AXIS, (y,x,y) matrix (REXI) and translation (TEXI)
        call SetAxisAngles(OBJ,OBJ%AX)
      ! Centers of the incident and exit frames are identical.
      ! No deflection by default,
      ! But this can be changed for deflecting components calling FRAME_DEFLECTION
        OBJ%TEXI = 0.D0  ! centre of the exit reference frame (incident frame)
        OBJ%TDEFL = 0.D0 ! beam centre at the exit (incident frame)
        call UNIMAT(OBJ%RDEFL,3,3) ! deflection matrix
      ! set GONIOMETER, Euler craddle, (y,x,y) matrix
        call SetGonioAngles(OBJ,OBJ%GON)
        Z=abs(OBJ%MLOC)
        Z=Z+OBJ%STA(1)**2+OBJ%STA(2)**2+OBJ%STA(3)**2
        Z=Z+OBJ%VEL(1)**2+OBJ%VEL(2)**2+OBJ%VEL(3)**2
        OBJ%ISMOVED=0
        if (Z>1.D-10) OBJ%ISMOVED=1
        !write(*,*) 'FRAME_INIT_MAT ',trim(OBJ%ID),' ',OBJ%MEXI,OBJ%REXI(1,3)
      end SUBROUTINE FRAME_INIT_MAT

!-------------------------------------------------
      SUBROUTINE FRAME_DEFLECTION(OBJ, RHO)
! Add deflection due to bending to transformation matrices
! Used by objects like curved guides or neutron benders
! RHO(2) is the horizontal and vertical curvature in 1/mm
! recalculates OBJ%REXI and OBJ%TEXI so that transformation from
! incident to exit frame has the form:
! R' = REXI.(R - TEXI)
!-------------------------------------------------
      TYPE(TFRAME),intent(out) :: OBJ
      REAL(KIND(1.D0)),intent(in) :: RHO(2)
      REAL(KIND(1.D0)) :: ADEF(3), DV(3),DX,DY,ALPHA,BETA
      ! recalculate REXI, including deflection angles
        call SetAxisAngles(OBJ,OBJ%AX)
      ! get deflection angles
        ALPHA=RHO(1)*OBJ%SIZE(3)
        BETA=RHO(2)*OBJ%SIZE(3)
        ADEF = (/ALPHA,-BETA,0.D0/)
      ! correct REXI matrix for deflection
        call getRotYXY(OBJ%AX+ADEF,OBJ%REXI)
        OBJ%MEXI = 1
      ! calculate deflection matrix
        call getRotYXY(ADEF,OBJ%RDEFL)
      ! calculate deflection vector =  beam centre at the exit (incident frame)
        DX=0.5*RHO(1)*OBJ%SIZE(3)**2
        DY=0.5*RHO(2)*OBJ%SIZE(3)**2
        OBJ%TDEFL = (/DX,DY,OBJ%SIZE(3)/)
      ! get origin of the exit reference frame, shifted due to deflection (incident frame)
      ! TEXI = exit centre (TDEFL), minus length along exit axis converted to incident frame
        call M3XV3(-1, OBJ%RDEFL,(/0.D0,0.D0,OBJ%SIZE(3)/),DV)
        OBJ%TEXI = OBJ%TDEFL - DV
      end SUBROUTINE FRAME_DEFLECTION

!---------------------------------------------------------
      SUBROUTINE FRAME_INP_INST(INST,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      integer,intent(in) :: INST
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      TYPE(TFRAME),POINTER :: OBJ
      if (.not.FRAME_isValid(INST)) return
      OBJ => AFRAMES(INST)%X
      call FRAME_INP_OBJ(OBJ,PNAME,ARG,NARG)
      END SUBROUTINE FRAME_INP_INST

!---------------------------------------------------------
      SUBROUTINE FRAME_INP_OBJ(OBJ,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      TYPE(TFRAME),intent(out) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      LR=0
      NARG=0
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call FIELD2ARRAY(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
          call FRAME_INP_R(OBJ,trim(ARG%ID),NUM,LR)
        case default
          write(*,*) 'FRAME_INP: undefined field type '//trim(ARG%DEF%TID)
      end select
      NARG=LR
      END SUBROUTINE FRAME_INP_OBJ

!---------------------------------------------------------
      SUBROUTINE FRAME_OUT_INST(INST,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      integer,intent(in) :: INST
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      TYPE(TFRAME),POINTER :: OBJ
      if (.not.FRAME_isValid(INST)) return
      OBJ => AFRAMES(INST)%X
      call FRAME_OUT(OBJ,PNAME,ARG,NARG)
      END SUBROUTINE FRAME_OUT_INST


!---------------------------------------------------------
      SUBROUTINE FRAME_OUT_OBJ(OBJ,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      TYPE(TFRAME),intent(in) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      LR=0
      NARG=0
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call FRAME_OUT_R(OBJ,trim(ARG%ID),NUM,LR)
          call ARRAY2FIELD(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
        case default
          write(*,*) 'FRAME_OUT: undefined field type '//trim(ARG%DEF%TID)
      end select
      NARG=LR
      END SUBROUTINE FRAME_OUT_OBJ

C---------------------------------------------------------
      SUBROUTINE FRAME_INP_R(OBJ,PNAME,ARG,NARG)
C input data from ARG to OBJ for parameter namespace PNAME
C INPUT
C    PNAME   .... parameter name
C    ARG     .... REAL*8 parameter values
C INPUT
C    OBJ     .... SLIT structure
C    NARG    .... number of items read from ARG
C---------------------------------------------------------
      TYPE(TFRAME),intent(out) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(in) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR
      LR=1
      SELECT CASE (trim(PNAME))
C distance from previous element [mm]
        CASE('DIST');  OBJ%DIST=ARG(1);
C shape: ellipsoid(0), cylinder(1), disc(2),rectangle(3)
        CASE('SHAPE'); OBJ%SHAPE=NINT(ARG(1))
C order stage/gonio
        CASE('ORDER'); OBJ%ORDER=NINT(ARG(1))
        !write(*,*) trim(OBJ%ID),' inp ORDER=',OBJ%ORDER
C dimensions [mm]
        CASE('SIZE');  OBJ%SIZE(1:3)=ARG(1:3);LR=3
        CASE('SIZEX'); OBJ%SIZE(1)=ARG(1)
        CASE('SIZEY'); OBJ%SIZE(2)=ARG(1)
        CASE('SIZEZ'); OBJ%SIZE(3)=ARG(1)
C arm axes [deg] - convert to [rad]
        CASE('AX');    OBJ%AX(1:3)=ARG(1:3)*deg;LR=3
        CASE('AXX');   OBJ%AX(2)=ARG(1)*deg
        CASE('AXY');   OBJ%AX(1)=ARG(1)*deg
        CASE('AXZ');   OBJ%AX(3)=ARG(1)*deg
C linear stage [mm]
        CASE('STA');   OBJ%STA(1:3)=ARG(1:3);LR=3
        CASE('STAX');  OBJ%STA(1)=ARG(1)
        CASE('STAY');  OBJ%STA(2)=ARG(1)
        CASE('STAZ');  OBJ%STA(3)=ARG(1)
C goniometer [deg] - convert to [rad]
        CASE('GON');   OBJ%GON(1:3)=ARG(1:3)*deg;LR=3
        CASE('GONX');  OBJ%GON(1)=ARG(1)*deg
        CASE('GONY');  OBJ%GON(2)=ARG(1)*deg
        CASE('GONZ');  OBJ%GON(3)=ARG(1)*deg
C velocity [m/s] - convert to wavevector
        CASE('VEL');   OBJ%VEL(1:3)=ARG(1:3)/hovm/1000;LR=3
        CASE('VELX');  OBJ%VEL(1)=ARG(1)/hovm/1000
        CASE('VELY');  OBJ%VEL(2)=ARG(1)/hovm/1000
        CASE('VELZ');  OBJ%VEL(3)=ARG(1)/hovm/1000
        CASE DEFAULT; LR=0
      END SELECT
      NARG=LR
      END SUBROUTINE FRAME_INP_R


C---------------------------------------------------------
      SUBROUTINE FRAME_OUT_R(OBJ,PNAME,ARG,NARG)
C output data from OBJ to ARG for parameter namespace PNAME
C INPUT
C    PNAME   .... parameter name
C    ARG     .... REAL*8 parameter values
C INPUT
C    OBJ     .... SLIT structure
C    NARG    .... number of items read from ARG
C---------------------------------------------------------
      TYPE(TFRAME),intent(in) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(out) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR
      LR=1
!      write(*,*) 'FRAME_OUT [',trim(PNAME),']'
      SELECT CASE (trim(PNAME))
C distance from previous element [mm]
          CASE('DIST');  ARG(1)=OBJ%DIST;
C shape: ellipsoid(0), cylinder(1), disc(2),rectangle(3)
          CASE('SHAPE')
            ARG(1)=OBJ%SHAPE
!
C order stage/gonio
          CASE('ORDER'); ARG(1)=OBJ%ORDER
          !write(*,*) trim(OBJ%ID),' out ORDER=',OBJ%ORDER
C dimensions [mm]
          CASE('SIZE');  ARG(1:3)=OBJ%SIZE(1:3);LR=3
          CASE('SIZEX'); ARG(1)=OBJ%SIZE(1)
          CASE('SIZEY'); ARG(1)=OBJ%SIZE(2)
          CASE('SIZEZ'); ARG(1)=OBJ%SIZE(3)
C arm axes [deg] - convert to [rad]
          CASE('AX');    ARG(1:3)=OBJ%AX(1:3)/deg;LR=3
          CASE('AXX');   ARG(1)=OBJ%AX(2)/deg
          CASE('AXY');   ARG(1)=OBJ%AX(1)/deg
          CASE('AXZ');   ARG(1)=OBJ%AX(3)/deg
C linear stage [mm]
          CASE('STA');   ARG(1:3)=OBJ%STA(1:3);LR=3
          CASE('STAX');  ARG(1)=OBJ%STA(1)
          CASE('STAY');  ARG(1)=OBJ%STA(2)
          CASE('STAZ');  ARG(1)=OBJ%STA(3)
C goniometer [deg] - convert to [rad]
          CASE('GON');   ARG(1:3)=OBJ%GON(1:3)/deg;LR=3
          CASE('GONX');  ARG(1)=OBJ%GON(1)/deg
          CASE('GONY');  ARG(1)=OBJ%GON(2)/deg
          CASE('GONZ');  ARG(1)=OBJ%GON(3)/deg
C velocity [m/s] - convert to wavevector
          CASE('VEL');   ARG(1:3)=OBJ%VEL(1:3)*hovm*1000;LR=3
          CASE('VELX');  ARG(1)=OBJ%VEL(1)*hovm*1000
          CASE('VELY');  ARG(1)=OBJ%VEL(2)*hovm*1000
          CASE('VELZ');  ARG(1)=OBJ%VEL(3)*hovm*1000
          CASE DEFAULT; LR=0
      END SELECT
!      if (LR.gt.0) write(*,*) 'FRAME_OUT [',trim(PNAME),'] = ',ARG(1:LR)
!      read(*,*)
      NARG=LR
      END SUBROUTINE FRAME_OUT_R


!------------------------------------------------------------
      SUBROUTINE FRAME_LOCAL(IT,OBJ,NEU)
! transform neutron coordinates between incident axis and local reference frames
! IT=-1 ... local -> axis
! else    ... axis -> local
! MUST remain in this module due to dependences
!------------------------------------------------------------
      use TRACINGDATA,only:NEUTRON,trace_cnt
      TYPE(NEUTRON) :: NEU
      TYPE(TFRAME) :: OBJ
      INTEGER,intent(in) :: IT
      REAL(KIND(1.D0)) :: V(3),VK(3)
      logical :: dbg=.false.
1     format(a,': ',7(1x,G10.4))
      !dbg=((trace_cnt>1990).and.(OBJ%ISMOVED>0))
      !if (OBJ%ISMOVED==0) return
      if (dbg) write(*,1) 'FRAME_LOCAL '//trim(OBJ%ID), IT
      if (dbg) write(*,1) '  in ',NEU%R, NEU%K, NEU%T

      select case(OBJ%ORDER)
      case(top_stage)
      ! stage on top of the gonio
        select case(IT)
        case(-1)
          CALL V3AV3(1,NEU%R,OBJ%STA+OBJ%DM,V)
          CALL V3AV3(1,NEU%K,OBJ%VEL,VK)  ! add m/h*velocity
          CALL M3XV3(-OBJ%MLOC,OBJ%RLOC,V,NEU%R)
          CALL M3XV3(-OBJ%MLOC,OBJ%RLOC,VK,NEU%K)
        case default
          CALL M3XV3(OBJ%MLOC,OBJ%RLOC,NEU%R,V)
          CALL M3XV3(OBJ%MLOC,OBJ%RLOC,NEU%K,VK)
          CALL V3AV3(-1,V,OBJ%STA-OBJ%DM,NEU%R)
          CALL V3AV3(-1,VK,OBJ%VEL,NEU%K) ! subtract m/h*velocity
        end select
      case default
      ! gonio on top of the stage
        select case(IT)
        case(-1)
          CALL M3XV3(-OBJ%MLOC,OBJ%RLOC,NEU%R,V)
          CALL M3XV3(-OBJ%MLOC,OBJ%RLOC,NEU%K,VK)
          CALL V3AV3(1,V,OBJ%STA+OBJ%DM,NEU%R)
          CALL V3AV3(1,VK,OBJ%VEL,NEU%K)  ! add m/h*velocity
        case default
          CALL V3AV3(-1,NEU%R,OBJ%STA-OBJ%DM,V)
          CALL V3AV3(-1,NEU%K,OBJ%VEL,VK) ! subtract m/h*velocity
          CALL M3XV3(OBJ%MLOC,OBJ%RLOC,V,NEU%R)
          CALL M3XV3(OBJ%MLOC,OBJ%RLOC,VK,NEU%K)
        end select
      end select
      if (dbg) write(*,1) '   out',NEU%R, NEU%K, NEU%T
      END SUBROUTINE FRAME_LOCAL

!---------------------------------------------------------------------
      subroutine Loc2Glob(OBJ,NEU,NGLOB)
! Convert an event from OBJ'th local to global coordinates:
!---------------------------------------------------------------------
      use TRACINGDATA,only:NEUTRON
      TYPE(TFRAME) :: OBJ
      TYPE(NEUTRON) :: NEU,NGLOB
      TYPE(NEUTRON) :: NEU1
        NEU1=NEU
      ! from local to incident axis
        call FRAME_LOCAL(-1,OBJ,NEU1)
      ! from OBJ'th incident axis to global coordinates
        call Ax2Glob(OBJ,NEU1,NGLOB)
      end subroutine Loc2Glob


!---------------------------------------------------------------------
      subroutine Ax2Glob(OBJ,NEU,NGLOB)
! Convert an event from OBJ'th incident axis to global coordinates:
!---------------------------------------------------------------------
      use TRACINGDATA,only:NEUTRON
      TYPE(TFRAME) :: OBJ
      TYPE(NEUTRON) :: NEU,NGLOB
        NGLOB=NEU
      ! rotate to lab. coord.
        call M3XV3(1,OBJ%RLAB,NEU%R,NGLOB%R)
        call M3XV3(1,OBJ%RLAB,NEU%K,NGLOB%K)
      ! transition to lab. coord.
        NGLOB%R = NGLOB%R + OBJ%TLAB
      end subroutine Ax2Glob

!---------------------------------------------------------------------
      subroutine Loc2GlobR(OBJ,R,RGL)
! Convert position from OBJ'th local to global coordinates:
!---------------------------------------------------------------------
      TYPE(TFRAME) :: OBJ
      real(KIND(1.D0)),intent(in) :: R(3)
      real(KIND(1.D0)),intent(out) :: RGL(3)
      real(KIND(1.D0)) :: V(3),V1(3)
      ! from local to incident axis
        select case(OBJ%ORDER)
        case(top_stage)
        ! stage on top of the gonio
          CALL V3AV3(1,R,OBJ%STA+OBJ%DM,V)
          CALL M3XV3(-OBJ%MLOC,OBJ%RLOC,V,V1)
        case default
        ! gonio on top of the stage
          CALL M3XV3(-OBJ%MLOC,OBJ%RLOC,R,V)
          CALL V3AV3(1,V,OBJ%STA+OBJ%DM,V1)
        end select
      ! rotate to lab. coord.
        call M3XV3(1,OBJ%RLAB,V1,V)
      ! transition to lab. coord.
        RGL=V+OBJ%TLAB
      end subroutine Loc2GlobR

!---------------------------------------------------------------------
      subroutine Acc2Loc(OBJ,G,GLOC)
! Convert acceleration vector from global to local coordinates
!---------------------------------------------------------------------
      TYPE(TFRAME) :: OBJ
      real(KIND(1.D0)) :: G(3),GLOC(3)
      real(KIND(1.D0)) :: G1(3)
      ! from laboratory to incident frame
        call M3XV3(-1,OBJ%RLAB,G,G1)
      ! from incident to local frame
        CALL M3XV3(OBJ%MLOC,OBJ%RLOC,G1,GLOC)
      end subroutine Acc2Loc


!------------------------------------------------------------
      SUBROUTINE FRAME_LOCAL_POS(IT,OBJ,R)
! transform neutron coordinates between incident axis and local reference frames
! IT=-1 ... local -> axis
! else    ... axis -> local
!------------------------------------------------------------
      REAL(KIND(1.D0)) :: R(3)
      TYPE(TFRAME) :: OBJ
      INTEGER,intent(in) :: IT
      REAL(KIND(1.D0)) :: V(3)
      select case(OBJ%ORDER)
      case(top_stage)
      ! stage on top of the gonio
        select case(IT)
        case(-1)
          CALL V3AV3(1,R,OBJ%STA+OBJ%DM,V)
          CALL M3XV3(-OBJ%MLOC,OBJ%RLOC,V,R)
        case default
          CALL M3XV3(OBJ%MLOC,OBJ%RLOC,R,V)
          CALL V3AV3(-1,V,OBJ%STA-OBJ%DM,R)
        end select
      case default
      ! gonio on top of the stage
        select case(IT)
        case(-1)
          CALL M3XV3(-OBJ%MLOC,OBJ%RLOC,R,V)
          CALL V3AV3(1,V,OBJ%STA+OBJ%DM,R)
        case default
          CALL V3AV3(-1,R,OBJ%STA-OBJ%DM,V)
        CALL M3XV3(OBJ%MLOC,OBJ%RLOC,V,R)
        end select
      end select
      END SUBROUTINE FRAME_LOCAL_POS


!-----------------------------------------------------
! Get the component center coordinates.
! in incident axis frame:
! X,Y,Z = STA(1),STA(2),DIST
! if ORDER=top_stage, then STA has to be first converted to the
! axis coordinates.
!-----------------------------------------------------
      subroutine FRAME_GET_CTR(OBJ,X,Y,DIST)
      TYPE(TFRAME),intent(inout) :: OBJ
      real(KIND(1.D0)),intent(out) :: X,Y,DIST
      real(KIND(1.D0)) :: R(3)
1     format(a,': ',6(1x,G12.5))
      select case (OBJ%ORDER)
      case(top_stage)
        R=OBJ%STA
        CALL M3XV3(-OBJ%MLOC,OBJ%RLOC,OBJ%STA,R)
        X=R(1)
        Y=R(2)
        DIST=OBJ%DIST+R(3)
      case default
        X=OBJ%STA(1)
        Y=OBJ%STA(2)
        DIST=OBJ%DIST+OBJ%STA(3)
      end select
      !write(*,1) 'GET_CTR_SHIFT '//trim(OBJ%ID),X,Y,DIST
      end subroutine FRAME_GET_CTR


!-------------------------------------------------------------
! Update transformation TLAB(3), RLAB(3,3) for conversion to laboratory system.
! Given position R in the OBJ's incident axis frame, conversion to lab. frame is
! RL = OBJ%RLAB . R + OBJ%TLAB
! This subr. must be called in sequence so that RLAB and TLAB contain
! values corresponding to the preceding component
! OBJ%TLAB: transition to OBJ entry centre (=incident frame origin) in laboratory system
! OBJ%RLAB: rotation from incident frame to laboratory frame
!-------------------------------------------------------------
      subroutine FRAME_UPDATE_LAB(OBJ,RLAB,TLAB)
      TYPE(TFRAME),intent(out) :: OBJ
      REAL(KIND(1.D0)),intent(inout) :: RLAB(3,3),TLAB(3)
      REAL(KIND(1.D0)) :: T(3),R(3,3)
      integer :: i,j,k

  ! Rotation matrix as calculated by preceding call
        OBJ%RLAB=RLAB
  ! for translation, add path to entry from preceding component
  ! TLAB = TLAB + RLAB^T . (0,0,DIST)
        T=(/0.D0,0.D0,OBJ%DIST/)
        do i=1,3
          do k=1,3
            TLAB(i)=TLAB(i)+RLAB(i,k)*T(k)
          enddo
        enddo
        OBJ%TLAB = TLAB

!1       format(a,': ',6(1x,G13.7))
        !write(*,1) 'FRAME_UPDATE_LAB '//trim(OBJ%ID), OBJ%TLAB, OBJ%RLAB(1,3)

  ! Calculate rotation matrix for the next component
  ! RLAB = RLAB.REXI^T
        if (OBJ%MEXI.NE.0) then
          do i=1,3
            do j=1,3
              RLAB(i,j)=0.D0
              do k=1,3
                RLAB(i,j) = RLAB(i,j) + OBJ%RLAB(i,k)*OBJ%REXI(j,k)
              enddo
            enddo
          enddo
        endif
  ! Add shift of the exit frame origin due to deflection
        if (OBJ%ISDEFL) then
      ! convert TEXI to the laboratory coordinates and add it to TLAB
          call M3XV3(1,OBJ%RLAB,OBJ%TEXI,T)
          TLAB = TLAB + T
          !write(*,1) '        TDEFL', OBJ%TDEFL
          !write(*,1) '   added TEXI', OBJ%TEXI
          !write(*,1) '     lab TEXI', T
          !write(*,1) '     new TLAB', TLAB

      ! Already included in REXI
      ! add deflection to the rotation matrix
      ! RLAB=RLAB.RDEFL^T
      !    R=RLAB
      !    do i=1,3
      !      do j=1,3
      !        RLAB(i,j)=0.D0
      !        do k=1,3
      !          RLAB(i,j)=RLAB(i,j)+R(i,k)*OBJ%RDEFL(j,k)
      !        enddo
      !      enddo
      !    enddo
        endif
      end subroutine FRAME_UPDATE_LAB

      end module FRAMES
