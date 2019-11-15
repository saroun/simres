!//////////////////////////////////////////////////////////////////////
!////  $Id: array3d.f,v 1.14 2019/08/15 15:02:08 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.14 $
!////     $Date: 2019/08/15 15:02:08 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Describes tracing through a 3D array of rectangular sectors
!////////////////////////////////////////////////////////////////////////
      MODULE ARRAY3D
      use TRACINGDATA
      use FRAMES
      use SIMATH
      USE TRACELIB
      implicit none
      private


! FOCARRAY
      TYPE TARRAY3D
        integer :: NS(3)            ! number of segments
        real(kind(1.D0)) :: SA(3)   ! array envelope size
        real(kind(1.D0)) :: SZ(3)   ! single slab size
        real(kind(1.D0)) :: DS(3)   ! gaps
        real(kind(1.D0)) :: RHO(3)  ! curvature
        real(kind(1.D0)) :: TGSTACK ! tangent of the stacking angle (stack shift DX = TGSTACK*DZ)
        logical :: STACKH,STACKV ! stacking is flat (0) or curved (1)
        logical :: BENTH,BENTV ! segments are flat (0) or bent (1)
      END TYPE TARRAY3D

      logical, save :: dbg=.false.
      integer,SAVE :: DBGIT

      public TARRAY3D,ARRAY3D_INIT,ARRAY3D_ENTRY,ARRAY3D_SEG_TRANS
      public NEUT_TO_SEGMENT,NEUT_TO_ARRAY, R_TO_SEGMENT, RK_TO_SEGMENT,R_TO_ARRAY,ARRAY3D_BORDERS_NEW
      public ARRAY3D_GET_INDEX
      ! public ARRAY3D_BORDERS_DBG
      public ARRAY3D_SETDBG
      contains


      subroutine ARRAY3D_SETDBG(db)
      logical :: db
        dbg=db
      end subroutine ARRAY3D_SETDBG

!---------------------------------------------------------------
      subroutine ARRAY3D_INIT(A3D,SZ,GAP,NS,RHO,TGSTACK,SH,SV,BH,BV)
!---------------------------------------------------------------
      TYPE(TARRAY3D),intent(out) :: A3D
      real(kind(1.D0)),intent(in) :: SZ(3),GAP(3),RHO(3),TGSTACK
      integer,intent(in) :: NS(3)
      logical,intent(in) :: SH,SV,BH,BV
      integer :: i
      real(kind(1.D0)) :: DZ
      A3D%NS=NS
      A3D%SZ=SZ
      A3D%DS=GAP
      A3D%RHO=RHO
      A3D%TGSTACK=TGSTACK
      A3D%STACKH=SH
      A3D%STACKV=SV
      A3D%BENTH=BH
      A3D%BENTV=BV
      ! calculate the bounding box
      do i=1,3
        A3D%SA(i)=NS(i)*(A3D%SZ(i)+A3D%DS(i))
      enddo
      DZ=0.D0
      ! add z-displacement due to stack curvature
      if (A3D%STACKH) then
        DZ=DZ+abs(0.5*A3D%RHO(1)*(0.5*A3D%SA(1))**2)
      else
        DZ=DZ+0.5*abs(A3D%RHO(1)*A3D%SA(1)*A3D%SZ(1))
      endif
      if (A3D%STACKV) then
        DZ=DZ+abs(0.5*A3D%RHO(2)*(0.5*A3D%SA(2))**2)
      else
        DZ=DZ+0.5*abs(A3D%RHO(2)*A3D%SA(2)*A3D%SZ(2))
      endif
      if (A3D%RHO(3).ne.0.D0) DZ=DZ+0.5*abs(A3D%RHO(3)*A3D%SA(3)*A3D%SZ(1))

      ! add z-displacement due to segment bending
      if (A3D%BENTH) then
        DZ=DZ+abs(0.5*A3D%RHO(1)*(0.5*A3D%SZ(1))**2)
      endif
      if (A3D%BENTV) then
        DZ=DZ+abs(0.5*A3D%RHO(2)*(0.5*A3D%SZ(2))**2)
      endif
      ! make the bounding box symmetric w.r.t. z=0 and large enough ...
      A3D%SA(3)=A3D%SA(3)+2.D0*DZ
      DBGIT=0
      end subroutine ARRAY3D_INIT

!-----------------------------------------------------------------------
      logical function ARRAY3D_ENTRY(A3D,R,K,T)
! Find position (RIN) and time (T) to the entry of the crystal array space.
! This space is defined as a bounding box to the crystal array, including gaps and bending.
!------------------------------------------------------------------------
      TYPE(TARRAY3D),intent(in) :: A3D
      REAL(kind(1.D0)),intent(in) :: R(3),K(3)
      REAL(kind(1.D0)),intent(out) :: T
      REAL(kind(1.D0)) :: tout
    ! find cross-section times with the bounding box A3D%SA
      call BORDER_BOX(R,K,A3D%SA,T,tout)
      ARRAY3D_ENTRY=(T.lt.tout)
      end function ARRAY3D_ENTRY

!---------------------------------------------------------------------------
      SUBROUTINE ARRAY3D_GET_INDEX(A3D,r,idx)
! Calculate segment indexes from given position.
! Account for array curvature and z-stacking angle
!---------------------------------------------------------------------------
      TYPE(TARRAY3D),intent(in) :: A3D
      REAL(kind(1.D0)),intent(in) :: r(3)
      integer,intent(out) :: idx(3)
      REAL(kind(1.D0)) :: curv(2)
      REAL(kind(1.D0)) :: dd(3),dz,dx
      integer :: i
1     format(a,': ',6(1x,G12.5))

      ! size per segment (dd) incl., gaps
      do i=1,3
        dd(i)=A3D%SZ(i)+A3D%DS(i)
      enddo
      curv=0.D0
      if (A3D%BENTH.or.A3D%STACKV) curv(1)=A3D%RHO(1)
      if (A3D%BENTV.or.A3D%STACKV) curv(2)=A3D%RHO(2)
      ! subtract position of the neutral surface
      dz=0.5*(curv(1)*r(1)**2+curv(2)*r(2)**2)
      idx(2)=NINT(r(2)/dd(2)+0.5D0*(A3D%NS(2)+1))
      idx(3)=NINT((r(3)-dz)/dd(3)+0.5D0*(A3D%NS(3)+1))
      dx=idx(3)*dd(3)*A3D%TGSTACK
      idx(1)=NINT((r(1)-dx)/dd(1)+0.5D0*(A3D%NS(1)+1))
      END SUBROUTINE ARRAY3D_GET_INDEX

!---------------------------------------------------------------------------
      SUBROUTINE ARRAY3D_INDEX_SEGMENTS(A3D,r,k,imi,ima)
! Find index ranges of the sectors crossed by the trajectory (r,k)
! No backtracing, segments before the current position are neglected.
! No return to the segments alread crossed:
! assumes monotonous sequence of crossed segments indexes
!---------------------------------------------------------------------------
      TYPE(TARRAY3D),intent(in) :: A3D
      REAL(kind(1.D0)),intent(in) :: r(3),k(3)
      integer,intent(out) :: imi(3),ima(3)
      REAL(kind(1.D0)) :: tin,tout
      REAL(kind(1.D0)) :: r1(3),r2(3),curv(2)
      REAL(kind(1.D0)) :: dd(3),da(3)
      integer :: i,idx1(3),idx2(3)
1     format(a,': ',6(1x,G12.5))

      ! size per segment (dd) and total (da) incl., gaps
      do i=1,3
        dd(i)=A3D%SZ(i)+A3D%DS(i)
        da(i)=A3D%NS(i)*dd(i)
      enddo
      !call BORDER_BOX(r,k,A3D%SA,tin,tout)
      curv=0.D0
      if (A3D%BENTH.or.A3D%STACKV) curv(1)=A3D%RHO(1)
      if (A3D%BENTV.or.A3D%STACKV) curv(2)=A3D%RHO(2)

      if (dbg) write(*,1) 'ARRAY3D_INDEX_SEGMENTS ',r,k
      ! get entry and exit times
      ! call BORDER_WAFER(r,k,da,rhoX,rhoY,tin,tout)
      call BORDER_SANDWICH(r,k,(/A3D%TGSTACK,0.D0/),curv,da,tin,tout,.false.)

      if (dbg) write(*,1) '              tin,tout ',tin,tout
      if (tin>tout) return
      ! get contact points for entry and exit surfaces
      do i=1,3
        r1(i)=r(i)+tin*k(i)
        r2(i)=r(i)+tout*k(i)
      enddo
      call ARRAY3D_GET_INDEX(A3D,r1,idx1)
      call ARRAY3D_GET_INDEX(A3D,r2,idx2)
      do i=1,3
        imi(i)=max(1,min(idx1(i),idx2(i)))
        ima(i)=min(A3D%NS(i),max(idx1(i),idx2(i)))
      enddo

      if (dbg) write(*,1) '      r1', r1, idx1
      if (dbg) write(*,1) '      r2', r2, idx2
      if (dbg) write(*,1) '      imi', imi
      if (dbg) write(*,1) '      ima', ima

      END SUBROUTINE ARRAY3D_INDEX_SEGMENTS


!---------------------------------------------------------------------
      SUBROUTINE ARRAY3D_BORDERS_NEW(A3D,R,K,TSEG0,TSEG,ISEG,NSEG)
! Traces through all segments along the neutron path (along k)
! Calculate the segments entry and exit times
! Exclude negative times (crystal area behind the current position is ignored)
! INPUT:
!   R,K starting position and K-vector
! OUTPUT
!   TSEG    times to the I-th segment exit, started at the assembly entry
!   TSEG0   times to the I-th segment entry, -"-
!   NSEG    number of crossed segments
!   ISEG    segments coordinates, dimension must be at least (3,NSEG)
!---------------------------------------------------------------------
      TYPE(TARRAY3D),intent(in) :: A3D
      REAL(KIND(1.D0)),intent(in) :: R(3),K(3)
      REAL(KIND(1.D0)),intent(out) :: TSEG(0:),TSEG0(0:)
      integer,intent(out) :: ISEG(:,:),NSEG
      INTEGER :: I,J,M,IS,MSEG,MLOC
      REAL(KIND(1.D0)) :: t1,t2,curv(2)
      REAL(KIND(1.D0)) :: R1(3),K1(3),TR(3),RLOC(3,3)
      REAL(KIND(1.D0)) :: TS1(256),TS2(256)
      integer :: ORDSEG(256),IDX(1:3,256),imi(3),ima(3)
1     format(a,6(1x,G12.5))
      IS=0
      NSEG=0
      TSEG(0)=0.D0
      TSEG0(0)=0.D0
      MSEG=SIZE(TSEG)-1
      curv=0.D0
      if (A3D%BENTH) curv(1)=A3D%RHO(1)
      if (A3D%BENTV) curv(2)=A3D%RHO(2)

      !DBGIT=DBGIT+1
      !dbg=(DBGIT<20)


      call ARRAY3D_INDEX_SEGMENTS(A3D,r,k,imi,ima)
      if (dbg) write(*,1) 'ARRAY3D_BORDERS_NEW imi',imi
      if (dbg) write(*,1) 'ARRAY3D_BORDERS_NEW ima',ima
      do i=imi(1),ima(1)
        do j=imi(2),ima(2)
          do m=imi(3),ima(3)
            call ARRAY3D_SEG_TRANS(A3D,(/i,j,m/),TR,RLOC,MLOC)
            call RK_TO_SEGMENT(TR,RLOC,MLOC,R,K,R1,K1)
            ! call BORDER_BOX(R1,K1,A3D%SZ,t1,t2)
            !if (dbg) then
            !  call BORDER_WAFER_DBG(R1,K1,curv,A3D%SZ,t1,t2)
            !else
              call BORDER_WAFER(R1,K1,curv,A3D%SZ,t1,t2)
            !endif
            if (dbg) write(*,1) '   i,j,m,is',i,j,m,is
            if (dbg) write(*,1) '   R,K loc',R1,K1
            if (dbg) write(*,1) '   t1,t2',t1,t2
            if ((IS<MSEG).and.(t2>=0.D0).and.(t2>t1)) then
              IS=IS+1
              TS1(IS)=max(0.D0,t1)
              TS2(IS)=t2
              IDX(1:3,IS)=(/i,j,m/)
            endif
          enddo
        enddo
      enddo
      NSEG=IS
      if (dbg) write(*,1) '     NSEG',NSEG
      if (NSEG>1) then
        call HEAPSORT(TS1,NSEG,1,ORDSEG)
        do J=1,NSEG
          TSEG0(ORDSEG(J))=TS1(J)
          TSEG(ORDSEG(J))=TS2(J)
          ISEG(1:3,ORDSEG(J))=IDX(1:3,J)
        enddo
      else if (NSEG==1) then
        TSEG0(1)=TS1(1)
        TSEG(1)=TS2(1)
        ISEG(1:3,1)=IDX(1:3,1)
      endif
      if (dbg) then
        do J=1,NSEG
          if (dbg) write(*,1) '   J,TIN,TOUT,ISEG',J,K(3)*TSEG0(J),K(3)*TSEG(J),ISEG(1:3,J)
        enddo
      endif

      END SUBROUTINE ARRAY3D_BORDERS_NEW

!---------------------------------------------------------------------
      SUBROUTINE ARRAY3D_SEG_TRANS(A3D,ISEG,TR,RLOC,MLOC)
! get transformation data (translation and rot. matrix) for given segment
! the matrix and translation convert neutron coordinates from array local
! to segment local system
!---------------------------------------------------------------------
      TYPE(TARRAY3D),intent(in) :: A3D
      integer,intent(in) :: ISEG(3)
      real(kind(1.D0)),intent(out) :: TR(3),RLOC(3,3)
      integer,intent(out) :: MLOC
      real(kind(1.D0)) :: A(3),DX,DY,DZ,Z
      integer :: i
1     format(a,' :',6(1x,G12.5))
! get coordinates of the segment center
      do i=1,3
        TR(i)=(A3D%SZ(i)+A3D%DS(i))*(ISEG(i)-0.5D0*(A3D%NS(i)+1.D0))
      enddo
      Z=ABS(A3D%RHO(1))+ABS(A3D%RHO(2))+ABS(A3D%RHO(3))
! add stacking angle
      TR(1)=TR(1)+TR(3)*A3D%TGSTACK
      if (Z>0.D0) then
! get tilt angles
        A(1)=-TR(1)*A3D%RHO(1)
        A(2)=TR(2)*A3D%RHO(2)
        A(3)=TR(3)*A3D%RHO(3)
        DX=TR(3)*A(1)
        TR(1)=TR(1)+DX
        DY=-TR(3)*A(2)
        TR(2)=TR(2)+DY
! add z-shift for curved segment stacking
        if ((A3D%STACKH).and.(abs(A3D%RHO(1))>0.D0)) then
          DZ=0.5D0*A3D%RHO(1)*TR(1)**2
          TR(3)=TR(3)+DZ
        endif
        if ((A3D%STACKV).and.(abs(A3D%RHO(2))>0.D0)) then
          DZ=0.5D0*A3D%RHO(2)*TR(2)**2
          TR(3)=TR(3)+DZ
        endif
        ! write(*,1) 'ARRAY3D_SEG_TRANS ',ISEG(2),A(2),TR(2:3)

        call CalRotMatrix(A,RLOC,MLOC)
      else
        MLOC=0
        call UNIMAT(RLOC,3,3)
      endif

      end SUBROUTINE ARRAY3D_SEG_TRANS


!---------------------------------------------------------------------
      SUBROUTINE RK_TO_SEGMENT(TR,RLOC,MLOC,R0,K0,R,K)
! convert neutron k-vector from array to segment coordinates
! RLOC is the transformation matrix
!---------------------------------------------------------------------
      real(kind(1.D0)),intent(in) :: TR(3),RLOC(3,3),R0(3),K0(3)
      integer,intent(in) :: MLOC
      real(kind(1.D0)),intent(out) :: R(3),K(3)
      real(kind(1.D0)) :: V(3)
      if (MLOC==0) then
        CALL V3AV3(-1,R0,TR,R)
        K=K0
      else
        CALL V3AV3(-1,R0,TR,V)
        CALL M3XV3(MLOC,RLOC,V,R)
        CALL M3XV3(MLOC,RLOC,K0,K)
      endif
      end SUBROUTINE RK_TO_SEGMENT

!---------------------------------------------------------------------
      SUBROUTINE R_TO_SEGMENT(TR,RLOC,MLOC,R0,R)
! convert neutron position from array to segment coordinates
! TR and RLOC are corresponding transformation arrays (translation and rotation)
!---------------------------------------------------------------------
      real(kind(1.D0)),intent(in) :: TR(3),RLOC(3,3),R0(3)
      integer,intent(in) :: MLOC
      real(kind(1.D0)),intent(out) :: R(3)
      real(kind(1.D0)) :: V(3)
      if (MLOC==0) then
        CALL V3AV3(-1,R0,TR,R)
      else
        CALL V3AV3(-1,R0,TR,V)
        CALL M3XV3(MLOC,RLOC,V,R)
      endif
      end SUBROUTINE R_TO_SEGMENT

!---------------------------------------------------------------------
      SUBROUTINE R_TO_ARRAY(TR,RLOC,MLOC,R0,R)
! convert neutron position from segment to array coordinates
! TR and RLOC are corresponding transformation arrays (translation and rotation)
!---------------------------------------------------------------------
      real(kind(1.D0)),intent(in) :: TR(3),RLOC(3,3),R0(3)
      integer,intent(in) :: MLOC
      real(kind(1.D0)),intent(out) :: R(3)
      real(kind(1.D0)) :: V(3)
      if (MLOC==0) then
        CALL V3AV3(1,R0,TR,R)
      else
        CALL M3XV3(-MLOC,RLOC,R0,V)
        CALL V3AV3(1,V,TR,R)
      endif
      end SUBROUTINE R_TO_ARRAY

!---------------------------------------------------------------------
      SUBROUTINE NEUT_TO_SEGMENT(TR,RLOC,MLOC,NEUT0,NEUT1)
! convert neutron coordinates from array local (NEUT0)
! to segment local (NEUT1) coordinates
! TR and RLOC are corresponding transformation data
!---------------------------------------------------------------------
      real(kind(1.D0)),intent(in) :: TR(3),RLOC(3,3)
      integer,intent(in) :: MLOC
      TYPE(NEUTRON),intent(in) :: NEUT0
      TYPE(NEUTRON),intent(out) :: NEUT1
      real(kind(1.D0)) :: V(3)
      if (MLOC==0) then
        CALL V3AV3(-1,NEUT0%R,TR,NEUT1%R)
        NEUT1%K=NEUT0%K
      else
        CALL V3AV3(-1,NEUT0%R,TR,V)
        CALL M3XV3(MLOC,RLOC,V,NEUT1%R)
        CALL M3XV3(MLOC,RLOC,NEUT0%K,NEUT1%K)
      endif
      end SUBROUTINE NEUT_TO_SEGMENT

!---------------------------------------------------------------------
      SUBROUTINE NEUT_TO_ARRAY(TR,RLOC,MLOC,NEUT0,NEUT1)
! inverse to NEUT_TO_SEGMENT
!---------------------------------------------------------------------
      real(kind(1.D0)),intent(in) :: TR(3),RLOC(3,3)
      integer,intent(in) :: MLOC
      TYPE(NEUTRON),intent(in) :: NEUT0
      TYPE(NEUTRON),intent(out) :: NEUT1
      real(kind(1.D0)) :: V(3)
      if (MLOC==0) then
        CALL V3AV3(1,NEUT0%R,TR,NEUT1%R)
        NEUT1%K=NEUT0%K
      else
        CALL M3XV3(-MLOC,RLOC,NEUT0%R,V)
        CALL M3XV3(-MLOC,RLOC,NEUT0%K,NEUT1%K)
        CALL V3AV3(1,V,TR,NEUT1%R)
      endif
      end SUBROUTINE NEUT_TO_ARRAY


!==================================
! obsolete
!===============================



!---------------------------------------------------------------------------
      SUBROUTINE ARRAY3D_SEGCOORD(A3D,R,I0,SEGT,SEGR,SEGM)
! return indices I0 of the sector, in which a particle at R resides
! return transformation arrays for this sector
! if R is outside the array boundaries, return I0=0
!---------------------------------------------------------------------------
      use CLASSES
      TYPE(TARRAY3D),intent(in) :: A3D
      REAL(kind(1.D0)),intent(in) :: R(3)
      REAL(kind(1.D0)),intent(out)  :: SEGT(3),SEGR(3,3)
      INTEGER,intent(out) :: I0(3),SEGM
      REAL(kind(1.D0)) :: R1(3),SZ(3),TR(3),RLOC(3,3)
      integer :: i,j,m,MLOC
      logical :: loop,LOG1

    ! find cross-section times with the array boundary
      SZ=A3D%SZ+A3D%DS
      loop=.true.
      i=0
      j=0
      m=0
      I0=0
      do while (loop.and.(i<A3D%NS(1)))
        i=i+1
        do while (loop.and.(j<A3D%NS(2)))
          j=j+1
          do while (loop.and.(m<A3D%NS(3)))
            m=m+1
            call ARRAY3D_SEG_TRANS(A3D,(/i,j,m/),TR,RLOC,MLOC)
            call R_TO_SEGMENT(TR,RLOC,MLOC,R,R1)
            LOG1=TLIB_INFRAME(R1,SZ,FRAME_SHAPE_BOX)
            if (LOG1) then
              I0=(/i,j,m/)
              SEGT=TR
              SEGR=RLOC
              SEGM=MLOC
              loop=.false.
            endif
          enddo
        enddo
      enddo
      END SUBROUTINE ARRAY3D_SEGCOORD

!---------------------------------------------------------------------
      SUBROUTINE ARRAY3D_BORDERS_OBSOLETE(A3D,R,K,TSEG0,TSEG,ISEG,NSEG)
! Traces through all segments along the neutron path
! Calculate the segments entry and exit times
! Assume that R is just at (inside) the crystal area
! INPUT:
!   R,K starting position and K-vector
! OUTPUT
!   TSEG    times to the I-th segment exit, started at the assembly entry
!   TSEG0   times to the I-th segment entry, -"-
!   NSEG    number of crossed segments
!   ISEG    segments coordinates, dimension must be at least (3,NSEG)
!---------------------------------------------------------------------
      TYPE(TARRAY3D),intent(in) :: A3D
      REAL(KIND(1.D0)),intent(in) :: R(3),K(3)
      REAL(KIND(1.D0)),intent(out) :: TSEG(0:),TSEG0(0:)
      integer,intent(out) :: ISEG(:,:),NSEG
      INTEGER :: I,J,I0(3),MSEG,MLOC
      REAL(KIND(1.D0)) :: T,DT,TIN1,TIN2,TOUT1,TOUT2
      REAL(KIND(1.D0)) :: V(3),R1(3),K1(3),TR(3),RLOC(3,3)
      logical :: isInside
1     format(a,6(1x,G12.5))
      J=0
      NSEG=0
      TSEG(0)=0.D0
      TSEG0(0)=0.D0
      MSEG=SIZE(TSEG)-1

    ! move to the crystal entry. T measures TOF.
    ! if R is inside, stay there
      isInside=ARRAY3D_ENTRY(A3D,R,K,T)
      if (.not.isInside) return
! trace while inside the volume of the crystal array
      DO WHILE (isInside.AND.(J.LT.MSEG))

! Get entry (TIN) and exit (TOUT) times of a neutron moving along K and starting at V
! (1) .. crossing segment sector
! (2) .. crossing the segment itself (differs by gaps between the segments, CR.DH ...)
      ! get segment coordinates
        call ARRAY3D_SEGCOORD(A3D,V,I0,TR,RLOC,MLOC)
!      write(*,1) 'ARRAY3D_BORDERS I0=',J+1, I0
        if (I0(1)==0) then
          isInside=.false.
        else
          call RK_TO_SEGMENT(TR,RLOC,MLOC,V,K,R1,K1)
          call BORDER_BOX(R1,K1,A3D%SZ+A3D%DS,TIN1,TOUT1)
          call BORDER_BOX(R1,K1,A3D%SZ,TIN2,TOUT2)
! Count only intersected segments
    ! must start before or inside the segment
    ! exclude corners and gaps
!        write(*,1) 'ARRAY3D_BORDERS TIN, TOUT=',TIN1,TOUT1,TIN2,TOUT2
          IF (((TOUT2-TIN2).GT.1.D-3).AND.(TOUT2.GT.1.D-10)) THEN
            J=J+1
! Move to the segment entry, if not already inside
            IF (TIN2.GT.0.) THEN
              DO I=1,3
                V(I)=V(I)+TIN2*K(I)
              ENDDO
              T=T+TIN2
!              write(*,1) ' to entry ',T
            ELSE
              TIN2=0.D0
            ENDIF
            TSEG0(J)=T     ! J-th segment starting time (with respect to the orig. position R)
            TSEG(J)=T+TOUT2-TIN2  ! J-th segment exit time
            ISEG(1:3,J)=I0(1:3)   ! coordinates of the segment
          ELSE
            TIN2=0.D0
          ENDIF
! Move to the entry of the next segment sector
    ! move slightly behind, we need to start INSIDE the next sector
          DT=1.000001D0*(TOUT1-TIN2)
          DO I=1,3
            V(I)=V(I)+DT*K(I)
          ENDDO
          T=T+DT
        endif
      ENDDO
!      write(*,1) 'ARRAY3D_BORDERS NSEG=',J
      NSEG=J
      END SUBROUTINE ARRAY3D_BORDERS_OBSOLETE


      end MODULE ARRAY3D