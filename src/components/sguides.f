!//////////////////////////////////////////////////////////////////////
!////  $Id: sguides.f,v 1.33 2019/08/16 17:16:25 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2012, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.33 $
!////     $Date: 2019/08/16 17:16:25 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Describes component class: SGUIDE
!////  Free-shape segmented guide
!////
!////////////////////////////////////////////////////////////////////////
      MODULE SGUIDES
      use CONSTANTS
      use CLASSES
      use FIELDDATA
      use FRAMES
      use MIRROR_TABLE
      use MIRRLOG
      use RNDGEN
      implicit none

      integer,parameter, PRIVATE :: MAX_SEG=128
      integer,parameter, PRIVATE :: sguide_prof_man=0
      integer,parameter, PRIVATE :: sguide_prof_flat=1
      integer,parameter, PRIVATE :: sguide_prof_para=2
      integer,parameter, PRIVATE :: sguide_prof_ell=3

      TYPE  SGUIDE; SEQUENCE
        TYPE(TFRAME) :: FRAME
        real(kind(1.D0)) :: W(MAX_SEG),H(MAX_SEG),L(MAX_SEG)   ! width,height,length of segments
        real(kind(1.D0)) :: MC(4,MAX_SEG)  ! critical angles of segments (left, right, top, bottom)
        real(kind(1.D0)) :: RHO(2)         ! curvatures (H&V)
        real(kind(1.D0)) :: GAP  ! spacing between segments [mm]
        real(kind(1.D0)) :: WAV  ! waviness (sigma, [rad])
        real(kind(1.D0)) :: EXI(2)   ! exit width and height
        real(kind(1.D0)) :: MISALIGN(2)   ! segments misalignment in mm (gaussian sigma)
        INTEGER :: NSEG      ! number of segments
        LOGICAL :: SMOOTH(4)    ! smooth or segmented wall
        LOGICAL :: ACTIVE(4)    ! false if the wall does not exist
        LOGICAL :: MONITOR(4)   ! monitor neutrons captured by the wall
        LOGICAL :: CLOSED       ! if true, no neutrons can path through
        LOGICAL :: FRONTMASK    ! if true, clip the beam by entrance window
        LOGICAL :: GAPABS       ! neutrons can't leak through the gaps
        integer :: IPROF(2)     ! profile (h,v): (0) manual, (1) straight, (2) parabolic, (3) elliptic
        integer :: LOGBNC       ! log bounces at reflecting walls - produce map in BNCREC
! calculated fields
        integer :: MLOG   ! reference to a mirror logger  (see MIRRLOG unit)
        ! nominal wavelength (used to calculate mean critical angle)
        real(kind(1.D0)) :: LAMBDA
        ! indexes pointing to the field containg reflectivity data
        INTEGER :: NR(4,MAX_SEG)
        ! positions of segments windows (center)
        real(kind(1.D0)) ::  WIN(3,0:MAX_SEG)
        ! cubic coefficients for interpolation (pos = D(1)*p**3+D(2)*p**2+D(3)*p+D(4))
        ! where p = z-z0, z0=position of the segment entry = WIN(3,iseg)
        ! 1st index is for the coefficients D(1) .. D(4)
        ! 2nd index is for right,left, top, bottom walls
        real(kind(1.D0)) ::  D(4,4,MAX_SEG)
        ! generated misalignments
        real(kind(1.D0)) ::  DX(2,0:MAX_SEG)
        real(kind(1.D0)) ::  DA(2,MAX_SEG)
! dbg
        REAL(kind(1.D0)) :: dbg_SUM(4)
        integer :: dbg_NSUM(4)
      END TYPE SGUIDE


      integer, parameter :: SGUIDES_DIM=16
! pointer type for SGUIDE
      type PSGUIDE
        integer :: IDX  ! instance index
        TYPE(SGUIDE),pointer :: X
      end type PSGUIDE
! instances of SGUIDE. ASGUIDES(0) is always unallocated
      integer :: SGUIDES_NC
      type(PSGUIDE) :: ASGUIDES(1:SGUIDES_DIM)

      contains

!---------------------------------------------------------
      SUBROUTINE SGUIDE_MISALIGN(OBJ)
! generate random misalignment for each segment
C---------------------------------------------------------
      TYPE(SGUIDE) :: OBJ
      real(kind(1.D0)) ::  DR0(2),L
      integer :: i,iseg
      DR0=(/0.D0,0.D0/)
      do i=1,2
        if (OBJ%MISALIGN(i).ne.0.D0) then
          OBJ%DX(i,0)=OBJ%MISALIGN(i)*GASDEV1(0.D0,3.D0)
          Do iseg=1,OBJ%NSEG
            L=OBJ%WIN(3,iseg)-OBJ%WIN(3,iseg-1)
            OBJ%DX(i,iseg)=OBJ%MISALIGN(i)*GASDEV1(0.D0,3.D0)
            ! modified 2018-08-28: consider only spatial misalignment
            ! OBJ%DA(i,iseg)=(OBJ%DX(i,iseg)-DR0(i))/L
            OBJ%DA(i,iseg)=0.D0

            DR0(i)=OBJ%DX(i,iseg)
          enddo
        endif
      enddo
      END SUBROUTINE SGUIDE_MISALIGN


!---------------------------------------------------------
      subroutine SGUIDE_REGMLOG(INST)
! register logger for mirror in MIRRLOG
!---------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(SGUIDE),POINTER :: OBJ
1     FORMAT(a,1x,6(1x,G11.5))
      if (.not.SGUIDE_isValid(INST)) return
      OBJ => ASGUIDES(INST)%X
      if (OBJ%LOGBNC>0) then
        OBJ%MLOG=MLOGS_REGISTER(OBJ%FRAME%ID,CCLS_SGUIDE,INST,OBJ%LOGBNC,OBJ%FRAME%SIZE(3))
      else
        OBJ%MLOG=0
      endif
      end subroutine SGUIDE_REGMLOG

!---------------------------------------------------------
      subroutine SGUIDE_PROFILE(INST,IU)
! write guide profile in the given file unit
!---------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(SGUIDE),POINTER :: OBJ
      integer,intent(in) :: IU
      integer :: i,j,ISEG,NP,k
      real(kind(1.D0)) :: Z,DZ,X(4),LN,R(3),R1(3),R2(3),M(4),gap1,gap2,rhomax
      character(LEN_LINE) :: LINE,LINE2
      character(32) :: CNUM
      logical :: curved
      integer,parameter :: ISIDE(4)=(/1,1,2,2/)

      if (.not.SGUIDE_isValid(INST)) return
      OBJ => ASGUIDES(INST)%X
! print this only to STDOUT
      if (IU==6) then
          write(IU,*) '# SGUIDE '//trim(OBJ%FRAME%ID),' NSEG=',OBJ%NSEG
          do i=1,OBJ%NSEG
          call ARRAY2STR(OBJ%D(1:4,1,i),4,LINE)
          call ARRAY2STR(OBJ%D(1:4,3,i),4,LINE2)
          write(IU,*) '# SEG=',i,', DH=',trim(LINE),' DV=',trim(LINE2)
        enddo
      endif
      ! generate segments misalignment if required
      call SGUIDE_MISALIGN(OBJ)
      curved=.false.
      rhomax=0.D0
      rhomax=max(rhomax,abs(OBJ%RHO(1)),abs(OBJ%RHO(2)))
      do j=1,4
        if (OBJ%ACTIVE(j).and.OBJ%SMOOTH(j).and.(rhomax.gt.0.D0)) curved=.true.
      enddo
      LN=0.D0
      ISEG=1
      Do ISEG=1,OBJ%NSEG
        if (curved) then
          NP=MAX(5,nint(OBJ%L(ISEG)*SQRT(rhomax/0.02)))
        else
          NP=1
        endif
        DZ=OBJ%L(ISEG)/NP
        if (ISEG==1) then
          gap1=0.D0
          gap2=0.5D0*OBJ%GAP
        else if (ISEG==OBJ%NSEG) then
          gap1=0.5D0*OBJ%GAP
          gap2=0.D0
        else
          gap1=0.5D0*OBJ%GAP
          gap2=0.5D0*OBJ%GAP
        endif
        do i=0,NP
          Z=i*DZ
          do j=1,4
            X(j)=SGUIDE_POS(0,OBJ%D(1:4,j,ISEG),Z+gap1)
            k=ISIDE(j)
            X(j)=X(j)+OBJ%DX(k,iseg-1)+OBJ%DA(k,iseg)*(Z+gap1)
            M(j)=OBJ%MC(j,ISEG)
          enddo
! convert to lab. coordinates:
          Z=Z+LN+gap1
          R=(/X(1),X(3),Z/)
          call Loc2GlobR(OBJ%FRAME,R,R1)
          R=(/X(2),X(4),Z/)
          call Loc2GlobR(OBJ%FRAME,R,R2)
          X=(/R1(1),R2(1),R1(2),R2(2)/)
          R=(/0.D0,0.D0,Z/)
          call Loc2GlobR(OBJ%FRAME,R,R1)
          Z=R1(3)
          call ARRAY2STR(X,4,LINE)
          call ARRAY2STR(M,4,LINE2)
          call FLOAT2STR(Z,CNUM)
          write(IU,*) iseg,' ',trim(CNUM),' ',trim(LINE),' ',trim(LINE2)
        enddo
        LN=LN+OBJ%L(ISEG)+gap1+gap2
      enddo
      end subroutine SGUIDE_PROFILE

!!---------------------------------------------------------
!      SUBROUTINE SGUIDE_ADJUST_SEGMENTS(INST)
! adjust number and size of segments automatically if NSEG<0
!---------------------------------------------------------
!      integer,intent(in) :: INST
!      real(kind(1.D0)) :: Ltot,L
!      integer :: i
!      TYPE(SGUIDE),POINTER :: OBJ
!      if (.not.SGUIDE_isValid(INST)) return
!      OBJ => ASGUIDES(INST)%X
!      if (OBJ%NSEG==0) OBJ%NSEG=-MIN(30,MAX(3,NINT(OBJ%FRAME%SIZE(3))))
!      if (OBJ%NSEG<0) then
!        OBJ%NSEG=MIN(ABS(OBJ%NSEG),MAX_SEG)
!        Ltot=OBJ%FRAME%SIZE(3)+(OBJ%NSEG-1)*OBJ%GAP
!        L=OBJ%FRAME%SIZE(3)/OBJ%NSEG
!        do i=1,OBJ%NSEG
!          OBJ%L(i)=L
!        enddo
!      endif
!      end SUBROUTINE SGUIDE_ADJUST_SEGMENTS


!---------------------------------------------------------
      SUBROUTINE SGUIDE_ADJUST_SEGMENTS(OBJ, NS)
! if (NS<=0), set automatically the number of segments and related parameters:
! NS=0: divide the guide in default 0.5 long segmens
! NS<0: change the number of segments to abs(NS)
! NS=NSEG: do nothing
! if NS<>NSEG: use the first segment to define the segment parameters.
! Adjust segment length so that the total length, OBJ%FRAME%SIZE(3) including gaps is preserved
!---------------------------------------------------------
      TYPE(SGUIDE) :: OBJ
      integer,intent(in) :: NS
      real(kind(1.D0)) :: L,W,H,M(4)
      integer :: i
      if ((NS>0).or.(abs(NS)==OBJ%NSEG)) return
      if (OBJ%NSEG>0) then
        W = OBJ%W(1)
        H = OBJ%H(1)
        M(1:4) = OBJ%MC(1:4,1)
      else
        W = OBJ%FRAME%SIZE(1)
        H = OBJ%FRAME%SIZE(2)
        M = (/2.,2.,2.,2./)
      endif
      if (NS==0) then
        OBJ%NSEG=-MIN(30,MAX(3,NINT(OBJ%FRAME%SIZE(3)/500)))
      else
        OBJ%NSEG=MIN(ABS(NS),MAX_SEG)
      endif
      !Ltot=OBJ%FRAME%SIZE(3)+(OBJ%NSEG-1)*OBJ%GAP
      !write(*,*) 'SGUIDE_ADJUST_SEGMENTS for SIZE(3) ',OBJ%FRAME%SIZE(3)
      L=(OBJ%FRAME%SIZE(3)-(OBJ%NSEG-1)*OBJ%GAP)/OBJ%NSEG
      do i=1,OBJ%NSEG
          OBJ%L(i)=L
          OBJ%W(i)=W
          OBJ%H(i)=H
          OBJ%MC(1:4,i)=M(1:4)
      enddo

      end SUBROUTINE SGUIDE_ADJUST_SEGMENTS

!---------------------------------------------------------
      SUBROUTINE SGUIDE_SET_PROFILE(INST,iside)
! set guide profile according to the user's input
! ISIDE=1 (horizontal) or 2 (vertical)
!---------------------------------------------------------
      integer,intent(in) :: INST,iside
      TYPE(SGUIDE),POINTER :: OBJ
      real(kind(1.D0)) :: xmin,dx,Ltot,f,foc,z,wmin,wmax,y2,y20,a2,eps
      integer :: i,k
1     format(a,': ',6(1x,G12.5))
      if (.not.SGUIDE_isValid(INST)) return
      OBJ => ASGUIDES(INST)%X

      ! guide length is the sum of segment lengths + gaps
      Ltot = 0.D0
      OBJ%WIN(3,0)=0.D0
      do i=1,OBJ%NSEG
        Ltot = Ltot + OBJ%L(i)
        if (i<OBJ%NSEG) Ltot = Ltot + 0.5*OBJ%GAP
        ! position in the middle of the gap, or guide exit
        OBJ%WIN(3,i)=Ltot
        if (i<OBJ%NSEG) Ltot = Ltot + 0.5*OBJ%GAP
      enddo

      select case(OBJ%IPROF(iside))
! manual profile setting
      case(sguide_prof_man)
        if (iside==1) then
          OBJ%EXI(iside)=OBJ%W(OBJ%NSEG)
        else
          OBJ%EXI(iside)=OBJ%H(OBJ%NSEG)
        endif
! flat walls
      case(sguide_prof_flat)
        xmin=OBJ%FRAME%SIZE(iside)
        dx=(OBJ%EXI(iside)-OBJ%FRAME%SIZE(iside))
        do i=1,OBJ%NSEG
          ! W, H are final dimensions of a segment in the middle of the gap (or guide exit)
          if (iside==1) then
            OBJ%W(i)=xmin+dx*OBJ%WIN(3,i)/Ltot
          else
            OBJ%H(i)=xmin+dx*OBJ%WIN(3,i)/Ltot
          endif
          !write(*,1) 'SGUIDE_SET_PROFILE', i, OBJ%WIN(3,i), OBJ%W(i), OBJ%H(i)
        enddo
! curved walls
      case(sguide_prof_para,sguide_prof_ell)
        wmin=OBJ%FRAME%SIZE(iside)
        wmax=OBJ%EXI(iside)
        if (wmin>wmax) then
        ! converging
          z=wmin
          wmin=wmax
          wmax=z
          k=0
        else
        ! expanding
          k=1
        endif
        eps=1.D0-(wmin/wmax)**2
      ! parallel walls
        if (eps<1.D-10) then
          do i=1,OBJ%NSEG
            if (iside==1) then
              OBJ%W(i)=wmin
            else
              OBJ%H(i)=wmin
            endif
          enddo
        else if (OBJ%IPROF(iside)==sguide_prof_para) then
  ! parabolic
          ! focal length
          f=wmax**2*eps/16.D0/Ltot
          ! distance to focal point
          foc=f*((wmin/4/f)**2-1.D0)
          y20=(1.D0-k)*Ltot+foc+f
          do i=1,OBJ%NSEG
            y2=y20+(2.D0*k-1.D0)*OBJ%WIN(3,i)
            if (iside==1) then
              OBJ%W(i)=4.D0*SQRT(abs(f*y2))
            else
              OBJ%H(i)=4.D0*SQRT(abs(f*y2))
            endif
          enddo
        else if (OBJ%IPROF(iside)==sguide_prof_ell) then
  ! semi ellipse
          a2=Ltot**2/EPS
          do i=1,OBJ%NSEG
            y2=wmax*SQRT(1.D0-((1-2.D0*k)*OBJ%WIN(3,i)+k*Ltot)**2/a2)
            if (iside==1) then
              OBJ%W(i)=y2
            else
              OBJ%H(i)=y2
            endif
          enddo
        endif
      end select
      end SUBROUTINE SGUIDE_SET_PROFILE


!---------------------------------------------------------
      real(kind(1.D0)) function SGUIDE_POS(ider,d,z)
! calculate wall position or derivative at given z coordinate (z=0 at the segment entry)
! d(4) ... wall parameters (= OBJ%D)
! ider=0..3 for 0-3nd derivative
!---------------------------------------------------------
      integer,intent(in) :: ider
      real(kind(1.D0)),intent(in) :: d(4),z
      real(kind(1.D0)) :: Y
      Y=0.D0
      if (ider==0) then
        Y=d(1)*z**3+d(2)*z**2+d(3)*z+d(4)
      else if (ider==1) then
        Y=3*d(1)*z**2+2*d(2)*z+d(3)
      else if (ider==2) then
        Y=6*d(1)*z+2*d(2)
      else if (ider==3) then
        Y=6*d(1)
      endif
      SGUIDE_POS=Y
      end function SGUIDE_POS


!---------------------------------------------------------
      SUBROUTINE SGUIDE_CALC_PARAM(INST)
! calculate parameters of all segments walls cubic spline,
! pos = D(1)*[x-x(i)]**3+D(2)*[x-x(i)]**2+D(3)*[x-x(i)]+D(4)
!---------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(SGUIDE),POINTER :: OBJ
      integer,parameter :: SIDE(4)=(/1,1,2,2/)
      integer,parameter :: SG(4)=(/1,-1,1,-1/)
      real(kind(1.D0)) :: x(MAX_SEG+2),y(MAX_SEG+2)
      real(kind(1.D0)) :: LN,h,d(3),space,delta
      integer :: i,j,k,n
1     format(a,6(1x,G12.5))
      if (.not.SGUIDE_isValid(INST)) return
      OBJ => ASGUIDES(INST)%X

      OBJ%D=0.D0
      do j=1,4
    !  if (OBJ%ACTIVE(j)) then
        k=SIDE(j)
! calculate edge positions
        x(1)=0.D0
        y(1)=SG(j)*0.5D0*OBJ%FRAME%SIZE(k)
        LN=0.D0
        OBJ%WIN(k,0)=0.D0
        OBJ%WIN(3,0)=0.D0
        do i=1,OBJ%NSEG
          space=0.5D0*OBJ%GAP
          if (i==OBJ%NSEG) space=0.D0
          LN=LN+OBJ%L(i)+space
          !write(*,*) 'SGUIDE_CALC_PARAM ',i,LN
          OBJ%WIN(k,i)=0.5D0*OBJ%RHO(k)*LN**2
          OBJ%WIN(3,i)=LN
          x(i+1)=LN
          if (k==1) then
            y(i+1)=OBJ%WIN(k,i)+SG(j)*0.5D0*OBJ%W(i)
          else
            y(i+1)=OBJ%WIN(k,i)+SG(j)*0.5D0*OBJ%H(i)
          endif
          OBJ%NR(j,i)=MIRROR_GETTABLE(OBJ%MC(j,i))
          LN=LN+space
        ! DBG
        !  if ((j==1).and.(i==1)) then
        !    write(*,1) 'SGUIDE_CALC_PARAM ',OBJ%NR(j,i),OBJ%MC(j,i)
        !    if (OBJ%NR(j,i)>0) then
        !      write(*,1) '  REF=',MIRROR_REF(OBJ%NR(j,i),TWOPI*1.5*GammaNi,1.D0)
        !    endif
        !  endif

        enddo

        n=OBJ%NSEG+1
        !if (j==1) then
        !  do i=1,n
        !    write(*,1) 'i,x,y= ',i,x(i),y(i)
        !  enddo
        !endif
        if (OBJ%SMOOTH(j).and.(n>2)) then
! get boundary conditions = 2nd derivatives at the ends
          !call quadinterp3(x(1:3),y(1:3),d)
          !alpha=d(1)
! get coefficients for interpolation - quadratic approximation with continuous 1st derivative
          do i=1,OBJ%NSEG
            if (i<OBJ%NSEG) then
              call quadinterp3(x(i:i+2),y(i:i+2),d)
              delta=0.D0
            else
              call quadinterp3(x(i-1:i+1),y(i-1:i+1),d)
              delta=x(i)-x(i-1)
            endif
            OBJ%D(1,j,i)=0.D0
            OBJ%D(2,j,i)=d(1)
            OBJ%D(3,j,i)=d(2)+2*d(1)*delta
            OBJ%D(4,j,i)=d(3)+d(2)*delta+d(1)*delta**2

          !  h=x(i+1)-x(i)
          !  OBJ%D(1,j,i)=0.D0
          !  OBJ%D(2,j,i)=(y(i+1)-y(i)-alpha*h)/h**2
          !  OBJ%D(3,j,i)=alpha
          !  OBJ%D(4,j,i)=y(i)
          !  alpha=alpha+2*OBJ%D(2,j,i)*h
          enddo
          !do i=1,4
          !  OBJ%D(i,j,OBJ%NSEG)=OBJ%D(i,j,OBJ%NSEG-1)
          !enddo
          ! if (j==1) then
              !write(*,*)
              !do i=1,n
              !  write(*,1) ' j,i,d(1:3)= ',i,OBJ%D(4,j,i)/1.D3,OBJ%D(3,j,i),OBJ%D(2,j,i)*1.D3
              !enddo
           ! endif

! linear interpolation
        else
          do i=1,OBJ%NSEG
            h=x(i+1)-x(i)
            OBJ%D(1,j,i)=0.D0
            OBJ%D(2,j,i)=0.D0
            OBJ%D(3,j,i)=(y(i+1)-y(i))/h
            OBJ%D(4,j,i)=y(i)
          enddo
        endif
        OBJ%FRAME%SIZE(3)=LN
        !write(*,*) 'SGUIDE_CALC_PARAM len=',LN
    !  endif
      enddo
! update total guide length
      end SUBROUTINE SGUIDE_CALC_PARAM


!-------------------------------------------------------------------------
! Creator for SGUIDE, return instance number
! Memory is allocated on the first free element of ASGUIDES array
!-------------------------------------------------------------------------
      integer function SGUIDE_CREATE(ID,NAMESTR)
      character(*) :: ID,NAMESTR
      integer :: ierr,i
      character(32) :: S
        ierr=1
        i=1
        do while ((i.lt.SGUIDES_DIM).and.(ASGUIDES(i)%IDX.gt.0))
          i=i+1
        enddo
        if (ASGUIDES(i)%IDX.le.0) then
          allocate(ASGUIDES(i)%X,STAT=ierr)
          if(ierr.eq.0) then
            SGUIDES_NC=SGUIDES_NC+1
            call SGUIDE_DEFAULT(i)
            ASGUIDES(i)%IDX=i
            ASGUIDES(i)%X%FRAME%ID=trim(ID)
            ASGUIDES(i)%X%FRAME%NAME=trim(NAMESTR)
          endif
        endif
        if (ierr.eq.0) then
          SGUIDE_CREATE=i
        else
          call INT2STR(SGUIDES_DIM,S)
          call MSG_ERROR('SGUIDE_CREATE','Reached maximum number or instances ('//trim(S)//')',1,1)
          SGUIDE_CREATE=0
        endif
      end function SGUIDE_CREATE

!---------------------------------------------------------
      subroutine SGUIDE_DISPOSE(INST)
!---------------------------------------------------------
        integer,intent(in) :: INST
        if ((inst.gt.0).and.(inst.le.SGUIDES_DIM)) then
          if (associated(ASGUIDES(inst)%X)) then
            DEallocate(ASGUIDES(inst)%X)
            nullify(ASGUIDES(inst)%X)
            SGUIDES_NC=SGUIDES_NC-1
          endif
          ASGUIDES(inst)%IDX=0
        endif
      end subroutine SGUIDE_DISPOSE

!---------------------------------------------------------
      subroutine SGUIDE_DISPOSE_ALL
!---------------------------------------------------------
        integer :: i
        do i=1,SGUIDES_DIM
          call SGUIDE_DISPOSE(i)
        enddo
      end subroutine SGUIDE_DISPOSE_ALL

!-------------------------------------------------------------------------
! Add a component
!-------------------------------------------------------------------------
      INTEGER function AddSGUIDE(INST)
      integer,intent(in) :: INST
      integer :: i
      i=0
      if (SGUIDE_isValid(INST)) then
        i=SGUIDE_CREATE(ASGUIDES(INST)%X%FRAME%ID,ASGUIDES(INST)%X%FRAME%NAME)
        if (i.gt.0) then
          ASGUIDES(i)%X=ASGUIDES(INST)%X
          ASGUIDES(i)%IDX=i
        !  write(*,*) 'AddSGUIDE ',trim(ASGUIDES(i)%X%FRAME%ID),i
        endif
      endif
      AddSGUIDE=i
      end function AddSGUIDE

!------------------------------------
      SUBROUTINE SGUIDE_PREPARE(INST,IERR)
! set default parameters
!------------------------------------
      integer,intent(in) :: INST
      integer,intent(out) :: IERR
      REAL(KIND(1.D0)) :: V(3), DV(3), R(3,3), ALPHA, BETA, DX, DY
      TYPE(SGUIDE),POINTER :: OBJ
        IERR=0
        if (.not.SGUIDE_isValid(INST)) return
        OBJ => ASGUIDES(INST)%X
      ! calculate wall profiles
        call SGUIDE_SET_PROFILE(INST,1)
        call SGUIDE_SET_PROFILE(INST,2)
      ! calculate curved wall parameters
        call SGUIDE_CALC_PARAM(INST)
      ! initialize rotation matrices
        call FRAME_INIT_MAT(OBJ%FRAME)
      ! for curved guides, add deflection
        OBJ%FRAME%ISDEFL=(ABS(OBJ%RHO(1))+ABS(OBJ%RHO(2))>1.D-9)
        if (OBJ%FRAME%ISDEFL) call FRAME_DEFLECTION(OBJ%FRAME, OBJ%RHO)
      END SUBROUTINE SGUIDE_PREPARE

!-------------------------------------------------------------
      logical function SGUIDE_isValid(INST)
! check index and association
!-------------------------------------------------------------
      integer,intent(in) :: INST
      SGUIDE_isValid= ((INST.gt.0).and.(INST.le.SGUIDES_DIM).and.associated(ASGUIDES(INST)%X))
      end function SGUIDE_isValid


!-------------------------------------------------------------
      SUBROUTINE SGUIDE_GET(INST,OBJ)
! return reference to given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(PSGUIDE),intent(out) :: OBJ
        if (SGUIDE_isValid(INST)) then
          OBJ%IDX=ASGUIDES(INST)%IDX
          OBJ%X => ASGUIDES(INST)%X
        else
          NULLIFY(OBJ%X)
          OBJ%IDX=0
        endif
      end SUBROUTINE SGUIDE_GET

!-------------------------------------------------------------
      SUBROUTINE SGUIDE_GET_FRAME(INST,OBJ)
! return reference to given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(PFRAME),intent(out) :: OBJ
        if (SGUIDE_isValid(INST)) then
          OBJ%IDX=ASGUIDES(INST)%IDX
          OBJ%X => ASGUIDES(INST)%X%FRAME
        else
          NULLIFY(OBJ%X)
          OBJ%IDX=0
        endif
      end SUBROUTINE SGUIDE_GET_FRAME

!------------------------------------
      SUBROUTINE SGUIDE_DEFAULT(INST)
! set default parameters
!------------------------------------
      integer,intent(in) :: INST
      TYPE(SGUIDE),POINTER :: OBJ
        if (.not.SGUIDE_isValid(INST)) return
        OBJ => ASGUIDES(INST)%X
        call FRAME_CLEAR(OBJ%FRAME)
        OBJ%FRAME%SIZE(3)=1.D3
        OBJ%FRAME%CLASS=CCLS_SGUIDE
        OBJ%FRAME%SHAPE=FRAME_SHAPE_BOX
        OBJ%SMOOTH=.false.
        OBJ%ACTIVE=.true.
        OBJ%MONITOR=.false.
        OBJ%NSEG=1
        OBJ%W=40.D0
        OBJ%H=100.D0
        OBJ%L=1.D3
        OBJ%RHO=0.D0
        OBJ%MC=1.D0
        OBJ%NR=-1
        OBJ%GAP=0.D0
        OBJ%WAV=2.D-4
        OBJ%FRONTMASK=.true.
        OBJ%CLOSED=.false.
        OBJ%GAPABS=.true.
        OBJ%EXI=(/40.D0,80.D0/)
        OBJ%IPROF=(/0,0/)
        OBJ%LOGBNC=0
        OBJ%MLOG=0
        OBJ%MISALIGN=(/0.D0,0.D0/)
      END SUBROUTINE SGUIDE_DEFAULT


!---------------------------------------------------------
      SUBROUTINE SGUIDE_INP(INST,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(SGUIDE),POINTER :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      if (.not.SGUIDE_isValid(INST)) return
      OBJ => ASGUIDES(INST)%X
      LR=0
      NARG=0
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call FIELD2ARRAY(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
          call SGUIDE_INP_R(OBJ,trim(ARG%ID),NUM,LR)
        case(ID_FIELD_TABLE)
          call SGUIDE_INP_TAB(OBJ,ARG,LR)
        case default
          write(*,*) 'SGUIDE_INP: undefined field type '//trim(ARG%DEF%TID)
      end select
      NARG=LR
      END SUBROUTINE SGUIDE_INP


!---------------------------------------------------------
      SUBROUTINE SGUIDE_OUT(INST,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(SGUIDE),POINTER :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      if (.not.SGUIDE_isValid(INST)) return
      OBJ => ASGUIDES(INST)%X
      LR=0
      NARG=0
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call SGUIDE_OUT_R(OBJ,trim(ARG%ID),NUM,LR)
          call ARRAY2FIELD(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
        case (ID_FIELD_TABLE)
          call SGUIDE_OUT_TAB(OBJ,ARG,LR)
        case default
          write(*,*) 'SGUIDE_OUT: undefined field type '//trim(ARG%DEF%TID)
      end select
      NARG=LR
      END SUBROUTINE SGUIDE_OUT

!---------------------------------------------------------
      SUBROUTINE SGUIDE_INP_R(OBJ,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... REAL*8 parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      TYPE(SGUIDE),intent(out) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(in) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR,i,N
      LR=1
      SELECT CASE (trim(PNAME))
! profile horizontal
          CASE('PROFH')
            OBJ%IPROF(1)=MAX(MIN(NINT(ARG(1)),3),0)
! profile vertical
          CASE('PROFV')
            OBJ%IPROF(2)=MAX(MIN(NINT(ARG(1)),3),0)
! exit size (H,V) [mm]
          CASE('EXIT');    OBJ%EXI(1)=ARG(1);OBJ%EXI(2)=ARG(2);LR=2
          CASE('EXITH');   OBJ%EXI(1)=ARG(1)
          CASE('EXITV');   OBJ%EXI(2)=ARG(1)
! no. of segments
          !CASE('NSEG'); OBJ%NSEG=MAX(1,MIN(NINT(ARG(1)),MAX_SEG))
          CASE('NSEG')
            N=MIN(NINT(ARG(1)),MAX_SEG)
            call SGUIDE_ADJUST_SEGMENTS(OBJ, N)
! curvatures [1/m] - convert to [1/mm]
          CASE('RHO'); OBJ%RHO(1)=ARG(1)/1000.; OBJ%RHO(2)=ARG(2)/1000.;LR=2
! reflectivity
        !  CASE('REF0'); OBJ%REF0=ARG(1)
        !  CASE('REFM'); OBJ%REF1=ARG(1)
        !  CASE('REFW'); OBJ%REFW=ARG(1)
! gap [mm]
          CASE('GAP'); OBJ%GAP=ARG(1)
! misalignment [mm]
          CASE('MISALIGN');    OBJ%MISALIGN(1)=ARG(1);OBJ%MISALIGN(2)=ARG(2);LR=2
! waviness [mrad]
          CASE('WAV'); OBJ%WAV=ARG(1)/1.D3
! smoothness (bit mask)
          case('SMOOTH')
            do i=1,4
              OBJ%SMOOTH(i)=(ARG(i).ne.0.D0)
            enddo
            LR=4
! active (bit mask)
          case('ACTIVE')
            do i=1,4
             OBJ%ACTIVE(i)=(ARG(i).ne.0.D0)
            enddo
            LR=4
! monitor (bit mask)
          case('MONITOR')
            do i=1,4
             OBJ%MONITOR(i)=(NINT(ARG(i)).eq.1)
            enddo
            LR=4
! closed guide
          case('CLOSED'); OBJ%CLOSED=(NINT(ARG(1)).eq.1)
! front mask
          case('FRONTMASK'); OBJ%FRONTMASK=(NINT(ARG(1)).eq.1)
! absorbing gap
          case('GAPABS'); OBJ%GAPABS=(NINT(ARG(1)).eq.1)
! log bounces
          case('LOGBNC'); OBJ%LOGBNC=MAX(MIN(NINT(ARG(1)),7),0)
! try FRAME parameters
          CASE DEFAULT; CALL FRAME_INP_R(OBJ%FRAME,trim(PNAME),ARG(1),LR)
      END SELECT
      NARG=LR
      END SUBROUTINE SGUIDE_INP_R

!---------------------------------------------------------
      SUBROUTINE SGUIDE_OUT_R(OBJ,PNAME,ARG,NARG)
! output data from OBJ to ARG for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... REAL*8 parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      TYPE(SGUIDE),intent(in) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(out) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR,i
      LR=1
      SELECT CASE (trim(PNAME))
! profile horizontal
          CASE('PROFH'); ARG(1)=OBJ%IPROF(1)
! profile vertical
          CASE('PROFV'); ARG(1)=OBJ%IPROF(2)
! exit size (H,V) [mm]
          CASE('EXIT');    ARG(1)=OBJ%EXI(1);ARG(2)=OBJ%EXI(2);LR=2
          CASE('EXITH');   ARG(1)=OBJ%EXI(1)
          CASE('EXITV');   ARG(1)=OBJ%EXI(2)
! no. of segments
          CASE('NSEG'); ARG(1)=OBJ%NSEG
! curvatures [1/m] - convert to [1/mm]
          CASE('RHO');  ARG(1)=OBJ%RHO(1)*1000.;ARG(1+1)=OBJ%RHO(2)*1000.;LR=2
! reflectivity
          !CASE('REF0'); ARG(1)=OBJ%REF0
          !CASE('REFM'); ARG(1)=OBJ%REF1
          !CASE('REFW'); ARG(1)=OBJ%REFW
! gap [mm]
          CASE('GAP'); ARG(1)=OBJ%GAP
! misalignment [mm]
          CASE('MISALIGN');    ARG(1)=OBJ%MISALIGN(1);ARG(2)=OBJ%MISALIGN(2);LR=2
! waviness [mrad]
          CASE('WAV'); ARG(1)=OBJ%WAV*1.D3
! smoothness (bit mask)
          case('SMOOTH')
            do i=1,4
              ARG(i)=LOG2INT(OBJ%SMOOTH(i))
            enddo
            LR=4
! active (bit mask)
          case('ACTIVE')
            do i=1,4
              ARG(i)=LOG2INT(OBJ%ACTIVE(i))
            enddo
            LR=4
! monitor (bit mask)
          case('MONITOR')
            do i=1,4
              ARG(i)=LOG2INT(OBJ%MONITOR(i))
            enddo
            LR=4
! closed guide
          case('CLOSED'); ARG(1)=LOG2INT(OBJ%CLOSED)
! front mask
          case('FRONTMASK'); ARG(1)=LOG2INT(OBJ%FRONTMASK)
! absorbing gap
          case('GAPABS'); ARG(1)=LOG2INT(OBJ%GAPABS)
! log bounces
          case('LOGBNC'); ARG(1)=OBJ%LOGBNC
! try FRAME parameters
          CASE DEFAULT; CALL FRAME_OUT_R(OBJ%FRAME,trim(PNAME),ARG(1),LR)
      END SELECT
      NARG=LR
      END SUBROUTINE SGUIDE_OUT_R


!---------------------------------------------------------
      SUBROUTINE SGUIDE_INP_TAB(OBJ,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD structure of SELECT type
! OUTPUT
!    OBJ     .... XTAL structure
!    NARG    .... number of items read from ARG
C---------------------------------------------------------
      TYPE(SGUIDE),intent(out) :: OBJ
      TYPE(TFIELD),intent(in) :: ARG
      integer,intent(out) :: NARG
      integer :: LR,ierr,I,IL
1     format(a,6(1x,G12.5))
      real(kind(1.D0)) :: NUM(7)
      LR=0
      ierr=0
      SELECT CASE (trim(ARG%ID))
        CASE('SEG')
          OBJ%NSEG=ARG%NP
          do i=1,ARG%NP
            call ROW2ARRAY(ARG,i,0,NUM,7,IL)
            ! write(*,1) 'SGUIDE_INP_TAB ',i,NUM(1:3)
            OBJ%W(i)=NUM(1)
            OBJ%H(i)=NUM(2)
            OBJ%L(i)=NUM(3)
            OBJ%MC(1:4,i)=NUM(4:7)
            LR=LR+IL
          enddo
          LR=MAX_SEG*7 ! return max size
      END SELECT
      NARG=LR
      END SUBROUTINE SGUIDE_INP_TAB


!---------------------------------------------------------
      SUBROUTINE SGUIDE_OUT_TAB(OBJ,ARG,NARG)
! input data from OBJ to SARG for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    OBJ     .... XTAL structure
! OUTPUT
!    ARG     .... TFIELD structure of SELECT type
!    NARG    .... number of items read from ARG
C---------------------------------------------------------
      TYPE(SGUIDE),intent(in) :: OBJ
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR,IL,i
1     format(a,6(1x,G12.5))
      real(kind(1.D0)) :: NUM(7)
      LR=0
      SELECT CASE (trim(ARG%ID))
        CASE('SEG')
          ARG%NP=OBJ%NSEG
          do i=1,OBJ%NSEG
            NUM=(/OBJ%W(i),OBJ%H(i),OBJ%L(i),OBJ%MC(1,i),OBJ%MC(2,i),OBJ%MC(3,i),OBJ%MC(4,i)/)
            !write(*,1) 'SGUIDE_OUT_TAB ',i,OBJ%W(i),OBJ%H(i),OBJ%L(i)
            call ARRAY2ROW(ARG,i,0,NUM,7,IL)
            LR=LR+IL
          enddo
          LR=MAX_SEG*7 ! return max size
      END SELECT
      NARG=LR
      END SUBROUTINE SGUIDE_OUT_TAB


      end MODULE SGUIDES