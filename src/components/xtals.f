!//////////////////////////////////////////////////////////////////////
!////  $Id: xtals.f,v 1.50 2018/10/29 17:50:48 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2009, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.50 $
!////     $Date: 2018/10/29 17:50:48 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Describes component class: XTAL
!////  bent single-crystal
!////////////////////////////////////////////////////////////////////////
      MODULE XTALS
      use CONSTANTS
      use ENUMDEF
      use CLASSES
      use FRAMES
      use LATTICE
      use TABLE_CRYSTALS
      use XTALS_REF
      ! use TRACINGDATA
      use XMLUTILS
      use FOCARRAYS
      use ARRAY3D
      implicit none

! static storage of XTAL instances
      integer, parameter :: XTALS_DIM=MAX_CR

! bent crystal
      TYPE XTAL
        TYPE(TFRAME) :: FRAME
        character*8 :: CRNAME ! crystal name = ID in TABLE_CRYSTALS
        integer :: SGN   ! sign of the take-off angle
        integer :: HKL(3)  ! selected node [hkl] -> direction of tracing
        real(kind(1.D0)) :: ORI(3)  ! [hkl] nominal reflection (used for positioning)
        real(kind(1.D0)) :: A(3),B(3)  ! [hkl] lattice vectors which define reference plane (together with CHI,PSI define the B-matrix)
        real(kind(1.D0)) :: CHI,PSI ! cutting angles (round y,x) [deg]
        real(kind(1.D0)) :: POI ! Poisson const.
        real(kind(1.D0)) :: T   ! temperature [K]
        real(kind(1.D0)) :: RH ! horizontal bending curvature [m^-1]
        real(kind(1.D0)) :: RV ! vertical bending curvature [m^-1]
        real(kind(1.D0)) :: SWGAP ! gap between sandwich slabs [mm]
        real(kind(1.D0)) :: DT(3) ! linear temperature gradient [K/mm]
        real(kind(1.D0)) :: DKK ! delta_k/k limit used in searching for reflections
        logical :: AUTOADJ   ! automatically adjust in reflection position
        integer :: UMWEG   ! trace Renninger reflections
        integer :: MAXHKL   ! max. reflection index of searched reflection
        integer :: NSW      ! number of sandwich slabs
      ! subclasses
        TYPE(TFOCARRAY) :: VFOC  ! vertically focusing array
        TYPE(TFOCARRAY) :: HFOC  ! horizontally focusing array

      ! calculated fields
        INTEGER :: ITAB        ! index to TABLE_CRYSTALS
        INTEGER :: ICR         ! index to XREF_TAB and NODES_TAB arrays (from XTALS_REF)
        INTEGER :: MAX_ITER_COUNT         ! counter of events exceeding max. iteration number
        logical :: FCONE       ! flat-cone adjustment
      !  integer :: NNODES      ! number of reflecting rec. lattice nodes
      !  TYPE(RLNODE) :: NODES(0:MAX_NODES) ! reflecting rec. lattice nodes
        real(kind(1.D0)) :: LAMBDA  ! nominal wavelength (see XTAL_SET_LAMBDA)
        INTEGER :: REFLIST(10,MAX_NODES)   ! list of branches leading to selected node
        ! <= 10 reflections per branch, <= MAX_NODES reflections
        INTEGER :: REFLIST_N(MAX_NODES) ! number of reflections per branch
        INTEGER :: REFLIST_NB  ! number branches

      !  real(kind(1.D0)) :: DHKL ! dhkl for the nominal reflection (ORI)
      !  real(kind(1.D0)) :: G(3) ! G-vector for the nominal reflection (ORI) in local coordinates [1/A]
        real(kind(1.D0)) :: BMAT(3,3)   ! B-matrix (conversion from hkl to local coordinates)
        real(kind(1.D0)) :: SIGABS      ! absorption cross-section [mm^-1.A^-1]
        real(kind(1.D0)) :: SIGINC      ! incoherent scattering cross-section [mm^-1]
        real(kind(1.D0)) :: SIGTDS(4)   ! coefficients for thermal scattering (see comments at GET_CR_SIGMA)
        TYPE(TARRAY3D) :: A3D           ! used for tracing in 3D segments array


!        REAL*8 HMOS,VMOS,VOL ! mosaicity (H,V), cell volume,
!        REAL*8 DLAM,DEXT,EXT1,DELTA         ! extinction parameters: lamellae thickness, ext. length, primary ext., Darwin box
!        REAL*8 THB,LAMBDA,STMCH,CTMCH,QHKL,REF,DETA,TANVOL ! Bragg angle, wavelength, sin(t), cos(t), kin. refl., ...
!        REAL*8 G(3),GAMA(3),GTOT,DGR,DGA ! deformation gradient, diff. vector, d-gradient, gradient orientation
!        LOGICAL*4 MAPG(3)
      END TYPE XTAL
! pointer type for XTAL
      type PXTAL
        integer :: IDX  ! instance index
        TYPE(XTAL),pointer :: X
      end type PXTAL
! instances of XTAL. AXTALS(0) is always unallocated
      integer :: XTALS_NC
      type(PXTAL) :: AXTALS(1:XTALS_DIM)
! private flag for errors in the recurrent calls to CAL_SREF
      integer :: NODE_ERR=0
      private NODE_ERR

      contains

!-------------------------------------------------------------------------
! Creator for XTAL, return instance number
! Memory is allocated on the first free element of AXTALS array
!-------------------------------------------------------------------------
      integer function XTAL_CREATE(ID,NAMESTR)
      character(*) :: ID,NAMESTR
      integer :: ierr,i
      character(32) :: S
        ierr=1
        i=1
        do while ((i.lt.XTALS_DIM).and.(AXTALS(i)%IDX.gt.0))
          i=i+1
        enddo
        if (AXTALS(i)%IDX.le.0) then
          allocate(AXTALS(i)%X,STAT=ierr)
          if(ierr.eq.0) then
            XTALS_NC=XTALS_NC+1
            call XTAL_DEFAULT(i)
            AXTALS(i)%IDX=i
            AXTALS(i)%X%FRAME%ID=trim(ID)
            AXTALS(i)%X%FRAME%NAME=trim(NAMESTR)
            AXTALS(i)%X%ICR=XREF_ADD_CR()
          endif
        endif
        if (ierr.eq.0) then
          XTAL_CREATE=i
        else
          call INT2STR(XTALS_DIM,S)
          call MSG_ERROR('XTAL_CREATE','Reached maximum number or instances ('//trim(S)//')',1,1)
          XTAL_CREATE=0
        endif
      end function XTAL_CREATE

!---------------------------------------------------------
      subroutine XTAL_DISPOSE(INST)
!---------------------------------------------------------
        integer,intent(in) :: INST
        if ((inst.gt.0).and.(inst.le.XTALS_DIM)) then
          if (associated(AXTALS(inst)%X)) then
            DEallocate(AXTALS(inst)%X)
            nullify(AXTALS(inst)%X)
            XTALS_NC=XTALS_NC-1
          endif
          AXTALS(inst)%IDX=0
        endif
      end subroutine XTAL_DISPOSE

!---------------------------------------------------------
      subroutine XTAL_DISPOSE_ALL
!---------------------------------------------------------
        integer :: i
        do i=1,XTALS_DIM
          call XTAL_DISPOSE(i)
        enddo
      end subroutine XTAL_DISPOSE_ALL

!-------------------------------------------------------------------------
! Add a component
!-------------------------------------------------------------------------
      INTEGER function AddXTAL(INST)
      integer,intent(in) :: INST
      integer :: i
      i=0
      if (XTAL_isValid(INST)) then
        i=XTAL_CREATE(AXTALS(INST)%X%FRAME%ID,AXTALS(INST)%X%FRAME%NAME)
        if (i.gt.0) then
          AXTALS(i)%X=AXTALS(INST)%X
          AXTALS(i)%IDX=i
        !  write(*,*) 'AddXTAL ',trim(AXTALS(i)%X%FRAME%ID),i
        endif
      endif
      AddXTAL=i
      end function AddXTAL

!------------------------------------
      SUBROUTINE XTAL_PREPARE(INST,IERR)
! prepare dependent fields after input
!------------------------------------
      integer,intent(in) :: INST
      integer,intent(out) :: IERR
      TYPE(XTAL),POINTER :: OBJ
        IERR=0
        if (.not.XTAL_isValid(INST)) return
        OBJ => AXTALS(INST)%X
        call FRAME_INIT_MAT(OBJ%FRAME)
      END SUBROUTINE XTAL_PREPARE

!-------------------------------------------------------------
      logical function XTAL_isValid(INST)
! check index and association
!-------------------------------------------------------------
      integer,intent(in) :: INST
      XTAL_isValid= ((INST.gt.0).and.(INST.le.XTALS_DIM).and.associated(AXTALS(INST)%X))
      end function XTAL_isValid

!------------------------------------
      SUBROUTINE XTAL_DEFAULT(INST)
! set default parameters
!------------------------------------
      integer,intent(in) :: INST
      TYPE(XTAL),POINTER :: OBJ
        if (.not.XTAL_isValid(INST)) return
        OBJ => AXTALS(INST)%X
        call FRAME_CLEAR(OBJ%FRAME)
        call FOCARRAY_DEFAULT(OBJ%VFOC)
        call FOCARRAY_DEFAULT(OBJ%HFOC)
        OBJ%FRAME%CLASS=CCLS_XTAL
        OBJ%FRAME%SHAPE=FRAME_SHAPE_BOX
        OBJ%FRAME%SIZE=(/120.D0,40.D0,5.D0/)
        OBJ%CRNAME="Si"
        OBJ%SGN=1
        OBJ%AUTOADJ=.true.
        OBJ%UMWEG=0
        OBJ%MAXHKL=10
        OBJ%HKL=(/0,0,0/)
        OBJ%ORI=(/2.D0,2.D0,0.D0/)
        OBJ%A=(/2.D0,2.D0,0.D0/)
        OBJ%B=(/0.D0,0.D0,1.D0/)
        OBJ%CHI=0.D0*deg
        OBJ%PSI=0.D0*deg
        OBJ%POI=0.3D0
        OBJ%T=298.D0
        OBJ%RH=1.D-1
        OBJ%RV=0.D0
        OBJ%DT=0.D0
        OBJ%DKK=1.D-2
        OBJ%LAMBDA=0.D0
        OBJ%FCONE=.false.
        OBJ%NSW=1
        OBJ%SWGAP=0.05
        OBJ%ITAB=GET_CR_INDEX(trim(OBJ%CRNAME))
      END SUBROUTINE XTAL_DEFAULT

!-------------------------------------------------------------
      SUBROUTINE XTAL_GET(INST,OBJ)
! return reference to given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(PXTAL),intent(out) :: OBJ
        if (XTAL_isValid(INST)) then
          OBJ%IDX=AXTALS(INST)%IDX
          OBJ%X => AXTALS(INST)%X
        else
          NULLIFY(OBJ%X)
          OBJ%IDX=0
        endif
      end SUBROUTINE XTAL_GET

!-------------------------------------------------------------
      SUBROUTINE XTAL_GET_FRAME(INST,OBJ)
! return reference to given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(PFRAME),intent(out) :: OBJ
        if (XTAL_isValid(INST)) then
          OBJ%IDX=AXTALS(INST)%IDX
          OBJ%X => AXTALS(INST)%X%FRAME
        else
          NULLIFY(OBJ%X)
          OBJ%IDX=0
        endif
      end SUBROUTINE XTAL_GET_FRAME


!------------------------------------------------------------
      SUBROUTINE XTAL_AUTOFOCUS(CR)
! Adjust crystal focal distances
!------------------------------------------------------------
      TYPE (XTAL) :: CR
      real(kind(1.D0)) :: K0,STH,LL,TH,TH1,TH2
      TYPE(TCR_ITEM) :: CDATA
      real(kind(1.D0)) :: DHKL
      REAL(KIND(1.D0)) :: DS(3),R(3)
      integer :: NS(3)
      call GET_CR_ITEM(CR%ITAB,CDATA)
      DHKL=GET_DHKL(CDATA%LAT,CR%ORI,CR%T)
      if ((CR%LAMBDA.gt.0.D0).and.(CR%LAMBDA.lt.2.D0*DHKL)) then
        K0=2.D0*PI/CR%LAMBDA
        STH=CR%LAMBDA/(2.D0*DHKL) ! sin(thetaB)
        ! if (ABS(STH).ge.1.D0) return
        TH=ASIN(STH)
        TH1=SIN(TH-CR%SGN*CR%CHI) ! angle on the side of FH1
        TH2=SIN(TH+CR%SGN*CR%CHI) ! angle on the side of FH2
    ! set horizontal curvatures
        if (CR%HFOC%AUTO) THEN
          LL=CR%HFOC%FOC1*CR%HFOC%FOC2
          if (LL.gt.0.D0) then
            CR%HFOC%RHO=(CR%HFOC%FOC1*TH2+CR%HFOC%FOC2*TH1)/LL*0.5D-3
          endif
        endif

    ! set vertical curvatures
        if (CR%VFOC%AUTO) THEN
          LL=CR%VFOC%FOC1*CR%VFOC%FOC2*COS(CR%CHI)
          if (LL.gt.0.D0) then
            CR%VFOC%RHO=(CR%VFOC%FOC1+CR%VFOC%FOC2)/LL/STH*0.5D-3
          endif
        endif
      endif

    ! set dependences
      if ((CR%HFOC%NSEG.eq.1).and.(CR%HFOC%RHO.NE.0.D0)) then
        CR%HFOC%BENT=.true.
      endif
      if (CR%HFOC%BENT) THEN
        CR%RH=CR%HFOC%RHO
      else
        CR%RH=0.D0
      endif
      if (CR%VFOC%BENT) THEN
        CR%RV=CR%VFOC%RHO
      else
        CR%RV=0.D0
      endif

    ! initialize A3D data
      ! get gaps
      DS=(/CR%HFOC%GAP,CR%VFOC%GAP,CR%SWGAP/)
      ! get segment number vector
      NS=(/CR%HFOC%NSEG,CR%VFOC%NSEG,CR%NSW/)
      ! get curvatures
      R=(/CR%HFOC%RHO,CR%VFOC%RHO,0.D0/)
      call ARRAY3D_INIT(CR%A3D,CR%FRAME%SIZE,DS,NS,R,0.D0,CR%HFOC%STACK,CR%VFOC%STACK,CR%HFOC%BENT,CR%VFOC%BENT)

      end SUBROUTINE XTAL_AUTOFOCUS


!----------------------------------------------------------------
      SUBROUTINE XTAL_ORIENT(CR,IERR)
! Adjust crystal in reflection position for nominal diffraction vector (ORI)
! NOTE: This procedure only aligns G-vector with respect to exit axis,
! but Bragg position requires also G^2 + 2.K.G = 0, which is not checked here
! NOTE2: Axis angle is adjusted only if lambda<2.dhkl
!-----------------------------------------------------------------
      TYPE (XTAL) :: CR
      integer,intent(out) :: IERR
      TYPE(TCR_ITEM) :: CDATA
      real(kind(1.D0)) :: KI(3),KF(3),Q(3),GON(3),THETA,G(3),DHKL
      ierr=1
      call GET_CR_ITEM(CR%ITAB,CDATA)
      DHKL=GET_DHKL(CDATA%LAT,CR%ORI,CR%T)
      if ((CR%LAMBDA.gt.0.D0).and.(CR%LAMBDA.le.2.D0*DHKL)) then
        THETA=2.D0*ASIN(CR%LAMBDA/2.D0/DHKL)*CR%SGN
      !  write(*,*) 'XTAL_ORIENT theta=',THETA*180/PI,' DHKL=',DHKL
    ! flat-cone rotation is applied when required
        if (CR%FCONE) then
          call SetAxisAngles(CR%FRAME,(/0.D0,-ABS(THETA),0.D0/))
        else
          call SetAxisAngles(CR%FRAME,(/THETA,0.D0,0.D0/))
        endif
        ierr=0
      else
        call SetAxisAngles(CR%FRAME,(/0.D0,0.D0,0.D0/))
        call SetGonioAngles(CR%FRAME,(/0.D0,0.D0,0.D0/))
        return
      endif
    ! adjust crystal in Bragg position
      KI=(/0.D0,0.D0,1.D0/)
      call M3xV3(-CR%FRAME%MEXI,CR%FRAME%REXI,KI,KF)
      Q=KF-KI
      call M3xV3(1,CR%BMAT,CR%ORI,G)
      call CalGonioAngles(Q,G,.false.,GON)
      call SetGonioAngles(CR%FRAME,GON)
      end SUBROUTINE XTAL_ORIENT


!-----------------------------------------------------------------------------
      SUBROUTINE XTAL_GLOC(CR,XR,R,G)
! Calculate G-vector at given position in the array
! Assumes R,G in local coordinates of the whole array
! INPUT:
!   XR  .. reflection data
!   R   .. local position
! OUTPUT:
!   G   .. local G-vector
!-----------------------------------------------------------------------------
      TYPE (XTAL),intent(in) :: CR
      TYPE (XREF),intent(in) :: XR
      real(KIND(1.D0)), intent(in)  :: R(3)
      real(KIND(1.D0)), intent(out) :: G(3)
      real(KIND(1.D0)) :: TR(3),RLOC(3,3),R1(3),G1(3)
      integer :: idx(3),MLOC
      call ARRAY3D_GET_INDEX(CR%A3D,R,idx)
      call ARRAY3D_SEG_TRANS(CR%A3D,idx,TR,RLOC,MLOC)
      call R_TO_SEGMENT(TR,RLOC,MLOC,R,R1)
      call XTAL_SEG_GLOC(XR,R1,G1)
      call R_TO_ARRAY((/0.D0,0.D0,0.D0/),RLOC,MLOC,G1,G)
      end SUBROUTINE XTAL_GLOC

!-----------------------------------------------------------------------------
      SUBROUTINE XTAL_SEG_GLOC(XR,R,G)
! Calculate local G-vector
! Assumes R,G in local coordinates of a single segment
! INPUT:
!   XR  .. reflection data
!   R   .. local position
! OUTPUT:
!   G   .. local G-vector
!-----------------------------------------------------------------------------
      TYPE (XREF),intent(in) :: XR
      real(KIND(1.D0)), intent(in)  :: R(3)
      real(KIND(1.D0)), intent(out) :: G(3)
      real(KIND(1.D0)) :: zs,dz,DG(3)
      ! neutral surface
        zs=0.5*(XR%RHO(1)*R(1)**2+XR%RHO(2)*R(2)**2)
      ! distance from the neutral surface
        dz=R(3)-zs
      ! local deviation from nominal G
        call M3xV3(1,XR%DG_DR,(/R(1),R(2),dz/),DG)
        G=XR%G+DG
      end SUBROUTINE XTAL_SEG_GLOC

! -------------------------------------------
      subroutine GET_DKK_LIMIT(CR,K,XR,DKMIN,DKMAX)
! get limits as dK/K for Bragg condition inside the array of crystals
! K ... Ki, incident wave-vector in local coordinates
! -------------------------------------------
      TYPE (XTAL),intent(in) :: CR
      real(kind(1.D0)),intent(in) :: K(3)
      TYPE (XREF),intent(in) :: XR
      real(kind(1.D0)),intent(out) :: DKMIN,DKMAX
      integer :: i,j,m
      real(kind(1.D0)) :: DELTA,R(3),G(3),GG,KdotG
      integer,parameter :: SG(2)=(/-1,1/)

      DKMIN=1.D30
      DKMAX=-1.D30
      ! check deviation from Bragg condition in the corners of the bounding box
      ! and choose max/min values
      do i=1,2
        R(1)=SG(i)*CR%A3D%SA(1)
        do j=1,2
          R(2)=SG(j)*CR%A3D%SA(2)
          do m=1,2
            R(3)=SG(m)*CR%A3D%SA(3)
            call XTAL_GLOC(CR,XR,R,G)
            GG = V3xV3(G,G)
            KdotG=V3xV3(K,G)
            if (abs(KdotG).gt.1.D-8) then
              DELTA = -0.5D0*GG/KdotG - 1.D0
              DKMIN=min(DKMIN,DELTA)
              DKMAX=max(DKMAX,DELTA)
            endif
          enddo
        enddo
      enddo
      if (DKMIN.GE.DKMAX) then
      ! no reflection is possible
        DKMIN=1.D30
        DKMAX=1.D30
      else
      ! bound by the pre-selected k-range
        DKMIN=max(DKMIN,-abs(CR%DKK))
        DKMAX=min(DKMAX,abs(CR%DKK))
      endif
      end subroutine GET_DKK_LIMIT

! -------------------------------------------------------
      real(kind(1.D0)) function GET_XREF(CR,HKL,K,applyLimits,XR)
! get XR = reflection data for given HKL and Ki
! K ... Ki, incident wave-vector in local coordinates
! return deviation from Bragg condition as delta_k/k
! if applyLimits=true, check that results fits within limits calculated
! by GET_DKK_LIMIT. Otherwise use CR%DKK as limits. If result does not fit
! within the limits, set XR%QML=0
! -------------------------------------------------------
      TYPE (XTAL),intent(in) :: CR
      integer,intent(in) :: HKL(3)
      real(kind(1.D0)),intent(in) :: K(3)
      logical, intent(in) :: applyLimits
      TYPE (XREF),intent(out) :: XR
      TYPE(TCR_ITEM) :: CDATA
      integer :: i,j
      real(kind(1.D0)) :: GG,DELTA,A(3),AA,DHKL0,KdotG,DKMIN,DKMAX,dif
      real(kind(1.D0)) :: K0,K1(3),KaG(3),DGK(3),KGK,RT
      GET_XREF=0.D0
1     format(a,': ',8(2x,G11.5))
      call XREF_CLEAR(XR)
! get crystal data
      if (CR%ITAB.le.0) return
      call GET_CR_ITEM(CR%ITAB,CDATA)
    ! assign HKL
      XR%HKL=HKL
    ! DHKL and reflectivity
      XR%DHKL=GET_DHKL(CDATA%LAT,HKL,CR%T)

      XR%QML=GET_QML(CDATA%LAT,HKL,CR%T)
    ! get diffracton vector in local coordinates
      do i=1,3
        A(i)=1.D0*HKL(i)
      enddo
      call M3xV3(1,CR%BMAT,A,XR%G)
    ! account for thermal expansion
      DHKL0=GET_DHKL(CDATA%LAT,HKL,298.D0)
      do i=1,3
        XR%G(i)=XR%G(i)*DHKL0/XR%DHKL
      enddo
    ! get vector GAMMA perpendicular to G in (K,G) plane, normalized
      call V3cV3(K,XR%G,A)
      call V3cV3(A,XR%G,XR%GAMMA)
      AA=SQRT(XR%GAMMA(1)**2+XR%GAMMA(2)**2+XR%GAMMA(3)**2)
      do i=1,3
        XR%GAMMA(i)=XR%GAMMA(i)/AA
      enddo

    ! G-gradient
      XR%DG_DR=0.D0
    ! cylindrical bending
      if (CR%RH.NE.0.D0) then
        XR%DG_DR(1,1)=-CR%RH*XR%G(3)
        XR%DG_DR(1,3)=CR%RH*XR%G(1)
        XR%DG_DR(3,1)=CR%RH*XR%G(1)
        XR%DG_DR(3,3)=-CR%POI*CR%RH*XR%G(3)
      endif
      if (CR%RV.NE.0.D0) then
       ! incl. vertical bending
        XR%DG_DR(2,2)=XR%DG_DR(2,2)-CR%RV*XR%G(3)
        XR%DG_DR(2,3)=CR%RV*XR%G(2)
        XR%DG_DR(3,2)=CR%RV*XR%G(2)
        XR%DG_DR(3,3)=XR%DG_DR(3,3)-CR%POI*CR%RV*XR%G(3)
      endif
    ! neutral surface curvature
      XR%RHO(1)=CR%RH
      XR%RHO(2)=CR%RV

! experimenal: torsion along x
      RT=0.D0
      ! RT=0.1D-3
      if (RT.NE.0.D0) then
        XR%DG_DR(1,2)=XR%DG_DR(1,2)-RT*XR%G(3)
        XR%DG_DR(1,3)=XR%DG_DR(1,3)+RT*XR%G(2)
        XR%DG_DR(2,1)=XR%DG_DR(1,2)
        XR%DG_DR(3,1)=XR%DG_DR(1,3)
      endif

    ! temperature gradient
      AA=0.D0
      do i=1,3
        A(i)=CDATA%LAT%ALPHA(i)
        AA=AA+ABS(A(i))
      enddo
      if (AA.gt.0.D0) then
        DO i=1,3
          do j=1,3
            XR%DG_DR(i,j) = XR%DG_DR(i,j) - A(j)*CR%DT(j)*XR%G(i)
          enddo
        enddo
      ENDIF

    ! get deviation from Bragg condition as d_k/k
      GG = V3xV3(XR%G,XR%G)
      KdotG=V3xV3(K,XR%G)
      if (abs(KdotG).gt.1.D-8*GG) then
      !  DELTA = (0.5D0*GG + V3xV3(K,XR%G))/SQRT(GG)
        DELTA = -0.5D0*GG/KdotG - 1.D0
      else
        DELTA=1.D30
      endif
      ! difference w.r.t. orientation hkl
      dif=0.D0
      do i=1,3
        dif=dif+abs(hkl(i)-CR%ORI(i))
      enddo

    ! check required limits for dK/K
    ! apply limits only if HKL <> ORI
    !  if (dif.gt.1.D-1) then
      if (applyLimits) then
        call GET_DKK_LIMIT(CR,K,XR,DKMIN,DKMAX)
        if ((DELTA.lt.DKMIN).or.(DELTA.GT.DKMAX)) XR%QML=0.D0
        if (DKMIN*DKMAX.gt.0.D0) XR%QML=0.D0
      else
        if (ABS(DELTA).gt.ABS(CR%DKK)) XR%QML=0.D0
      endif
    !  endif
      ! if ((XR%QML>0.D0).and.applyLimits) write(*,1) '   GET_XREF k-limits ',HKL,DELTA,DKMIN,DKMAX,XR%QML

      ! get gradient of dk/k
      ! GRADKK= 2*(K+G).DG_DR.K/|G|^2/|k|
      ! where K is such that |K+G|^2=|K|^2
      if ((DELTA>-1.D0).and.(DELTA<1.D30)) then
        K0=SQRT(V3xV3(K,K))*(1.D0+DELTA)
        K1=K*(1.D0+DELTA)
        DO I=1,3
          KaG(I)=K1(I)+XR%G(I)
        ENDDO
        call M3xV3(1,XR%DG_DR,K1,DGK)
        KGK=V3xV3(KaG,DGK)
        XR%GRADKK=2.D0*KGK/GG/K0
      else
        XR%GRADKK=0.D0
      endif


    !  if ((HKL(1).eq.1).and.(HKL(2).eq.1).and.(HKL(3).eq.1)) then
    !    write(*,*) 'GET_XREF hkl=',HKL,' DELTA=',DELTA,' DKK=',CR%DKK
    !    write(*,1) '    G=',XR%G, 'K=',K,' DELTA=',DELTA
    !  endif
      GET_XREF=DELTA
      end function GET_XREF


!--------------------------------------------------------
      real(kind(1.D0)) FUNCTION XTAL_GETMI(CR,K0)
! Calculates total attenuation cross-section for given K0=2.PI/LAMBDA
! assume that SIGABS,SIGINC and SIGTDS have previously been calculated
! for CR%LAMBDA and CR%T
C--------------------------------------------------------
      TYPE (XTAL),intent(in) :: CR
      real(kind(1.D0)),intent(in) :: K0
      real(kind(1.D0)) :: Z,lambda,lam2
      Z=1.D0
      lambda=2.D0*PI/K0
    ! correction for lambda <> CR%LAMBDA (approximation only)
      lam2=lambda**2
      Z=CR%SIGTDS(3)/lam2 + CR%SIGTDS(4)/lam2**2
      XTAL_GETMI=CR%SIGABS*lambda + CR%SIGINC + CR%SIGTDS(1)*lambda + CR%SIGTDS(2)*(1.D0-exp(-Z))
      END FUNCTION XTAL_GETMI

!-------------------------------------------------------------------------
      subroutine XTAL_GET_MAXHKL(CR,HKL)
! return HKL(3) = limits of h,k,l indices for reflections search
!--------------------------------------------------------------------------
      TYPE (XTAL),intent(in) :: CR

      integer,intent(out) :: HKL(3)
      TYPE(TCR_ITEM) :: CDATA
      real(kind(1.D0)):: LAMBDA
      integer :: i
      HKL=0
      if (CR%ITAB.le.0) return
      call GET_CR_ITEM(CR%ITAB,CDATA)
    ! shortest wavelength from the search range
      LAMBDA=MAX(1.D-1,CR%LAMBDA*(1.D0-CR%DKK))
      HKL(1)=INT(2.D0*GET_DHKL(CDATA%LAT,(/1,0,0/),CR%T)/LAMBDA)
      HKL(2)=INT(2.D0*GET_DHKL(CDATA%LAT,(/0,1,0/),CR%T)/LAMBDA)
      HKL(3)=INT(2.D0*GET_DHKL(CDATA%LAT,(/0,0,1/),CR%T)/LAMBDA)
    ! apply preset limits
      do i=1,3
        HKL(i)=MIN(CR%MAXHKL,HKL(i))
      enddo
      end subroutine XTAL_GET_MAXHKL

!---------------------------------------------------------
      subroutine XTAL_SORT_NODES(CR)
! sorting of all nodes for given crystal ICR
! sort by DHKL (descending) using Heapsort algorithm
!---------------------------------------------------------
      TYPE (XTAL),intent(in) :: CR
      TYPE(PRLNODE) :: PNODE
      TYPE(TCR_ITEM) :: CDATA
      real(kind(1.D0)) :: D(MAX_NODES)
      integer :: i,N
    !  write(*,*) 'XTAL_SORT_NODES icr=',CR%ICR,' ID=',trim(CR%FRAME%ID),' NCR=',XREF_NCR,' NN=',RLNODE_GET_LAST(CR%ICR)
  ! check index validity
      if (CR%ICR.gt.XREF_NCR) return
    !  write(*,*) 'XTAL_SORT_NODES icr=',CR%ICR,' ID=',trim(CR%FRAME%ID),' ITAB=',CR%ITAB
      if (CR%ITAB.le.0) return
  ! clear ord arrays
      RLNODE_ORD(1:MAX_NODES,CR%ICR)=0
      D=0.D0
  ! get lattice data
      call GET_CR_ITEM(CR%ITAB,CDATA)
  ! get number of associated nodes
      N=min(RLNODE_GET_LAST(CR%ICR),MAX_NODES)
    !  write(*,*) 'XTAL_SORT_NODES icr=',CR%ICR,' ID=',trim(CR%FRAME%ID),' NN=',N
      if (N.le.0) return
  ! fill array D with DHKL values
      do i=1,N
        call RLNODE_GET_ITEM(CR%ICR,i,PNODE)
        if (associated(PNODE%N)) then
          D(i)=GET_DHKL(CDATA%LAT,PNODE%N%HKL,CR%T)
        endif
      enddo
  ! 1st node should be (0,0,0) with dhkl=0
  ! ensure that the (0,0,0) node is the first in the order
      if (D(1).le.0.D0) D(1)=1.D30
  ! sort according to DHKL
      call HEAPSORT(D,N,-1,RLNODE_ORD(1,CR%ICR))
      end subroutine XTAL_SORT_NODES

!---------------------------------------
      SUBROUTINE CAL_PREF(CR,KI)
! find primary reflections
!---------------------------------------
      TYPE (XTAL) :: CR
      real(kind(1.D0)),intent(in) :: KI(3)
      integer :: NDEF,NNODES
      TYPE(RLNODE) :: NEWNODE
      CR%REFLIST=0
      CR%REFLIST_N=0
      CR%REFLIST_NB=0
      call XREF_CLEAR_CR(CR%ICR)
      call RLNODE_CLEAR_CR(CR%ICR)
      call RLNODE_CLEAR(NEWNODE)
      call RLNODE_ADD(NEWNODE,CR%ICR,NODE_ERR)
      NDEF=0
      NNODES=RLNODE_GET_LAST(CR%ICR)
1     format('Search reflection for ',a,' k0=',G12.6,' ',$)
    !  write(*,*) 'CAL_PREF crystal=',CR%ICR,' NNODES=',NNODES
      write(*,1) trim(CR%FRAME%ID),SQRT(KI(1)**2+KI(2)**2+KI(3)**2)
      write(*,*)
      do while ((NODE_ERR.eq.0).and.(NDEF.lt.NNODES))
        NDEF=NDEF+1
        call CAL_SREF(CR,KI,NDEF)
        NNODES=RLNODE_GET_LAST(CR%ICR)
      !  write(*,*) 'CAL_PREF NDEF=',NDEF,' NNODES=',NNODES,' err=',NODE_ERR
      enddo
      write(*,*) ' done'
      call XTAL_SORT_NODES(CR)
      end SUBROUTINE CAL_PREF

!---------------------------------------
      SUBROUTINE CAL_SREF(CR,KI,IPR)
! find secondary reflections for given primary reflection
!---------------------------------------
      TYPE (XTAL) :: CR
      real(kind(1.D0)),intent(in) :: KI(3)
      integer,intent(in) :: IPR
      TYPE (XREF) :: XR
      integer :: h,k,l,IREF
      real(kind(1.D0)) :: KG(3),KK0,DKK,KI0 !,DKMIN,DKMAX
      ! TYPE(RLNODE) :: NEWNODE
      TYPE(PRLNODE) :: PNODE
      integer :: HMIN(3),HMAX(3),HKL(3),MAXHKL(3),cnt
      logical :: skip
      integer,parameter :: progfrq=5000  ! frequency for showig progress dots
  ! get pointer to the primary node
    !  write(*,*) 'CAL_SREF primary node=',IPR
      call RLNODE_GET_ITEM(CR%ICR,IPR,PNODE)
      if (.not.associated(PNODE%N)) return

    !  write(*,*) 'level=',PNODE%N%LEVEL,' hkl=','node=',PNODE%N%HKL
      if (PNODE%N%LEVEL.le.CR%UMWEG) then
      ! scan h,k,l and look for reflections
        call XTAL_GET_MAXHKL(CR,MAXHKL)
        HMIN=-MAXHKL
        HMAX=+MAXHKL
      else
        return
      endif
      KI0=SQRT(KI(1)**2+KI(2)**2+KI(3)**2)
    ! get KG = incident K for this node
      if (PNODE%N%LEVEL.eq.0) then
        KK0=KI0
        KG=KI
      else
        ! get |K| for this branch
        KK0=2.D0*PI/PNODE%N%LAMBDA
        ! get diffracton vector leading to this node in local coordinates
        DKK=GET_XREF(CR,PNODE%N%HKL,KI*KK0/KI0,.FALSE.,XR)
        KG=KI*KK0/KI0+XR%G
      ! write(*,*) 'CAL_SREF (',PNODE%N%LEVEL,') ',PNODE%N%HKL,KK0,(sqrt(V3xV3(KG,KG))-KK0)/KI0,DKK
      endif
      cnt=0
      do h=HMIN(1),HMAX(1)
        do k=HMIN(2),HMAX(2)
          do l=HMIN(3),HMAX(3)
            if (NODE_ERR.ne.0) exit
          ! exclude (0,0,0)
            if ((abs(h)+abs(k)+abs(l)).gt.0) then
          ! only backward reflections are allowed if umweg level is reached
              skip=((PNODE%N%LEVEL.gt.0).and.(PNODE%N%LEVEL.ge.CR%UMWEG))

              if (skip) then
                HKL=(/h,k,l/)+PNODE%N%HKL
                skip=((abs(HKL(1))+abs(HKL(2))+abs(HKL(3))).gt.0)
              endif
              if (.not.skip) then
                cnt=cnt+1
                if (MOD(cnt,progfrq).eq.0)  write(*,"('.',$)")
                DKK=GET_XREF(CR,(/h,k,l/),KG,(PNODE%N%LEVEL.gt.0),XR)
                if (XR%QML.gt.1.D-10) then
                  ! for primary reflection, select only those within defined k-shell
                  if (PNODE%N%LEVEL.eq.0) then
                  !  write(*,*) 'CAL_SREF ',DKK,CR%DKK
                  !  if (abs(DKK).gt.abs(CR%DKK)) exit
                    KK0=KI0*(1.D0+DKK)
                  !  write(*,*) 'CAL_SREF accepted ',(/h,k,l/)+PNODE%N%HKL,KK0
                  endif

                 ! call GET_DKK_LIMIT(CR,KG,XR,DKMIN,DKMAX)
                 ! write (*,*) '   CAL_SREF (',(/h,k,l/)+PNODE%N%HKL,') limits= ',DKMIN,DKMAX,DKK

                  IREF=RLNODE_ASSIGN_XREF(PNODE%N,CR%ICR,XR,KK0,NODE_ERR)
                !  if (IREF.gt.0) then
                 !   HKL=(/-h,-k,-l/)+CR%HKL
                   ! selected=((abs(HKL(1))+abs(HKL(2))+abs(HKL(3))).gt.0)
                 ! endif
                endif
              endif
            endif
          enddo
        enddo
      enddo
      call RLNODE_SORT(PNODE%N)
      end SUBROUTINE CAL_SREF


!---------------------------------------------------------
      logical function XTAL_ADJUST_LAMBDA(INST)
! Set OBJ%LAMBDA for selected reflection if
! a) selected reflection <> (0,0,0)
! b) selected reflection is primary one (it is on the list of the 1st node)
!---------------------------------------------------------
      integer,intent(in) :: INST
      integer :: i,j
      TYPE(XTAL),POINTER :: OBJ
      TYPE(PRLNODE) :: NODE
      logical :: RES
        RES=.false.
        if (XTAL_isValid(INST)) then
        OBJ => AXTALS(INST)%X
        I=RLNODE_GET(OBJ%ICR,OBJ%HKL)
        if (I.gt.1) then
          call RLNODE_GET_ITEM(OBJ%ICR,1,NODE) ! primary node (0,0,0)
          do j=1,NODE%N%NP
            if (NODE%N%NEXT(j).eq.I) then
              OBJ%LAMBDA=NODE%N%LAM(j)
              RES=.true.
              exit
            endif
          enddo
        endif
        endif
        XTAL_ADJUST_LAMBDA=RES
      end function XTAL_ADJUST_LAMBDA

!---------------------------------------------------------
      SUBROUTINE XTAL_INP(INST,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... XTAL structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(XTAL),POINTER :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      ! character(256) :: SARG
      integer :: LR,i
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      if (.not.XTAL_isValid(INST)) return
      OBJ => AXTALS(INST)%X
      LR=0
      NARG=0
       ! call subclass _INP ?
      i=INDEX(PNAME,'.')
      ! if (INDEX(PNAME,'VFOC').gt.0) write(*,*) '  XTAL_INP [',trim(PNAME),']'
      if (i.gt.1) then
        select case (PNAME(:i-1))
        case('VFOC')
          call FOCARRAY_INP(OBJ%VFOC,PNAME(i+1:),ARG,NARG)
        case('HFOC')
          call FOCARRAY_INP(OBJ%HFOC,PNAME(i+1:),ARG,NARG)
          select case(PNAME(i+1:))
          case('RHO')
            OBJ%RH=OBJ%HFOC%RHO
          end select
        end select
        return
      endif
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call FIELD2ARRAY(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
          call XTAL_INP_R(OBJ,trim(ARG%ID),NUM,LR)
        case(ID_FIELD_SELECT)
          call XTAL_INP_SEL(OBJ,trim(ARG%ID),ARG,LR)
        case(ID_FIELD_CLASS)
          ! accept, but do nothing
        case default
          write(*,*) 'XTAL_INP: undefined field type '//trim(ARG%DEF%TID)
      end select
      NARG=LR
      END SUBROUTINE XTAL_INP


!---------------------------------------------------------
      SUBROUTINE XTAL_OUT(INST,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(XTAL),POINTER :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR,i
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      ! character(256) :: SARG
      if (.not.XTAL_isValid(INST)) return
      OBJ => AXTALS(INST)%X
      LR=0
      NARG=0
    ! call subclass _OUT ?
      i=INDEX(PNAME,'.')
      ! if (INDEX(PNAME,'VFOC').gt.0) write(*,*) '  XTAL_OUT [',trim(PNAME),']'
      if (i.gt.1) then
        select case (PNAME(:i-1))
        case('VFOC')
          call FOCARRAY_OUT(OBJ%VFOC,PNAME(i+1:),ARG,NARG)
        case('HFOC')
          call FOCARRAY_OUT(OBJ%HFOC,PNAME(i+1:),ARG,NARG)
        end select
        return
      endif
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call XTAL_OUT_R(OBJ,trim(ARG%ID),NUM,LR)
          call ARRAY2FIELD(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
        case(ID_FIELD_SELECT)
          call XTAL_OUT_SEL(OBJ,trim(ARG%ID),ARG,LR)
        case(ID_FIELD_CLASS)
          ! accept, but do nothing
        case default
          write(*,*) 'XTAL_OUT: undefined field type '//trim(ARG%DEF%TID)
      end select
      NARG=LR
      END SUBROUTINE XTAL_OUT


!---------------------------------------------------------
      SUBROUTINE XTAL_INP_R(OBJ,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... REAL*8 parameter values
! OUTPUT
!    OBJ     .... XTAL structure
!    NARG    .... number of items read from ARG
C---------------------------------------------------------
      TYPE(XTAL),intent(out) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(in) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR,I
      LR=1
      SELECT CASE (trim(PNAME))
! crystal ID and index to TABLE_CRYSTALS
          CASE('CRID')
            call GetEnumString('CRNAMES',NINT(ARG(1)),OBJ%CRNAME)
            OBJ%ITAB=GET_CR_INDEX(OBJ%CRNAME)
! sign of the take-off angle
! NOTE: ARG contains ordinal value, i.e. ARG=0 for SGN=-1
          CASE('SGN'); OBJ%SGN=INT2SGN(NINT(ARG(1)))
! main reflection (used for orienting) [hkl]
          CASE('ORI');OBJ%ORI(1:3)=ARG(1:3);LR=3
! reference vectors [hkl]
          CASE('A')
            do i=1,3
              OBJ%A(i)=NINT(ARG(i))
            enddo
            LR=3
          CASE('B')
            do i=1,3
              OBJ%B(i)=NINT(ARG(i))
            enddo
            LR=3
! cutting angle [deg]
          CASE('CHI');   OBJ%CHI=ARG(1)*deg
! vertical cutting angle [deg]
          CASE('PSI');   OBJ%PSI=ARG(1)*deg
! Poisson number
          CASE('POISS');   OBJ%POI=ARG(1)
! temperature [K]
          CASE('T');   OBJ%T=ARG(1)
! bending curvature [m^-1]
          CASE('RHO');   OBJ%RH=ARG(1)*1.D-3;
          CASE('RHV');   OBJ%RV=ARG(1)*1.D-3;
! gap between sandwich slabs
          CASE('SWGAP ');   OBJ%SWGAP=ARG(1)
! number of sandwich slabs
          CASE('NSW ');   OBJ%NSW=max(1,min(128,NINT(ARG(1))))
! linear temperature gradient [K/mm]
          CASE('DT');OBJ%DT(1:3)=ARG(1:3);LR=3
! delta_k/k limit used in searching for reflections
          CASE('DKK');   OBJ%DKK=ARG(1)
! automatically adjust in reflection position
          CASE('AUTOADJ');   OBJ%AUTOADJ=(ARG(1).eq.1.D0)
! trace Renninger reflections
          CASE('UMWEG');   OBJ%UMWEG=max(0,min(6,NINT(ARG(1))))
! max. reflection index of searched reflection
          CASE('MAXHKL');   OBJ%MAXHKL=max(1,min(20,NINT(ARG(1))))
          CASE DEFAULT; CALL FRAME_INP_R(OBJ%FRAME,trim(PNAME),ARG(1),LR)
      END SELECT
      NARG=LR
      END SUBROUTINE XTAL_INP_R

!---------------------------------------------------------
      SUBROUTINE XTAL_INP_SEL(OBJ,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD structure of SELECT type
! OUTPUT
!    OBJ     .... XTAL structure
!    NARG    .... number of items read from ARG
C---------------------------------------------------------
      TYPE(XTAL),intent(out) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD),intent(in) :: ARG
      integer,intent(out) :: NARG
      integer :: LR,ierr
      ! type(PRLNODE) :: PNODE
      LR=1
      ierr=0
      SELECT CASE (trim(PNAME))
        CASE('REF')
          if (ARG%NP.ge.0) then
            if (ARG%EORD.ge.0) then
              READ(ARG%SFIELD(ARG%EORD+1),*,IOSTAT=ierr) OBJ%HKL
            else
              ierr=1
            endif
          !  write(*,*) 'XTAL_INP_SEL err=',ierr,' ord=',ARG%EORD,' HKL=',OBJ%HKL
            if (ierr.ne.0) OBJ%HKL=(/0,0,0/)
          else
            OBJ%HKL=(/0,0,0/)
          endif
      END SELECT
      NARG=LR
      END SUBROUTINE XTAL_INP_SEL

C---------------------------------------------------------
      SUBROUTINE XTAL_OUT_R(OBJ,PNAME,ARG,NARG)
! output data from OBJ to ARG for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    OBJ     .... XTAL structure
! INPUT
!    SARG    .... string parameter values
!    NARG    .... number of items read from ARG
C---------------------------------------------------------
      TYPE(XTAL),intent(in) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(out) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR,ITYPE,IORD
      LR=1
      SELECT CASE (trim(PNAME))
! crystal ID
          CASE('CRID')
         ! get enum type corresponding to this variable
            call GetEnumOrd('CRNAMES',trim(OBJ%CRNAME),ITYPE,IORD)
            ARG(1)=IORD
! sign of the take-off angle
! NOTE: ARG contains ordinal value, i.e. ARG=0 for SGN=-1
          CASE('SGN'); ARG(1)=SGN2INT(OBJ%SGN)
! main reflection (used for orienting) [hkl]
          CASE('ORI');ARG(1:3)=OBJ%ORI(1:3);LR=3
! reference vectors [hkl]
          CASE('A');ARG(1:3)=OBJ%A(1:3);LR=3
          CASE('B');ARG(1:3)=OBJ%B(1:3);LR=3
! cutting angle [deg]
          CASE('CHI');   ARG(1)=OBJ%CHI/deg
! vertical cutting angle [deg]
          CASE('PSI');   ARG(1)=OBJ%PSI/deg
! Poisson number
          CASE('POISS');   ARG(1)=OBJ%POI
! temperature [K]
          CASE('T');   ARG(1)=OBJ%T
! bending curvature [m^-1]
          CASE('RHO');   ARG(1)=OBJ%RH*1.D3
          CASE('RHV');   ARG(1)=OBJ%RV*1.D3
! gap between sandwich slabs
          CASE('SWGAP ');   ARG(1)=OBJ%SWGAP
! number of sandwich slabs
          CASE('NSW ');   ARG(1)=OBJ%NSW
! linear temperature gradient [K/mm]
          CASE('DT');ARG(1:3)=OBJ%DT(1:3);LR=3
! delta_k/k limit used in searching for reflections
          CASE('DKK');   ARG(1)=OBJ%DKK
! automatically adjust in reflection position
          CASE('AUTOADJ');   ARG(1)=LOG2INT(OBJ%AUTOADJ)
! trace Renninger reflections
          CASE('UMWEG');   ARG(1)=OBJ%UMWEG
! max. reflection index of searched reflection
          CASE('MAXHKL');   ARG(1)=OBJ%MAXHKL
          CASE DEFAULT; CALL FRAME_OUT_R(OBJ%FRAME,trim(PNAME),ARG(1),LR)
      END SELECT
      NARG=LR
      END SUBROUTINE XTAL_OUT_R


!---------------------------------------------------------
      SUBROUTINE XTAL_OUT_SEL(OBJ,PNAME,ARG,NARG)
! input data from OBJ to SARG for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    OBJ     .... XTAL structure
! OUTPUT
!    ARG     .... TFIELD structure of SELECT type
!    NARG    .... number of items read from ARG
C---------------------------------------------------------
      TYPE(XTAL),intent(in) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR,ierr
      ! character(512) :: SNODES
      ! character(64) :: S
      LR=1
      SELECT CASE (trim(PNAME))
        CASE('REF')
          call XREF_LIST_NODES(OBJ%ICR,OBJ%HKL,ARG,IERR)
        !  write(*,*) 'XTAL_OUT_SEL err=',ierr,' ord=',ARG%EORD,' HKL=',OBJ%HKL
          if (ierr.ne.0) LR=0
      END SELECT
      NARG=LR
      END SUBROUTINE XTAL_OUT_SEL

      end MODULE XTALS