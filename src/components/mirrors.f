!//////////////////////////////////////////////////////////////////////
!////  $Id: mirrors.f,v 1.6 2019/08/15 15:02:07 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2016, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.6 $
!////     $Date: 2019/08/15 15:02:07 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Ray-tracing subroutines for a double coated supermirror
!////  Handles reflection and transmission through both the coating and substrate
!////  Default: NiTi coating on Si wafer
!////
!////  MIRROR is not a standalone component, it is used as a part of other components
!////
!////////////////////////////////////////////////////////////////////////
      MODULE MIRRORS
      use TRACELIB
      use MIRROR_TABLE
      use TABLE_CRYSTALS
      use TRACINGDATA
      use RNDGEN
      use SIMATH
      use FRAMES
      implicit none
      private

      TYPE  MIRROR; SEQUENCE
        real(kind(1.D0)) :: MC  ! m-value
        real(kind(1.D0)) :: MU  ! default absorption coefficient [1/mm/Ang]
        real(kind(1.D0)) :: W,H1,H2,L    ! width, height entry, height exit, length
        real(kind(1.D0)) :: RLOC(3,3)    ! orientation matrix
        real(kind(1.D0)) :: STA(3)       ! position of the origin
        real(kind(1.D0)) :: CSIG1,CSIG2  ! coating absorption cross-sections (SIGMA=CSIG1=CSIG2*lambda) [1/mm]
        real(kind(1.D0)) :: CT  ! coating thickness [mm]
        INTEGER :: MATER        ! substrate material: (0) absorbing (1) MU*lambda (2) Si (3) Al2O3
        integer :: IC           ! index to reflectivity tables
        INTEGER :: TRFRONT      ! transparent front surface of lamellae, if >0
        integer :: LSIDE,RSIDE  ! > 0 if the right/left side is coated
        integer :: MLOC         ! tag indicating the rotation
! calculated fields
        real(kind(1.D0)) :: MU2PI,CMU2PI  ! transmission probability = EXP(-MU2PI*DT)
        real(kind(1.D0)) :: CRHO,SRHO  ! 4*pi*rho in coating and substrate
        real(kind(1.D0)) :: LAMBDA  ! wavelength
        real(kind(1.D0)) :: DH  ! slope of top/bottom edges
        integer :: SIDE     ! left/right = 1, top/bottom = 2
      END TYPE MIRROR

      ! TODO
      ! Publish after testing, the module is not yet in use ...
      !public MIRROR
      !public MIRROR_INIT,MIRROR_GO,MIRROR_GO_INSIDE,MIRROR_GETREF
      !public MIRROR_SET_COORD,MIRROR_SET_MATER,MIRROR_SET_SIZE

      contains

!---------------------------------------------------------------------
      subroutine MIRROR_SET_MATER(OBJ,IC,MATER,MC,MU)
! set material parameters
! IC ... coatig lookup table index
! MATER ... substrate material index
! MC ... default m-value for reflectivity edge
! MU ... default substrate absorption coefficient
! MC and MU can be overriden with the lookup table data
!---------------------------------------------------------------------
      TYPE(MIRROR) :: OBJ
      integer, intent(in) :: IC,MATER
      real(kind(1.D0)), intent(in) :: MC,MU
      real(kind(1.D0)) :: par(7)
      OBJ%IC=IC
      OBJ%MATER=MATER
      OBJ%MC=MC
      OBJ%MU=MU
      ! get 4*pi*rho for the substrate
      select case (OBJ%MATER)
      case(3) ! Al2O3
        OBJ%SRHO=4.D0*PI*5.715D-6
      case default ! Si
        OBJ%SRHO=4.D0*PI*2.073D-6
      end select
      ! get coating parameters
      if (OBJ%IC>0) then
        call MIRROR_GET_PARAM(OBJ%IC,par,OBJ%MC)
        OBJ%CT=par(1)! thickness [mm]
        OBJ%CRHO=4.D0*PI*par(2) ! 4PI*average scattering density [Ang^-2]
        ! att. cross-sections S1 and S2 (Stot=S0+lambda*S1) [mm-1,mm-1Ang-1]
        OBJ%CSIG1=par(3)
        OBJ%CSIG2=par(4)
        OBJ%MC=par(5) ! edge m-value (needed to calculate transmission)
      else
      ! set default
        OBJ%CT=0.003
        OBJ%CRHO=4.D0*PI*4.13D-6
        OBJ%CSIG1=0.1011
        OBJ%CSIG2=0.0211
      endif

      end subroutine MIRROR_SET_MATER

!---------------------------------------------------------------------
      subroutine MIRROR_SET_SIZE(OBJ,W,H1,H2,L)
! set size blade dimensions
!---------------------------------------------------------------------
      TYPE(MIRROR) :: OBJ
      real(kind(1.D0)), intent(in) :: W,H1,H2,L
      OBJ%W=W
      OBJ%H1=H1
      OBJ%H2=H2
      OBJ%L=L
      OBJ%DH=0.5D0*(OBJ%H2-OBJ%H1)/OBJ%L
      end subroutine MIRROR_SET_SIZE

!---------------------------------------------------------------------
      subroutine MIRROR_SET_COORD(OBJ,x,angle,side)
! set transformation arrays
! x ... position
! a ... tilt angle
! side ... left/right (1), top/bottom (2)
! called at the begining of each trace
!---------------------------------------------------------------------
      TYPE(MIRROR) :: OBJ
      real(kind(1.D0)), intent(in) :: x,angle
      integer, intent(in) :: side
      real(kind(1.D0)) :: GON(3)
      real(kind(1.D0)),parameter :: PIHALF=0.5D0*PI
      OBJ%SIDE=side
      select case(side)
      case(2)
        GON=(/-PIHALF,PIHALF,PIHALF+angle/)
        OBJ%STA=(/0.D0,x,0.D0/)
      case default
        GON=(/angle,0.D0,0.D0/)
        OBJ%STA=(/x,0.D0,0.D0/)
      end select
      call CalRotMatrix(GON,OBJ%RLOC,OBJ%MLOC)
      end subroutine MIRROR_SET_COORD

!---------------------------------------------------------------------
      subroutine MIRROR_INIT(OBJ,K0)
! initialize wavelength dependent parameters  = absorption coefficient
! called at the begining of each trace
! K0 = 2 PI / lambda
!---------------------------------------------------------------------
      TYPE(MIRROR) :: OBJ
      real(kind(1.D0)), intent(in) :: K0

      OBJ%LAMBDA=TWOPI/K0
! absorption coefficient in mm^-1.A^-1)
! transmission probability = EXP(-MU2PI*DT), where DT=path/k
      select case (OBJ%MATER)
      case(1)
        OBJ%MU2PI=OBJ%MU*TWOPI
      case(2)
        OBJ%MU2PI=GET_CR_ABS(MU_Si,K0)
      case(3)
        OBJ%MU2PI=GET_CR_ABS(MU_Al2O3,K0)
      case default
        OBJ%MU2PI=1.D10
      end select
      OBJ%CMU2PI=(OBJ%CSIG1+OBJ%CSIG1*OBJ%LAMBDA)*K0
      end subroutine MIRROR_INIT

!---------------------------------------------------------------------
      real(kind(1.D0)) function MIRROR_GETREF(qT,IC,MC,S)
!---------------------------------------------------------------------
      real(kind(1.D0)), intent(in) :: qT,MC,S
      integer, intent(in) :: IC
      real(kind(1.D0)) :: PR
      PR=MIRROR_REF(IC,MC,qT,S)
      MIRROR_GETREF=PR
      end function MIRROR_GETREF


!---------------------------------------------------------------------
      subroutine MIRROR_GO_INSIDE(OBJ,IRES,TC)
! tracing inside a double-coated mirror blade
! neglects transport through coating and refraction
! TC: apply coordinate transformation
! IRES = exit code:
! -1  ... absorbed
!  0  ... no cross-section
!  1  ... exit left
!  2  ... exit right
!  3  ... exit rear edge
!---------------------------------------------------------------------
      TYPE(MIRROR) :: OBJ
      integer,intent(out) :: ires
      logical, intent(in) :: TC
      real(kind(1.D0)) :: V(3),RN,PR,PT
      real(kind(1.D0)) :: tin,tout,qT0,PATH,SGX,SGZ,Z1,S
      integer :: iloop
      logical :: CONT
      TYPE(NEUTRON) :: NEU
      integer,parameter :: SPOL(2)=(/2,1/)

      IRES=-1

! transform neutron coordinates R0,K0 to the mirror's local frame, R,K
      if (TC) then
! transform neutron coordinates R0,K0 to the mirror's local frame, R,K
        NEU=NEUT
        CALL V3AV3(-1,NEU%R,OBJ%STA,V)
        CALL M3XV3(OBJ%MLOC,OBJ%RLOC,V,NEUT%R)
        CALL M3XV3(OBJ%MLOC,OBJ%RLOC,NEU%K,NEUT%K)
      endif
  ! distinguish forward and backward tracing
      ! directions
      SGZ=SIGN(1.D0,NEUT%K(3))
      ! front and exit edge
      if (SGZ>0) then
        Z1=OBJ%L
      else
        Z1=0.D0
      endif
      iloop=0
      CONT=.true.
      PATH=0.D0
      S=NEUT%S(SPOL(OBJ%SIDE))
      do while (CONT)
        iloop=iloop+1
        SGX=SIGN(1.D0,NEUT%K(1))
      ! check cross-times
        call BORDER_WALL(NEUT%R(1),NEUT%K(1),OBJ%W,tin,tout)
      ! something is wrong = no cross-section
        if ((tin>tout).or.(tout<=0.D0)) then
          IRES=0
          goto 10
        endif
        call Transport(tout)
        if (abs(NEUT%R(2))>0.5*OBJ%H1+OBJ%DH*NEUT%R(3)) then
        ! exit through top/bottom edges = absorb
          goto 10
        endif
        PATH=PATH+tout
        RN=1.D0*RAN1()
      ! exit at the end?
        if (abs(NEUT%R(3)-Z1)<=1.D-5) then
        ! decide between absorption and transmission
          PT=exp(-OBJ%MU2PI*PATH)
          IF (RN>PT) IRES=3 ! transmitted
          CONT=.false.
        else
        ! calculate qT
          qT0=abs(NEUT%K(1))
        ! calculate reflectivity
          PR=MIRROR_GETREF(qT0,OBJ%IC,OBJ%MC,S)
          IF (RN<=PR) THEN
        ! reflection and continue
            NEUT%K(1)=-NEUT%K(1)
          else
        ! exit through the coating
            if (SGX>0.D0) then
              IRES=1 ! left
            else
              IRES=2 ! right
            endif
            CONT=.false.
          endif
        endif
        CONT=(iloop<50)
      enddo

! final action
! transform back to the original frame
10    if (TC) then
        CALL M3XV3(-OBJ%MLOC,OBJ%RLOC,NEUT%R,V)
        CALL V3AV3(1,V,OBJ%STA,NEU%R)
        CALL M3XV3(-OBJ%MLOC,OBJ%RLOC,NEUT%K,NEU%K)
        NEUT=NEU
      endif
      if (IRES<=0) then
        NEUT%P=0.D0
      endif

      end subroutine MIRROR_GO_INSIDE

!---------------------------------------------------------------------
      subroutine MIRROR_GO(OBJ,IRES)
! tracing through a double-coated mirror blade
! IRES = exit code:
! -1  ... absorbed
!  0  ... no cross-section
!  1  ... exit left
!  2  ... exit right
!  3  ... exit rear edge
!---------------------------------------------------------------------
      TYPE(MIRROR) :: OBJ
      integer,intent(out) :: ires
      real(kind(1.D0)) :: V(3),RN,PR,PTC(0:3),PATH(3),Z0,Z1
      real(kind(1.D0)) :: tin,tout,qT0,qT02,qTC,qTS,DTC,DTS,DX,DZ,DT,SGX,SGZ,S
      integer :: i
      logical :: edge
      TYPE(NEUTRON) :: NEU
      integer,parameter :: SPOL(2)=(/2,1/)

      IRES=-1
! transform neutron coordinates R0,K0 to the mirror's local frame, R,K
      NEU=NEUT
      CALL V3AV3(-1,NEU%R,OBJ%STA,V)
      CALL M3XV3(OBJ%MLOC,OBJ%RLOC,V,NEUT%R)
      CALL M3XV3(OBJ%MLOC,OBJ%RLOC,NEU%K,NEUT%K)

! proceed to the entry surface
      call BORDER_BLADE(NEUT%R,NEUT%K,OBJ%W,OBJ%H1,OBJ%H2,OBJ%L,tin,tout)
! no cross-section => exit with IRES=0
      if (tin>=tout) then
        IRES=0
        goto 10
      endif
      call Transport(tin)

  ! distinguish forward and backward tracing
      ! directions
      SGZ=SIGN(1.D0,NEUT%K(3))
      SGX=SIGN(1.D0,NEUT%K(1))
      ! front and exit edge
      if (SGZ>0) then
        Z0=0.D0
        Z1=OBJ%L
      else
        Z0=OBJ%L
        Z1=0.D0
      endif

      ! entry through the front edge ?
      if (abs(NEUT%R(3)-Z0)<=1.D-5) then
       ! front edge
        if (OBJ%TRFRONT>0) then
          call MIRROR_GO_INSIDE(OBJ,IRES,.false.)
        endif
        goto 10
      endif
! calculate qT in vacuum
      qT0=abs(NEUT%K(1))
! calculate reflectivity
      S=NEUT%S(SPOL(OBJ%SIDE))
      PR=MIRROR_GETREF(qT0,OBJ%IC,OBJ%MC,S)
! decide on reflection or transmission
! if reflected => invert K(1) and exit
      RN=1.D0*RAN1()
      IF (RN<=PR) THEN
        IRES=1
        NEUT%K(1)=-NEUT%K(1)
        if (SGX<0.D0) then
          IRES=1 ! left
        else
          IRES=2 ! right
        endif
        goto 10
      endif

! calculate probability of transmission or absorption
  ! for qT<m*qTNi, coating transmission prob. is calculated from the reflection prob.
      if (qT0<=OBJ%MC*qTNi) then
        ! approximate absorption above the reflectivity curve for m<MC
        ! approximate reflectivity + absorption by extrapolating with r'(MC-m)=1-r(MC+m)
        PR=1.D0-MIRROR_GETREF(2.D0*OBJ%MC*qTNi-qT0,OBJ%IC,OBJ%MC,S)
        ! absorb everything below the new PR
        if (RN<=PR) goto 10
      endif
! get absorption probabilities in each layer
  ! calculate qT for coating and substrate
      qT02=qT0**2
      ! not reflected below critical angle = absorb
      if ((qT02<=OBJ%CRHO).or.(qT02<=OBJ%SRHO)) goto 10
      qTS=SQRT(qT02-OBJ%SRHO)
      qTC=SQRT(qT02-OBJ%CRHO)
  ! cross times in coating and substrate
      DTC=OBJ%CT/qTC
      DTS=(OBJ%W-2.D0*OBJ%CT)/qTS
  ! cal. position and time at the exit side
      DX=OBJ%W
      DT=2.D0*DTC+DTS
! handle exit through the front edge
  ! z-path to the exit
      DZ=DT/NEUT%K(3)
  ! exit through the edge ?
      edge=.false.
      if (SGZ*(NEUT%R(3)+DZ)>SGZ*Z1) then
        edge=.true.
        DT=(Z1-NEUT%R(3))/NEUT%K(3)
        if (DT<DTC) then
          DX=DT*qTC
          PATH=(/DT,0.D0,0.D0/)
        else if (DT<DTC+DTS) then
          DX=OBJ%CT+(DT-DTC)*qTS
          PATH=(/DTC,DT,0.D0/)
        else if (DT<2.D0*DTC+DTS) then
          DX=OBJ%W-OBJ%CT+(DT-DTC-DTS)*qTC
          PATH=(/DTC,DTS,DT/)
        endif
      else
        PATH=(/DTC,DTS,DTC/)
      endif
      ! cummulative probabilities of interaction: reflection, 1st coating, substrate, 2nd coating
      PTC(0)=PR
      PTC(1)=PTC(0)+(1.D0-PTC(0))*exp(-PATH(1)*OBJ%CMU2PI)
      PTC(2)=PTC(1)+(1.D0-PTC(1))*exp(-PATH(2)*OBJ%MU2PI)
      PTC(3)=PTC(2)+(1.D0-PTC(2))*exp(-PATH(3)*OBJ%CMU2PI)
      ! decide what happens
      i=1
      do while((i<4).and.(RN>PTC(i)))
        i=i+1
      enddo
      ! effective k-vector for transport
      V(1)=SGX*DX/DT
      V(2:3)=NEUT%K(2:3)
      if (i==4) then ! transmission
      ! move to the other side and exit
        if (edge) then
          if (OBJ%TRFRONT>0) IRES=3 ! exit through edge side
        else if (SGX>0.D0) then
          IRES=1 ! left
        else
          IRES=2 ! right
        endif
        call TransportK(V,DT)
      else ! absorption
        select case(i)
        case(1) ! absorption in the 1st coating
          RN=(PTC(1)-PTC(0))*RAN1()
          DT=-LOG(1.D0-RN)/OBJ%CMU2PI
          DX=DT*qTC
        case(2) ! absorption in the substrate
          RN=(PTC(2)-PTC(1))*RAN1()
          DT=-LOG(1.D0-RN)/OBJ%MU2PI
          DX=OBJ%CT+DT*qTS
          DT=DT+DTC
        case(3) ! absorption in the 2nd coating
          RN=(PTC(3)-PTC(2))*RAN1()
          DT=-LOG(1.D0-RN)/OBJ%CMU2PI
          DX=OBJ%W-OBJ%CT+DT
          DT=DT+DTC+DTS
        end select
        V(1)=SGX*DX/DT
        call TransportK(V,DT)
      endif

! final action
! transform back to the original frame
10    CALL M3XV3(-OBJ%MLOC,OBJ%RLOC,NEUT%R,V)
      CALL V3AV3(1,V,OBJ%STA,NEU%R)
      CALL M3XV3(-OBJ%MLOC,OBJ%RLOC,NEUT%K,NEU%K)
      NEUT=NEU
      if (IRES<0) then
        NEUT%P=0.D0
      endif
      end subroutine MIRROR_GO


      end module MIRRORS