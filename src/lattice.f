!//////////////////////////////////////////////////////////////////////
!////  $Id: lattice.f,v 1.15 2015/06/24 22:32:45 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.15 $
!////     $Date: 2015/06/24 22:32:45 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Describes crystal lattice and related transformations
!////
!////////////////////////////////////////////////////////////////////////
      MODULE LATTICE
      use CONSTANTS
      use SIMATH
      use XMLINFO
      use XMLUTILS
      use TABLE_ATOMS
      implicit none
      private

      integer,parameter :: MAX_CELL_ATOMS = 32

      TYPE TCR_LATTICE; SEQUENCE
        integer :: NAT                   !   number of atoms in unit cell
        INTEGER ::  AT(MAX_CELL_ATOMS)   !   proton numbers of compound atoms (indices to TABLE_ATOMS)
        real(kind(1.D0)) ::  AX(3)       !   axes [A]
        real(kind(1.D0)) ::  ANG(3)      !   angles [deg]
        integer :: NEP                   !   number of equivalent positions
        real(kind(1.D0)) ::  EPOS(3,4)   !   equivalent positions
        real(kind(1.D0)) ::  APOS(3,MAX_CELL_ATOMS)   !  atom positions
        real(kind(1.D0)) ::  BM(4,MAX_CELL_ATOMS)   !   magnetic scattering lengths [fm] associated with atoms
        real(kind(1.D0)) ::  ALPHA(3)    ! linear thermal expansion coefficient [K^-1]
        real(kind(1.D0)) ::  thetaD      ! Debye temperature (K)
      ! calculated fields
        real(kind(1.D0)) ::  VOL         !   cell volume [A^3]
        real(kind(1.D0)) ::  G(3,3)      !   metric tensor for rec. lattice
      end TYPE TCR_LATTICE

      interface GET_DHKL
        MODULE procedure GET_DHKL_R
        MODULE procedure GET_DHKL_I
      end interface GET_DHKL

      public TCR_LATTICE,LAT_INIT,LAT_CLEAR
      public GET_DHKL,GET_DHKL_R,GET_DHKL_I,CalBMatrix,GET_FHKL2,GET_QML,AdotB
      public startLATTICE,dataLATTICE,endLATTICE,getLATTICE,LAT_DW

      ! local instance for XML handlers
      TYPE(TCR_LATTICE) :: myLAT
      contains

!---------------------------------------------------------
      subroutine LAT_CLEAR(LAT)
!---------------------------------------------------------
      TYPE(TCR_LATTICE),intent(out) :: LAT
        LAT%AT=0
        LAT%AX=0.D0
        LAT%ANG=0.D0
        LAT%VOL=0.D0
        LAT%NEP=0
        LAT%NAT=0
        LAT%EPOS=0.D0
        LAT%APOS=0.D0
        LAT%BM=0.D0
        LAT%ALPHA=0.D0
        LAT%thetaD=0.D0
      end subroutine LAT_CLEAR

!-------------------------------------------------------------
      subroutine LAT_INIT(LAT)
! calculate dependent fields in LAT (metric tensor, cell volume, ...)
!-------------------------------------------------------------
      TYPE (TCR_LATTICE) :: LAT
      call GET_METRIC(LAT)
      end subroutine LAT_INIT

!--------------------------------------------------------
      subroutine LAT_DW(LAT,ATOM,T,B0,BT)
! return  B0 and BT terms for calculation of the Debye-Waller factor
! Wd = (B0 + BT)(sin(thetaB)/lambda)^2 = 1/4*(B0 + BT)*d^-2
! for given atom
! Our definition is from [1]
! ASSUMPTIONS:
!  - cubic crystals
!  - Debye model with single Debye frequency for all atoms
! [1] A.Freund, Nucl.Inst.Meth. 213 (1983) 495-501.
!--------------------------------------------------------
      TYPE (TCR_LATTICE),intent(in)  :: LAT
      TYPE(TTAB_ATOMS),intent(in)  :: ATOM
      real(kind(1.D0)),intent(in) :: T
      real(kind(1.D0)),intent(out) :: B0,BT
      real(kind(1.D0)) :: x
    ! B0=3h^2/(2.k_B.thetaD.M)  ATOM%M is in atomic mass units
      B0=2872.558/ATOM%M/LAT%thetaD
    ! unphysical case T>>thetaD, but let's handle it
      if (LAT%thetaD.lt.1.D-3*T) then
        BT=1.D3
      else
        x=LAT%thetaD/T
        BT=4.D0*B0*DW_PHI(x)/x**2
      endif
      end subroutine LAT_DW

!-------------------------------------------------------------
      real(kind(1.D0)) function GET_FHKL2(LAT,HKL,T)
! calculate |Fhkl|^2 [fm^2] for given reflection inidices
! includes Debye-Waller factor for given temperature
!-------------------------------------------------------------
      TYPE (TCR_LATTICE),intent(in)  :: LAT
      integer,intent(in) :: HKL(3)
      real(kind(1.D0)),intent(in) :: T
      integer :: i,j,k
      TYPE(TTAB_ATOMS) :: ATOM
      real(kind(1.D0)) :: ImF,ReF,phi,FHKL2,Wd,DHKL,B0,BT
      real(kind(1.D0)),parameter :: EPS=1.D-15
!      logical :: dbg
!      dbg=((HKL(1).eq.1).and.(HKL(2).eq.1).and.(HKL(3).eq.1))
      FHKL2=0.D0
      ImF=0.D0
      ReF=0.D0
      DHKL=GET_DHKL(LAT,HKL,T)
    ! sum over atoms
      do i=1,LAT%NAT
        call GET_ATOM_DATA(LAT%AT(i),ATOM)
        call LAT_DW(LAT,ATOM,T,B0,BT)
        Wd=exp(-(B0+BT)/(2.D0*DHKL)**2)
      ! sum over equivalent positions
        do j=1,LAT%NEP
          phi=0.D0
          do k=1,3
            phi=phi+HKL(k)*(LAT%APOS(k,i)+LAT%EPOS(k,j))
          enddo
          phi=2.D0*PI*phi
          ImF=ImF + ATOM%BC*sin(phi)*Wd
          ReF=ReF + ATOM%BC*cos(phi)*Wd
        enddo
      enddo
      FHKL2=ImF**2 + ReF**2
      if (FHKL2.lt.EPS) FHKL2=0.D0
      GET_FHKL2=FHKL2
!      if (dbg) then
!        write(*,*) 'GET_FHKL2 ',Wd,B0,BT,FHKL2
!      endif
      end function GET_FHKL2

!-----------------------------------------------------------------
      real(kind(1.D0)) function GET_QML(LAT,HKL,T)
! calculate 4*PI*(F*dhkl/V0)**2 [ A^-1 cm^-1] (Maier-Leibnitz reflectivity)
! includes thermal expansion and DW factor
!-----------------------------------------------------------------
      TYPE (TCR_LATTICE)  :: LAT
      integer,intent(in) :: HKL(3)
      real(kind(1.D0)),intent(in) :: T
      if (LAT%VOL.gt.0.D0) then
        GET_QML=4.0D-2*pi*GET_FHKL2(LAT,HKL,T)*(GET_DHKL(LAT,HKL,T)/LAT%VOL)**2
      else
        GET_QML=0.D0
      endif
      end function GET_QML

!-------------------------------------------------------------
      real(kind(1.D0)) function GET_DHKL_R(LAT,HKL,T)
! calculate dhkl [A] for given plane and temperature
! use metric tensor calculated previously by GET_METRIC
!-------------------------------------------------------------
      TYPE (TCR_LATTICE)  :: LAT
      real(kind(1.D0)),intent(in) :: HKL(3)
      real(kind(1.D0)),intent(in) ::T
      integer :: i,j,k
      real(kind(1.D0)) :: DHKL,A(3)
      integer,parameter :: IT(6)=(/2,3,3,1,1,2/)
      DHKL=0.D0
    ! use symmetry of G(i,j)=G(j,i)
      do i=1,3
        A(i)=(1.D0+LAT%ALPHA(i)*(T-298.D0))
      enddo
      do i=1,3
        DHKL=DHKL+HKL(i)**2*LAT%G(i,i)/A(i)**2
        j=IT(2*i-1)
        k=IT(2*i)
        DHKL=DHKL+2.D0*HKL(j)*HKL(k)*LAT%G(j,k)/(A(j)*A(k))
      enddo
      if (DHKL.gt.0.D0) DHKL=1.D0/SQRT(DHKL)
      GET_DHKL_R=DHKL
      end function GET_DHKL_R

!-------------------------------------------------------------
      real(kind(1.D0)) function GET_DHKL_I(LAT,HKL,T)
! calculate dhkl [A] for given plane
! use metric tensor calculated previously by GET_METRIC
!-------------------------------------------------------------
      TYPE (TCR_LATTICE)  :: LAT
      integer,intent(in) :: HKL(3)
      real(kind(1.D0)),intent(in) ::T
      GET_DHKL_I=GET_DHKL_R(LAT,1.D0*HKL,T)
      end function GET_DHKL_I

!-------------------------------------------------------------
      real(kind(1.D0)) function GET_CELLVOL(LAT,T)
! calculate unit cell volume in [A^3]
!-------------------------------------------------------------
      TYPE (TCR_LATTICE),intent(in)  :: LAT
      real(kind(1.D0)),intent(in) ::T
      integer :: k
      real(kind(1.D0)) :: COSA(3),VOL,CC
      CC=0.D0
      VOL=1.D0
      do k=1,3
        COSA(k)=cos(LAT%ANG(k)*deg)
        CC=CC+COSA(k)**2
        VOL=VOL*LAT%AX(k)*(1.D0+LAT%ALPHA(k)*(T-298.D0))
      enddo
      if (CC.GT.0.D0) then
        VOL=VOL*SQRT(1.D0-CC+2.D0*COSA(1)*COSA(2)*COSA(3))
      else
        VOL=0.D0
      endif
      GET_CELLVOL=VOL
      end function GET_CELLVOL

!-------------------------------------------------------------
      subroutine GET_METRIC(LAT)
! calculate metric tensor for reciprocal lattice at T=298 K
!-------------------------------------------------------------
      TYPE (TCR_LATTICE) :: LAT
      real(kind(1.D0)) :: G(3,3)
      integer :: k
      real(kind(1.D0)) :: AX(3),SINA(3),COSA(3),VOL2
      LAT%VOL=GET_CELLVOL(LAT,298.D0)
      if (LAT%VOL.gt.0.D0) then
        VOL2 = LAT%VOL**2
        do k=1,3
          COSA(k)=cos(LAT%ANG(k)*deg)
          SINA(k)=sin(LAT%ANG(k)*deg)
          AX(k)=LAT%AX(k)
        enddo
        G(1,1)=(AX(2)*AX(3)*SINA(1))**2/VOL2
        G(2,2)=(AX(1)*AX(3)*SINA(2))**2/VOL2
        G(3,3)=(AX(1)*AX(2)*SINA(3))**2/VOL2
        G(1,2)=AX(1)*AX(2)*AX(3)**2*(COSA(1)*COSA(2)-COSA(3))/VOL2
        G(2,1)=G(1,2)
        G(1,3)=AX(1)*AX(3)*AX(2)**2*(COSA(1)*COSA(3)-COSA(2))/VOL2
        G(3,1)=G(1,3)
        G(2,3)=AX(2)*AX(3)*AX(1)**2*(COSA(2)*COSA(3)-COSA(1))/VOL2
        G(3,2)=G(2,3)
      else
        G=0.D0
      endif
      LAT%G = G
      end subroutine GET_METRIC


!--------------------------------------------------------------------
      SUBROUTINE CalBMatrix(LAT,A,B,CHI,PSI,T,BMAT)
!  Calculate B-matrix for given two reference hkl-vectors
! B-matrix transforms lattice vectors [hkl] to local coordinates in [A^-1]
!  A,B are lattice vectors [hkl] in the reference plane
!  CHI,PSI are cutting angles [rad] around y and x axes
!--------------------------------------------------------------------
      TYPE (TCR_LATTICE) :: LAT
      real(KIND(1.D0)),intent(in) :: A(3),B(3),CHI,PSI,T
      real(KIND(1.D0)),intent(out) :: BMAT(3,3)
      real(KIND(1.D0)), PARAMETER :: EPS=1.0D-10
      real(KIND(1.D0)) :: C1,C2,C3,BB(3),COSB(3),SINB(3),XBB(3,3),RCHI(3,3),RPSI(3,3),R(3,3)
      REAL*8 U(3,3),V1(3),V2(3),V3(3),W(3,3)
      INTEGER :: I,J,K
      integer,parameter :: IT(6)=(/2,3,3,1,1,2/)

      BMAT=0.D0
      IF(LAT%VOL.LE.0.) return
    ! get rec. lattice axes and angles
      DO i=1,3
        BB(i)=2.D0*PI*SQRT(LAT%G(i,i))/(1.D0+LAT%ALPHA(i)*(T-298.D0))
        j=IT(2*i-1)
        k=IT(2*i)
        COSB(i)=LAT%G(j,k)/SQRT(LAT%G(j,j)*LAT%G(k,k))
        SINB(i)=SQRT(1.-COSB(I)**2)
      ENDDO

C XBB(i,j) are projections of rec. lat. base vectors on the orthonormal base:
C Let (a,b,c) and (a*,b*,c*) are direct and reciprocal lattice base vectors
C assume a parallel to a*
C XBB(i,1) ... parallel to a*
C XBB(i,2) ... parallel to c x a*
C XBB(i,3) ... parallel to c
C i.e. the columns are a*,b*,c* in the new orthonormal base
      XBB(1,1)=BB(1)
      XBB(2,1)=0.D0
      XBB(3,1)=0.D0
      XBB(1,2)=BB(2)*COSB(3)
      XBB(2,2)=BB(2)*SINB(3)
      XBB(3,2)=0.D0
      XBB(1,3)=BB(3)*COSB(2)
      XBB(2,3)=-BB(3)*SINB(2)*COS(LAT%ANG(1)*deg)
      XBB(3,3)=2.D0*PI/LAT%AX(3)/(1.D0+LAT%ALPHA(3)*(T-298.D0))

! convert A,B to the orthonormal system XBB:
      call M3xV3(1,XBB,A,V1)
      call M3xV3(1,XBB,B,V2)
!      write(*,1) 'A  ',A
!      write(*,1) 'V1 ',V1
! get V3 perpendicular to V1,V2
      call V3cV3(V1,V2,V3)
! get V2 perpendicular to V3,V1
      call V3cV3(V3,V1,V2)
! get norms of V1,V2,V3
      C1=SQRT(V1(1)**2+V1(2)**2+V1(3)**2)
      C2=SQRT(V2(1)**2+V2(2)**2+V2(3)**2)
      C3=SQRT(V3(1)**2+V3(2)**2+V3(3)**2)
      IF (C1*C2*C3.LT.EPS) return ! check scattering plane
! rows of U(i,j) are normalized vectors V1..V3 in XBB coordinates
! => U is transformation matrix from the XBB base to the orthonormal base V1,V2,V3
! NOTE: V1//z, V2//x and V3//y for CHI=0,PSI=0, hence the following permutaion of rows
      DO I=1,3
        U(1,I)=V2(I)/C2
        U(2,I)=V3(I)/C3
        U(3,I)=V1(I)/C1
      ENDDO
!1     format(a,3(1x,G10.4))
!      write(*,1) 'B(1,i) ',U(1,1:3)
!      write(*,1) 'B(2,i) ',U(2,1:3)
!      write(*,1) 'B(3,i) ',U(3,1:3)
! get rotation matrix for CHI,PSI rotation
    ! RCHI rotates vectors by CHI around y-axis
      call MK_ROT3(2,-CHI,RCHI)
    ! RPSI rotates vectors by PSI around x-axis
      call MK_ROT3(1,-PSI,RPSI)
    ! R transforms vectors from V1,V2,V3 coordinates to x,y,z coordinates (coordinates fixed to the crystal body)
      call M3xM3(1,RPSI,RCHI,R)
! construct B-matrix: sequence of transformations (hkl) -> (XBB) -> (V1,V2,V3) -> (x,y,z)
      call M3xM3(1,R,U,W)
      call M3xM3(1,W,XBB,BMAT)

      END SUBROUTINE CalBMatrix


!----------------------------------------------------
      REAL(KIND(1.D0)) FUNCTION AdotB(LAT,A,B)
! dot product of hkl vectors (result in [A^-2])
!----------------------------------------------------
      TYPE (TCR_LATTICE) :: LAT
      REAL(KIND(1.D0)),intent(in) :: A(3),B(3)
      REAL(KIND(1.D0)) :: C(3)
      call M3xV3(-1,LAT%G,A,C)
      AdotB=C(1)*B(1)+C(2)*B(2)+C(3)*B(3)
      end FUNCTION AdotB


!----------------------------------------------------------------------
      subroutine GET_POSITION(LINE,XPOS,error)
! read a line with atom position, which can be expressed by fractions
! e.g. LINE="0  1/3  1/3"
!----------------------------------------------------------------------
      character(*),intent(in) :: LINE
      real(KIND(1.D0)),intent(out) :: XPOS(3)
      logical,intent(out) :: error
      character(256) :: STRARRAY
      character*32 :: S
      integer :: NA,i,j,IS,IL
      logical :: isFloat
      real(KIND(1.D0)) :: A1,A2

      call READ_STRARRAY(LINE,'|',STRARRAY,NA)
      if (NA.ge.3) then
        do i=1,3
          call FINDSTRPAR(STRARRAY,'|',i,IS,IL)
          if (IL.gt.0) then
            S=STRARRAY(IS:IS+IL-1)
            j=INDEX(S,'/')
            if (j.gt.0) then
              if ((j.lt.2).or.(.not.IsFloat(S(1:j-1),A1))) A1=0.D0
              if ((j.gt.IL-1).or.(.not.IsFloat(S(j+1:IL),A2))) A2=0.D0
              XPOS(i)=A1/A2
            else if (IsFloat(S,A1)) then
              XPOS(i)=A1
            else
              XPOS(i)=0.D0
            endif
          endif
        enddo
        error=.false.
      else
        error=.true.
      endif
      end subroutine GET_POSITION


!***************************************
! XML handlers
!***************************************

!--------------------------------------------------------
      subroutine startLATTICE(tag,attribs,error)
!--------------------------------------------------------
        character(len=*)                  :: tag
        character(len=*), dimension(:,:)  :: attribs
        logical                           :: error
        character(32) :: S
        integer :: IA,NA
        logical :: b
        real(kind(1.D0)) :: ARG(4)
    ! new table
        if (error) return
        select case(trim(tag))
        case('LATTICE')
          call LAT_CLEAR(myLAT)
        case('ATOM')
          call readAttrib(attribs,'id',S,error)
          if (.not.error) then
            IA=GET_ATOM_INDEX(trim(S))
            if (IA.le.0) then
              error=.true.
              call MSG_ERROR('LATTICE','Unknown atom symbol: '//trim(S),0,1)
            else if (myLAT%NAT.lt.MAX_CELL_ATOMS) then
              myLAT%NAT=myLAT%NAT+1
              myLAT%AT(myLAT%NAT)=IA
              myLAT%APOS(1:3,myLAT%NAT)=(/0.D0,0.D0,0.D0/)
            ! check "mag" attribute
              call readAttrib(attribs,'mag',S,b)
              if (.not.b) then
                call GETLINPARG(S,ARG,4,NA)
                if (NA.ge.4) then
                  myLAT%BM(1:4,myLAT%NAT)=ARG(1:4)
                endif
              endif
            else
              error=.true.
              call INT2STR(MAX_CELL_ATOMS,S)
              call MSG_ERROR('LATTICE','Too many atom positions, max='//trim(S),0,1)
            endif
          else
            call MSG_ERROR('LATTICE','invalid attributes in '//trim(tag),0,1)
          endif

        case('EQPOS')
          if (myLAT%NEP.lt.4) then
            myLAT%NEP=myLAT%NEP+1
            myLAT%EPOS(1:3,myLAT%NEP)=(/0.D0,0.D0,0.D0/)
          else
            error=.true.
            call INT2STR(MAX_CELL_ATOMS,S)
            call MSG_ERROR('LATTICE','Too many equivalent positions, max=4',0,1)
          endif
        end select
      end subroutine startLATTICE

!--------------------------------------------------------
      subroutine dataLATTICE(tag, xmldata, error )
! XML handler for CRTABLE
!--------------------------------------------------------
        character(len=*)                  :: tag
        character(len=*)                  :: xmldata
        logical                           :: error
        real(KIND(1.D0)) :: ARG(3)
        integer :: NA
      ! process only non-empty data inside class entries
        if (error) return
        if (len_trim(xmldata).le.0) return
        select case(trim(tag))
          case('AX')
            call GETLINPARG(trim(xmldata),myLAT%AX,3,NA)
            error=(NA.lt.3)
          case('ANG')
            call GETLINPARG(trim(xmldata),myLAT%ANG,3,NA)
            error=(NA.lt.3)
          ! equivalent positions
          case('EQPOS')
            call GET_POSITION(trim(xmldata),ARG,error)
            if (.not.error) myLAT%EPOS(1:3,myLAT%NEP)=ARG(1:3)
          case('ATOM')
            call GET_POSITION(trim(xmldata),ARG,error)
            if (.not.error) myLAT%APOS(1:3,myLAT%NAT)=ARG(1:3)
          case('ALPHA')
            call GETLINPARG(xmldata,ARG,3,NA)
            if (NA.ge.3) myLAT%ALPHA(1:3)=ARG(1:3)*1.D-6
          case('THETAD')
            read(xmldata,*) myLAT%thetaD
        end select
      end subroutine dataLATTICE

!--------------------------------------------------------
      subroutine endLATTICE(tag,error)
!--------------------------------------------------------
      character(len=*)                  :: tag
      logical                           :: error
      if (error) return
      select case(trim(tag))
      case('LATTICE')
        call LAT_INIT(myLAT)
      end select
      end subroutine endLATTICE

!--------------------------------------------------------
      type(TCR_LATTICE) function getLATTICE()
!--------------------------------------------------------
      getLATTICE=myLAT
      end function getLATTICE

      end MODULE LATTICE