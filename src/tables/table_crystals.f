!//////////////////////////////////////////////////////////////////////
!////  $Id: table_crystals.f,v 1.15 2019/08/15 15:02:09 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.15 $
!////     $Date: 2019/08/15 15:02:09 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Dynamic table with crystal properties
!////  + related subroutines (e.g. structure factors)
!//////////////////////////////////////////////////////////////////////
      MODULE TABLE_CRYSTALS
      use CONSTANTS
      use TABLE_ATOMS
      use LATTICE
      use XMLINFO
      use XMLUTILS
      use SIMATH
      implicit none
      private

      real(kind(1.D0)),parameter :: MU_Si(4)=(/0.1018D0, 6.054D-3, 0.38D0, 0.D0/)
      real(kind(1.D0)),parameter :: MU_Al2O3(4)=(/0.2120D0, 8.11D-3, 0.16D0, 0.129D0/)

      integer,parameter :: MAX_CELL_ATOMS = 32
      integer,parameter :: MAX_NCR_TAB = 100

      TYPE TCR_ITEM; SEQUENCE
        character*8 ::  ID               ! id string
        character*32 :: name             ! name (descriptive)
        TYPE(TCR_LATTICE) :: LAT         ! lattice data
        real(kind(1.D0)) ::  C2          ! constant from the Freund's paper  [A^-2 meV^-1]
        real(kind(1.D0)) ::  ISO(3)      ! isotropic axis (h,k,l)
        real(kind(1.D0)) ::  XHAR        ! harmonic term
        real(kind(1.D0)) ::  XANH        ! anharmonic term
      end TYPE TCR_ITEM

! encapsulated pointer to TCR_ITEM
      type PCR_ITEM
        TYPE(TCR_ITEM),pointer :: C
      end TYPE PCR_ITEM

      TYPE(PCR_ITEM) :: CR_TAB(MAX_NCR_TAB)
      integer :: NCR_TAB=0

      ! local variable for XML import
      TYPE(TCR_ITEM) :: myITEM
      integer :: XML_LEVEL=-1
      integer,parameter :: XML_LEVEL_TABLE=0
      integer,parameter :: XML_LEVEL_DATA=1
      integer,parameter :: XML_LEVEL_LATTICE=2

      public startCRTABLE,dataCRTABLE,endCRTABLE,CLEAR_TABLE_CRYSTALS,GET_CR_INDEX
      public GET_CR_ITEM,TCR_ITEM,GET_CR_SIGMA,GET_CR_ABS
      public MU_Si,MU_Al2O3

      contains


!---------------------------------------------------------
      subroutine ADD_CRITEM(C)
!---------------------------------------------------------
      TYPE(TCR_ITEM),intent(in) :: C
        if (NCR_TAB.lt.MAX_NCR_TAB) then
          if (GET_CR_INDEX(C%ID).gt.0) then
            call MSG_ERROR('ADD_CRITEM','Duplicite crystal definition for '//trim(C%ID),0,1)
          else if (ALLOC_CRITEM(CR_TAB(NCR_TAB+1))) then
            NCR_TAB=NCR_TAB+1
            CR_TAB(NCR_TAB)%C=C
           ! lattice is initialized in its XML parser
           !  call LAT_INIT(CR_TAB(NCR_TAB)%C%LAT)
            if (CR_TAB(NCR_TAB)%C%LAT%VOL.le.0.D0) then
              call MSG_ERROR('ADD_CRITEM','Wrong unit cell definition for '//trim(C%ID),0,1)
            endif
          endif
        endif
      end subroutine ADD_CRITEM

!---------------------------------------------------------
      subroutine CLEAR_TABLE_CRYSTALS
!---------------------------------------------------------
        do while (NCR_TAB.gt.0)
          call DISPOSE_CRITEM(CR_TAB(NCR_TAB))
          NCR_TAB=NCR_TAB-1
        enddo
      end subroutine CLEAR_TABLE_CRYSTALS

!---------------------------------------------------------
      subroutine CLEAR_CRITEM(C)
!---------------------------------------------------------
      TYPE(TCR_ITEM),intent(out) :: C
        C%ID=' '
        C%NAME=' '
        C%C2=0.D0
        C%ISO=0.D0
        C%XHAR=1.D0
        C%XANH=0.D0
        call LAT_CLEAR(C%LAT)
      end subroutine CLEAR_CRITEM

!---------------------------------------------------------
      logical function ALLOC_CRITEM(C)
!---------------------------------------------------------
      TYPE(PCR_ITEM) :: C
      integer :: ierr
        allocate(C%C,STAT=ierr)
        if(ierr.eq.0) then
          call CLEAR_CRITEM(C%C)
        endif
        ALLOC_CRITEM=(ierr.eq.0)
      end function ALLOC_CRITEM

!---------------------------------------------------------
      subroutine DISPOSE_CRITEM(C)
!---------------------------------------------------------
      TYPE(PCR_ITEM) :: C
        if (associated(C%C)) then
          DEallocate(C%C)
          nullify(C%C)
        endif
      end subroutine DISPOSE_CRITEM

!---------------------------------------------------------
      integer function GET_CR_INDEX(ID)
! return structure with atom data from SC_TABLE by index
!---------------------------------------------------------
        character(*),intent(in) :: ID
        integer :: i,ia
        i=0
        ia=0
        do while ((i.lt.NCR_TAB).and.(ia.eq.0))
          i=i+1
          if (trim(ID).EQ.trim(CR_TAB(i)%C%ID)) ia=i
        enddo
        GET_CR_INDEX=ia
      end function GET_CR_INDEX

!---------------------------------------------------------
      subroutine GET_CR_ITEM(IDX,CR)
! return structure with atom data from SC_TABLE by index
!---------------------------------------------------------
        type(TCR_ITEM),intent(out) :: CR
        integer,intent(in) :: IDX
        if ((IDX.gt.0).and.(IDX.le.NCR_TAB)) then
          CR=CR_TAB(IDX)%C
        else
          call CLEAR_CRITEM(CR)
        endif
      end subroutine GET_CR_ITEM

!----------------------------------------------------
      REAL(kind(1.D0)) function GET_CR_ABS(P,K)
! Returns absorption coefficient for a single crystal
! for given material parameters (P) and wave-vector magnitude (K)
! Absorption formula:
!	mu = P(2)*lambda + P(1)*(1 - exp(-P(3)/lambda^2 - P(4)/lambda^4)
!	Parameters, Units:
!	P(1) = s_free [1/cm]
!	P(2) [1/cm/A]
!	P(3) [A^2]
!	P(4) [A^4]
!----------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: P(4),K
      REAL(KIND(1.D0)) :: lam,lam2,lam4,ex,mu
        lam=TWOPI/K
        lam2=lam**2
        lam4=lam2**2
        ex=P(3)/lam2 + P(4)/lam4
        ! mu in 1/cm
        mu=P(2)*lam + P(1)*(1 - exp(-ex))
        ! convert to K*MU in 1/mm/A
        GET_CR_ABS=mu*K*1.D-1
      end function GET_CR_ABS


!--------------------------------------------------------------
      subroutine GET_CR_SIGMA(CR,LAMBDA,T,SIGABS,SIGINC,SIGTDS)
! calculate macroscopic absorption and incoherent cross-sections
! INPUT:
!   CR    .. crystal data data
!   LAMBDA .. wavelength
!   T      .. temperature
! OUTOUT:
!   SIGABS ... absorption [mm^-1.A^-1]
!   SIGINC ... incoherent scattering [mm^-1]
!   SIGTDS ... coefficients for thermal scattering (see below)
!NOTE:
! SIGTDS has three elements:
! SIGTDS(1) = single phonon contribution (A_sph from [2]), [mm^-1.A^-1]
! SIGTDS(2) = multi-phonon contribution in low energy limit [mm^-1] = SIGMA_fa from [2]
! SIGTDS(3) = (B_0+B_T)<sin(theta/2)^2> = B from [2] = (B_0+B_T)*C2*E from [1] [A^2]
! SIGTDS(4) = D from [2] = [A^4]
! NOTE: we use additional weighting terms for harmonic and anharmonic contributions
! therefore SIGTDS(3) = XHAR*B, SIGTDS(4)=XANH*B^2, where B=C2*h^2/2/m_n (= B from [2])
! in theory, C2 should be = 4*m_n/h^2*(1/2 - 1/(3*M))
! [1] A.Freund, Nucl.Inst.Meth. 213 (1983) 495-501.
! [2] J. G. Barker et al., J. Appl. Cryst. 41 (2008) 1003-1008.
!--------------------------------------------------------------
      type(TCR_ITEM),intent(in)  :: CR
      real(kind(1.D0)),intent(in) :: LAMBDA,T
      real(kind(1.D0)),intent(out) :: SIGABS,SIGINC,SIGTDS(4)
      integer i,j
      TYPE(TTAB_ATOMS) :: ATOM
      real(kind(1.D0)) :: X,Xn,R,B(0:30),E,IFact,SIGB,SII,B0,BT,ksi

      call BERNOULLI(30,B)
      X=CR%LAT%thetaD/max(1.D-1,T)
      R=0.D0
      Ifact=1.D0
      Xn=1.D0/X
      DO I=0,30
        R=R+B(I)*Xn/(Ifact*(I+2.5D0))
        Xn=Xn*X
        Ifact=Ifact*(I+1)
      ENDDO
      SIGINC=0.D0
      SIGABS=0.D0
      SIGTDS=0.D0
      E=HSQOV2M*(2*PI/LAMBDA)**2
! sum over atoms
      do i=1,CR%LAT%NAT
        SIGB=0.D0
        SII=0.D0
        call GET_ATOM_DATA(CR%LAT%AT(i),ATOM)
        call LAT_DW(CR%LAT,ATOM,T,B0,BT)
      ! sum over equivalent positions
        do j=1,CR%LAT%NEP
          SIGABS = SIGABS + ATOM%SIGA
          SIGINC = SIGINC + ATOM%SIGI
          SIGB = SIGB + ATOM%BC**2
          SII = SII + ATOM%SIGI
        enddo
      ! sig_coh + sig_inc in [barn]
        SIGB=4.D-2*PI*SIGB+SII
        if (x.le.6) then
          SIGTDS(1)=SIGTDS(1)+3.D0*SIGB/ATOM%M*SQRT(k_b*CR%LAT%thetaD/E)*R
        else
          SIGTDS(1)=SIGTDS(1)+3.D0*SIGB/ATOM%M*SQRT(k_b*CR%LAT%thetaD/E)*3.3D0/SQRT(x**7)
        endif
        SIGTDS(2)=SIGTDS(2)+SIGB*(ATOM%M/(ATOM%M+1.008D0))**2
        SIGTDS(3)=SIGTDS(3)+SIGB*(ATOM%M/(ATOM%M+1.008D0))**2*exp(-(B0+BT)*CR%C2*E)
      enddo
! Definne  SIGTDS(2),SIGTDS(3) so that multiphonon cross-section is
! SIGTDS(2)*(1 - exp(-SIGTDS(3))
      if ((SIGTDS(2).gt.0).and.(SIGTDS(3).gt.0)) then
        ksi=-log(SIGTDS(3)/SIGTDS(2))
      else
        ksi=0.D0
      endif
      SIGTDS(1)=SIGTDS(1)/LAMBDA
      SIGTDS(3)=CR%XHAR*ksi*LAMBDA**2
      SIGTDS(4)=CR%XANH*(ksi*LAMBDA**2)**2
    ! per unit volume, convert to mm^-1 [barn/A^3 = 1/cm]
      SIGABS=SIGABS/CR%LAT%VOL/10.D0
      SIGINC=SIGINC/CR%LAT%VOL/10.D0
      do i=1,2
        SIGTDS(i)=SIGTDS(i)/CR%LAT%VOL/10.D0
      enddo
      end subroutine GET_CR_SIGMA

!*******************************
! XML handler for table input
!********************************

!--------------------------------------------------------
      subroutine startCRTABLE(tag,attribs,error)
!--------------------------------------------------------
        character(len=*)                  :: tag
        character(len=*), dimension(:,:)  :: attribs
        logical                           :: error
        character(32) :: S,ID,CNAME
    ! new table
        ! write(*,*) 'startCRTABLE ' //trim(tag)
        if (error) return
        select case(XML_LEVEL)
        case(-1)
          if (trim(tag).eq.'CRTABLE') then
            call CLEAR_TABLE_CRYSTALS
            XML_LEVEL=XML_LEVEL_TABLE
            return
          endif
        case(XML_LEVEL_TABLE)
          if (trim(tag).eq.'CRDATA') then
            if (NCR_TAB.lt.MAX_NCR_TAB) then
              call readAttrib(attribs,'id',ID,error)
              if (.not.error) call readAttrib(attribs,'name',CNAME,error)
              if (.not.error) then
                call CLEAR_CRITEM(myITEM)
                myITEM%ID=trim(adjustl(ID))
                myITEM%NAME=trim(adjustl(CNAME))
                XML_LEVEL=XML_LEVEL_DATA
              else
                call MSG_ERROR('CRTABLE','invalid attributes in '//trim(tag),0,1)
              endif
            else
              XML_LEVEL=-1
              call INT2STR(MAX_NCR_TAB,S)
              call MSG_WARN('Too many items in crystal table, max='//trim(S),1)
            endif
          endif
        case(XML_LEVEL_DATA)
          if (trim(tag).eq.'LATTICE') then
            XML_LEVEL=XML_LEVEL_LATTICE
            call startLATTICE(tag,attribs,error)
          endif
        case(XML_LEVEL_LATTICE)
          call startLATTICE(tag,attribs,error)
        end select
      end subroutine startCRTABLE

!--------------------------------------------------------
      subroutine dataCRTABLE(tag, xmldata, error )
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
        select case(XML_LEVEL)
        case(XML_LEVEL_DATA)
          select case(trim(tag))
          case('C2')
            read(xmldata,*) myITEM%C2
            myITEM%C2=myITEM%C2/1.D3 ! convert from ev^-1 to meV^-1
          case('XHAR')
            read(xmldata,*) myITEM%XHAR
          case('XANH')
            read(xmldata,*) myITEM%XANH
          case('ISO')
            call GETLINPARG(xmldata,ARG,3,NA)
            if (NA.ge.3) myITEM%ISO=ARG
          end select
        case(XML_LEVEL_LATTICE)
          call dataLATTICE(tag, xmldata, error)
        end select
      end subroutine dataCRTABLE

!--------------------------------------------------------
      subroutine endCRTABLE( tag, error )
!--------------------------------------------------------
      character(len=*)                  :: tag
      logical                           :: error
      if (error) return
      select case(XML_LEVEL)
      case(XML_LEVEL_DATA)
        select case(trim(tag))
        case('CRDATA')
        !  write(*,*) 'endCRTABLE ',myITEM%ID
          call ADD_CRITEM(myITEM)
          XML_LEVEL=XML_LEVEL_TABLE
        end select
      case(XML_LEVEL_LATTICE)
        select case(trim(tag))
        case('LATTICE')
          call endLATTICE(tag, error)
          myITEM%LAT=getLATTICE()
          XML_LEVEL=XML_LEVEL_DATA
        end select
      case(XML_LEVEL_TABLE)
        select case(trim(tag))
        case('CRTABLE')
          XML_LEVEL=-1
        end select
      end select
      end subroutine endCRTABLE

      end module TABLE_CRYSTALS
