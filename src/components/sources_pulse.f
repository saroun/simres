
!//////////////////////////////////////////////////////////////////////
!////  $Id: sources_pulse.f,v 1.13 2016/11/13 15:20:31 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2011, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.13 $
!////     $Date: 2016/11/13 15:20:31 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Calculation of spallation source pulse source
!////  Using F. Mezei, "ESS reference moderator characteristics  ..." report, Dec/2000
!////  updated for the ESS moderator (Dec 2013)
!////  types:
!////  (1) short coupled
!////  (2) short decoupled
!////  (3) long
!////  (4) ESS cold, models 2013-2014 (TDR and pancake)
!////  (5) ESS thermal, models 2013-2014 (TDR and pancake)
!////  (6) ESS cold, 2015 (butterfly)
!////  (7) ESS thermal, 2015 (butterfly)
!////////////////////////////////////////////////////////////////////////////////////
! The ESS source is modelled according to:
! "Analytical  fits  to  the  ESS  cold  and  thermal  moderator  spectra",  T.  Schönfeldt, 23/5/2013
! typ=3   T>273, thermal source, eq. 3
! typ=4   cold source, eq. 4
! typ=5   thermal source, eq. 3
! typ=5,6 as in McStas model, uses new time structure with a single
! relaxation time depending on wavelength, tau = A(1 + B*lambda^2) exp(-C(lambda+D)^-E)
      MODULE SOURCES_PULSE
      use SIMATH
      implicit none
      save
      private

      integer,parameter :: pulse_spc=1 ! short, coupled
      integer,parameter :: pulse_spd=2 ! short, decoupled
      integer,parameter :: pulse_lp=3  ! long, coupled
      integer,parameter :: pulse_ess_cold=4  ! ESS cold (para-H2), 2013 (Schonfeldt, eq. 4)
      integer,parameter :: pulse_ess_therm=5  ! ESS thermal, 2013 (Schonfeldt, eq. 3)
      integer,parameter :: pulse_ess_2015C=6  ! ESS 2015 (butterfly), cold
      integer,parameter :: pulse_ess_2015T=7  ! ESS 2015 (butterfly), thermal

      ! defined names
      character(25),parameter :: NAMDEF123='F1:TEMP:ISD:ALPHASD:LAMSD'
      character(54),parameter :: NAMDEF4='I1:ALPHA1:I2:ALPHA2:ALPHAC:LAMC:EXPO:ISD:ALPHASD:LAMSD'
      character(40),parameter :: NAMDEF5='F1:TEMP1:F2:TEMP2:ISD:ALPHASD:LAMSD:EXPO'
      character(54),parameter :: NAMDEF6='I1:ALPHA1:I2:ALPHA2:ALPHAC:LAMC:EXPO:ISD:ALPHASD:LAMSD'
      character(32),parameter :: NAMDEF7='F1:TEMP1:DELTA:ISD:ALPHASD:LAMSD'
      character(41),parameter :: NAMDEFTIME='TAU1:TAUA:N1:TAU2:N2:WIDTH:TMAX:LMIN:LMAX'
      character(30),parameter :: NAMDEFTIME2='TAU1:TAUP:WIDTH:TMAX:LMIN:LMAX'
      integer,parameter :: NPAR_TIME=9
      integer,parameter :: NPAR_TIME2=6

      TYPE SRCPULSE; SEQUENCE
  ! Maxwell spectrum
        real(kind(1.D0)) :: TEMP ! temperature [K]
        real(kind(1.D0)) :: F1 ! flux [1e10/s/cm2/ster]
        real(kind(1.D0)) :: TEMP2 ! temperature [K]
        real(kind(1.D0)) :: F2 ! flux [1e10/s/cm2/ster]
        real(kind(1.D0)) :: DELTA ! lambda exponent, difference from 5
  ! slowing down
        real(kind(1.D0)) :: ISD ! flux [1e10/s/cm2/ster]
        real(kind(1.D0)) :: ALPHASD ! cut-off parameter [A^-1]
        real(kind(1.D0)) :: LAMSD ! cut-off wavelength [A]
  ! para-H2, first term (type 4)
        real(kind(1.D0)) :: I1 ! flux 1 [1e10/s/cm2/ster]
        real(kind(1.D0)) :: I2 ! flux 2 [1e10/s/cm2/ster]
        real(kind(1.D0)) :: ALPHA1 ! exp. decay coefficient 1 [A^-1]
        real(kind(1.D0)) :: ALPHA2 ! exp. decay coefficient 2 [A^-1]
  ! para-H2, cut-off (type 4)
        real(kind(1.D0)) :: ALPHAC ! exp. decay coefficient [A^-1]
        real(kind(1.D0)) :: LAMC  ! cut-off wavelength [A]
        real(kind(1.D0)) :: EXPO ! exponent for the cut-off term
  ! pulse shape
        real(kind(1.D0)) :: TAU1 ! decay time 1 [us]
        real(kind(1.D0)) :: TAU2 ! decay time 2 [us]
        real(kind(1.D0)) :: D    ! proton pulse duration [us]
        real(kind(1.D0)) :: TAUA ! decay time for ambiente H2O coupled [us]
        real(kind(1.D0)) :: TAUP(4) ! parameters for decay time wavelength dependence
        integer :: N1 ! pulse shape parameter 1
        integer :: N2 ! pulse shape parameter 2
        integer :: typ ! source type
  ! table range
        real(kind(1.D0)) :: TMAX,LMIN,LMAX ! time [us] and wavelength [A} limits
      END TYPE SRCPULSE
      real(kind(1.D0)),parameter :: TAMB = 273.D0 ! temperature limit for using ambiente version

      integer,parameter :: SPSMAX=2
      TYPE(SRCPULSE) :: SPS(SPSMAX)

      public SRCPULSE,SPS,PULSE_SET_DEFAULT
      public PULSE_GET_PARAM,PULSE_GET_NAMES,PULSE_READ_PARAM
      public PULSE_FLUX
      public pulse_spc,pulse_spd,pulse_lp,pulse_ess_cold,pulse_ess_therm,pulse_ess_2015T,pulse_ess_2015C
      public LPSSA_DEBUG,PULSE_FLUX_DBG


      contains


!-------------------------------------------------------------------------
      subroutine PULSE_SET_DEFAULT(it)
! set default pulse data - 1 MW short pulse ambiente decoupled unposioned
!-------------------------------------------------------------------------
      integer, intent(in) :: it
      select case(SPS(it)%TYP)
      case(pulse_spc,pulse_spd,pulse_lp)
      ! Maxwell spectrum
        SPS(it)%TEMP=325.D0    ! temperature [K]
        SPS(it)%F1=160.7       ! flux 1 [10^10/s/cm2/ster]
      ! slowing down
        SPS(it)%ISD=32.84      ! flux 2 [10^10/s/cm2/ster]
        SPS(it)%ALPHASD=2.5D0  ! cut-off parameter [A^-1]
        SPS(it)%LAMSD=8.8D-1   ! cut-off wavelength [A]
      case(pulse_ess_cold)
      ! para-H2, first term
        SPS(it)%I1=361.3       ! flux, exp. decay 1 [n/s/cm2/ster]
        SPS(it)%I2=5.774       ! flux, exp. decay 2 [n/s/cm2/ster]
        SPS(it)%ALPHA1=0.6175  ! exp. decay coefficient 1 [A^-1]
        SPS(it)%ALPHA2=0.223   ! exp. decay coefficient 2 [A^-1]
      ! para-H2, cut-off
        SPS(it)%ALPHAC=-6.714  ! exp. decay coefficient [A^-1]
        SPS(it)%LAMC=2.421     ! cut-off wavelength [A]
        SPS(it)%EXPO=0.5D0     ! exponent for the cut-off term
      ! slowing down
        SPS(it)%ISD=32.84      ! flux 2 [10^10/s/cm2/ster]
        SPS(it)%ALPHASD=9.0D-1 ! cut-off parameter [A^-1]
        SPS(it)%LAMSD=2.44D0   ! cut-off wavelength [A]
      case(pulse_ess_2015C)
      ! para-H2, first term
        SPS(it)%I1=1680.8114       ! flux, exp. decay 1 [n/s/cm2/ster]
        SPS(it)%I2=74.864       ! flux, exp. decay 2 [n/s/cm2/ster]
        SPS(it)%ALPHA1=0.7975  ! exp. decay coefficient 1 [A^-1]
        SPS(it)%ALPHA2=0.325   ! exp. decay coefficient 2 [A^-1]
      ! para-H2, cut-off
        SPS(it)%ALPHAC=-13.5  ! exp. decay coefficient [A^-1]
        SPS(it)%LAMC=2.53     ! cut-off wavelength [A]
        SPS(it)%EXPO=0.1496     ! exponent for the cut-off term
      ! slowing down
        SPS(it)%ISD=24.1143      ! flux 2 [10^10/s/cm2/ster]
        SPS(it)%ALPHASD=2.5 ! cut-off parameter [A^-1]
        SPS(it)%LAMSD=2.2   ! cut-off wavelength [A]
      case(pulse_ess_therm)
      ! Maxwell spectrum
        SPS(it)%TEMP=325.D0  ! temperature [K]
        SPS(it)%F1=124.667   ! flux 1 [10^10/s/cm2/ster]
        SPS(it)%TEMP2=20.D0  ! 2nd Maxw. component, temperature [K]
        SPS(it)%F2=0.0       ! 2nd Maxw. component, flux [10^10/s/cm2/ster]
      ! slowing down
        SPS(it)%ISD=21.393    ! flux 2 [10^10/s/cm2/ster]
        SPS(it)%ALPHASD=2.5D0 ! cut-off parameter [A^-1]
        SPS(it)%LAMSD=8.8D-1  ! cut-off wavelength [A]
        SPS(it)%EXPO=1.0D0    ! exponent for the slow-down term
      case(pulse_ess_2015T)
      ! Maxwell spectrum
        SPS(it)%F1=304.6536   ! flux 1 [10^10/s/cm2/ster]
        SPS(it)%TEMP=325  ! temperature [K]
        SPS(it)%DELTA=-0.289
      ! slowing down
        SPS(it)%ISD=43.9971    ! flux 2 [10^10/s/cm2/ster]
        SPS(it)%ALPHASD=2.5D0 ! cut-off parameter [A^-1]
        SPS(it)%LAMSD=8.8D-1  ! cut-off wavelength [A]
      end select
  ! pulse shape
      SPS(it)%TAU1=35.D0  ! decay time 1 [us]
      SPS(it)%TAU2=12.D0  ! decay time 2 [us]
      SPS(it)%D=2.86      ! proton pulse duration [us]
      SPS(it)%TAUA=400.D0 ! decay time for ambiente H2O coupled [us]
      SPS(it)%N1=20       ! pulse shape parameter 1
      SPS(it)%N2=5        ! pulse shape parameter 2
   ! coefficients for calculation of wavelength - dependent relaxation time
      SPS(it)%TAUP=(/1.195D-2,39.109D0,0.98799D0,7.65675D0/)
      if (SPS(it)%TYP==pulse_ess_2015T) then
        SPS(it)%TAU1=308.81
      else if (SPS(it)%TYP==pulse_ess_2015C) then
        SPS(it)%TAU1=285.372
        SPS(it)%TAUP=(/4.37125D-3,66.3537D0,0.9D0,8.64455D0/)
      endif
  ! table limits
      SPS(it)%TMAX=10.D0
      SPS(it)%LMIN=0.3D0
      SPS(it)%LMAX=20.D0
      end subroutine PULSE_SET_DEFAULT

!---------------------------------------------------------------
      subroutine PULSE_GET_NAMES(it,NAMESPC)
! return parameter name space for given type
!---------------------------------------------------------------
      integer, intent(in) :: it
      character(*),intent(out) :: NAMESPC
      select case(SPS(it)%TYP)
      case(pulse_spc,pulse_spd,pulse_lp)
        NAMESPC=NAMDEFTIME//':'//NAMDEF123
      case(pulse_ess_cold)
        NAMESPC=NAMDEFTIME//':'//NAMDEF4
      case(pulse_ess_therm)
        NAMESPC=NAMDEFTIME//':'//NAMDEF5
      case(pulse_ess_2015C)
        NAMESPC=NAMDEFTIME2//':'//NAMDEF6
      case(pulse_ess_2015T)
        NAMESPC=NAMDEFTIME2//':'//NAMDEF7
      case default
        NAMESPC=''
      end select
      end subroutine PULSE_GET_NAMES

!---------------------------------------------------------------
      subroutine PULSE_READ_PARAM(it,LINE,ID,IDX,IERR)
! Read source parameter value from a text line, assumed format is ID=value
!---------------------------------------------------------------
      integer, intent(in) :: it
      character(*),intent(in) :: LINE
      character(*),intent(out) :: ID
      integer, intent(out) :: IDX,IERR
      integer :: IE
1     format(a,':',6(1x,G12.5))
      select case(SPS(it)%TYP)
      case(pulse_spc,pulse_spd,pulse_lp)
        call READ_PARAM123(it,LINE,ID,IDX,IE)
      case(pulse_ess_cold,pulse_ess_2015C)
        ! call READ_PARAM4(LINE,ID,IDX,IE)
        call READ_PARAM_COLD(it,LINE,SPS(it)%TYP,ID,IDX,IE)
      case(pulse_ess_therm,pulse_ess_2015T)
        ! call READ_PARAM5(LINE,ID,IDX,IE)
        call READ_PARAM_THERM(it,LINE,SPS(it)%TYP,ID,IDX,IE)
        !write(*,1) 'READ_PARAM_THERM '//trim(ID)//' IDX',idx
      case default
        IE=-1
      end select
      if ((IE==0).AND.(IDX>0)) then
        select case(SPS(it)%TYP)
        case(pulse_ess_2015C,pulse_ess_2015T)
          IDX=IDX+NPAR_TIME2
        case default
          IDX=IDX+NPAR_TIME
        end select
      else if (IE==1) then
        call READ_PARAM_TIME(it,LINE,SPS(it)%TYP,ID,IDX,IE)
      endif
      !write(*,1) 'PULSE_READ_PARAM '//trim(ID)//' IDX',idx
      IERR=IE
      end subroutine PULSE_READ_PARAM


!---------------------------------------------------------------
      subroutine VALID_ID(it,S,ID)
! Convert given string to valid ID:
! convert to upper case
! change obsolete ID's the current namespace
!---------------------------------------------------------------
      integer, intent(in) :: it
      character(*),intent(in) :: S
      character(*),intent(out) :: ID

      character(16) :: STR
      STR=trim(S)
      call MKUPCASE(STR)
      ID=trim(STR)
      select case(SPS(it)%TYP)
      case(pulse_spc,pulse_spd,pulse_lp)
        select case(trim(STR))
          case('F2')
            ID='ISD'
          case('K2')
            ID='ALPHASD'
        end select
      case(pulse_ess_therm)
          select case(trim(STR))
          case('TEMP')
            ID='TEMP1'
          end select
      end select
      end subroutine VALID_ID

!---------------------------------------------------------------
      subroutine PULSE_GET_PARAM(it,ID,PAR,NP,IDX)
! Return parameter values for given ID
! IDX=0 if not found
!---------------------------------------------------------------
      integer, intent(in) :: it
      character(*),intent(in) :: ID
      real(kind(1.D0)),intent(out) :: PAR(:)
      integer, intent(out) :: NP,IDX
      character(256) :: NAM
      integer :: i
      real(kind(1.D0)) :: Z(64)
      character(16) :: S

      IDX=0
      Z=0.D0
      NP=0
      call PULSE_GET_NAMES(it,NAM)
      call VALID_ID(it,ID,S)
      call GETPARINDX(NAM,':',trim(S),i)
      if (i>0) then
        NP=1
        IDX=i
        select case(trim(S))
! Maxwell spectrum
            case('TEMP')
              Z(1)=SPS(it)%TEMP
            case('F1')
              Z(1)=SPS(it)%F1
            case('TEMP2')
              Z(1)=SPS(it)%TEMP2
            case('F2')
              Z(1)=SPS(it)%F2
! slowing down
            case('ISD')
              Z(1)=SPS(it)%ISD
            case('ALPHASD')
              Z(1)=SPS(it)%ALPHASD
            case('LAMSD')
              Z(1)=SPS(it)%LAMSD
 ! para-H2, first term
            case('I1')
              Z(1)=SPS(it)%I1
            case('I2')
              Z(1)=SPS(it)%I2
            case('ALPHA1')
              Z(1)=SPS(it)%ALPHA1
            case('ALPHA2')
              Z(1)=SPS(it)%ALPHA2
! para-H2, cut-off
            case('ALPHAC')
              Z(1)=SPS(it)%ALPHAC
            case('LAMC')
              Z(1)=SPS(it)%LAMC
            case('EXPO')
              Z(1)=SPS(it)%EXPO
! pulse shape
            case('TAU1')
              Z(1)=SPS(it)%TAU1
            case('TAUA')
              Z(1)=SPS(it)%TAUA
            case('N1')
              Z(1)=SPS(it)%N1
            case('TAU2')
              Z(1)=SPS(it)%TAU2
            case('N2')
              Z(1)=SPS(it)%N2
            case('WIDTH')
              Z(1)=SPS(it)%D
            case('TAUP')
              Z(1:4)=SPS(it)%TAUP(1:4)
              NP=4
! time range
            case('TMAX')
              Z(1)=SPS(it)%TMAX
! wavelength range
            case('LMIN')
              Z(1)=SPS(it)%LMIN
            case('LMAX')
              Z(1)=SPS(it)%LMAX
            case default
              IDX=0
        end select
        if ((IDX>0).AND.(NP>0).AND.(SIZE(PAR)>=NP)) PAR(1:NP)=Z(1:NP)
      endif

      end subroutine PULSE_GET_PARAM

!---------------------------------------------------------------
      subroutine READ_PARAM123(it,LINE,ID,IDX,IERR)
! For source type = 1..3
! Read source parameter value from a text line, assumed format is ID=value
! ID      parameter name
! LINE    input text
! return: ID and index of the parameter from NAMDEF123
! IERR=-1  ... wrong syntax
! IERR=-2  ... wrong number format
! IERR=1   ... unknown ID
!---------------------------------------------------------------
      integer, intent(in) :: it
      character(*),intent(in) :: LINE
      character(*),intent(out) :: ID
      integer, intent(out) :: IDX,IERR
      character*(16) :: S,S2
      real(kind(1.D0)) :: Z(1)
      integer :: IE,NA
1     format(a,'=',6(1x,G12.5))
      ID=''
      IDX=0
      call READ_ID_R8(LINE,S,Z,1,NA,IE)
      if (IE==0) then
        call VALID_ID(it,S,S2)
        select case(trim(S2))
    ! Maxwell spectrum
            case('TEMP')
              SPS(it)%TEMP=Z(1)
            case('F1')
              SPS(it)%F1=Z(1)
    ! slowing down
            case('ISD')
              SPS(it)%ISD=Z(1)
            case('ALPHASD')
              SPS(it)%ALPHASD=Z(1)
            case('LAMSD')
              SPS(it)%LAMSD=Z(1)
            case default
              IE=1 ! unkown ID
        end select
        if (IE==0) then
          ID=trim(S2)
          call GETPARINDX(NAMDEF123,':',ID,IDX)
          !write(*,1) 'READ_PARAM123: '//trim(ID),Z,IDX
        endif
      else if (IE==-1) then
        ID=trim(S)
      endif
      IERR=IE
      end subroutine READ_PARAM123


!---------------------------------------------------------------
      subroutine READ_PARAM4(it,LINE,ID,IDX,IERR)
! For source type = 4
! Read source parameter value from a text line, assumed format is ID=value
! ID      parameter name
! LINE    input text
! return: ID and index of the parameter from NAMDEF123
! IERR=-1  ... wring syntax
! IERR=-2  ... wrong number format
! IERR=1   ... unknown ID
!---------------------------------------------------------------
      integer, intent(in) :: it
      character(*),intent(in) :: LINE
      character(*),intent(out) :: ID
      integer, intent(out) :: IDX,IERR
      character*(16) :: S,S2
      real(kind(1.D0)) :: Z(1)
      integer :: IE,NA
1     format(a,'=',6(1x,G12.5))
      ID=''
      IDX=0
      call READ_ID_R8(LINE,S,Z,1,NA,IE)
      if (IE==0) then
        call VALID_ID(it,S,S2)
        select case(trim(S2))
  ! slowing down
              case('ISD')
                SPS(it)%ISD=Z(1)
              case('ALPHASD')
                SPS(it)%ALPHASD=Z(1)
              case('LAMSD')
                SPS(it)%LAMSD=Z(1)
  ! para-H2, first term
              case('I1')
                SPS(it)%I1=Z(1)
              case('I2')
                SPS(it)%I2=Z(1)
              case('ALPHA1')
                SPS(it)%ALPHA1=Z(1)
              case('ALPHA2')
                SPS(it)%ALPHA2=Z(1)
  ! para-H2, cut-off
              case('ALPHAC')
                SPS(it)%ALPHAC=Z(1)
              case('LAMC')
                SPS(it)%LAMC=Z(1)
              case('EXPO')
                SPS(it)%EXPO=Z(1)
            case default
              IE=1 ! unkown ID
        end select
        if (IE==0) then
          ID=trim(S2)
          call GETPARINDX(NAMDEF4,':',ID,IDX)
          !write(*,1) 'READ_PARAM4: '//trim(ID),Z,IDX
        endif
      else if (IE==-1) then
        ID=trim(S)
      endif
      IERR=IE
      end subroutine READ_PARAM4


!---------------------------------------------------------------
      subroutine READ_PARAM5(it,LINE,ID,IDX,IERR)
! For source type = 5
! Read source parameter value from a text line, assumed format is ID=value
! ID      parameter name
! LINE    input text
! return: ID and index of the parameter from NAMDEF123
! IERR=-1  ... wring syntax
! IERR=-2  ... wrong number format
! IERR=1   ... unknown ID
!---------------------------------------------------------------
      integer, intent(in) :: it
      character(*),intent(in) :: LINE
      character(*),intent(out) :: ID
      integer, intent(out) :: IDX,IERR
      character*(16) :: S,S2
      real(kind(1.D0)) :: Z(1)
      integer :: IE,NA
1     format(a,'=',6(1x,G12.5))
      ID=''
      IDX=0
      call READ_ID_R8(LINE,S,Z,1,NA,IE)
      if (IE==0) then
        call VALID_ID(it,S,S2)
        select case(trim(S2))
  ! Maxwell spectrum
              case('TEMP1')
                SPS(it)%TEMP=Z(1)
              case('F1')
                SPS(it)%F1=Z(1)
              case('TEMP2')
                SPS(it)%TEMP2=Z(1)
              case('F2')
                SPS(it)%F2=Z(1)
  ! slowing down
              case('ISD')
                SPS(it)%ISD=Z(1)
              case('ALPHASD')
                SPS(it)%ALPHASD=Z(1)
              case('LAMSD')
                SPS(it)%LAMSD=Z(1)
              case('EXPO')
                SPS(it)%EXPO=Z(1)
              case default
                IE=1 ! unkown ID
        end select
        if (IE==0) then
          ID=trim(S2)
          call GETPARINDX(NAMDEF5,':',ID,IDX)
          !write(*,1) 'READ_PARAM5: '//trim(ID),Z,IDX
        endif
      else if (IE==-1) then
        ID=trim(S)
      endif
      IERR=IE
      end subroutine READ_PARAM5

!---------------------------------------------------------------
      subroutine READ_PARAM_COLD(it,LINE,TYP,ID,IDX,IERR)
! For cold source type = 4,6
! Read source parameter value from a text line, assumed format is ID=value
! TYP     source type - reference to the name space (4 or 6)
! ID      parameter name
! LINE    input text
! return: ID and index of the parameter from NAMDEF123
! IERR=-1  ... wrong syntax
! IERR=-2  ... wrong number format
! IERR=1   ... unknown ID
!---------------------------------------------------------------
      integer, intent(in) :: it
      character(*),intent(in) :: LINE
      integer,intent(in) :: TYP
      character(*),intent(out) :: ID
      integer, intent(out) :: IDX,IERR
      character*(16) :: S,S2
      real(kind(1.D0)) :: Z(5)
      integer :: IE,NA
1     format(a,'=',6(1x,G12.5))
      ID=''
      IDX=0
      call READ_ID_R8(LINE,S,Z,5,NA,IE)
      if (IE==0) then
        call VALID_ID(it,S,S2)
        select case(trim(S2))
  ! slowing down
              case('ISD')
                SPS(it)%ISD=Z(1)
              case('ALPHASD')
                SPS(it)%ALPHASD=Z(1)
              case('LAMSD')
                SPS(it)%LAMSD=Z(1)
  ! para-H2, first term
              case('I1')
                SPS(it)%I1=Z(1)
              case('I2')
                SPS(it)%I2=Z(1)
              case('ALPHA1')
                SPS(it)%ALPHA1=Z(1)
              case('ALPHA2')
                SPS(it)%ALPHA2=Z(1)
  ! para-H2, cut-off
              case('ALPHAC')
                SPS(it)%ALPHAC=Z(1)
              case('LAMC')
                SPS(it)%LAMC=Z(1)
              case('EXPO')
                SPS(it)%EXPO=Z(1)
            case default
              IE=1 ! unkown ID
        end select
        if (IE==0) then
          ID=trim(S2)
          select case (TYP)
          case(pulse_ess_cold)
            call GETPARINDX(NAMDEF4,':',ID,IDX)
          case(pulse_ess_2015C)
            call GETPARINDX(NAMDEF6,':',ID,IDX)
          end select
        endif
      else if (IE==-1) then
        ID=trim(S)
      endif
      IERR=IE
      end subroutine READ_PARAM_COLD


!---------------------------------------------------------------
      subroutine READ_PARAM_THERM(it,LINE,TYP,ID,IDX,IERR)
! For thermal source type = 5,7
! Read source parameter value from a text line, assumed format is ID=value
! TYP     source type - reference to the name space (5 or 7)
! ID      parameter name
! LINE    input text
! return: ID and index of the parameter from NAMDEF123
! IERR=-1  ... wrong syntax
! IERR=-2  ... wrong number format
! IERR=1   ... unknown ID
!---------------------------------------------------------------
      integer, intent(in) :: it
      character(*),intent(in) :: LINE
      integer,intent(in) :: TYP
      character(*),intent(out) :: ID
      integer, intent(out) :: IDX,IERR
      character*(16) :: S,S2
      real(kind(1.D0)) :: Z(5)
      integer :: IE,NA
1     format(a,'=',6(1x,G12.5))
      ID=''
      IDX=0
      call READ_ID_R8(LINE,S,Z,5,NA,IE)
      if (IE==0) then
        call VALID_ID(it,S,S2)
        select case(trim(S2))
  ! Maxwell spectrum
              case('TEMP1')
                SPS(it)%TEMP=Z(1)
              case('F1')
                SPS(it)%F1=Z(1)
              case('TEMP2')
                SPS(it)%TEMP2=Z(1)
              case('F2')
                SPS(it)%F2=Z(1)
              case('DELTA')
                SPS(it)%DELTA=Z(1)
  ! slowing down
              case('ISD')
                SPS(it)%ISD=Z(1)
              case('ALPHASD')
                SPS(it)%ALPHASD=Z(1)
              case('LAMSD')
                SPS(it)%LAMSD=Z(1)
              case('EXPO')
                SPS(it)%EXPO=Z(1)
              case default
                IE=1 ! unkown ID
        end select
        if (IE==0) then
          ID=trim(S2)
          select case (TYP)
          case(pulse_ess_therm)
            call GETPARINDX(NAMDEF5,':',ID,IDX)
          case(pulse_ess_2015T)
            call GETPARINDX(NAMDEF7,':',ID,IDX)
            !write(*,1) 'READ_PARAM_THERM '//trim(ID)//' IDX=',IDX
          end select
        endif
      else if (IE==-1) then
        ID=trim(S)
      endif
      IERR=IE
      end subroutine READ_PARAM_THERM


!---------------------------------------------------------------
      subroutine READ_PARAM_TIME(it,LINE,TYP,ID,IDX,IERR)
! For all source types, time parameters and table range
! Read source parameter value from a text line, assumed format is ID=value
! ID      parameter name
! LINE    input text
! return: ID and index of the parameter from NAMDEF123
! IERR=-1  ... wring syntax
! IERR=-2  ... wrong number format
! IERR=1   ... unknown ID
!---------------------------------------------------------------
      integer, intent(in) :: it
      character(*),intent(in) :: LINE
      integer,intent(in) :: TYP
      character(*),intent(out) :: ID
      integer, intent(out) :: IDX,IERR
      character*(16) :: S,S2
      real(kind(1.D0)) :: Z(5)
      integer :: IE,NA
1     format(a,'=',6(1x,G12.5))
      ID=''
      IDX=0
      ! read scalars
      call READ_ID_R8(LINE,S,Z,5,NA,IE)
      if (IE==0) then
        call VALID_ID(it,S,S2)
        if (NA==1) then
          select case(trim(S2))
  ! pulse shape
              case('TAU1')
                SPS(it)%TAU1=Z(1)
              case('TAUA')
                SPS(it)%TAUA=Z(1)
              case('N1')
                SPS(it)%N1=NINT(Z(1))
              case('TAU2')
                SPS(it)%TAU2=Z(1)
              case('N2')
                SPS(it)%N2=NINT(Z(1))
              case('WIDTH')
                SPS(it)%D=Z(1)
  ! time range
              case('TMAX')
                SPS(it)%TMAX=Z(1)
  ! wavelength range
              case('LMIN')
                SPS(it)%LMIN=Z(1)
              case('LMAX')
                SPS(it)%LMAX=Z(1)
              case default
                IE=1 ! unkown ID
          end select
        else if (NA>0) then
      ! try array parameters
          select case(trim(S2))
            case('TAUP')
            if (NA==4) then
              SPS(it)%TAUP(1:4)=Z(1:4)
            else
              IE=-2
            endif
          case default
            IE=1 ! unkown ID
          end select
        endif
      endif
      if (IE==0) then
        ID=trim(S2)
        select case (TYP)
          case(pulse_ess_2015T,pulse_ess_2015C)
            call GETPARINDX(NAMDEFTIME2,':',ID,IDX)
          case default
            call GETPARINDX(NAMDEFTIME,':',ID,IDX)
        end select
      else if (IE==-1) then
        ID=trim(S)
      endif
      IERR=IE
      end subroutine READ_PARAM_TIME


!-------------------------------------------------------------------------
      real(kind(1.D0)) function PULSE_FLUX_DBG(it,t,lam,dbg)
! Instant flux for spallation source in [10^10/cm^2/us/A]
!-------------------------------------------------------------------------
      integer, intent(in) :: it
      real(kind(1.D0)),intent(in) :: t,lam
      real(kind(1.D0)) :: RES
      logical,intent(in) :: dbg
      RES=0.D0
      select case(SPS(it)%typ)
      case(pulse_spc)
        if (SPS(it)%TEMP.gt.TAMB) then
          RES=SPSSA_FLUX(it,t,lam)
        else
          RES=SPSS_FLUX(it,t,lam)
        endif
      case(pulse_spd)
        RES=SPSS_FLUX(it,t,lam)
      case(pulse_lp)
        if (SPS(it)%TEMP.gt.TAMB) then
          RES=LPSSA_FLUX(it,t,lam)
        else
          RES=LPSS_FLUX(it,t,lam)
        endif
      case(pulse_ess_cold)
        RES=LPSS_ESS_COLD(it,t,lam)
      case(pulse_ess_therm)
        RES=LPSS_ESS_THERM(it,t,lam)
      case(pulse_ess_2015C)
        RES=ESS_2015_COLD(it,t,lam)
      case(pulse_ess_2015T)
        ! write(*,*) 'PULSE_FLUX ',t,lam
        ! RES=0
        RES=ESS_2015_THERM_DBG(it,t,lam,dbg)
      end select
! convert to proper units: multiply by 4PI
! note: RES is in 10^10/cm^2/us/ster]
      PULSE_FLUX_DBG=2.D0*TWOPI*RES
      end function PULSE_FLUX_DBG

!-------------------------------------------------------------------------
      real(kind(1.D0)) function PULSE_FLUX(it,t,lam)
! Instant flux for spallation source in [10^10/cm^2/us/A]
!-------------------------------------------------------------------------
      integer, intent(in) :: it
      real(kind(1.D0)),intent(in) :: t,lam
      real(kind(1.D0)) :: RES
      RES=0.D0
      select case(SPS(it)%typ)
      case(pulse_spc)
        if (SPS(it)%TEMP.gt.TAMB) then
          RES=SPSSA_FLUX(it,t,lam)
        else
          RES=SPSS_FLUX(it,t,lam)
        endif
      case(pulse_spd)
        RES=SPSS_FLUX(it,t,lam)
      case(pulse_lp)
        if (SPS(it)%TEMP.gt.TAMB) then
          RES=LPSSA_FLUX(it,t,lam)
        else
          RES=LPSS_FLUX(it,t,lam)
        endif
      case(pulse_ess_cold)
        RES=LPSS_ESS_COLD(it,t,lam)
      case(pulse_ess_therm)
        RES=LPSS_ESS_THERM(it,t,lam)
      case(pulse_ess_2015C)
        RES=ESS_2015_COLD(it,t,lam)
      case(pulse_ess_2015T)
        ! write(*,*) 'PULSE_FLUX ',t,lam
        ! RES=0
        RES=ESS_2015_THERM(it,t,lam)
      end select
! convert to proper units: multiply by 4PI
! note: RES is in 10^10/cm^2/us/ster]
      PULSE_FLUX=2.D0*TWOPI*RES
      end function PULSE_FLUX

!-------------------------------------------------------------------------
      real(kind(1.D0)) function SPSS_F(t,tau,n)
! pulse shape function for short pulse
!-------------------------------------------------------------------------
      real(kind(1.D0)),intent(in) :: t,tau
      integer,intent(in) :: n
      real(kind(1.D0)) :: Z,RES
      RES=0.D0
      if ((tau>0.D0).and.(n>1)) then
        Z=t/tau
        RES=(exp(-Z)-exp(-n*Z))*n/((n-1)*tau)
      endif
      SPSS_F=RES
      end function SPSS_F

!-------------------------------------------------------------------------
      real(kind(1.D0)) function IEXP(t,tau,d)
! pulse decay function for long pulse
!-------------------------------------------------------------------------
      real(kind(1.D0)),intent(in) :: t,tau,d
      real(kind(1.D0)) :: RES
      if (t<=0.D0) then
        RES=0.D0
      else if (t<d) then
        RES=tau*(1.D0-exp(-t/tau))
      else
        RES=tau*(exp((d-t)/tau) - exp(-t/tau))
      endif
      IEXP=RES
      end function IEXP

!-------------------------------------------------------------------------
      real(kind(1.D0)) function LPSS_F(t,tau,n,d)
! pulse shape function for long pulse
!-------------------------------------------------------------------------
      real(kind(1.D0)),intent(in) :: t,tau,d
      integer,intent(in) :: n
      real(kind(1.D0)) :: RES
      RES=0.D0
      if ((tau>0.D0).and.(n>1)) then
        RES=(IEXP(t,tau,d)-IEXP(t,tau/n,d))*n/((n-1)*tau*d)
      endif
      LPSS_F=RES
      end function LPSS_F

!-------------------------------------------------------------------------
      real(kind(1.D0)) function FCO(lam,alpha_c,lam_c)
! Cut off function
!-------------------------------------------------------------------------
      real(kind(1.D0)),intent(in) :: lam,alpha_c,lam_c
      real(kind(1.D0)) :: RES
      RES=1.D0+exp(alpha_c*(lam-lam_c))
      FCO=1.D0/RES
      end function FCO


!-------------------------------------------------------------------------
      real(kind(1.D0)) function MAXW_LAM(lam,T)
! Normalized Maxwell distribution for wavelengths
!-------------------------------------------------------------------------
      real(kind(1.D0)),intent(in) :: lam,T
      real(kind(1.D0)) :: Z,RES,LT2
      RES=0.D0
      if (lam*T>0.D0) then
        LT2=HSQOV2M*TWOPI**2/(k_B*T)
        Z=LT2/lam**2
        RES=2.D0*(Z**2)/lam*exp(-Z)
      endif
      MAXW_LAM=RES
      end function MAXW_LAM

!-------------------------------------------------------------------------
      real(kind(1.D0)) function SPSS_FLUX(it,t,lam)
! Short pulse flux at given time and wavelength
! Single Maxwellian with single time profile + slow down term
!-------------------------------------------------------------------------
      integer, intent(in) :: it
      real(kind(1.D0)),intent(in) :: lam,T
      real(kind(1.D0)) :: RES
      RES=0.D0
      if (t>0.D0) then
        RES=SPS(it)%F1*MAXW_LAM(lam,SPS(it)%TEMP)*SPSS_F(t,SPS(it)%TAU1,SPS(it)%N1)
        RES=RES+SPS(it)%ISD/lam*FCO(lam,SPS(it)%ALPHASD,SPS(it)%LAMSD)*SPSS_F(t,SPS(it)%TAU2*lam,SPS(it)%N2)
      endif
      SPSS_FLUX=RES
      end function SPSS_FLUX


!-------------------------------------------------------------------------
      real(kind(1.D0)) function SPSSA_FLUX(it,t,lam)
! Short pulse flux at given time and wavelength
! 2nd version with 2 shape components for ambiente coupled moderator
! Single Maxwellian with double time profile + slow down term
!-------------------------------------------------------------------------
      integer, intent(in) :: it
      real(kind(1.D0)),intent(in) :: lam,T
      real(kind(1.D0)) :: RES
      RES=0.D0
      if (t>0.D0) then
        RES=SPS(it)%F1*MAXW_LAM(lam,SPS(it)%TEMP)*(SPSS_F(t,SPS(it)%TAU1,SPS(it)%N1)+SPSS_F(t,SPS(it)%TAUA,SPS(it)%N1))
        RES=RES+SPS(it)%ISD/lam*FCO(lam,SPS(it)%ALPHASD,SPS(it)%LAMSD)*SPSS_F(t,SPS(it)%TAU2*lam,SPS(it)%N2)
      endif
      SPSSA_FLUX=RES
      end function SPSSA_FLUX

!-------------------------------------------------------------------------
      real(kind(1.D0)) function LPSS_FLUX(it,t,lam)
! Long pulse flux at given time and wavelength
! Single Maxwellian with single time profile + slow down term
!-------------------------------------------------------------------------
      integer, intent(in) :: it
      real(kind(1.D0)),intent(in) :: lam,T
      real(kind(1.D0)) :: RES
      RES=0.D0
      if (t>0.D0) then
        RES=SPS(it)%F1*MAXW_LAM(lam,SPS(it)%TEMP)*LPSS_F(t,SPS(it)%TAU1,SPS(it)%N1,SPS(it)%d)
        RES=RES+SPS(it)%ISD/lam*FCO(lam,SPS(it)%ALPHASD,SPS(it)%LAMSD)*LPSS_F(t,SPS(it)%TAU2*lam,SPS(it)%N2,SPS(it)%d)
      endif
      LPSS_FLUX=RES
      end function LPSS_FLUX

!-------------------------------------------------------------------------
      real(kind(1.D0)) function LPSSA_FLUX(it,t,lam)
! Long pulse flux at given time and wavelength
! 2nd version with 2 shape components for ambiente coupled moderator
! Single Maxwellian with double time profile + slow down term
!-------------------------------------------------------------------------
      integer, intent(in) :: it
      real(kind(1.D0)),intent(in) :: lam,T
      real(kind(1.D0)) :: RES,J1,J2
      RES=0.D0
      if (t>0.D0) then
        J1=LPSS_F(t,SPS(it)%TAU1,SPS(it)%N1,SPS(it)%d)
        J2=LPSS_F(t,SPS(it)%TAUA,SPS(it)%N1,SPS(it)%d)
        RES=SPS(it)%F1*MAXW_LAM(lam,SPS(it)%TEMP)*(J1+J2)
        RES=RES+SPS(it)%ISD/lam*FCO(lam,SPS(it)%ALPHASD,SPS(it)%LAMSD)*LPSS_F(t,SPS(it)%TAU2*lam,SPS(it)%N2,SPS(it)%d)
      endif
      LPSSA_FLUX=RES
      end function LPSSA_FLUX

!-------------------------------------------------------------------------
      real(kind(1.D0)) function LPSS_ESS_THERM(it,t,lam)
! Long pulse, ESS thermal moderator,  2013 (Schonfeldt, eq. 3)
! Single Maxwellian with double time profile + slow down term
!-------------------------------------------------------------------------
      integer, intent(in) :: it
      real(kind(1.D0)),intent(in) :: lam,T
      real(kind(1.D0)) :: RES,I1,I2,SD,Z
2     format(a,' ',6(1x,G12.5))
      RES=0.D0
      if ((t>0.D0).and.(t<SPS(it)%TMAX)) then
      ! pulse shape 1
        I1=0.5*(LPSS_F(t,SPS(it)%TAU1,SPS(it)%N1,SPS(it)%d)+LPSS_F(t,SPS(it)%TAUA,SPS(it)%N1,SPS(it)%d))
      ! pulse shape 2
        I2=LPSS_F(t,SPS(it)%TAU2*lam,SPS(it)%N2,SPS(it)%d)
      ! Maxwell term
        RES=I1*SPS(it)%F1*MAXW_LAM(lam,SPS(it)%TEMP)
      ! add 2nd term if required
        if (SPS(it)%F2>0.D0) RES=RES+I2*SPS(it)%F2*MAXW_LAM(lam,SPS(it)%TEMP2)
      ! slowing down term
        Z=FCO(lam,SPS(it)%ALPHASD,SPS(it)%LAMSD)
        if (SPS(it)%EXPO.ne.1.D0) then
          SD=exp(SPS(it)%EXPO*log(Z))
        else
          SD=Z
        endif
        RES=RES+I2*SPS(it)%ISD/lam*SD
      endif
      LPSS_ESS_THERM=RES
      end function LPSS_ESS_THERM


!-------------------------------------------------------------------------
      real(kind(1.D0)) function ESS_2015_THERM(it,t,lam)
! Long pulse, ESS thermal moderator,  2015 (butterfly)
!-------------------------------------------------------------------------
      integer, intent(in) :: it
      real(kind(1.D0)),intent(in) :: lam,T
      real(kind(1.D0)) :: RES,I1,SD,MXW,Z,EZ,tau
2     format(a,' ',6(1x,G12.5))
      RES=0.D0
      if ((t>0.D0).and.(t<SPS(it)%TMAX)) then
      ! pulse shape
        Z=lam+SPS(it)%TAUP(3)
        EZ=exp(-SPS(it)%TAUP(4)*LOG(Z))
        tau=SPS(it)%TAU1*(1.D0+SPS(it)%TAUP(1)*lam**2)*exp(-SPS(it)%TAUP(2)*EZ)
        I1=IEXP(t,tau,SPS(it)%d)/(tau*SPS(it)%d)
      ! Maxwell term
        MXW=SPS(it)%F1*MAXW_LAM(lam,SPS(it)%TEMP)
        if (SPS(it)%DELTA.ne.0) then
          EZ=exp(-SPS(it)%DELTA*LOG(lam))
          MXW=MXW*EZ
        endif
      ! slowing down term
        SD=SPS(it)%ISD/lam*FCO(lam,SPS(it)%ALPHASD,SPS(it)%LAMSD)
        RES=(MXW+SD)*I1
      endif
      ESS_2015_THERM=RES
      end function ESS_2015_THERM

!-------------------------------------------------------------------------
      real(kind(1.D0)) function ESS_2015_THERM_DBG(it,t,lam,dbg)
! Long pulse, ESS thermal moderator,  2015 (butterfly)
!-------------------------------------------------------------------------
      integer, intent(in) :: it
      real(kind(1.D0)),intent(in) :: lam,T
      logical,intent(in) :: dbg
      real(kind(1.D0)) :: RES,I1,SD,MXW,Z,EZ,tau
2     format(a,' ',6(1x,G12.5))
      RES=0.D0
      if ((t>0.D0).and.(t<SPS(it)%TMAX)) then
      ! pulse shape
        Z=lam+SPS(it)%TAUP(3)
        EZ=exp(-SPS(it)%TAUP(4)*LOG(Z))
        tau=SPS(it)%TAU1*(1.D0+SPS(it)%TAUP(1)*lam**2)*exp(-SPS(it)%TAUP(2)*EZ)
        I1=IEXP(t,tau,SPS(it)%d)/(tau*SPS(it)%d)
      ! Maxwell term
        MXW=SPS(it)%F1*MAXW_LAM(lam,SPS(it)%TEMP)
        if (SPS(it)%DELTA.ne.0) then
          EZ=exp(-SPS(it)%DELTA*LOG(lam))
          MXW=MXW*EZ
        endif
      ! slowing down term
        SD=SPS(it)%ISD/lam*FCO(lam,SPS(it)%ALPHASD,SPS(it)%LAMSD)
        if (dbg) write(*,2) 'ESS_2015_THERM ',t,lam,tau,MXW,SD,I1
        RES=(MXW+SD)*I1
      endif
      ESS_2015_THERM_DBG=RES
      end function ESS_2015_THERM_DBG

!-------------------------------------------------------------------------
      real(kind(1.D0)) function LPSS_ESS_COLD(it,t,lam)
! Long pulse, ESS thermal moderator,  2013 (Schonfeldt, eq. 4)
! Double exponential x cut-off function  + slow down term
!-------------------------------------------------------------------------
      integer, intent(in) :: it
      real(kind(1.D0)),intent(in) :: lam,T
      real(kind(1.D0)) :: Z,RES,I1,I2,CO1
      RES=0.D0
      if ((t>0.D0).and.(t<SPS(it)%TMAX)) then
      ! pulse shape 1
        I1=LPSS_F(t,SPS(it)%TAU1,SPS(it)%N1,SPS(it)%d)
      ! pulse shape 2
        I2=LPSS_F(t,SPS(it)%TAU2*lam,SPS(it)%N2,SPS(it)%d)
      ! para-H2, cut off term 1
        Z=FCO(lam,SPS(it)%ALPHAC,SPS(it)%LAMC)
        CO1=exp(SPS(it)%EXPO*log(Z))
      ! para-H2 moderator, term 1
        RES=I1*CO1*(SPS(it)%I1*exp(-SPS(it)%ALPHA1*lam)+ SPS(it)%I2*exp(-SPS(it)%ALPHA2*lam))
      ! para-H2 moderator, slowing down term
        RES=RES+I2*SPS(it)%ISD/lam*FCO(lam,SPS(it)%ALPHASD,SPS(it)%LAMSD)
      endif
      LPSS_ESS_COLD=RES
      end function LPSS_ESS_COLD


!-------------------------------------------------------------------------
      real(kind(1.D0)) function ESS_2015_COLD(it,t,lam)
! Long pulse, ESS cold moderator,  2015 (butterfly)
!-------------------------------------------------------------------------
      integer, intent(in) :: it
      real(kind(1.D0)),intent(in) :: lam,T
      real(kind(1.D0)) :: Z,RES,I1,PH2,CO1,SD,EZ,tau
      RES=0.D0
      if ((t>0.D0).and.(t<SPS(it)%TMAX)) then
      ! pulse shape 1
        Z=lam+SPS(it)%TAUP(3)
        EZ=exp(-SPS(it)%TAUP(4)*LOG(Z))
        tau=SPS(it)%TAU1*(1.D0+SPS(it)%TAUP(1)*lam**2)*exp(-SPS(it)%TAUP(2)*EZ)
        I1=IEXP(t,tau,SPS(it)%d)/(tau*SPS(it)%d)
      ! para-H2, cut off term 1
        Z=FCO(lam,SPS(it)%ALPHAC,SPS(it)%LAMC)
        CO1=exp(SPS(it)%EXPO*log(Z))
      ! para-H2 moderator, term 1
        PH2=SPS(it)%I1*exp(-SPS(it)%ALPHA1*lam)+ SPS(it)%I2*exp(-SPS(it)%ALPHA2*lam)
      ! slowing down term
        SD=SPS(it)%ISD/lam*FCO(lam,SPS(it)%ALPHASD,SPS(it)%LAMSD)
        RES=I1*(CO1*PH2 + SD)
      endif
      ESS_2015_COLD=RES
      end function ESS_2015_COLD


!-------------------------------------------------------------------------
      SUBROUTINE LPSSA_DEBUG(it,t,lam)
! Long pulse flux at given time and wavelength
! 2nd version with 2 shape components for ambiente coupled moderator
!-------------------------------------------------------------------------
      integer, intent(in) :: it
      real(kind(1.D0)),intent(in) :: lam,T
      real(kind(1.D0)) :: RES,J1,J2
1     format(a,' = ',3(G12.5,2x))
      RES=0.D0
      if (t>0.D0) then
        write(*,1) 'LPSSA_DEBUG t,lam ',t,lam
        write(*,1) 'LPSSA_DEBUG F1,ISD ',SPS(it)%F1,SPS(it)%ISD
        write(*,1) 'LPSSA_DEBUG M ',MAXW_LAM(lam,SPS(it)%TEMP)
        write(*,1) 'LPSSA_DEBUG I1a',LPSS_F(t,SPS(it)%TAU1,SPS(it)%N1,SPS(it)%d)
        write(*,1) 'LPSSA_DEBUG I1b',LPSS_F(t,SPS(it)%TAUA,SPS(it)%N1,SPS(it)%d)
        write(*,1) 'LPSSA_DEBUG ksi',1.D0/lam*FCO(lam,SPS(it)%ALPHASD,SPS(it)%LAMSD)
        write(*,1) 'LPSSA_DEBUG I2',LPSS_F(t,SPS(it)%TAU2*lam,SPS(it)%N2,SPS(it)%d)
        J1=LPSS_F(t,SPS(it)%TAU1,SPS(it)%N1,SPS(it)%d)
        J2=LPSS_F(t,SPS(it)%TAUA,SPS(it)%N1,SPS(it)%d)
        RES=SPS(it)%F1*MAXW_LAM(lam,SPS(it)%TEMP)*(J1+J2)
        RES=RES+SPS(it)%ISD/lam*FCO(lam,SPS(it)%ALPHASD,SPS(it)%LAMSD)*LPSS_F(t,SPS(it)%TAU2*lam,SPS(it)%N2,SPS(it)%d)
        write(*,1) 'LPSSA_DEBUG RES',RES
      endif
      end SUBROUTINE LPSSA_DEBUG


      end MODULE SOURCES_PULSE
