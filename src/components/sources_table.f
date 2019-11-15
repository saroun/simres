!//////////////////////////////////////////////////////////////////////
!////  $Id: sources_table.f,v 1.27 2019/08/16 17:16:26 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2011, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.27 $
!////     $Date: 2019/08/16 17:16:26 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Lookup tables for source flux distribution
!////
!////////////////////////////////////////////////////////////////////////
      MODULE SOURCES_TABLE
      use SIMATH
      use SOURCES_PULSE
      use SOURCES_ISIS
      USE MESSAGES
      use TRACELIB,only:SURFACE_XZ,TLIB_INF
      use FILETOOLS
      implicit none
      save
      private

      ! names of all tables separated by :
      character*(256) :: FLUX_TABLE_NAMES
      integer,parameter :: SRF_MAX=16
      integer,parameter :: PLS_NT=128
      integer,parameter :: PLS_NLAM=256

      TYPE SOURCE_TABLE; SEQUENCE
! flux lookup table, lambda dependence
        REAL(KIND(1.D0)) :: FLXLAM(256),FLXDIST(256),FLXDLAM ! dPhi/dLambda [1e12 n/s/cm2/Ang]
        INTEGER :: FLXN,FLXLOG ! number of points, log-scale
        character*(128) :: TABNAME

! pulse spallation source lookup table (time x lambda) (col x row)
! PLS_FLUX is dPhi/dLam/dTime in [10^10/s/cm^2/A/us]


        logical :: PLS_valid ! set true if table is ready
        REAL(KIND(1.D0)) :: PLS_FLUX(PLS_NLAM,PLS_NT),PLS_FLUX_AVE(PLS_NLAM),PLS_FLUX_INT
        REAL(KIND(1.D0)) :: PLS_DT,PLS_DLAM,PLS_L0 ! time-step, lambda-step, lambda_min

! moderator surface shape (z-profile)
! composed of flat surfaces (polygonal) defined by pairs of linear coefficients
! zs = SRF_A + SRF_B*x

        integer :: SRF_NS
        REAL(KIND(1.D0)) :: SRF_A(SRF_MAX),SRF_B(SRF_MAX)

! horizontal distribution (x, alpha)
        REAL(KIND(1.D0)) :: FLXHP(64,64),FLXHX,FLXHA,FLXHX0,FLXHA0 ! array, spatial and angular half-widths, centres
        INTEGER :: FLXHNX,FLXHNA ! actual table size

! vertical distribution (y, phi)
        REAL(KIND(1.D0)) :: FLXVP(64,64),FLXVX,FLXVA,FLXVX0,FLXVA0 ! array, spatial and angular half-widths, centres
        INTEGER :: FLXVNX,FLXVNA ! actual table size

! flux lookup table, pulse shape
        REAL(KIND(1.D0)) :: FLXTIME(256),FLXPULSE(256),FLXDTIME, FLXPULSE_AVE ! time [us], P(time), step, integral [us]
        INTEGER :: FLXNP ! number of points

! radial distribution
        REAL(KIND(1.D0)) :: FLXRAD(256),FLXRDIST(256),FLXRSTEP ! radius [mm], weight, step [mm]
        INTEGER :: FLXNR ! number of points

      END TYPE SOURCE_TABLE

      integer :: STN
      integer,parameter :: STMAX=2
      TYPE(SOURCE_TABLE) :: ST(STMAX)

! clip range of ki
      REAL(KIND(1.D0)) :: CLIP_LMIN,CLIP_LMAX

! 2*mn*kB/h_bar^2 [Angstrom^-2*K^-1]
      REAL(KIND(1.D0)),parameter :: CMAXW=4.15871716D-2

      logical,private :: dbg=.false.

      ! info about table limits
      public SOURCE_TABLE_NUM, FLUX_RANGE,PULSE_RANGE,FLUX_RANGE_DIV,FLUX_SET_RANGE
      ! read table from file
      public READ_FLUX_TABLE,REPORT_FLUX_TABLE,FLUX_TABLE_NAMES,CLEAR_FLUX_TABLE,GET_FLXPULSE_AVE
      ! get value functions
      public HOMOG_TAB,PULSE_TAB,FLUX_TAB,SPS_TAB_AVE,SPS_TAB,FLUX_SURFACE_TIME
      ! is table defined?
      public FLUX_TABLE_DEF,PULSE_TABLE_DEF,SPS_TABLE_DEF,HOMOG_TABLE_DEF,SURFACE_TABLE_DEF
      public SOURCE_TABLE_IS_PULSED
      contains

!---------------------------------------------------------
      SUBROUTINE CLEAR_FLUX_TABLE
!---------------------------------------------------------
      integer :: k
      STN=0
      FLUX_TABLE_NAMES='none'
      do k=1,STMAX
        ST(k)%flxn=0
        ST(k)%flxnr=0
        ST(k)%flxnp=0
        ST(k)%flxhnx=0
        ST(k)%flxvnx=0
        ST(k)%flxhna=0
        ST(k)%flxvna=0
        ST(k)%FLXHX0=0.D0
        ST(k)%FLXHA0=0.D0
        ST(k)%FLXVX0=0.D0
        ST(k)%FLXVA0=0.D0
        ST(k)%PLS_valid=.false.
        ST(k)%SRF_NS=0
        ST(k)%TABNAME='none'
      enddo
      end SUBROUTINE CLEAR_FLUX_TABLE



!---------------------------------------------------------
      logical function SOURCE_TABLE_IS_PULSED()
!---------------------------------------------------------
      logical :: LG1
      integer :: it
      LG1=.false.
      do it=1,STN
        LG1=(LG1.or.SPS_TABLE_DEF(it).or.PULSE_TABLE_DEF(it))
      enddo
      SOURCE_TABLE_IS_PULSED =LG1
      end function SOURCE_TABLE_IS_PULSED

!---------------------------------------------------------
      integer function SOURCE_TABLE_NUM()
!---------------------------------------------------------
      SOURCE_TABLE_NUM=STN
      end function SOURCE_TABLE_NUM

!---------------------------------------------------------
      logical function SURFACE_TABLE_DEF(k)
!---------------------------------------------------------
      integer,intent(in) :: k
      if ((k<=0).or.(k>STN)) then
        SURFACE_TABLE_DEF=.false.
        return
      endif
      SURFACE_TABLE_DEF=(ST(k)%SRF_NS.gt.0)
      end function SURFACE_TABLE_DEF

!---------------------------------------------------------
      logical function FLUX_TABLE_DEF(k)
!---------------------------------------------------------
      integer,intent(in) :: k
      if ((k<=0).or.(k>STN)) then
        FLUX_TABLE_DEF=.false.
        return
      endif
      FLUX_TABLE_DEF=(ST(k)%FLXN.gt.1)
      end function FLUX_TABLE_DEF

!---------------------------------------------------------
      logical function PULSE_TABLE_DEF(k)
!---------------------------------------------------------
      integer,intent(in) :: k
      if ((k<=0).or.(k>STN)) then
        PULSE_TABLE_DEF=.false.
        return
      endif
      PULSE_TABLE_DEF=(ST(k)%FLXNP.gt.1)
      end function PULSE_TABLE_DEF

!---------------------------------------------------------
      logical function SPS_TABLE_DEF(k)
!---------------------------------------------------------
      integer,intent(in) :: k
      if ((k<=0).or.(k>STN)) then
        SPS_TABLE_DEF=.false.
        return
      endif
      SPS_TABLE_DEF=((ST(k)%PLS_valid).and.(ST(k)%PLS_L0>0.D0))
      end function SPS_TABLE_DEF

!---------------------------------------------------------
      logical function HOMOG_TABLE_DEF(k)
!---------------------------------------------------------
      integer,intent(in) :: k
      if ((k<=0).or.(k>STN)) then
        HOMOG_TABLE_DEF=.false.
        return
      endif
      HOMOG_TABLE_DEF=((ST(k)%FLXNR>1).or.(ST(k)%FLXHNX>1).or.(ST(k)%FLXVNX>1))
      end function HOMOG_TABLE_DEF

!---------------------------------------------------------
      SUBROUTINE FLUX_SET_RANGE(k0,dkk)
!  Set clip range for ki if dk>0
!---------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: k0,dkk
      if (dkk.le.0.D0) then
        CLIP_LMIN=0.D0
        CLIP_LMAX=0.D0
      else
        CLIP_LMIN=max(0.D0,TWOPI/k0*(1.D0-0.5D0*dkk))
        CLIP_LMAX=TWOPI/k0*(1.D0+0.5D0*dkk)
      endif
      end SUBROUTINE FLUX_SET_RANGE

!---------------------------------------------------------
      logical FUNCTION FLUX_SURFACE_TIME(it,R,K,DT,only_positive)
! get shortest time>0 to the surface defined in the source table
! only_positive > 0: allow only positive transition times
! R,K current neutron coordinates
!---------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: R(3),K(3)
      integer,intent(in) :: it,only_positive
      REAL(KIND(1.D0)),intent(out) :: DT
      REAL(KIND(1.D0)) :: t
      integer :: i
      logical :: LOG1
      LOG1=.false.
      DT=TLIB_INF
      do i=1,ST(it)%SRF_NS
        call SURFACE_XZ(R,K,ST(it)%SRF_A(i),ST(it)%SRF_B(i),1,t)
        if ((abs(t)<DT).and.((only_positive.le.0).or.(t>0.D0))) then
          DT=t
          LOG1=.true.
        endif
      enddo
      FLUX_SURFACE_TIME=LOG1
      END FUNCTION FLUX_SURFACE_TIME




!---------------------------------------------------------
      REAL(KIND(1.D0)) FUNCTION GET_FLXPULSE_AVE(it)
      integer,intent(in) :: it
        GET_FLXPULSE_AVE=ST(it)%FLXPULSE_AVE
      end FUNCTION GET_FLXPULSE_AVE


!---------------------------------------------------------
      REAL(KIND(1.D0)) FUNCTION SPS_TAB(it,FLX,time,k)
! Returns instanteneous flux dPhi/d^3K/dTime in [1e14/us/cm^2*Ang^3]
! time is in sec.hbar/m
! Using SPS table (SOURCE_PULSE module)
! FLX is normalizing factor
! Note: Use FLX=1 to respect the table scaling.
!---------------------------------------------------------
      integer,intent(in) :: it
      REAL(KIND(1.D0)),intent(in) :: FLX,time,k
      REAL(KIND(1.D0)) :: V(2),RES,K02,T,T0
      integer :: NL,NT
      integer,parameter :: meth=1
  ! center table at t0=middle of the range
      T0=0.5D0*ST(it)%PLS_DT*(PLS_NT-1)
      T=time+T0
      select case(meth)
      case(1)
  ! search 2D-table
        V=(/TWOPI/k,T/)
        NL=PLS_NLAM
        NT=PLS_NT
        RES=LOGLINTERP2D(ST(it)%PLS_FLUX,NL,NT,NL,ST(it)%PLS_L0,0.D0,ST(it)%PLS_DLAM,ST(it)%PLS_DT,V)
      case(2)
  ! alternative: calculate analytically
        RES=PULSE_FLUX(it,T,TWOPI/k)
      end select
      ! table is dPhi/dLam/dTime in [10^10/us/cm^2/A]
      ! convert to dPhi/d^3K/dTime in [1e14/us/cm^2*Ang^3]
	  ! NOTE; program integrates over time x hbar/m
      K02=k**2
      SPS_TAB=FLX/2.D0/K02**2*RES*1.D-4
      end FUNCTION SPS_TAB

!---------------------------------------------------------
      REAL(KIND(1.D0)) FUNCTION SPS_TAB_AVE(it,FLX,k)
! Returns time-integrated flux dPhi/d^3K in [1e14/pulse/cm^2*Ang^3]
! Integration is over one pulse only. Multiply by source frequency in Hz
! to get the steady-state averaged flux.
! Using SPS table (see SOURCE_PULSE module).
! FLX is normalizing factor.
! Note: Use FLX=1 to respect the table scaling.
!---------------------------------------------------------
      integer,intent(in) :: it
      REAL(KIND(1.D0)),intent(in) :: FLX,k
      REAL(KIND(1.D0)) :: RES,K02
      RES=LOGINTERP(ST(it)%PLS_FLUX_AVE(1),PLS_NLAM,ST(it)%PLS_L0,ST(it)%PLS_DLAM,TWOPI/k)
      ! table is dPhi/dLam in [10^10/pulse/cm^2/A]
      ! convert to dPhi/d^3K in [1e14/pulse/cm^2*Ang^3]
      K02=k**2
      SPS_TAB_AVE=FLX/2.D0/K02**2*RES*1.D-4
      end FUNCTION SPS_TAB_AVE

!---------------------------------------------------------
      REAL(KIND(1.D0)) FUNCTION FLUX_TAB(it,FLX,K0)
! Returns dPhi/d^3K in [1e14/s/cm^2*Ang^3] from the table
! FLX is the integral neutron flux [1e14/s/cm^2]
! Note: Use FLX=1 to respect the table scaling.
!---------------------------------------------------------
      integer,intent(in) :: it
      REAL(KIND(1.D0)),intent(in) :: FLX,K0
      REAL(KIND(1.D0)) :: RES,LAM,K02
! Lookup table
! FLXDIST = dPhi/dLambda in [1e12/s/cm^2/Ang]
      LAM=2*PI/K0
! search wavelength distribution
      if (ST(it)%FLXLOG.GT.0) then
  ! logarithmic scale
        RES=LOGINTERP(ST(it)%FLXDIST(1),ST(it)%FLXN,ST(it)%FLXLAM(1),ST(it)%FLXDLAM,LAM)
      else
  ! linear scale
        RES=LINTERP(ST(it)%FLXDIST(1),ST(it)%FLXN,ST(it)%FLXLAM(1),ST(it)%FLXDLAM,LAM)
      endif
! convert to k3-distribution
! FLXLAM is dPhi/dLam in [10^12/s/cm^2/A]
! FLUX_TAB returns dPhi/d^3K in [1e14/s/cm^2*Ang^3]
      K02=K0**2
      FLUX_TAB=FLX/2.D0/K02**2*RES*1.D-2
      END FUNCTION FLUX_TAB

!---------------------------------------------------------
      REAL(KIND(1.D0)) FUNCTION HOMOG_TAB(it,R,K,K0)
! Returns common weight for homogeneity distributions
!---------------------------------------------------------
      integer,intent(in) :: it
      REAL(KIND(1.D0)),intent(in) :: R(3),K(3),K0
      REAL(KIND(1.D0)) :: RES
      RES=1.D0
      if (ST(it)%FLXNR>0) then
        if(dbg) write(*,*) 'R= ',R
        RES=RES*RADIAL_TAB(it,SQRT(R(1)**2+R(2)**2+R(3)**2))
      else
        RES=RES*DIVH_TAB(it,R(1),K(1),K0)
        RES=RES*DIVV_TAB(it,R(2),K(2),K0)
      endif
      HOMOG_TAB=RES
      END FUNCTION HOMOG_TAB

!---------------------------------------------------------
      REAL(KIND(1.D0)) FUNCTION DIVH_TAB(it,X,KX,K0)
! Returns weight for horizontal space x divergence distribution
!---------------------------------------------------------
      integer,intent(in) :: it
      REAL(KIND(1.D0)),intent(in) :: X,KX,K0
      REAL(KIND(1.D0)) :: RES,V(2),DI,DJ,X0,A0
      RES=1.D0
      if (ST(it)%FLXHNX.GT.0) THEN
        V(1)=X
        V(2)=KX/K0
        X0=-ST(it)%FLXHX+ST(it)%FLXHX0
        A0=-ST(it)%FLXHA+ST(it)%FLXHA0
        DI=2.D0*ST(it)%FLXHX/(ST(it)%FLXHNX-1)
        DJ=2.D0*ST(it)%FLXHA/(ST(it)%FLXHNA-1)
        RES=LINTERP2D(ST(it)%FLXHP,ST(it)%FLXHNX,ST(it)%FLXHNA,64,X0,A0,DI,DJ,V)
      endif
      DIVH_TAB=RES
      END FUNCTION DIVH_TAB

!---------------------------------------------------------
      REAL(KIND(1.D0)) FUNCTION DIVV_TAB(it,Y,KY,K0)
! Returns weight for horizontal space x divergence distribution
!---------------------------------------------------------
      integer,intent(in) :: it
      REAL(KIND(1.D0)),intent(in) :: Y,KY,K0
      REAL(KIND(1.D0)) :: RES,V(2),DI,DJ,X0,A0
      RES=1.D0
      if (ST(it)%FLXVNX.GT.0) THEN
        V(1)=Y
        V(2)=KY/K0
        X0=-ST(it)%FLXVX+ST(it)%FLXVX0
        A0=-ST(it)%FLXVA+ST(it)%FLXVA0
        DI=2.D0*ST(it)%FLXVX/(ST(it)%FLXVNX-1)
        DJ=2.D0*ST(it)%FLXVA/(ST(it)%FLXVNA-1)
        RES=LINTERP2D(ST(it)%FLXVP,ST(it)%FLXVNX,ST(it)%FLXVNA,64,X0,A0,DI,DJ,V)
      endif
      DIVV_TAB=RES
      END FUNCTION DIVV_TAB

!---------------------------------------------------------
      REAL(KIND(1.D0)) FUNCTION PULSE_TAB(it,T)
! Returns pulse amplitude for given time in [us]
!---------------------------------------------------------
      integer,intent(in) :: it
      REAL(KIND(1.D0)),intent(in) :: T
      REAL(KIND(1.D0)) :: RES
      RES=1.D0
      if (ST(it)%FLXNP>0) then
        RES=LINTERP(ST(it)%FLXPULSE(1),ST(it)%FLXNP,ST(it)%FLXTIME(1),ST(it)%FLXDTIME,T)
      endif
      PULSE_TAB=RES
      END FUNCTION PULSE_TAB

!---------------------------------------------------------
      REAL(KIND(1.D0)) FUNCTION RADIAL_TAB(it,R)
! Returns weight for radial distribution
!---------------------------------------------------------
      integer,intent(in) :: it
      REAL(KIND(1.D0)),intent(in) :: R
      REAL(KIND(1.D0)) :: RES
      RES=1.D0
      if (ST(it)%FLXNR>0) then
        RES=LINTERP(ST(it)%FLXRDIST(1),ST(it)%FLXNR,ST(it)%FLXRAD(1),ST(it)%FLXRSTEP,R)
      endif
      RADIAL_TAB=RES
      END FUNCTION RADIAL_TAB

!---------------------------------------------------------
      SUBROUTINE FLUX_RANGE(temp,kmin,kmax)
!  Get range of k-vector value, either from table or Maxwell
!---------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: temp
      REAL(KIND(1.D0)),intent(out) :: kmin,kmax
      REAL(KIND(1.D0)) :: VK,lmax,kmi,kma
      integer :: it
1     format(a,': ',6(G12.5,1x))
      kmax=-1e30
      kmin=1e30
      do it=1,STN
  ! lookup table
        !write(*,1) 'FLUX_RANGE L0, valid',ST(it)%PLS_L0,ST(it)%PLS_valid
        if (SPS_TABLE_DEF(it)) then
          lmax=ST(it)%PLS_L0*exp((PLS_NLAM-1)*ST(it)%PLS_DLAM)
          kma=TWOPI/ST(it)%PLS_L0
          kmi=TWOPI/lmax
          kmin=min(kmi,kmin)
          kmax=max(kma,kmax)
          !write(*,1) 'FLUX_RANGE SPS L0,DLAM',ST(it)%PLS_L0,ST(it)%PLS_DLAM
          !write(*,1) 'FLUX_RANGE SPS it, kmi,kma',it,kmi,kma,kmin,kmax
        else if (FLUX_TABLE_DEF(it)) then
          kma=TWOPI/ST(it)%FLXLAM(1)
          kmi=TWOPI/ST(it)%FLXLAM(ST(it)%FLXN)
          kmin=min(kmi,kmin)
          kmax=max(kma,kmax)
          !write(*,1) 'FLUX_RANGE TAB',it,kmi,kma,kmin,kmax
        endif
      enddo
   ! Maxwell
      if (kmax<kmin) then
        VK=SQRT(CMAXW*TEMP)
        kmin=5.D-2*VK
        kmax=5.D0*VK
        !write(*,1) 'FLUX_RANGE Maxwell',kmin,kmax
      endif
  ! apply clip range
      if (CLIP_LMAX.gt.0.D0) then
        kmin=max(kmin,TWOPI/CLIP_LMAX)
        kmax=min(kmax,TWOPI/CLIP_LMIN)
        !write(*,1) 'FLUX_RANGE CLIP',kmin,kmax
        !write(*,*) 'FLUX_RANGE CLIP ',CLIP_LMIN,CLIP_LMAX
      endif
      if (kmax<kmin) then
        kmin=0.D0
        kmax=0.D0
      endif
      end SUBROUTINE FLUX_RANGE


!---------------------------------------------------------
      SUBROUTINE PULSE_RANGE(tmin,tmax)
! Get time range [m/h units] in pulse mode
!---------------------------------------------------------
      REAL(KIND(1.D0)),intent(out) :: tmin,tmax
      REAL(KIND(1.D0)) :: Z,tmi,tma
      integer :: it
      tmax=-1e30
      tmin=1e30
      tmi=0.D0
      tma=0.D0
      do it=1,STN
        if (SPS_TABLE_DEF(it)) then
          Z=0.5D0*ST(it)%PLS_DT*(PLS_NT-1)*HOVM
          tmi=-Z
          tma=+Z
        else if (PULSE_TABLE_DEF(it)) then
          tmi=ST(it)%FLXTIME(1)*HOVM
          tma=ST(it)%FLXTIME(ST(it)%FLXNP)*HOVM
        endif
        tmin=min(tmi,tmin)
        tmax=max(tma,tmax)
      enddo
      if (tmax<tmin) then
        tmin=0.D0
        tmax=0.D0
      endif
      end SUBROUTINE PULSE_RANGE


!---------------------------------------------------------
      SUBROUTINE FLUX_RANGE_DIV(divx,divy)
!  Get range of k-vector value
!---------------------------------------------------------
      REAL(KIND(1.D0)),intent(out) :: divx,divy
      integer :: it
      divx=0.D0
      divy=0.D0
      do it=1,STN
        if (ST(it)%FLXHNX.GT.0) divx=max(divx,2*ST(it)%FLXHA)
        if (ST(it)%FLXVNX.GT.0) divy=max(divy,2*ST(it)%FLXVA)
      enddo
      end SUBROUTINE FLUX_RANGE_DIV

!---------------------------------------------------------------
      SUBROUTINE READ_FLUX_SPS(it,IU,ilin,LINE,MSG,ierr)
! read spallation source datadefined in the SOURCE_PULSE module
!---------------------------------------------------------------
      INTEGER,intent(in) :: IU,it
      INTEGER,intent(out) :: ilin,ierr
      character(*),intent(out) :: LINE
      CHARACTER*(1024) S
      integer :: I,J,IE,IS,IL
      integer,parameter :: ND=19 ! max. number of parameters
      integer :: NPAR,IDX,NP
      character*(96) :: NAM
      logical:: IRD(ND),LOG1
      character*(*) :: MSG
      character*(16) :: ID
      character*256 :: CNUM1
      real(kind(1.D0)) :: Z(64),t,lam,y(PLS_NT)

      ierr=i_ERR
      MSG=' '
3     format(a,': ',6(1x,G12.5))

! evaluate header: look for type flag
! set the name space according to the type
      call MKUPCASE(LINE)
      J=INDEX(LINE,'TYPE')
      if (J.gt.0) then
        call READ_I4('TYPE', LINE(J:), I, IERR)
        if (IERR.eq.0) then
          if ((I.gt.0).and.(I.le.7)) then
            SPS(it)%TYP=I
            ! write(*,3) 'READ_FLUX_SPS type=',SPS(it)%TYP
            ierr=0
            call PULSE_SET_DEFAULT(it)
          else
            MSG='Unknown TYPE.'
          endif
        else
          MSG='Missing TYPE.'
        endif
        call PULSE_GET_NAMES(it,NAM)
        call COUNTPAR(NAM,':',NPAR)
      endif
      if (ierr.ne.0) return
! read parameters in format NAME=value
      S=' '
      ierr=0
      IRD=.false.
      MSG=' '
      !write(*,3) 'READ_FLUX_SPS type=',SPS(it)%TYP,' NPAR=',NPAR
      !write(*,*) trim(NAM)
      do while (IERR.eq.0)
        call READ_TABLE_ROW(IU,ilin,S,ierr)
        call MKUPCASE(S)
      ! end tag - stop parsing and read next line
        if (INDEX(S,'END SPALLATION_SOURCE').eq.1) then
          if (ierr.ne.i_EOF) then
            call READ_TABLE_ROW(IU,ilin,S,IE)
            if (IE.eq.0) ierr=i_VALID
          endif
        endif
        if ((ierr.eq.0).and.len_trim(S)>0) then
          call PULSE_READ_PARAM(it,S,ID,IDX,IE)
          !write(*,3) trim(ID), IDX
          if ((IE==0).and.(IDX<=NPAR).and.(IDX>0)) then
            IRD(IDX)=.true.
            !write(*,3) trim(S)//' | '//trim(ID),IDX,IE,IRD(IDX)
          else
            !write(*,3) ' unknnown ',IDX,IE,IRD(IDX)
          endif
        endif
      enddo
      !write(*,*) 'END SPALLATION_SOURCE ierr=',ierr
      ! return last line
      select case(ierr)
      case(i_EOF,i_VALID)
        LINE=trim(S)
      case default
        LINE=' '
      end select
! check that all required parameters were found
! if not, print warning
      LOG1=.true.
      !write(*,*) 'READ_FLUX_SPS read'
      do i=1,NPAR
        LOG1=(LOG1.and.IRD(i))
        !call FINDSTRPAR(NAM,':',i,IS,IL)
        !write(*,*) i,' ',NAM(IS:IS+IL-1),' ',IRD(i)
      enddo
      if (.not.LOG1) then
        !write(*,*) 'undefined parameters:'
        S=' '
        do i=1,NPAR
          if (.not.IRD(i)) then
            call FINDSTRPAR(NAM,':',i,IS,IL)
            if (IL>0) then
              ID=NAM(IS:IS+IL-1)
              call PULSE_GET_PARAM(it,ID,Z,NP,IDX)
              call ARRAY2STR(Z(1),NP,CNUM1)
              call FLOAT2STR(Z,CNUM1)
              !write(*,3) NAM(IS:IS+IL-1)//'='//trim(CNUM1),i
              S=trim(S)//', '//NAM(IS:IS+IL-1)//'='//trim(CNUM1)
            endif
          endif
        enddo
        if (len_trim(S)>0) MSG=trim(MSG)//' Missig parameters, using defaults'//trim(S)
      endif
      select case(ierr)
      case(i_EOF,i_VALID)
! valid input, calculate lookup  table
1     format(a,$)
2     format(a,' ',6(1x,G12.5))
        ST(it)%PLS_DLAM=LOG(SPS(it)%LMAX/SPS(it)%LMIN)/(PLS_NLAM-1)
        ST(it)%PLS_L0=SPS(it)%LMIN
        ST(it)%PLS_DT=SPS(it)%TMAX/(PLS_NT-1)
        write(*,2) 'calculating pulse table ',it
        !write(*,2) 'lmin,lmax',SPS(it)%LMIN,SPS(it)%LMAX
        !write(*,2) 'L0,DLAM,DT',ST(it)%PLS_L0,ST(it)%PLS_DLAM,ST(it)%PLS_DT
        !write(*,*) 'dim=',PLS_NT,PLS_NLAM
        do j=1,PLS_NT
          t=(j-1)*ST(it)%PLS_DT
          !! write(*,2) 'READ_FLUX_SPS ',j,t
          do i=1,PLS_NLAM
            lam=SPS(it)%LMIN*exp((i-1)*ST(it)%PLS_DLAM)
            !PLS_FLUX(i,j)=0
            !dbg=((t<2500).and.(t>2000).and.(lam<6).and.(lam>1))
            !dbg=((j==PLS_NT).and.(i>PLS_NLAM-2))
            !PLS_FLUX(i,j)=PULSE_FLUX_DBG(t,lam,dbg)
            ST(it)%PLS_FLUX(i,j)=PULSE_FLUX(it,t,lam)
          enddo
        enddo
        !write(*,2) 'calc average '
        ! get flux per pulse
        do i=1,PLS_NLAM
          do j=1,PLS_NT
            y(j)=ST(it)%PLS_FLUX(i,j)
          enddo
          ST(it)%PLS_FLUX_AVE(i)=INTEG1D(y,ST(it)%PLS_DT,PLS_NT)
        enddo
        ST(it)%PLS_FLUX_INT=INTEG1DLOG(ST(it)%PLS_FLUX_AVE,ST(it)%PLS_L0,ST(it)%PLS_DLAM,PLS_NLAM)
        write(*,*) ' done.'
        ST(it)%PLS_valid=.true.
      case default
        ST(it)%PLS_valid=.false.
      end select
      end SUBROUTINE READ_FLUX_SPS


!---------------------------------------------------------------
      SUBROUTINE READ_FLUX_ISIS(it,IU,ilin,LINE,ierr)
! read spallation source data defined in the SOURCE_PULSE module
!---------------------------------------------------------------
      INTEGER,intent(in) :: it
      INTEGER,intent(in) :: IU
      INTEGER,intent(out) :: ilin,ierr
      character(*),intent(out) :: LINE
      CHARACTER*(1024) S
      integer :: IE,NL,NT
      character*(16) :: ID
      real(kind(1.D0)) :: DL,L0

      call ISIS_DEFAULT

3     format(a,': ',6(1x,G12.5))
! read parameters in format NAME=value
      S=' '
      ierr=0
      do while (IERR.eq.0)
        call READ_TABLE_ROW(IU,ilin,S,ierr)
        !write(*,3) trim(S) , ilin, ierr
      ! end tag - stop parsing and read next line
        if (INDEX(S,'END SPALLATION_SOURCE_ISIS').eq.1) then
          if (ierr.ne.i_EOF) then
            call READ_TABLE_ROW(IU,ilin,S,ierr)
            !write(*,3) 'END: '//trim(S) , ilin, ierr
            if (ierr.eq.i_OK) ierr=i_VALID
          endif
        endif
        if ((ierr.eq.0).and.(len_trim(S)>0)) then
          call READ_PARAM_ISIS(S,ID,IE)
        endif
      enddo
      !write(*,*) 'END SPALLATION_SOURCE ierr=',ierr
      ! return last line
      select case(ierr)
      case(i_EOF,i_VALID)
        LINE=trim(S)
      case default
        LINE=' '
      end select
      select case(ierr)
      case(i_EOF,i_VALID)
! valid input, calculate lookup  table
2       format(a)
        write(*,2) 'reading ISIS source table '//trim(ISIS%FNAME)
        NL=PLS_NLAM
        NT=PLS_NT
        call PARSE_ISIS_TABLE(ST(it)%PLS_FLUX,ST(it)%PLS_FLUX_AVE,L0,DL,ST(it)%PLS_DT,NL,NT,IE)
        !write(*,*) 'ISIS source table done, ierr=',IE,' ISIS%VALID=',ISIS%VALID, NL, NT
        select case(IE)
        case(i_EOF,i_OK)
          ST(it)%PLS_FLUX_INT=INTEG1DLOG(ST(it)%PLS_FLUX_AVE,L0,DL,NL)
          ST(it)%PLS_L0=L0
          ST(it)%PLS_DLAM=DL
          ST(it)%PLS_valid=ISIS%VALID
        case default
          ST(it)%PLS_valid=.false.
        end select
      case default
        ST(it)%PLS_valid=.false.
      end select
      !write(*,*) 'READ_FLUX_ISIS ierr=',ierr,', ilin=',ilin,', ISIS%VALID=',ISIS%VALID
      end SUBROUTINE READ_FLUX_ISIS

!---------------------------------------------------------------
      SUBROUTINE READ_FLUX_LAMBDA(it,IU,ilin,LINE,ierr)
! read wavelength distribution
! 1 line header + 2 columns: Lambda, dPhi/dLambda
! units are [A], [1e12/s/cm^2/A] (isotropic)
!---------------------------------------------------------------
      INTEGER,intent(in) :: it
      INTEGER,intent(in) :: IU
      INTEGER,intent(out) :: ilin,ierr
      character(*),intent(out) :: LINE
      CHARACTER*(1024) S
      integer :: J
! evalueate header: look for LOGSCALE flag
      ST(it)%FLXLOG=INDEX(trim(LINE),'LOGSCALE')
      if (INDEX(trim(LINE),'LOGSCALE=no').gt.0) ST(it)%FLXLOG=0
! read dPhi/dLambda
      ST(it)%flxn=0
      S=' '
      ierr=0
      do while (IERR.eq.0)
        call READ_TABLE_ROW(IU,ilin,S,ierr)
        J=LEN_TRIM(S)
        !write(*,*) 'READ_FLUX_LAMBDA ',trim(S),' ierr=',ierr
        if ((ierr.eq.0).and.(ST(it)%flxn.ge.256)) ierr=i_VALID
        if (ierr.eq.0) then
          Read(S(1:J),*,iostat=IERR,ERR=50) ST(it)%flxlam(ST(it)%flxn+1), ST(it)%flxdist(ST(it)%flxn+1)
50          if (ierr.ne.0) then
            ierr=i_VALID
          else
            ST(it)%flxn=ST(it)%flxn+1
          endif
        endif
      enddo
      if (ST(it)%flxn.le.1) ierr=i_ERR
      select case(ierr)
      case(i_EOF,i_VALID)
        LINE=trim(S)
      ! set step in lambda for log or linear scale
        IF (ST(it)%FLXLOG.GT.0) THEN
          ST(it)%FLXDLAM=LOG(ST(it)%FLXLAM(ST(it)%FLXN)/ST(it)%FLXLAM(1))/(ST(it)%FLXN-1)
        ELSE
          ST(it)%FLXDLAM=(ST(it)%FLXLAM(ST(it)%FLXN)-ST(it)%FLXLAM(1))/(ST(it)%FLXN-1)
        endif
      case default
        LINE=' '
      end select
      end SUBROUTINE READ_FLUX_LAMBDA


!---------------------------------------------------------------
      SUBROUTINE READ_FLUX_HOR(it,IU,ilin,LINE,ierr)
! read horizontal position x divergence distribution
! 1 line header contains: HORIZONTAL nx,na,dx,da
! na,da  .. divergence: number of colums and range (-da;+da) in [rad]
! nx,dx  .. position: number of rows and range (-dx;+dx) in [cm]
! angle is actually interpretted as sinus of the direction angle
! values: multiplies the flux obtained from a lambda table (or others with abs. flux)
!---------------------------------------------------------------
      INTEGER,intent(in) :: it
      INTEGER,intent(in) :: IU
      INTEGER,intent(out) :: ilin,ierr
      character(*),intent(out) :: LINE
      CHARACTER*(1024) S
      integer :: J,i,ic,NA
      REAL(KIND(1.D0)) :: A(6)
      ic=0
      ST(it)%FLXHNX=0
      ST(it)%FLXHNA=0
      ST(it)%FLXHX0=0.D0
      ST(it)%FLXHA0=0.D0
      S=' '
      ierr=0
      J=LEN_TRIM(LINE)
! evaluate header: read table ranges
      call READ_ARRAY8(LINE(11:J),6,A,NA)
      if (NA<4) then
        ierr=i_ERR
        return
      else
        ST(it)%FLXHNX=NINT(A(1))
        ST(it)%FLXHNA=NINT(A(2))
        ST(it)%FLXHX=A(3)
        ST(it)%FLXHA=A(4)
        if (NA>4) ST(it)%FLXHX0=A(5)
        if (NA>5) ST(it)%FLXHA0=A(6)
      endif
      ! Read(LINE(11:J),*,iostat=IERR,err=40) FLXHNX,FLXHNA,FLXHX,FLXHA
      ! call READ_ID_R8(LINE(11:J), ID, PARAM, MA, NA,IERR)
!40    if (ierr.ne.0) then
!        ierr=i_ERR
!        return
!      endif
      ST(it)%FLXHX=ST(it)%FLXHX*10.D0
      ST(it)%FLXHX0=ST(it)%FLXHX0*10.D0
      if ((ST(it)%FLXHNX.GT.64).OR.(ST(it)%FLXHNA.GT.64))  then
        ierr=i_ERR
        return
      endif
      do while (IERR.eq.0)
        call READ_TABLE_ROW(IU,ilin,S,ierr)
        J=LEN_TRIM(S)
        if ((ierr.eq.0).and.(ic.ge.ST(it)%FLXHNX)) ierr=i_VALID
        if ((J>0).and.(ierr.eq.0)) then
          Read(S(1:J),*,iostat=IERR,ERR=50) (ST(it)%FLXHP(ic+1,i),i=1,ST(it)%FLXHNA)
50        if (ierr.ne.0) then
            ierr=i_VALID
          else
            ic=ic+1
          endif
        endif
      enddo
      !write(*,*) 'READ_FLUX_HOR ic=',ic,' ierr=',ierr
      if (ic.ne.ST(it)%FLXHNX) ierr=i_ERR
      select case(ierr)
      case(i_EOF,i_VALID)
        LINE=trim(S)
      case default
        LINE=' '
      end select
      end SUBROUTINE READ_FLUX_HOR


!---------------------------------------------------------------
      SUBROUTINE READ_FLUX_VER(it,IU,ilin,LINE,ierr)
! read vertical position x divergence distribution
! 1 line header contains: VERTICAL ny,na,dy,da
! na,da  .. divergence: number of colums and range (-da;+da) in [rad]
! ny,dy  .. position: number of rows and range (-dx;+dx) in [cm]
! angle is actually interpretted as sinus of the direction angle
! values: multiplies the flux obtained from a lambda table (or others with abs. flux)
!---------------------------------------------------------------
      INTEGER,intent(in) :: it
      INTEGER,intent(in) :: IU
      INTEGER,intent(out) :: ilin,ierr
      character(*),intent(out) :: LINE
      CHARACTER*(1024) S
      integer :: J,i,ic,NA
      REAL(KIND(1.D0)) :: A(6)
      ic=0
      ST(it)%FLXVNX=0
      ST(it)%FLXVNA=0
      ST(it)%FLXVX0=0.D0
      ST(it)%FLXVA0=0.D0
      S=' '
      ierr=0
      J=LEN_TRIM(LINE)
! evalueate header: read table ranges
      call READ_ARRAY8(LINE(9:J),6,A,NA)
      if (NA<4) then
        ierr=i_ERR
        return
      else
        ST(it)%FLXVNX=NINT(A(1))
        ST(it)%FLXVNA=NINT(A(2))
        ST(it)%FLXVX=A(3)
        ST(it)%FLXVA=A(4)
        if (NA>4) ST(it)%FLXVX0=A(5)
        if (NA>5) ST(it)%FLXVA0=A(6)
      endif
!      Read(LINE(9:J),*,iostat=IERR,err=40) FLXVNX,FLXVNA,FLXVX,FLXVA
!40    if (ierr.ne.0) then
!        ierr=i_ERR
!        return
!      endif
      ST(it)%FLXVX=ST(it)%FLXVX*10.D0
      ST(it)%FLXVX0=ST(it)%FLXVX0*10.D0
      if ((ST(it)%FLXVNX.GT.64).OR.(ST(it)%FLXVNA.GT.64))  then
        ierr=i_ERR
        return
      endif
      do while (IERR.eq.0)
        call READ_TABLE_ROW(IU,ilin,S,ierr)
        J=LEN_TRIM(S)
        if ((ierr.eq.0).and.(ic.ge.ST(it)%FLXVNX)) ierr=i_VALID
        if (ierr.eq.0) then
          Read(S(1:J),*,iostat=IERR,ERR=50) (ST(it)%FLXVP(ic+1,i),i=1,ST(it)%FLXVNA)
50          if (ierr.ne.0) then
            ierr=i_VALID
          else
            ic=ic+1
          endif
        endif
      enddo
      !write(*,*) 'READ_FLUX_VER ic=',ic,' ierr=',ierr
      if (ic.ne.ST(it)%FLXVNX) ierr=i_ERR
      select case(ierr)
      case(i_EOF,i_VALID)
        LINE=trim(S)
      case default
        LINE=' '
      end select
      end SUBROUTINE READ_FLUX_VER


!---------------------------------------------------------------
      SUBROUTINE READ_FLUX_RAD(it,IU,ilin,LINE,ierr)
! read radial source distribution
! two columns with R [mm] and I [normalized to 1 at R=0)
!---------------------------------------------------------------
      INTEGER,intent(in) :: it
      INTEGER,intent(in) :: IU
      INTEGER,intent(out) :: ilin,ierr
      character(*),intent(out) :: LINE
      CHARACTER*(1024) S
      integer :: j
! read dPhi/dLambda
      ST(it)%FLXNR=0
      S=' '
      ierr=0
      do while (IERR.eq.0)
        call READ_TABLE_ROW(IU,ilin,S,ierr)
        J=LEN_TRIM(S)
        if ((ierr.eq.0).and.(ST(it)%flxnr.ge.256)) ierr=i_VALID
        if (ierr.eq.0) then
          Read(S(1:J),*,iostat=IERR,ERR=50) ST(it)%FLXRAD(ST(it)%FLXNR+1), ST(it)%FLXRDIST(ST(it)%FLXNR+1)
          !if (dbg) write(*,*) FLXNR+1,FLXRAD(FLXNR+1), FLXRDIST(FLXNR+1)
50          if (ierr.ne.0) then
            ierr=i_VALID
          else
            ST(it)%FLXNR=ST(it)%FLXNR+1
          endif
        endif
      enddo
      if (ST(it)%flxnr.le.1) ierr=i_ERR
      select case(ierr)
      case(i_EOF,i_VALID)
        LINE=trim(S)
      ! set step for radial distribution
        if (ST(it)%FLXNR.gt.1) then
          ST(it)%FLXRSTEP=(ST(it)%FLXRAD(ST(it)%FLXNR)-ST(it)%FLXRAD(1))/(ST(it)%FLXNR-1)
        ENDIF
      case default
        LINE=' '
      end select
      end SUBROUTINE READ_FLUX_RAD

!---------------------------------------------------------------
      SUBROUTINE READ_FLUX_PULSE(it,IU,ilin,LINE,ierr)
! read time  distribution
! two columns with time [us] and intensity (will be normlized to unit integral)
! modified Jan 2014: pulse profile is NOT normalized on unit integral
! instant brilliance  = flux(lambda) x pulse(time)
!---------------------------------------------------------------
      INTEGER,intent(in) :: it
      INTEGER,intent(in) :: IU
      INTEGER,intent(out) :: ilin,ierr
      character(*),intent(out) :: LINE
      CHARACTER*(1024) S
      logical :: LOG1
      real(kind(1.D0)) :: Z
      integer :: i,j
      real(kind(1.D0)) ,parameter :: EPS=2.D-2
! read dPhi/dLambda
      ST(it)%flxnp=0
      S=' '
      ierr=0
      do while (IERR.eq.0)
        call READ_TABLE_ROW(IU,ilin,S,ierr)
        J=LEN_TRIM(S)
        if ((ierr.eq.0).and.(ST(it)%flxnp.ge.256)) ierr=i_VALID
        if (ierr.eq.0) then
          Read(S(1:J),*,iostat=IERR,ERR=50) ST(it)%FLXTIME(ST(it)%FLXNP+1), ST(it)%FLXPULSE(ST(it)%FLXNP+1)
50        if (ierr.ne.0) then
            ierr=i_VALID
          else
            ST(it)%flxnp=ST(it)%flxnp+1
          endif
        endif
      enddo
      if (ST(it)%flxnp.le.1) ierr=i_ERR
      select case(ierr)
      case(i_EOF,i_VALID)
        LINE=trim(S)
      ! table step for interpolation
        ST(it)%FLXDTIME=(ST(it)%FLXTIME(ST(it)%FLXNP)-ST(it)%FLXTIME(1))/(ST(it)%FLXNP-1)
      ! get pulse integral is sec
        ST(it)%FLXPULSE_AVE=INTEG1D(ST(it)%FLXPULSE,ST(it)%FLXDTIME,ST(it)%FLXNP)/1.D6
! check that data are equidistant
        LOG1=.true.
        do i=1,ST(it)%FLXNP
          !FLXPULSE(i)=FLXPULSE(i)/Z
          if (i.gt.1) then
            LOG1=(LOG1.and.(abs(1.D0-(ST(it)%FLXTIME(i)-ST(it)%FLXTIME(i-1))/ST(it)%FLXDTIME).lt.EPS))
			!write(*,*) i, LOG1, FLXDTIME, FLXTIME(i)-FLXTIME(i-1)
          endif
        enddo
! set time=0 in the middle of the table range
        Z=(ST(it)%FLXTIME(ST(it)%FLXNP)+ST(it)%FLXTIME(1))/2.D0
        do i=1,ST(it)%FLXNP
          ST(it)%FLXTIME(i)=ST(it)%FLXTIME(i)-Z
        enddo
        if (.not.LOG1) ierr=i_ERR
      case default
        LINE=' '
      end select
      end SUBROUTINE READ_FLUX_PULSE


!---------------------------------------------------------------
      SUBROUTINE READ_FLUX_SURFACE(it,IU,ilin,LINE,ierr)
! read parameters of moderator surface
! header: SURFACE
! format: LINE= A B
! footer: END SURFACE
! each line defines one segment of a polygon, zs=A + B*x
!---------------------------------------------------------------
      INTEGER,intent(in) :: it
      INTEGER,intent(in) :: IU
      INTEGER,intent(out) :: ilin,ierr
      character(*),intent(out) :: LINE
      CHARACTER*(128) :: S
      CHARACTER*(16) :: ID
      integer :: NA,IE
      real(kind(1.D0)) :: Z(3)
      ST(it)%SRF_NS=0
      S=' '
      ierr=0
      do while (IERR.eq.0)
        call READ_TABLE_ROW(IU,ilin,S,ierr)
        call MKUPCASE(S)
      ! end tag - stop parsing and read next line
        if (INDEX(S,'END SURFACE').eq.1) then
          if (ierr.ne.i_EOF) then
            call READ_TABLE_ROW(IU,ilin,S,IE)
            if (IE.eq.0) ierr=i_VALID
          endif
        endif
        if (ST(it)%SRF_NS.ge.SRF_MAX) ierr=i_VALID
        if ((ierr.eq.0).and.len_trim(S)>0) then
          call READ_ID_R8(S,ID,Z,5,NA,IE)
          if ((IE==0).and.(NA>1)) then
            select case(trim(ID))
            case('LINE')
              ST(it)%SRF_NS=ST(it)%SRF_NS+1
              ST(it)%SRF_A(ST(it)%SRF_NS)=Z(1)
              ST(it)%SRF_B(ST(it)%SRF_NS)=Z(2)
            end select
          endif
        !else
        !  ierr=i_ERR
        endif
      enddo
      select case(ierr)
      case(i_EOF,i_VALID)
        LINE=trim(S)
      case default
        LINE=' '
      end select
      end SUBROUTINE READ_FLUX_SURFACE

!---------------------------------------------------------------
      SUBROUTINE READ_FLUX_NEW(it,IU,LINE,ilin,IERR,MSG)
! Read contents of flux table in new format
!---------------------------------------------------------------
      INTEGER,intent(in) :: it
      integer,intent(in) :: IU
      character(*),intent(out) :: LINE
      integer,intent(out) :: IERR,ilin
      character(*),intent(out) :: MSG
      integer :: IS,IL
      MSG=' '
      ierr=i_VALID
      !write(*,*) 'READ_FLUX_NEW, LINE=',trim(LINE)
      do while (ierr.eq.i_VALID)
        call FINDSTRPAR(LINE,' ',1,IS,IL)
        !write(*,*) 'READ_FLUX_NEW, LINE=',trim(LINE)
        if ((LEN_TRIM(LINE).gt.0).and.(IL.gt.0)) then
        select case(trim(LINE(IS:IS+IL-1)))
! wavelebgth distribution
        case('LAMBDA')
          call READ_FLUX_LAMBDA(it,IU,ilin,LINE,ierr)
          !write(*,*) 'READ_FLUX_NEW, lambda, flxn=',flxn,ierr
! time distribution
        case('PULSE')
          call READ_FLUX_PULSE(it,IU,ilin,LINE,ierr)
          !write(*,*) 'READ_FLUX_NEW, pulse, flxnp=',flxnp,ierr
! horizontal position x divergence
        case('HORIZONTAL')
          call READ_FLUX_HOR(it,IU,ilin,LINE,ierr)
          !write(*,*) 'READ_FLUX_NEW, horizontal, flxhn=',flxhnx,flxhna,ierr
          !write(*,*) 'LINE=[',trim(LINE),']'
          !ierr=i_EOF
! vertical position x divergence
        case('VERTICAL')
          call READ_FLUX_VER(it,IU,ilin,LINE,ierr)
          !write(*,*) 'READ_FLUX_NEW, vertical, flxvn=',flxvnx,flxvna,ierr
! radial distribution
        case('RADIAL')
          call READ_FLUX_RAD(it,IU,ilin,LINE,ierr)
          !write(*,*) 'READ_FLUX_NEW, radial, flxnr=',flxnr,ierr
! pulsed spallation source
        case('SPALLATION_SOURCE')
          !write(*,*) 'READ_FLUX_NEW, going to READ_FLUX_SPS'
          call READ_FLUX_SPS(it,IU,ilin,LINE,MSG,ierr)
          !write(*,*) 'READ_FLUX_NEW, spallation source ',ierr
! pulsed spallation source
        case('SPALLATION_SOURCE_ISIS')
          !write(*,*) 'READ_FLUX_NEW, going to READ_FLUX_ISIS'
          call READ_FLUX_ISIS(it,IU,ilin,LINE,ierr)
          !write(*,*) 'READ_FLUX_NEW, spallation source isis',ierr
! moderator surface
        case('SURFACE')
          call READ_FLUX_SURFACE(it,IU,ilin,LINE,ierr)
        case default
          call READ_TABLE_ROW(IU,ilin,LINE,ierr)
        end select
        !write(*,*) 'READ_FLUX_NEW, flxn=',flxn,' tag=[',trim(LINE(IS:IS+IL-1)),'], IL=',IL
        endif
      enddo
      select case (ierr)
      case(0,i_EOF)
        if (len_trim(MSG)==0) MSG='OK'
      case(i_ERR)
        MSG='Wrong flux table structure. '//trim(MSG)
      case default
        MSG='Error while reading flux table'//trim(MSG)
      end select
      end SUBROUTINE READ_FLUX_NEW

!---------------------------------------------------------------
      SUBROUTINE READ_FLUX_TABLE(FNAME,ilin)
! Read flux table from a file
!---------------------------------------------------------------

!      USE IO
      CHARACTER*(*),intent(in) :: FNAME
      integer,intent(out) :: ilin
      character(MAX_FNAME_LENGTH) :: FRES,FRESPATH
      INTEGER :: IU,IDCOMPARE
      INTEGER :: IERR,IS,IL,NTABS,it
      CHARACTER*512 :: S,MSG,FN
      CHARACTER*1024 :: LINE

      !write(*,*) 'READ_FLUX_TABLE [',trim(FNAME),']'

      IERR=0
      ilin=0
! empty string => clear table and exit => default = Maxwell
      if (LEN_TRIM(FNAME).eq.0) then
        call CLEAR_FLUX_TABLE
        STN=1
        ST(1)%TABNAME='none'
        Return
      endif
! nothing changed
      if (IDCOMPARE(FLUX_TABLE_NAMES,FNAME).eq.0) then
        return
      endif
      IERR=0
      ilin=1 ! to raise report message later
! from now, assume that a new table will be read
      call CLEAR_FLUX_TABLE

! no table => exit with clear table, Maxwell spectrum
      if (IDCOMPARE('none',trim(FNAME)).eq.0) then
        STN=1
        ST(1)%TABNAME='none'
        Return
      endif
! read all tables defined in FNAME
! FNBAME should be a string array of table filenames separated by :
      S=trim(FNAME)
      call COUNTPAR(S,':',NTABS) ! number of tables
      NTABS=MIN(NTABS,STMAX)
      do it=1,NTABS
        call FINDSTRPAR(S,':',it,IS,IL)
        if (IL>0) then
          FN=S(IS:IS+IL-1)
          if (trim(FN).ne.'none') then
! open file
            call OPENRESFILE(trim(FN),' ',0,FRES,FRESPATH,IU,IERR)
            IF(IERR.NE.0) GOTO 100
! assume 1-line header
            call READ_TABLE_ROW(IU,ilin,LINE,ierr)
            if (ierr.eq.i_EOF) goto 30
            if (ierr.ne.0) goto 40
! new table format
            if (trim(LINE).eq.'SIMRES FLUX TABLE') then
              call READ_TABLE_ROW(IU,ilin,LINE,ierr)
              if (ierr.eq.i_EOF) goto 30
              if (ierr.ne.0) goto 40
              call READ_FLUX_NEW(it,IU,LINE,ilin,IERR,MSG)
            else
! old table format
              call READ_FLUX_NEW(it,IU,LINE,ilin,IERR,MSG)
        !write(*,*) 'READ_FLUX_TABLE ierr=',ierr
            endif
            select case (ierr)
            case(0,i_EOF)
              goto 30
            case(i_ERR)
              goto 50
            case default
              goto 40
            end select
! close file, all is OK
30          CLOSE(IU)
! consolidate input
  ! if SPS table is read, wavelength and pulse tables are ignored
            if (ST(it)%PLS_valid) then
              ST(it)%FLXN=0
              ST(it)%FLXNP=0
            endif
          endif
          ST(it)%TABNAME=trim(FN)
          STN=STN+1
        endif
      enddo
      FLUX_TABLE_NAMES=trim(FNAME)
      if ((len_trim(MSG)>0).and.(trim(MSG).ne.'OK')) call MSG_INFO(MSG,1)
      return

! error while reading
4     FORMAT('Error ',I5,' while reading flux table, line ',I5)
40    CLOSE(IU)
      write(S,4) IERR,ilin
      call MSG_ERROR('READ_FLUX_TABLE',S,0,1)
      call CLEAR_FLUX_TABLE
      return

! structure error, e.g. non-equidistant data
5     FORMAT('Error in table structure, ',a,': ',a)
50    CLOSE(IU)
      write(S,5) trim(FNAME),trim(MSG)
      call MSG_ERROR('READ_FLUX_TABLE',S,0,1)
      call CLEAR_FLUX_TABLE
      return

! error on file open
100   call MSG_ERROR('READ_FLUX_TABLE','Cannot open flux table: ['//trim(FNAME)//'].',0,1)
      call CLEAR_FLUX_TABLE
      return

      end SUBROUTINE READ_FLUX_TABLE


!---------------------------------------------------------------
      SUBROUTINE REPORT_FLUX_TABLE(ilin)
! Send info about current source definition
!---------------------------------------------------------------
      use IO
      use XMLINFO
      integer,intent(in) :: ilin
      character*(512) :: MSG
      character*(64) :: CNUM1,CNUM2,CNUM3,CNUM4,CNUM5,CNUM6
      integer :: im,im1,it
      IF (XMLOUT.EQ.0) RETURN

      do it=1,STN

      MSG=' '
      !call XML_RSXDUMP(SMES,' ',1)
      !write(SMES,*) '<INFO priority="low">'
      im= ADD_MESSAGE(m_info,'REPORT_FLUX_TABLE','',m_low)
      im1=0
100   format(a)

      if (trim(ST(it)%TABNAME).ne.'none') then
1       FORMAT('Source table: ',a,', ',I4,' lines. ')
        write(MSG,1) trim(ST(it)%TABNAME),ilin
        !write(SMES,100) trim(MSG)
        call APPEND_TEXT(im,MSG)
        !call MSG_INFO(trim(MSG),0)
      endif

      if (ST(it)%FLXN.gt.1) then
2       FORMAT('Wavelength distribution: range=',a,' .. ',a,' [A], ',a,' rows. ')
        call FLOAT2STR(ST(it)%FLXLAM(1),CNUM1)
        call FLOAT2STR(ST(it)%FLXLAM(ST(it)%FLXN),CNUM2)
        call INT2STR(ST(it)%FLXN,CNUM3)
        write(MSG,2) trim(CNUM1),trim(CNUM2),trim(CNUM3)
      else if (ST(it)%PLS_valid.and.ISIS%VALID) then
21      FORMAT('Spallation source, ISIS table ',a,', flux=',a)
        call FLOAT2STR(ST(it)%PLS_FLUX_INT*1e10,CNUM1)
        write(MSG,21) trim(ISIS%FNAME),trim(CNUM1)
      else if (ST(it)%PLS_valid) then
22      FORMAT('Spallation source,',a,', T=',a,', flux=',a)
        select case(SPS(it)%TYP)
        case(pulse_spc)
          CNUM1='short pulse, coupled'
        case(pulse_spd)
          CNUM1='short pulse, decoupled'
        case(pulse_lp)
          call FLOAT2STR(SPS(it)%D/1.D3,CNUM3)
          CNUM1='long pulse, ver. 2012, d='//trim(CNUM3)//' [ms]'
        case(pulse_ess_cold)
          call FLOAT2STR(SPS(it)%D/1.D3,CNUM3)
          CNUM1='long pulse, para-H2, d='//trim(CNUM3)//' [ms]'
        case(pulse_ess_therm)
          call FLOAT2STR(SPS(it)%D/1.D3,CNUM3)
          CNUM1='long pulse, thermal, d='//trim(CNUM3)//' [ms]'
        case(pulse_ess_2015T)
          call FLOAT2STR(SPS(it)%D/1.D3,CNUM3)
          CNUM1='ESS 2015, thermal, d='//trim(CNUM3)//' [ms]'
        case(pulse_ess_2015C)
          call FLOAT2STR(SPS(it)%D/1.D3,CNUM3)
          CNUM1='ESS 2015, cold, d='//trim(CNUM3)//' [ms]'
        case default
          call INT2STR(SPS(it)%TYP,CNUM3)
          CNUM1='unknown type '//trim(CNUM3)
        end select
        call FLOAT2STR(SPS(it)%TEMP,CNUM2)
        call FLOAT2STR(ST(it)%PLS_FLUX_INT*1e10,CNUM3)
        write(MSG,22) trim(CNUM1),trim(CNUM2),trim(CNUM3)
      else
        write(MSG,100) 'No flux table, using Maxwell spectrum. '
      endif
      ! write(SMES,100) trim(MSG)
      call APPEND_TEXT(im,MSG)
      !call MSG_INFO(trim(MSG),0)

      if (ST(it)%FLXNR.gt.1) then
3       FORMAT('Radial distribution: range=',a,' .. ',a,' [mm], ',a,' rows. ')
        call FLOAT2STR(ST(it)%FLXRAD(1),CNUM1)
        call FLOAT2STR(ST(it)%FLXRAD(ST(it)%FLXNR),CNUM2)
        call INT2STR(ST(it)%FLXNR,CNUM3)
        write(MSG,3) trim(CNUM1),trim(CNUM2),trim(CNUM3)
        if (im1==0) im1=ADD_MESSAGE(m_info,'REPORT_FLUX_TABLE','',m_low)
        call APPEND_TEXT(im1,MSG)
        ! write(SMES,100) trim(MSG)
        !call MSG_INFO(trim(MSG),0)
      endif

      if (ST(it)%FLXHNX.gt.1) then
4       FORMAT('Horizontal(',a,'x',a,'), range = (',a,',',a,')+-(',a,',',a,') [mm,rad]. ')
        call INT2STR(ST(it)%FLXHNX,CNUM1)
        call INT2STR(ST(it)%FLXHNA,CNUM2)
        call FLOAT2STR(ST(it)%FLXHX0,CNUM3)
        call FLOAT2STR(ST(it)%FLXHA0,CNUM4)
        call FLOAT2STR(ST(it)%FLXHX,CNUM5)
        call FLOAT2STR(ST(it)%FLXHA,CNUM6)
        write(MSG,4) trim(CNUM1),trim(CNUM2),trim(CNUM3),trim(CNUM4),trim(CNUM5),trim(CNUM6)
        if (im1==0) im1=ADD_MESSAGE(m_info,'REPORT_FLUX_TABLE','',m_low)
        call APPEND_TEXT(im1,MSG)
        ! write(SMES,100) trim(MSG)
        !call MSG_INFO(trim(MSG),0)
      endif
      if (ST(it)%FLXVNX.gt.1) then
5       FORMAT('Vertical(',a,'x',a,'), range = (',a,',',a,')+-(',a,',',a,') [mm,rad]. ')
        call INT2STR(ST(it)%FLXVNX,CNUM1)
        call INT2STR(ST(it)%FLXVNA,CNUM2)
        call FLOAT2STR(ST(it)%FLXVX0,CNUM3)
        call FLOAT2STR(ST(it)%FLXVA0,CNUM4)
        call FLOAT2STR(ST(it)%FLXVX,CNUM5)
        call FLOAT2STR(ST(it)%FLXVA,CNUM6)
        write(MSG,5) trim(CNUM1),trim(CNUM2),trim(CNUM3),trim(CNUM4),trim(CNUM5),trim(CNUM6)
        if (im1==0) im1=ADD_MESSAGE(m_info,'REPORT_FLUX_TABLE','',m_low)
        call APPEND_TEXT(im1,MSG)
        !write(SMES,100) trim(MSG)
        !call MSG_INFO(trim(MSG),0)
      endif

      if (ST(it)%FLXNP.gt.1) then
6       FORMAT('Pulse time distribution: range=',a,' .. ',a,' [us], ',a,' rows. ')
        call FLOAT2STR(ST(it)%FLXTIME(1),CNUM1)
        call FLOAT2STR(ST(it)%FLXTIME(ST(it)%FLXNP),CNUM2)
        call INT2STR(ST(it)%FLXNP,CNUM3)
        write(MSG,6) trim(CNUM1),trim(CNUM2),trim(CNUM3)
        if (im1==0) im1=ADD_MESSAGE(m_info,'REPORT_FLUX_TABLE','',m_low)
        call APPEND_TEXT(im1,MSG)
        !write(SMES,100) trim(MSG)
        !call MSG_INFO(trim(MSG),0)
      endif
      if (ST(it)%SRF_NS.gt.0) then
7       FORMAT('Polygonal surface, ',a,' elements. ')
        call INT2STR(ST(it)%SRF_NS,CNUM1)
        write(MSG,7) trim(CNUM1)
        if (im1==0) im1=ADD_MESSAGE(m_info,'REPORT_FLUX_TABLE','',m_low)
        call APPEND_TEXT(im1,MSG)
        !write(SMES,100) trim(MSG)
        !call MSG_INFO(trim(MSG),0)
      endif
      !write(SMES,*)'</INFO>'
      !call XML_RSXDUMP(SMES,' ',0)

      enddo

      end SUBROUTINE REPORT_FLUX_TABLE

      end MODULE SOURCES_TABLE
