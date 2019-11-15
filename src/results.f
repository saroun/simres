!//////////////////////////////////////////////////////////////////////
!////  $Id: results.f,v 1.55 2019/08/15 15:02:06 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2008, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.55 $
!////     $Date: 2019/08/15 15:02:06 $
!//////////////////////////////////////////////////////////////////////
!////
!////  data structures with basic results of SIMRES
!////  some basic procedures for data collection
!//////////////////////////////////////////////////////////////////////
      MODULE RESULTS
      use SIMATH
      use TRACINGDATA
      use NSTORE
      use COMMANDS
      use SPECTROMETER
      implicit none
      SAVE

    ! these constants should be consistent with those in GRFDATA !!
      integer,parameter :: SPCM=4096 ! max dimension of 1D spectra
      integer,parameter :: MAX2DIM=512 ! max dimension of 2D spectra
      integer,parameter :: MAXPKS=64 ! max number of peaks in one spectrum
      integer,parameter :: MAXSCAN1D=2048 ! max steps of a scan
      integer,parameter :: MAXSCAN2D=256 ! max steps of a scan
      integer,parameter :: MAXVARS=10 ! max variables in a 1D scan

  ! ID of phase-space variables for graphics output
      integer,parameter :: RK_MAX=9
      character(55),parameter :: RK_NAMES='X:Y:Z:\gDk\dx\u/k:\gDk\dy\u/k:\gDk\dz\u/k:time:\gDE:\gl'
      character(24),parameter :: RK_UNITS='mm:mm:mm: : : :us:meV:\A'

  ! ID of resolution variables for graphics output
      integer,parameter :: QE_MAX=6
      character(36),parameter :: QE_NAMES='\gdQ/Q:theta:EN:Q\dx\u:Q\dy\u:Q\dz\u'
      character(35),parameter :: QE_UNITS=':deg:meV:\A\u-1\d:\A\u-1\d:\A\u-1\d'

  ! ID of detector variables for graphics output
      integer,parameter :: DT_MAX=4
      character(27),parameter :: DT_NAMES='theta:d\dhkl\u:d\dhkl\u:TOF'
      character(10),parameter :: DT_UNITS='deg:A:A:us'

! type of results set
      integer,parameter :: rset_basic=0
      integer,parameter :: rset_bpar=1
      integer,parameter :: rset_rpar=2
      integer,parameter :: rset_dpar=3
      integer,parameter :: rset_fmerit=4
      integer,parameter :: rset_bpeak=5
      integer,parameter :: rset_rpeak=6
      integer,parameter :: rset_dpeak=7
! type of output data (beam monitor, detector data, or gauge map)
      integer, parameter :: show_monitor=0
      integer, parameter :: show_detector=1
      integer, parameter :: show_gauge=2
! type of FM formula
      integer,parameter :: fmerit_int=0    ! intensity
      integer,parameter :: fmerit_w=1      ! 1/width
      integer,parameter :: fmerit_intw=2   ! intensity/width
      integer,parameter :: fmerit_intw2=3  ! intensity/width^2


! data structure for evaluated variable
! its meaning depends on rset_xxx value, can be IT,dINT,dE, or e.g. peak parameters
! HDR defines evaluated quatity, like intensity, peak width etc.
! IV defines the independent variable of intensity profile evaluated, e.g. x-coordinate or wavelength
      type TEVAL ;SEQUENCE
        character(LEN_LINE) :: HDR  ! header string (| separated strings), e.g. INT|dE|TOF
        character(2*LEN_NAME) :: CAP  ! caption
        INTEGER :: N        ! number of items defined in HDR
        integer :: TYP      ! type of the evaluated value (=rset_xxx)
        INTEGER :: IV       ! e.g. for beam profiles: index of phase-space variable
        REAL :: V0,DV       ! e.g. for beam profiles: center and range of phase-space variable
      end type TEVAL


! Array for simulated 1-dim spectra
! N .. nuber of points
! NPK .. number of peaks
! PKPOS(NP) ... peak position (estimates)
      type TRES1D ;SEQUENCE
        INTEGER :: N,NPK
        logical :: APPEND
        integer :: FILT
        REAL :: X(SPCM),Y(SPCM),DY(SPCM)
        REAL :: PKPOS(MAXPKS),PKINT(MAXPKS)
        integer :: PKMIN,PKMAX
      end type TRES1D
      type(TRES1D) :: RES1D
      integer, parameter :: RES1D_BINDIM=90
      real :: RES1D_BINS(SPCM,0:RES1D_BINDIM)
      real :: RES1D_BINS_STEP
      character(256) :: RES1D_BINS_COM
      integer :: RES1D_N
      integer :: RES1D_BINMAX
      integer :: RES1D_BINMIN

! Array for simulated 2-dim spectra
      type TRES2D ;SEQUENCE
        INTEGER :: NX,NY
        logical :: APPEND
        integer :: FILT
        REAL :: DX,DY
        REAL :: Z(MAX2DIM,MAX2DIM)
      end type TRES2D
      type(TRES2D) :: RES2D

! 1D scan results, permits multiple correlated variables and multiple values per step
      type TSCAN1D ;SEQUENCE
        INTEGER :: NVAR,NS  ! numer of variables and steps
        integer :: IX   ! index of X variables to plot
        integer :: IV   ! value index to plot
        integer :: MVAR, MSCAN, MPNT ! allocated array dimensions
        character(LEN_LINE) :: HDRVAR ! | separated list of labels for step variables
        TYPE(TEVAL) :: EVAL  ! definition of evaluated quatity
        real, allocatable :: VARIABLES(:,:) ! alocate to  MVAR:MSCAN
        REAL, allocatable :: VALUES(:,:) ! alocate to  MPNT:MSCAN
        REAL, allocatable :: ERRORS(:,:) ! alocate to MPNT:MSCAN
      end type TSCAN1D
      type(TSCAN1D) :: SCAN1D

! 2D scan results, permits 2 indepenent variables and a single value per step
      type TSCAN2D ;SEQUENCE
        INTEGER :: NX,NY    ! number of steps in each direction
        INTEGER :: IV   ! value index (depends on TYPE)
        character(LEN_LINE) :: HXVAR,HYVAR ! step variables
        real :: XMIN,XMAX  ! range of the X-variable
        real :: YMIN,YMAX  ! range of the Y-variable
        TYPE(TEVAL) :: EVAL  ! definition of evaluated quatity
        REAL :: VALUES(MAXSCAN2D,MAXSCAN2D)
      end type TSCAN2D
      type(TSCAN2D) :: SCAN2D

! bunch of peaks parameters
      type TPKBUNCH ;SEQUENCE
        INTEGER :: NP    ! number of peaks
        INTEGER :: REF(MAXPKS)   ! peak indices
        real :: XMIN,XMAX  ! range of the X-variable
        real :: YMIN,YMAX  ! range of the Y-variable
        REAL :: IMAX ! maximum intensity in the bunch
        REAL :: MAXPOS ! position of the bin with maximum intensity
        REAL :: ISUM ! sum of intensities in the bunch
        REAL :: XSCALE
      end type TPKBUNCH

! Array with basic results (intensity etc.) for each dataset and subset
      real(KIND(1.D0)) :: MCRESULT(6,MSETS,MDATS)
      real(KIND(1.D0)) :: MCNORM(2,MSETS,MDATS)
      real(KIND(1.D0)) :: MCTRIALS(MSETS,MDATS)
      real(KIND(1.D0)) :: MCVOLI(2,MSETS,MDATS)
      real(KIND(1.D0)) :: MCVOLF(2,MSETS,MDATS)

      contains



!-------------------------------------------------------------------
      SUBROUTINE DEALLOC_SCAN1D(ierr)
!-------------------------------------------------------------------
      integer, intent(out) :: ierr
      ierr=0
      SCAN1D%MVAR=0
      SCAN1D%MSCAN=0
      SCAN1D%MPNT=0
      if (allocated(SCAN1D%VARIABLES)) deallocate(SCAN1D%VARIABLES,STAT=ierr)
      if (allocated(SCAN1D%VALUES)) deallocate(SCAN1D%VALUES,STAT=ierr)
      if (allocated(SCAN1D%ERRORS)) deallocate(SCAN1D%ERRORS,STAT=ierr)
      end SUBROUTINE DEALLOC_SCAN1D

!-------------------------------------------------------------------
      SUBROUTINE ALLOC_SCAN1D(NVAR, NSCAN, NP, ierr)
!-------------------------------------------------------------------
      integer, intent(in) :: NVAR, NSCAN, NP
      integer, intent(out) :: ierr
      call DEALLOC_SCAN1D(ierr)
      if ((ierr==0).and.NVAR*NSCAN>0) then
        allocate(SCAN1D%VARIABLES(1:NVAR,1:NSCAN),STAT=ierr)
        if (ierr==0) then
          allocate(SCAN1D%VALUES(1:NP,1:NSCAN),STAT=ierr)
        endif
        if (ierr==0) then
          allocate(SCAN1D%ERRORS(1:NP,1:NSCAN),STAT=ierr)
        endif
        if (ierr==0) then
          SCAN1D%MVAR=NVAR
          SCAN1D%MSCAN=NSCAN
          SCAN1D%MPNT=NP
        else
          call DEALLOC_SCAN1D(ierr)
          call MSG_ERROR('ALLOC_SCAN1D','Cannot allocate SCAN1D',1,1)
        endif
      endif
      end SUBROUTINE ALLOC_SCAN1D


!-------------------------------------------------------------------
      SUBROUTINE GET_EVAL_RANGE(E,XLIM)
! return variable range for given EVAL structure as double array
!-------------------------------------------------------------------
      TYPE(TEVAL),intent(in) :: E
      real(KIND(1.D0)),intent(out) :: XLIM(2)
      XLIM = (/1.D0*E%V0 - 0.5D0*E%DV, 1.D0*E%V0 + 0.5D0*E%DV/)
      end SUBROUTINE GET_EVAL_RANGE

!-------------------------------------------------------------------
      SUBROUTINE VAR_GETNAME(IVAR,UNITS,VARNAME)
! Get the name and (if UNITS=true) physical units of a variable as a string.
! The names are in PGPLOT format for showing as graph labels.
! IVAR .. index of the variabe in RK_NAMES string array
!-------------------------------------------------------------------
      integer,intent(in) :: IVAR
      logical,intent(in) :: UNITS
      character(*) :: VARNAME
      integer :: IS,IL,IS1,IL1
      VARNAME='undefined'
      if ((IVAR>0).and.(IVAR<=RK_MAX)) then
        call FINDSTRPAR(RK_NAMES,':',IVAR,IS,IL)
        call FINDSTRPAR(RK_UNITS,':',IVAR,IS1,IL1)
        if (UNITS.and.(IL1>0)) then
          VARNAME=RK_NAMES(IS:IS+IL-1)//' ['//RK_UNITS(IS1:IS1+IL1-1)//']'
        else
          VARNAME=RK_NAMES(IS:IS+IL-1)
        endif
      endif
      end SUBROUTINE VAR_GETNAME

!-------------------------------------------------------------------
      SUBROUTINE QVAR_GETNAME(IVAR,UNITS,VARNAME)
! Get the name and (if UNITS=true) physical units of a variable as a string.
! The names are in PGPLOT format for showing as graph labels.
! IVAR .. index of the variabe in QE_NAMES string array
!-------------------------------------------------------------------
      integer,intent(in) :: IVAR
      logical,intent(in) :: UNITS
      character(*) :: VARNAME
      integer :: IS,IL,IS1,IL1
      VARNAME='undefined'
      if ((IVAR>0).and.(IVAR<=QE_MAX)) then
        call FINDSTRPAR(QE_NAMES,':',IVAR,IS,IL)
        call FINDSTRPAR(QE_UNITS,':',IVAR,IS1,IL1)
        if (UNITS.and.(IL1>0)) then
          VARNAME=QE_NAMES(IS:IS+IL-1)//' ['//QE_UNITS(IS1:IS1+IL1-1)//']'
        else
          VARNAME=QE_NAMES(IS:IS+IL-1)
        endif
      endif
      end SUBROUTINE QVAR_GETNAME


!-------------------------------------------------------------------
      SUBROUTINE DVAR_GETNAME(IVAR,UNITS,VARNAME)
! Get the name and (if UNITS=true) physical units of a variable as a string.
! The names are in PGPLOT format for showing as graph labels.
! IVAR .. index of the variabe in DT_NAMES string array
!-------------------------------------------------------------------
      integer,intent(in) :: IVAR
      logical,intent(in) :: UNITS
      character(*) :: VARNAME
      integer :: IS,IL,IS1,IL1
      VARNAME='undefined'
      if ((IVAR>0).and.(IVAR<=DT_MAX)) then
        call FINDSTRPAR(DT_NAMES,':',IVAR,IS,IL)
        call FINDSTRPAR(DT_UNITS,':',IVAR,IS1,IL1)
        if (UNITS.and.(IL1>0)) then
          VARNAME=DT_NAMES(IS:IS+IL-1)//' ['//DT_UNITS(IS1:IS1+IL1-1)//']'
        else
          VARNAME=DT_NAMES(IS:IS+IL-1)
        endif
      endif
      end SUBROUTINE DVAR_GETNAME


!-------------------------------------------------------------------
      SUBROUTINE CLEAR_RES1D(RES)
! clear data in RES
!-------------------------------------------------------------------
      TYPE(TRES1D) :: RES
      INTEGER :: I
      DO I=1,SPCM
        RES%X(I)=0.D0
        RES%Y(I)=0.D0
        RES%DY(I)=0.D0
        RES1D_BINS(I,0:RES1D_BINDIM)=0.0
      enddo
      RES1D_BINMAX=0
      RES1D_BINMIN=RES1D_BINDIM+1
      RES1D_N=0
      RES%N=0
      RES%NPK=0
      RES%PKPOS=0.0
      end SUBROUTINE CLEAR_RES1D

!-------------------------------------------------------------------
      SUBROUTINE CLEAR_RES2D(RES)
! clear data in RES
!-------------------------------------------------------------------
      TYPE(TRES2D) :: RES
      RES%DX=0.D0
      RES%DY=0.D0
      RES%Z=0.0
      end SUBROUTINE CLEAR_RES2D

C-------------------------------------------------------------------
      REAL(KIND(1.D0)) function PEAKFWHM(RES)
C Get true FWHM of a peak stored in RES
C-------------------------------------------------------------------
      TYPE(TRES1D),intent(in) :: RES
      INTEGER*4 I,IMAX,I1,I2
      REAL*4 z1,z2
      REAL*8 YMAX,X1,X2
      YMAX=RES%Y(1)
      IMAX=1
      DO I=1,RES%N
        IF(YMAX.LE.RES%Y(I)) THEN
          IMAX=I
          YMAX=RES%Y(I)
        ENDIF
      ENDDO
      I1=1
      DO WHILE ((RES%Y(I1).LT.YMAX/2.D0).AND.(I1.LT.RES%N))
        I1=I1+1
      ENDDO
      I2=RES%N
      DO WHILE ((RES%Y(I2).LT.YMAX/2.D0).AND.(I2.GT.0))
        I2=I2-1
      ENDDO
      IF((I1.GT.1).AND.(I1.LT.RES%N).AND.(I2.GT.1).AND.(I2.LT.RES%N).AND.(I1.LE.I2)) then
        z1= RES%Y(I1)-RES%Y(I1-1)
        z2=RES%Y(I2+1)-RES%Y(I2)
        if (z1.ne.0.and.z2.ne.0) then
          X1=RES%X(I1-1)+(YMAX/2.D0-RES%Y(I1-1))/(RES%Y(I1)-RES%Y(I1-1))*(RES%X(I1)-RES%X(I1-1))
          X2=RES%X(I2)+(YMAX/2.D0-RES%Y(I2))/(RES%Y(I2+1)-RES%Y(I2))*(RES%X(I2+1)-RES%X(I2))
          PEAKFWHM=x2-x1
        else
          PEAKFWHM=0.
        endif
      ELSE
       PEAKFWHM=0
      ENDIF
      END function PEAKFWHM

C-------------------------------------------------------------------
      SUBROUTINE PEAKPARAM(RES,suma,center,fwhm,wspread)
C Get parameters of a peak stored in RES
C-------------------------------------------------------------------
      TYPE(TRES1D),intent(in) :: RES
      REAL(KIND(1.D0)),intent(out) :: suma,center,fwhm,wspread
      REAL(KIND(1.D0)) :: S0,S1,S2,Z
      INTEGER*4 I
      suma=0.D0
      center=0.D0
      wspread=0.D0
      fwhm=0.D0
      suma=0.
      S0=0
      S1=0
      S2=0
      DO I=1,RES%N
        Z=RES%X(I)
        S0=S0+RES%Y(I)
        S1=S1+RES%Y(I)*Z
        S2=S2+RES%Y(I)*Z**2
        suma=suma+RES%Y(I)
      ENDDO
      IF(S0.GT.0) then
        center=S1/S0
        wspread=R8LN2*SQRT(ABS(S2/S0-center**2))
        fwhm=PEAKFWHM(RES)
      ENDIF
      end SUBROUTINE PEAKPARAM


!-------------------------------------------------------------------
      SUBROUTINE MPEAK_GETBUNCHES(RES,BS,NBMAX,NB)
! Scan data stored in RES and return parameters of all peak bunches found
! Initial positions must already be set in RES%PKPOS
! A bunch of peaks is defined by a continuous range of non-zero RES%Y values
! BS .. array of TPKBUNCH data structures
! NBMAX .. max. number of bunches
! NB .. actual number of bunches
!-------------------------------------------------------------------
      TYPE(TRES1D) :: RES
      TYPE(TPKBUNCH) :: BS(:)
      integer,intent(in) :: NBMAX
      integer,intent(out) :: NB
      integer :: i,n,nmax,j1,j2,k,NVAL,NVAL0,NPK,PKMIN,PKMAX
      REAL :: S0,S1,DX,XSTEP,LIM1,LIM2,SCA
      logical :: dbg=.false.
1     format(a,2x,6(G13.5,2x))

      XSTEP=(RES%X(RES%N)-RES%X(1))/(RES%N-1)
    ! clean records
      do i=1,NBMAX
        BS(i)%NP=0
        BS(i)%XMIN=RES%Y(RES%N)
        BS(i)%XMAX=RES%Y(1)
        BS(i)%YMIN=1.E10
        BS(i)%YMAX=0.0
        BS(i)%IMAX=0.0
        BS(i)%MAXPOS=0.0
        BS(i)%ISUM=0.0
        BS(i)%XSCALE=1.0
      enddo
      i=RES%N
      n=0
      DO while (i>0)
      !dbg=(n>17)
      ! skip a range with no data
        j1=0
        j2=0
        S0=0.D0
        S1=0.D0
        nmax=0
      ! find start and stop positions of the next bunch
      ! the start is defined by at least 2 consecutive RES%Y=0 values or i=RES%N
      ! a bunch must be a sequence of at least 2 nonzero points surrounded byat least 2 zeros
        do while ((i>=0).and.((j1==0).or.(j2==0)))
          ! scan a sequence of zeros
          NVAL0=0
          do while ((i>0).and.(RES%Y(i)==0.D0))
            i=i-1
            NVAL0=NVAL0+1
          enddo
          if ((NVAL0>1).or.(i>RES%N-2).or.(i==0)) then
            if (j1==0) then
              j1=min(i+1,RES%N)
            !else if (NVAL>1) then
            else
              j2=max(i+NVAL0,1)
              i=i+NVAL0
          !  else
          !    ! skip isolated non-zero bins
          !    j1=0
          !    j2=0
            endif
          endif
          ! scan a sequence of non-zeros
          if (j2==0) then
            NVAL=0
            do while ((i>0).and.(RES%Y(i)>0.D0))
              i=i-1
              NVAL=NVAL+1
            enddo
          endif
        enddo
        if (j1*j2==0) EXIT
        if (j2+2<j1) then
          S0=0.D0
          S1=0.D0
          nmax=0
          do k=j2,j1
            if (S0<RES%Y(k)) then
              S0=RES%Y(k)
              nmax=k
            endif
            S1=S1+RES%Y(k)
          enddo
        endif
        ! define new bunch parameters
        if ((n<NBMAX).and.(j2+2<j1).and.(S1>0.0)) then
          NVAL=0
          n=n+1
        ! set J2 and J1 as the nearest channels with zero intensity
          BS(n)%XMAX=RES%X(J1)
          BS(n)%XMIN=RES%X(J2)
          BS(n)%IMAX=S0
          BS(n)%MAXPOS=RES%X(nmax)
          BS(n)%ISUM=S1
          if (dbg) then
            if (n>0) write(*,1) '  BUNCH ',n,BS(n)%XMIN,BS(n)%XMAX,BS(n)%MAXPOS,BS(n)%IMAX,BS(n)%ISUM
          endif
        endif
      enddo
      NB=n
      ! default width value
      DX=2*XSTEP
      SCA=1.0
      PKMIN=MAXPKS
      PKMAX=0
      do I=1,NB
        S0=0.D0
        S1=0.D0
      ! search boundary is in the middle between the bunch limits
        if (NB==1) then
          !LIM1=RES%X(1)-DX
          !LIM2=RES%X(RES%N)+DX
          LIM1=BS(i)%XMIN-DX
          LIM2=BS(i)%XMAX+DX
        else if (I==1) then
          LIM1=BS(i)%XMIN
          !LIM2=RES%X(RES%N)+DX
          LIM2=BS(i)%XMAX+DX
        else if (I==NB) then
          !LIM1=RES%X(1)-DX
          LIM1=BS(i)%XMIN-DX
          LIM2=BS(i)%XMAX
        else
          LIM1=BS(i)%XMIN
          LIM2=BS(i)%XMAX
        endif
      ! search for peaks within the bunch limits
      ! apply peak shift factor SCA on the limits !!
        NPK=0
        DO k=1,RES%NPK
          if ((RES%PKPOS(k)*SCA>LIM1) .and. (RES%PKPOS(k)*SCA<=LIM2)) then
            NPK=NPK+1
            BS(i)%REF(NPK)=k
            if (dbg) write(*,1) '  peak: ',I,k,RES%PKPOS(k)*SCA
          endif
        enddo
        if (NPK>0) then
          BS(i)%NP=NPK
          if (NPK==1) then
            SCA=1.0+(BS(i)%MAXPOS-RES%PKPOS(BS(i)%REF(1)))/RES%PKPOS(BS(i)%REF(1))
            !if (dbg) write(*,1) '   x-scale: ',BS(i)%REF(1),BS(i)%MAXPOS,RES%PKPOS(BS(i)%REF(1)),SCA
          endif
          BS(i)%XSCALE=SCA
          PKMIN=min(PKMIN,BS(i)%REF(1),BS(i)%REF(NPK))
          PKMAX=max(PKMAX,BS(i)%REF(1),BS(i)%REF(NPK))
        endif
      enddo
      RES%PKMIN=PKMIN
      RES%PKMAX=PKMAX
      end SUBROUTINE MPEAK_GETBUNCHES

!---------------------------------------------------------------------------
      SUBROUTINE MPEAK_GETPARAM(BUNCH,INTS,PAR)
! Get estimated parameters of multiple peaks in a bunch
! INTS array contains nominal peak intensities for all peaks in the data set
! PAR(1) ... x-scale factor
! PAR(2*n) .. intensity of the n-th peak
! PARR(2*n+1) .. width of the n-th peak in steps
!----------------------------------------------------------------------------
      TYPE(TPKBUNCH) :: BUNCH
      REAL,intent(in) :: INTS(:)
      REAL,intent(out) :: PAR(:)
      REAL :: FAC,DX,S0,S1
      integer :: k,m
      S0=0.0
      S1=0.0
      DO k=1,BUNCH%NP
        S0=MAX(S0,INTS(BUNCH%REF(k)))
        S1=S1+INTS(BUNCH%REF(k))
      enddo
      if (S0<=0.0) return
    ! intensity scale factor = ratio Imax/F2max
      FAC=BUNCH%IMAX/S0
    ! estimated peak width in steps =
      DX=BUNCH%ISUM/(S1*FAC)/sqrt(TWOPI)*R8LN2
      do k=1,BUNCH%NP
        m=2*(k-1)+1
        PAR(m+1)=BUNCH%ISUM/DX*INTS(BUNCH%REF(k))/S1
        PAR(m+2)=DX
        !write(*,*) 'MPEAK_GETPARAM ',BUNCH%REF(k),PAR(m+1),PAR(m+2)
      enddo
      PAR(1)=BUNCH%XSCALE
      end SUBROUTINE MPEAK_GETPARAM

      !---------------------------------------------------------------------------
      SUBROUTINE MPEAK_GETPARAM2(BUNCH,INTS,PAR)
! Get estimated parameters of multiple peaks in a bunch
! INTS array contains nominal peak intensities for all peaks in the data set
! PAR(1) ... x-scale factor
! PAR(2) .. intensity of the n-th peak
! PAR(n+2) .. width of the n-th peak in steps
!----------------------------------------------------------------------------
      TYPE(TPKBUNCH) :: BUNCH
      REAL,intent(in) :: INTS(:)
      REAL,intent(out) :: PAR(:)
      REAL :: FAC,DX,S0,S1
      integer :: k
      S0=0.0
      S1=0.0
      DO k=1,BUNCH%NP
        S0=MAX(S0,INTS(BUNCH%REF(k)))
        S1=S1+INTS(BUNCH%REF(k))
      enddo
      if (S0<=0.0) return
    ! intensity scale factor = ratio Imax/F2max
      FAC=BUNCH%IMAX/S0
    ! estimated peak width in steps =
      DX=BUNCH%ISUM/(S1*FAC)/sqrt(TWOPI)*R8LN2
      do k=1,BUNCH%NP
        PAR(k+2)=BUNCH%ISUM/DX*INTS(BUNCH%REF(k))/S1
      enddo
      PAR(1)=BUNCH%XSCALE
      PAR(2)=DX
      end SUBROUTINE MPEAK_GETPARAM2

!-------------------------------------------------------------------
      SUBROUTINE MPEAKPARAM(RES,PAR,PARX)
! Get estimated parameters of multiple peaks
! Initial positions must already be set in RES%PKPOS
!-------------------------------------------------------------------
      !USE OPTIMIZATION
      TYPE(TRES1D) :: RES
      REAL,intent(out) :: PAR(:),PARX
      REAL(KIND(1.D0)) :: DX,XSTEP,FAC,SCA
      INTEGER :: I,J,K,M,IX,NPK,NPK0,PKMIN,PKMAX,NVAL
      INTEGER :: nmax,n1,n2,j1,j2,iref(MAXPKS),NB
      REAL(KIND(1.D0)) :: PBIMAX(1024),PBKMAX(1024),PBSUM(1024)
      REAL(KIND(1.D0)) :: PBMI(1024),PBMA(1024),LIM1,LIM2
      REAL :: S0,S1,S2,I0,Z,dZ,NCN
      logical :: dbg=.false.
1     format(a,2x,6(G13.5,2x))
      PAR=0.0
      PARX=1.0
      if (RES%N<2) return
    ! get peak positions and amplitudes within data range
      XSTEP=(RES%X(RES%N)-RES%X(1))/(RES%N-1)
      S0=0.0
      dZ=0.0
      Z=0.0
    ! sum over bins
      DO I=1,RES%N
        S0=S0+RES%Y(I)
        Z=max(Z,RES%Y(I))
        if (Z==RES%Y(I)) then
          dZ=RES%DY(I)
        endif
      enddo
      if (S0<=0.D0) return
    ! redefine error bars so that they correspond to counting statistics
      if (dbg) write(*,1) 'MPEAKPARAM eff. counts :',Z,dZ
      NCN=(Z/dZ)**2
      FAC=NCN/Z
      if (dbg) write(*,1) 'MPEAKPARAM eff. counts :',NCN
      do I=1,RES%N
        if (RES%Y(I)>0.D0) then
          n1=NINT(RES%Y(I)*FAC)+1
          RES%DY(I)=sqrt(1.D0*max(1,n1))/FAC
        endif
      enddo

    ! average peak intensity
      S1=0.0
      N1=0
      do J=1,RES%NPK
        K=3*(J-1)
        IX=NINT((RES%PKPOS(J)-RES%X(1))/XSTEP)
        if ((IX>0).and.(IX<=RES%N)) then
          S1=S1+RES%Y(IX)
          N1=N1+1
        endif
      enddo
      if (N1==0) return
      I0=S1/N1
      if (dbg) write(*,1) 'MPEAKPARAM I0:',I0,N1,XSTEP



    ! find all isolated bunches of peaks in the data range
      I=RES%N
      NB=0
      N2=0
      S2=0.0
      NVAL=0
      SCA=1.D0
      DO while (i>0)
      ! skip a range with no data
        do while ((I>0).and.(RES%Y(I)==0.D0))
          I=I-1
        enddo
        if (I<=0) exit
      ! the first nonzero point of the current bunch
        J1=I
        S0=0.D0
        S1=0.D0
        nmax=0
        do while ((I>0).and.(RES%Y(I)>0.D0))
          if (S0<RES%Y(I)) then
            S0=RES%Y(I)
            nmax=I
          endif
          S1=S1+RES%Y(I)
          I=I-1
          NVAL=NVAL+1
        enddo
        ! any bunch must be a sequence of at least 3 nonzero points
        if (NVAL>2) then
        NVAL=0
        J2=i+1 ! last non-zero channel
        NB=NB+1
        PBMA(NB)=RES%X(J1)
        PBMI(NB)=RES%X(J2)
        PBIMAX(NB)=S0
        PBKMAX(NB)=RES%X(nmax)
        PBSUM(NB)=S1
        if (dbg) then
          if (NB>0) write(*,1) '  BUNCH ',NB,PBMI(NB),PBMA(NB),PBKMAX(NB),PBIMAX(NB),PBSUM(NB)
        endif
        endif
      enddo
      !write(*,1) '  BUNCH ',NB,PB(NB),PBKMAX(NB),PBIMAX(NB),PBSUM(NB)
      ! for each bunch, find maximum intensity and assign peak parameters
      NPK=0
      PKMIN=RES%NPK
      PKMAX=0
    ! default width value
      DX=5*XSTEP
      do I=1,NB
        S0=0.D0
        S1=0.D0
      ! search boundary is in the middle between the bunch limits
        if (NB==1) then
          LIM1=RES%X(1)
          LIM2=RES%X(RES%N)
        else if (I==1) then
          LIM1=(PBMI(I)+PBMA(I+1))/2.D0
          LIM2=RES%X(RES%N)
        else if (I==NB) then
          LIM1=RES%X(1)
          LIM2=(PBMA(I)+PBMI(I-1))/2.D0
        else
          LIM1=(PBMI(I)+PBMA(I+1))/2.D0
          LIM2=(PBMA(I)+PBMI(I-1))/2.D0
        endif
      ! search for peaks within the bunch limits
      ! apply peak shift factor SCA on the limits !!
        NPK0=NPK
        DO k=1,RES%NPK
          if ((RES%PKPOS(k)*SCA>LIM1-DX) .and. (RES%PKPOS(k)*SCA<LIM2+DX)) then
            NPK=NPK+1
            IREF(NPK)=k
            S0=MAX(S0,RES%PKINT(k))
            S1=S1+RES%PKINT(k)
            PKMIN=MIN(PKMIN,k)
            PKMAX=MAX(PKMAX,k)
          endif
        enddo
        N1=NPK-NPK0
        if (N1>0) then
          FAC=PBIMAX(I)/S0
          if (dbg) write(*,1) 'peaks in bunch ',i,N1,LIM2,LIM1,PBSUM(i)
      ! estimated peak width
          DX=PBSUM(I)*XSTEP/(S1*FAC)/sqrt(TWOPI)*R8LN2
          DO k=NPK0+1,NPK
            m=3*(k-1)
            PAR(m+1)=FAC*RES%PKINT(IREF(k))
            PAR(m+2)=RES%PKPOS(IREF(k))
            PAR(m+3)=DX
            if (dbg) write(*,1) '   k=',k,iref(k),PAR(m+1),PAR(m+2),PAR(m+3)
          enddo
        endif
      ! if there is only 1 peak in the bunch, try to estimate correcting factor for position
        if (N1==1) then
          N2=N2+1
          SCA=1.0+(PBKMAX(I)-RES%PKPOS(IREF(NPK)))/RES%PKPOS(IREF(NPK))
          S2=S2+SCA
          if (dbg) write(*,1) '   x-scale: ',IREF(NPK),PBKMAX(I),RES%PKPOS(IREF(NPK)),SCA
        endif
      enddo
      if (N2>0) PARX=S2/N2
      RES%PKMIN=PKMIN
      RES%PKMAX=PKMAX

      do i=RES%PKMIN,RES%PKMAX
        k=3*(i-RES%PKMIN)
        PAR(k+2)=PAR(k+2)*PARX
        !write(*,1) 'MPEAKPARAM:',i,PAR(k+1),PAR(k+2),PAR(k+3)
      enddo

      end SUBROUTINE MPEAKPARAM

!----------------------------------------------------------------------
      SUBROUTINE Scan1DAssignX
! define SCAN1D using X-scan defined in SCAN2D
!----------------------------------------------------------------------
        integer :: I
        REAL :: DX
        SCAN1D%NVAR=1
        SCAN1D%NS=SCAN2D%NX
        SCAN1D%IX=1
        SCAN1D%IV=SCAN2D%IV
        SCAN1D%HDRVAR=SCAN2D%HXVAR
        SCAN1D%EVAL=SCAN2D%EVAL
        DX=(SCAN2D%XMAX-SCAN2D%XMIN)/(SCAN2D%NX-1)
        if (SCAN1D%MSCAN>=SCAN1D%NS) then
          do i=1,SCAN2D%NX
            SCAN1D%VARIABLES(1,i)=SCAN2D%XMIN+(i-1)*DX
          enddo
        endif
      end SUBROUTINE Scan1DAssignX

!----------------------------------------------------------------------
      SUBROUTINE WriteScan1D(iunit)
! Write SCAN1D data to given output unit
!----------------------------------------------------------------------
      integer,intent(in) :: iunit
      integer :: i,istep,IS,IL
      character*32 :: CNUM1,CNUM2,CNUM3
1     FORMAT('"',a,'"  ',$)
2     FORMAT(G15.7,2x,$)
3     FORMAT('SCAN_1D= "',a,'"')
4     FORMAT(a,'  ',$)

11     format('EVAL_VARIABLE="',a,'"')
13     format('V0=',a,' VRANGE=',a,' NP=',a)
14     format('VAR(',a,')_RANGE=(',a,',',a,')')
      if (iunit.le.0) return
      if (SCAN1D%MVAR<SCAN1D%NVAR) return
      if (SCAN1D%MSCAN<SCAN1D%NS) return
! write header line
      WRITE(iunit,11) trim(SCAN1D%EVAL%CAP)
      select case(SCAN1D%EVAL%TYP)
      case(rset_basic)
        write(iunit,3) 'basic parameters'
      case(rset_bpar,rset_dpar)
        write(iunit,3) 'peak parameters'
      case(rset_rpar)
        write(iunit,3) 'resolution parameters'
      case(rset_bpeak,rset_dpeak)
        write(iunit,3) 'beam profile'
        call float2STR(1.D0*SCAN1D%EVAL%V0,CNUM1)
        call float2STR(1.D0*SCAN1D%EVAL%DV,CNUM2)
        call INT2STR(SCAN1D%EVAL%N,CNUM3)
        write(iunit,13) trim(CNUM1),trim(CNUM2),trim(CNUM3)
        do i=1,SCAN1D%NVAR
          call FINDSTRPAR(SCAN1D%HDRVAR,':',i,IS,IL)
          call FLOAT2STR(1.D0*SCAN1D%VARIABLES(i,1),CNUM1)
          call FLOAT2STR(1.D0*SCAN1D%VARIABLES(i,SCAN1D%NS),CNUM2)
          write(iunit,14) SCAN1D%HDRVAR(IS:IS+IL-1),trim(CNUM1),trim(CNUM2)
        enddo
      end select

      select case(SCAN1D%EVAL%TYP)
      case(rset_bpeak,rset_rpeak,rset_dpeak)
      !  write(iunit,1) trim(SCAN1D%XCAP)
        do i=1,SCAN1D%NS
          call INT2STR(i,CNUM1)
          write(iunit,4) 'C'//trim(CNUM1)
        enddo
        write(iunit,*)
! write table
        do i=1,SCAN1D%EVAL%N
          do istep=1,SCAN1D%NS
            write(iunit,2) SCAN1D%VALUES(i,istep)
          enddo
          write(iunit,*)
        enddo
      case default
        do i=1,SCAN1D%NVAR
          call FINDSTRPAR(SCAN1D%HDRVAR,':',i,IS,IL)
          if (IL.gt.0) write(iunit,1) SCAN1D%HDRVAR(IS:IS+IL-1)
        enddo
        do i=1,SCAN1D%EVAL%N
          call FINDSTRPAR(SCAN1D%EVAL%HDR,'|',i,IS,IL)
          if (IL.gt.0) then
            write(iunit,1) SCAN1D%EVAL%HDR(IS:IS+IL-1)
            write(iunit,1) 'd'//SCAN1D%EVAL%HDR(IS:IS+IL-1)
          endif
        enddo
        write(iunit,*)
! write table
        do istep=1,SCAN1D%NS
          do i=1,SCAN1D%NVAR
            write(iunit,2) SCAN1D%VARIABLES(i,istep)
          enddo
          do i=1,SCAN1D%EVAL%N
            write(iunit,2) SCAN1D%VALUES(i,istep)
            write(iunit,2) SCAN1D%ERRORS(i,istep)
          enddo
          write(iunit,*)
        enddo
      end select
      end SUBROUTINE WriteScan1D


!----------------------------------------------------------------------
      SUBROUTINE WriteScan2D(iunit)
! Write SCAN2D data to given output unit
!----------------------------------------------------------------------
      integer,intent(in) :: iunit
      integer :: i,j
      character*32 :: CNUM1,CNUM2,CNUM3,CNUM4
1     FORMAT(a)
2     FORMAT(G13.5,2x,$)
3     FORMAT('SCAN_2D= "',a,'"')

20    format('XVAR (col): "',a,'"')
21    format('YVAR (row): "',a,'"')
22    format('dimension (col,row): (',a,',',a,')')
23    format('scale (col,row): ((',a,',',a,'),(',a,',',a,'))')

      if (iunit.le.0) return
! write headers
      WRITE(iunit,2) trim(SCAN2D%EVAL%CAP)
      select case(SCAN2D%EVAL%TYP)
      case(rset_basic)
        write(iunit,3) 'basic parameters'
      case(rset_bpar,rset_dpar)
        write(iunit,3) 'peak parameters'
      case(rset_rpar)
        write(iunit,3) 'resolution parameters'
      case default
        write(iunit,*) 'peak profiles can''t be recorded in 2D scan '
      end select
      write(iunit,20) trim(SCAN2D%HXVAR)
      write(iunit,21) trim(SCAN2D%HYVAR)
      call INT2STR(SCAN2D%NX,CNUM1)
      call INT2STR(SCAN2D%NY,CNUM2)
      write(iunit,22) trim(CNUM1),trim(CNUM2)
      call float2STR(1.D0*SCAN2D%XMIN,CNUM1)
      call float2STR(1.D0*SCAN2D%XMAX,CNUM2)
      call float2STR(1.D0*SCAN2D%YMIN,CNUM3)
      call float2STR(1.D0*SCAN2D%YMAX,CNUM4)
      write(iunit,23) trim(CNUM1),trim(CNUM2),trim(CNUM3),trim(CNUM4)
      write(iunit,1) 'DATA_:'
! write table
      do j=1,SCAN2D%NY
        do i=1,SCAN2D%NX
          call FLOAT2STR(1.D0*SCAN2D%VALUES(i,j),CNUM1)
          WRITE(iunit,2) trim(CNUM1)
        enddo
        WRITE(iunit,*)
      enddo
      end SUBROUTINE WriteScan2D

!-------------------------------------------------------------------
      SUBROUTINE FITGAUSS(PAR,X,Y,N)
! Fit GAUSSIAN to a peak in RES1D%X,Y
! PAR(3) ... amplitude, center and fwhm of the gaussian
!            input: initial guess, output: fitted values
! X(N),Y(N),N ... fitted curve
!-------------------------------------------------------------------
      USE OPTIMIZATION
      INTEGER*4 N,I
      REAL :: X(N),Y(N)
      REAL :: PAR(3),DPAR(3),PMI(3),PMA(3),TOL
1     format(a,6(2x,G12.5))
    !  EXTERNAL CHI2SPC
! default limits for gaussian fit
      DATA PMI/1.E-10,-1.E10,1.E-10/
      DATA PMA/1.E20,1.E10,1.E10/
      DO I=1,3
          DPAR(I)=0.
      ENDDO
      TOL=0.01
      !write(*,1) 'FITGAUSS INI PAR=',PAR(1:3)
      CALL LMOPT(CHI2GAUSS,PAR,PMI,PMA,DPAR,3,TOL,2)
      !write(*,1) 'FITGAUSS PAR=',PAR(1:3)
      IF(ABS(PAR(3)).GT.1E-10) THEN
        DO I=1,N
          Y(I)=PAR(1)*EXP(-0.5*(X(I)-PAR(2))**2/(PAR(3)/R8LN2)**2)
        ENDDO
      ELSE
        DO I=1,N
          Y(I)=0.D0
        ENDDO
      ENDIF
      END SUBROUTINE FITGAUSS

!-------------------------------------------------------------------
      REAL FUNCTION CHI2GAUSS(PAR)
! returns value to be minimized when fitting gaussian to the data in RES1D%X,Y,DY array.
!-------------------------------------------------------------------
      REAL,intent(inout) :: PAR(:)
      REAL(KIND(1.D0)) :: Z,Z1,Z2
      REAL :: CHI
      INTEGER :: I,J
1     format(a,6(2x,G12.5))
      CHI=0
      PAR(3)=ABS(PAR(3))
      J=0
      IF (ABS(PAR(3)).LT.1E-10) PAR(3)=1E-10
      DO I=1,RES1D%N
        IF(RES1D%DY(I).GT.0) THEN
          Z1=(RES1D%X(I)-PAR(2))**2
          Z2=(PAR(3)/R8LN2)**2
          Z=PAR(1)*EXP(-0.5*Z1/Z2)
          CHI=CHI+((RES1D%Y(I)-Z)/RES1D%DY(I))**2
          J=J+1
        ENDIF
      ENDDO

      if (RES1D%N>0) CHI=CHI/RES1D%N
      !write(*,1) 'CHI2GAUSS CHI=',CHI,PAR(1:3)
      CHI2GAUSS=CHI
      END FUNCTION CHI2GAUSS

!-------------------------------------------------------------------
      REAL FUNCTION CHI2MGAUSS(PAR)
! CHI2 for multi-peak Gaussian fitting
! N .. number of peaks
! PAR .. array with dimension 3xN
! Initial values should already be set in PAR
!-------------------------------------------------------------------
      USE OPTIMIZATION
      REAL :: PAR(:)
      REAL(KIND(1.D0)) :: Z
      REAL :: CHI
      INTEGER :: I,np,nd
      CHI=0
      nd=0
      np=RES1D%PKMAX-RES1D%PKMIN+1
      if (np<0) then
        CHI2MGAUSS=0.0
        return
      endif
      DO I=1,RES1D%N
        IF(RES1D%DY(I).gt.0) THEN
          nd=nd+1
          Z=MGAUSS(PAR,NP,RES1D%X(I))
          CHI=CHI+((RES1D%Y(I)-Z)/RES1D%DY(I))**2
        ENDIF
      ENDDO
      if (nd>3*np) CHI=CHI/(nd-3*np)
      CHI2MGAUSS=CHI
      END FUNCTION CHI2MGAUSS

!-------------------------------------------------------------------
      subroutine WriteBinData(IU)
!-------------------------------------------------------------------
      integer,intent(in) :: IU
      character(1024) :: LINE
      character(32) :: CNUM1,CNUM2
      integer :: i,j
1     format(a)
      call INT2STR(RES1D_BINMIN,CNUM1)
      call INT2STR(RES1D_BINMAX,CNUM2)
      ! comment line
      write(IU,1) 'COMMENT: '//trim(RES1D_BINS_COM)
      ! column headers
      LINE='X'
      do j=RES1D_BINMIN,RES1D_BINMAX
        call INT2STR(j,CNUM1)
        LINE=trim(LINE)//' Y'//trim(CNUM1)
      enddo
      write(IU,1) trim(LINE)
      ! column labels
      LINE='BINS'
      do j=RES1D_BINMIN,RES1D_BINMAX
        call FLOAT2STR((j-0.5D0)*RES1D_BINS_STEP,CNUM1)
        LINE=trim(LINE)//' Y='//trim(CNUM1)
      enddo
      write(IU,1) trim(LINE)
      ! data
      do i=1,RES1D_N
        call FLOAT2STR(1.D0*RES1D_BINS(i,0),CNUM1)
        LINE=trim(CNUM1)
        do j=RES1D_BINMIN,RES1D_BINMAX
          call FLOAT2STR(1.D0*RES1D_BINS(i,j),CNUM1)
          LINE=trim(LINE)//' '//trim(CNUM1)
        enddo
        write(IU,1) trim(LINE)
      enddo
      end subroutine WriteBinData


      end MODULE RESULTS
