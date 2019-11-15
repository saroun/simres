!//////////////////////////////////////////////////////////////////////
!////  $Id: beam1d.f,v 1.58 2019/08/16 17:16:25 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2009, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.58 $
!////     $Date: 2019/08/16 17:16:25 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Module for 1D beam profiles
!////
!//////////////////////////////////////////////////////////////////////
      MODULE BEAM1D
      use RESULTS
      use NSTORE
      use REPORTSOPT
      USE TRACINGDATA
      USE TRACINGOPT
      use COMPONENTS
      use EVENTMONITOR
      use DATASETS
      implicit none
      SAVE
      PRIVATE

! temporary variables, used during image evaluation

      real(KIND(1.D0)) :: tmp_theta
      logical :: dbg=.false.
      integer :: idbg=0

      public DETVAR,GETRANGE1D,DETRANGE1D,DOBEAM1D,BEAMPAR
      public BEAMPAR_MPK,HAS_MULTIPEAK_DATA,FITMGAUSS2

      contains

!---------------------------------------------------
      subroutine DETVARSTR(IVAR,SARG)
! get variable name by index
!---------------------------------------------------
      integer,intent(in) :: IVAR
      character(*),intent(out) :: SARG
      if (len(SARG).ge.22) then
        select case(IVAR)
        case(1)
          SARG='scattering angle [deg]'
        case(2)
          SARG='dhkl [A]'
        case(3)
          SARG='dhkl-TOF [A]'
        case(4)
          SARG='TOF [us]'
        case default
          SARG='undefined'
        end select
      endif
      end subroutine DETVARSTR

!------------------------------------------------------------
      real(KIND(1.D0)) function DETVAR(IVAR,NEU,DET,VALID)
! Get value of a detector variable defined by index=IVAR
! IVAR = TOF
! 1 ... THETA, scattering angle, radial, polar coord.
! 2 ... DHKL, steady state, powder peak d-spacing
! 3 ... DHKL, TOF, powder peak d-spacing
! 4 ... TOF
! NEU .. neutron state at the point of detection, in incident coordinates
! DET .. detector component
!-------------------------------------------------------------
      integer,intent(in) :: IVAR
      type(NEUTRON),intent(in) :: NEU
      TYPE(DETECTOR),intent(in) :: DET
      logical,intent(out) :: VALID
      real(KIND(1.D0)) :: XVAL,a,lambda,sinthb ! ,cosa,RN,R(3),t0,dt
1     format(a,': ',6(G12.5,1x))
        XVAL=0.D0
        VALID=.true.
        call FRAME_LOCAL(1,DET%FRAME,NEU)
        select case(IVAR)
        case(1,2,3,4)
          call DETECTOR_ANGLE(DET,NEU%R,NEU%T-NEU%T0,a,lambda)
          !R=NEU%R-(/0.D0,0.D0,-DET%RAD/)
          !cosa=DET%KI(1)*R(1)+DET%KI(2)*R(2)+DET%KI(3)*R(3)
          !RN=SQRT(R(1)**2+R(2)**2+R(3)**2)
          !a=acos(cosa/(RN*DET%KI0))
          VALID = (abs(a)>abs(DET%THMIN))
          VALID = (VALID.and.(abs(a)<abs(DET%THMAX)))
          !tmp_theta=a/deg
          if (VALID) then
          select case(IVAR)
      ! THETA, scattering angle, radial, polar coord.
          case(1)
            XVAL=a/deg
      ! DHKL, steady state, powder peak d-spacing
          case(2)
            !write(*,1) 'DETVAR theta, lambda', a/deg, DET%LAMBDA
            XVAL=0.5D0*DET%LAMBDA/sin(a/2.D0)
      ! DHKL, TOF, powder peak d-spacing
          case(3)
            sinthb=sin(a/2.D0)
            !if (modulation_d0>0.D0) then
            !  t0=2.D0*modulation_d0*sinthb*(DET%TOF-TOF_ZERO)/DET%LAMBDA
            !  dt=NEU%T-NEU%T0-t0
            !  ndt=NINT(dt/modulation_dt/HOVM)
            !  if (dbg) then
            !    write(*,1) 'DETVAR ',t0,NEU%T-NEU%T0,dt,ndt
            !    write(*,1) 'DETVAR ',DET%TOF,TOF_ZERO,NEU%T0
            !    write(*,1) ' L0 ',(DET%TOF-TOF_ZERO)/DET%LAMBDA*TWOPI
            !  endif
            !  dbg=.false.
            !  dt=NEU%T-NEU%T0-ndt*modulation_dt*HOVM

            !else
            !  dt=NEU%T-NEU%T0
            !endif
            !lambda=dt/(DET%TOF-TOF_ZERO)*DET%LAMBDA
            XVAL=0.5D0*lambda/sinthb
      ! TOF
          case(4)
            XVAL=NEU%T/HOVM
          end select
          endif
        end select
      ! PHI, azimuthal angle, polar coord.
      !  case(3)
      ! THETA_H, horizontal scatering angle
      !  case(4)
      ! THETA_V, vertical scatering angle
      !  case(5)
      ! QH, TOF, horizontal Q
      !  case(6)
      ! QV, TOF, vertical Q
       ! case(7)
      ! QRAD, TOF, radial Q-component
      !  case(8)
      ! EN, TOF, energy transfer
      !  case(9)
      !  end select
        DETVAR=XVAL
      end function DETVAR

!---------------------------------------------------
      subroutine GETRANGE1D(IVAR,FILT,X0,DX)
! Get range of a phase-space variable given by index
!---------------------------------------------------
      integer,intent(in) :: IVAR,FILT
      real(kind(1.D0)),intent(out) :: X0,DX
      real(kind(1.D0)) :: xmin,xmax,xval,K0
      integer :: i,NCNT
      type(NEUTRON) :: NEU
      integer,parameter :: MC=10
      integer :: IMON,ISTORE
      IMON=BMONITOR_REF
      ISTORE=BMONITORS(IMON)%ISTORE
      NCNT=NSTORE_MAXN(ISTORE) ! get number of events NCNT
      xmin=1.D20
      xmax=-1.D20
      X0=0.D0
      DX=10.D0
      K0=NSTORE_GETK0(ISTORE)
      IF (NCNT.GT.MC) THEN
        DO I=1,NCNT
          call BMONITORS_GETEVENT(IMON,i,FILT,NEU)
          if (NEU%P>0.D0) then
            XVAL=GETVAR(IVAR,NEU,K0)
            xmin=min(xmin,XVAL)
            xmax=max(xmax,XVAL)
          endif
        enddo
        X0=(xmin+xmax)/2.D0
        DX=(xmax-xmin)
      endif
      end subroutine GETRANGE1D


!---------------------------------------------------
      subroutine DETRANGE1D(IVAR,DET,X0,DX)
! Get range of a detector variable given by index
!---------------------------------------------------
      integer,intent(in) :: IVAR
      TYPE(DETECTOR),intent(in) :: DET
      real(kind(1.D0)),intent(out) :: X0,DX
      real(kind(1.D0)) :: xmin,xmax,xval
      integer :: i,NCNT,ISTORE
      type(NEUTRON) :: NEU
      integer,parameter :: MC=10
      logical :: VALID
      X0=0.D0
      DX=0.D0
      ISTORE=DET%FRAME%REGISTRY
      if (TROPT%IMODE==tr_primary) return
      NCNT=NSTORE_MAXN(ISTORE) ! get number of events NCNT
      xmin=1.D20
      xmax=-1.D20
      X0=0.D0
      DX=10.D0
      dbg=.false.
      idbg=0
      IF (NCNT.GT.MC) THEN
        DO I=1,NCNT
          call NSTORE_GETEVENT(ISTORE,i,NEU)
          XVAL=DETVAR(IVAR,NEU,DET,VALID)
          if (VALID.and.(NEU%P>0.D0)) then
            xmin=min(xmin,XVAL)
            xmax=max(xmax,XVAL)
          endif
        enddo
        X0=(xmin+xmax)/2.D0
        DX=(xmax-xmin)
      endif
      end subroutine DETRANGE1D


!----------------------------------------------------
      SUBROUTINE FILLBEAM1D(IVAR,what,XMIN,XMAX,N,RES)
! Get 1D beam profile for given variable
! IVAR = variable index (depends on the value of what)
! what=show_monitor  ... beam profile
! what=show_detector ... detector data
! XMIN, XMAX ... clip margins
! N = number of bins
! RES = data structure for output
!----------------------------------------------------
      integer,intent(in) :: IVAR,what,N
      REAL(KIND(1.D0)),intent(in) :: XMIN,XMAX
      TYPE(TRES1D) :: RES
      integer,parameter :: MC=10
      real :: EFY(SPCM,MC),Y,DY
      REAL(KIND(1.D0)) ::  dx,xval,CNORM,K0
      integer :: I,J,IP,NCNT,NN,ITH,IMON,ISTORE
      type(NEUTRON) :: NEU
      TYPE(PDETECTOR) :: DET
      logical VALID
      XVAL=0.D0
      idbg=0
1     format(a,6(G13.5,2x))
      select case(what)
      case(show_monitor)
        IMON=BMONITOR_REF
        ISTORE=BMONITORS(IMON)%ISTORE
      case(show_detector)
        if (GetDetector(DET)) then
          ISTORE=DET%X%FRAME%REGISTRY
        else
          return
        endif
      end select
      NN=MIN(N,SPCM)
      if ((what==show_detector).and.(IVAR==3)) then
        RES1D_N=NN
        RES1D_BINS_STEP=2.0
        RES1D_BINS_COM='binning by scattering angle [deg]'
      endif
      !ND=DSET_ISEL
      !NS=CSUBSET(IBEAM)
      !CNORM=MC*MCNORM(1,NS,ND)
      CNORM=MC*NSTORE_GETNORM(ISTORE)
      !write(*,1) 'FILLBEAM1D  WHAT,ISTORE,CNORM ',what,ISTORE,CNORM
      NCNT=NSTORE_MAXN(ISTORE) ! get number of events NCNT
      !write(*,1) 'FILLBEAM1D NCNT ',NCNT
      K0=NSTORE_GETK0(ISTORE)
      !write(*,1) 'FILLBEAM1D ',MCNORM(1,NS,ND),NSTORE_GETNORM(ISTORE)
      IF (NCNT.GT.MC) THEN
        EFY=0.0

        dx=(xmax-xmin)/NN
        DO I=1,NN
          RES%X(I)=xmin+(I-0.5D0)*dx
          DO J=1,MC
            EFY(I,J)=0.0
          enddo
          RES1D_BINS(I,0)=RES%X(I)
          !if (I==1) write(*,1) 'FILLBEAM1D ',RES%X(I),RES1D_BINS(I,0)
        ENDDO
        VALID=.true.
        DO I=1,NCNT
          ith=0
          if (what==show_monitor) then
            call BMONITORS_GETEVENT(IMON,i,RES1D%FILT,NEU)
            XVAL=GETVAR(IVAR,NEU,K0)
          else if (what==show_detector) then
            call NSTORE_GETEVENT(ISTORE,i,NEU)
            XVAL=DETVAR(IVAR,NEU,DET%X,VALID)
            if (IVAR==3) then
              ith=NINT(abs(tmp_theta)/RES1D_BINS_STEP+0.5)
            endif
            !dbg=(i.lt.10)
            !XVAL=DETVAR_dbg(dbg,IVAR,NEU,DET%X,VALID)
          endif
          if (VALID) then
            if ((XVAL.gt.xmin).and.(XVAL.lt.xmax)) then
              J=INT((XVAL-xmin)/dx)+1
              IP=MIN(MC,INT(1.0*MC*(1.0*I)/NCNT)+1)
              IF ((J.GE.1).AND.(J.LE.NN)) THEN
                EFY(J,IP)=EFY(J,IP)+NEU%P
                if ((ith>0).and.(ith<RES1D_BINDIM+1)) then
                  RES1D_BINMIN=MIN(RES1D_BINMIN,ith)
                  RES1D_BINMAX=MAX(RES1D_BINMAX,ith)
                  RES1D_BINS(J,ith)=RES1D_BINS(J,ith)+NEU%P
                endif
              ENDIF
            endif
          endif
        enddo
! renormalize
        DO I=1,NN
            Y=0.D0
            DY=0.D0
            do IP=1,MC
              Y=Y+EFY(I,IP)
              DY=DY+EFY(I,IP)**2
            enddo
            Y=Y/MC
            DY=DY/MC
            RES%DY(I)=0.D0
            if(DY.gt.Y**2) RES%DY(I)=SQRT(RES%DY(I)**2+CNORM**2*(DY-Y**2))
            RES%Y(I)=RES%Y(I)+CNORM*Y
            do J=RES1D_BINMIN,RES1D_BINMAX
              RES1D_BINS(I,J)=RES1D_BINS(I,J)*CNORM/MC
              !if (J==RES1D_BINMIN) write(*,1) 'FILLBEAM1D ',I,J,RES1D_BINS(I,J),RES%Y(I),CNORM
            enddo
        ENDDO
        RES%N=NN
      else
        RES%N=0
      ENDIF
      end SUBROUTINE FILLBEAM1D

!----------------------------------------------------
      SUBROUTINE DOBEAM1D(IVAR,what,XLIM,N,RES)
! Get 1D beam profile for given variable
! IVAR = variable index (depends on the value of what)
! what=show_monitor  ... beam profile
! what=show_detector ... detector data
! XLIM = limits of the independent variable
! N = number of bins
! RES = data structure for output
!----------------------------------------------------
      integer,intent(in) :: IVAR,what,N
      REAL(KIND(1.D0)),intent(in) :: XLIM(2)
      TYPE(TRES1D) :: RES
      REAL(KIND(1.D0)) :: PAR(3),dPAR(3)
      character(32) :: SARG
      integer :: i
1     format(a, ': ',6(G12.5,1x))
      if (.not.RES%APPEND) call CLEAR_RES1D(RES)
      select case(what)
      case(show_detector)
      ! for multiple peak powder model, get assumed peak positions and intensities
        if (HAS_MULTIPEAK_DATA(IVAR)) then
          RES%NPK=min(PCRYST_NREF,MAXPKS)
          do i=1,RES%NPK
            RES%PKPOS(i)=REFTAB(1,i)
            RES%PKINT(i)=REFTAB(2,i)*REFTAB(3,i)
          enddo
        endif
      end select
      call FILLBEAM1D(IVAR,what,XLIM(1),XLIM(2),N,RES)
      if ((RES%N>0).and.(REPOPT%RES)) then
        !  if (RES%NPK==0) then
            !write(*,1) 'DOBEAM1D ', XLIM
            call BEAMPAR(IVAR,what,RES%FILT, XLIM, PAR,dPAR)
            call XML_RSXDUMP(SOUT,' ',1)
            if (what==show_monitor) then
              call GETVARSTR(IVAR,SARG)
              call XML_TITLE(SOUT,'Beam profile for '//trim(SARG))
            else if (what==show_detector) then
              call DETVARSTR(IVAR,SARG)
              call XML_TITLE(SOUT,'Detector data for '//trim(SARG))
            endif
            call XML_FVALUE(SOUT,'Integrated intensity',' ',PAR(1),dPAR(1))
            call XML_FVALUE(SOUT,'Peak width',' ',PAR(2),dPAR(2))
            call XML_FVALUE(SOUT,'Peak center',' ',PAR(3),dPAR(3))
            call XML_RSXDUMP(SOUT,' ',0)
        !  endif
      endif
      END SUBROUTINE DOBEAM1D

!----------------------------------------------------
      SUBROUTINE BEAMPAR(IVAR,what,FILT,XLIM,PAR,dPAR)
! get beam profile parameters for the IVAR variable
! IVAR ... index of independent variable (the same as in BEAM1D and DET1D plots)
! what=0 ... beam profile
! what=1 ... detector data
! FILT ... fliter flag to filter monitor data
! XLIM ... min and max of the evaluated variable range
!          if min>= max, no limits are applied
! PAR ... Sum,Width,Center
! dPAR ...  Std errors of PAR
!----------------------------------------------------
      integer,intent(in) :: IVAR,what,FILT
      REAL(KIND(1.D0)),intent(in) :: XLIM(2)
      REAL(KIND(1.D0)),intent(out) :: PAR(3),dPAR(3)
      integer,parameter :: MC=10
      REAL(KIND(1.D0)) ::  XVAL,CNORM,S0,S1,S2,C,W,Z,K0
      integer :: I,IM,M(0:MC),NCNT,ISTORE,IMON, NSUM
      type(NEUTRON) :: NEU
      TYPE(PDETECTOR) :: DET
      logical :: VALID
1     format(a, 6(1x, G12.5))
      XVAL=0.D0
      PAR=0.D0
      dPAR=0.D0
      NSUM=0
      select case(what)
      case(show_monitor)
        IMON=BMONITOR_REF
        ISTORE=BMONITORS(IMON)%ISTORE
      case(show_detector)
        if (GetDetector(DET)) then
          ISTORE=DET%X%FRAME%REGISTRY
        else
          write(*,*) 'BEAMPAR: No detector registry found.'
          return
        endif
      end select

      CNORM=MC*NSTORE_GETNORM(ISTORE)
      NCNT=NSTORE_MAXN(ISTORE) ! get number of events NCNT
      if (NCNT.le.MC) return
      M(0)=0
      do IM=1,MC-1
        M(IM)=M(IM-1)+NCNT/MC
      enddo
      M(MC)=NCNT
      !write(*,*) 'BEAM1D, M: ',M
      K0=NSTORE_GETK0(ISTORE)
      VALID=.true.
      ! write(*,*) 'BEAMPAR NCNT, MC, IVAR: ',NCNT, MC, IVAR, XLIM
      IF (NCNT.GT.MC) THEN
        do IM=1,MC
          S0=0.D0
          S1=0.D0
          S2=0.D0
          DO I=M(IM-1)+1,M(IM)
            if (what==show_monitor) then
              call BMONITORS_GETEVENT(IMON,i,FILT,NEU)
              XVAL=GETVAR(IVAR,NEU,K0)
            else if (what==show_detector) then
              call NSTORE_GETEVENT(ISTORE,i,NEU)
              XVAL=DETVAR(IVAR,NEU,DET%X,VALID)

            endif
            ! apply range
            if (XLIM(2)>XLIM(1)) then
              VALID = ((XVAL>=XLIM(1)).and.(XVAL<=XLIM(2)))
            endif
            ! if (IM==1)  write(*,*) 'BEAMPAR XVAL, VALID: ',XVAL, VALID
            if (VALID) then
              Z=CNORM*NEU%P
              S0=S0+Z
              S1=S1+Z*XVAL
              S2=S2+Z*XVAL**2
              NSUM=NSUM+1
            endif
          enddo
          IF(S0.GT.0.D0) then
            C=S1/S0
            W=R8LN2*SQRT(ABS(S2/S0-C**2))
            PAR(1)=PAR(1)+S0
            dPAR(1)=dPAR(1)+S0**2
            PAR(2)=PAR(2)+W
            dPAR(2)=dPAR(2)+W**2
            PAR(3)=PAR(3)+C
            dPAR(3)=dPAR(3)+C**2
          endif
        enddo

        do I=1,3
          PAR(i)=PAR(i)/MC
          dPAR(i)=dPAR(i)/MC
          if (dPAR(i).gt.PAR(i)**2) dPAR(i)=SQRT(dPAR(i)-PAR(i)**2)
        enddo

        !write(*,1) 'BEAM1D, filter, SUM, dSUM, N: ', FILT, PAR(1), dPAR(1), NSUM
      ENDIF
      END SUBROUTINE BEAMPAR

!----------------------------------------------------
      logical function HAS_MULTIPEAK_DATA(IVAR)
! true only for DHKL plots (IVAR=2,3) and a powder sample
! with multiple peak data table
!----------------------------------------------------
      integer,intent(in) :: IVAR
      logical :: LOG1
      LOG1=.false.
      if ((IVAR==2).or.(IVAR==3)) then
        if ((SAMOBJ%ICLS>0).and.(SAMOBJ%SCLS==SCLS_POLY).and.(SAMOBJ%P_SAMPLE%TYP==isam_POWDER)) then
          if (PCRYST_NREF>0) LOG1=.true.
        endif
      endif
      HAS_MULTIPEAK_DATA=LOG1
      end function HAS_MULTIPEAK_DATA

!----------------------------------------------------
      REAL(KIND(1.D0)) function BEAMPAR_MPK(IVAR,MPAR,dMPAR,NP)
! get multiple peak parameters for the IVAR variable
! NP .. number of peaks found in plot range
!----------------------------------------------------
      integer,intent(in) :: IVAR
      REAL,intent(inout) :: MPAR(:),dMPAR(:)
      integer,intent(out) :: NP
      TYPE(PDETECTOR) :: DET
      logical :: auto,autoctr
      integer :: N
      REAL(KIND(1.D0)) :: X0,DX,ADX,CHISQR
      real :: YFIT(SPCM)
      CHISQR=0.D0
      NP=0
      BEAMPAR_MPK=0.D0
      if (.not.HAS_MULTIPEAK_DATA(IVAR)) return
      if (.not.GetDetector(DET)) return
! get plot range from DET1D command
  ! center
      X0=getCmdParam('DET1D','X0')
  ! range
      DX=getCmdParam('DET1D','DX')
  ! autoscale
      auto=(getCmdParam('DET1D','XAUTO').gt.0.D0)
      autoctr=(getCmdParam('DET1D','CAUTO').gt.0.D0)
      if (DX.le.0) auto=.true.
      if (auto.or.autoctr) then
        call DETRANGE1D(IVAR,DET%X,X0,ADX)
        if (auto) DX=ADX
      endif
      N=MIN(SPCM,NINT(getCmdParam('DET1D','NP')))
      RES1D%APPEND=.false.
      call DOBEAM1D(IVAR,show_detector,(/X0-0.5*DX, X0+0.5*DX/),N,RES1D)
      if (RES1D%NPK>0) then
        CHISQR= FITMGAUSS2(IVAR,MPAR,dMPAR,RES1D%X,YFIT,RES1D%N,NP)
      else
        CHISQR=0.0
      endif
    ! allow only peaks with centers within the data range
    !  i=0
    !  do while(i<NP)
    !    i=i+1
    !    k=3*(i-1)
    !    if ((MPAR(k+2)<X0-DX/2.D0).or.(MPAR(k+2)>X0+DX/2.D0)) then
    !      do j=i,NP-1
    !        do k1=3*(j-1)+1,3*(j-1)+3
    !          MPAR(k1)=MPAR(k1+3)
    !        enddo
    !      enddo
    !      NP=NP-1
    !    endif
    !  enddo
      BEAMPAR_MPK=CHISQR
      END function BEAMPAR_MPK


!-------------------------------------------------------------------
      REAL(KIND(1.D0)) function FITMGAUSS(IVAR,PAR,ePAR,X,Y,N,NPK)
! Fit Gaussian curves to multiple peaks in RES1D%X,Y
! PAR(3*NP) ... amplitude, center and fwhm of the gaussian
!            input: initial guess, output: fitted values
! X(N),Y(N),N ... fitted curve
!-------------------------------------------------------------------
      USE OPTIMIZATION
      INTEGER,intent(in) :: IVAR,N
      INTEGER,intent(out) :: NPK
      REAL :: PAR(:),ePAR(:)
      REAL :: X(:),Y(:)
      INTEGER :: I,J,K,M,ISIL,IERR,ib,KMIN,KMAX,NX,NP
      REAL :: DPAR(2*MAXPKS+1),PMI(2*MAXPKS+1),PMA(2*MAXPKS+1),TOL,DY0
      REAL :: MPAR(2*MAXPKS+1),eMPAR(2*MAXPKS+1),CTR(MAXPKS),XSTEP,XSTEP0
      REAL(kind(1.D0)) :: CHISQR,CHISQR1

      TYPE(TRES1D) :: RES0
      TYPE(TPKBUNCH) :: BS(MAXPKS)
      integer :: NB

      logical :: dbg=.false.
      REAL :: DBG_XYF(SPCM,3*MAXPKS)  ! values of dhkl, INT, FIT for all peaks
      integer :: DBG_NX(MAXPKS),DBG_NY,MAXROW
      ISIL=4
      if (REPOPT%MUTE) ISIL=4
      if (REPOPT%PROG) ISIL=4

1     format(a,1x,6(G13.5,2x))
    ! estimate peak parameters
      call MPEAK_GETBUNCHES(RES1D,BS,MAXPKS,NB)
      !write(*,1) 'FITMGAUSS  bunches, peaks:',NB,RES1D%PKMIN,RES1D%PKMAX
    ! store temporarily RES1D
      RES0=RES1D
      CHISQR1=0.0
      NPK=0
      RES1D%APPEND=.false.
      XSTEP0=(RES1D%X(RES1D%N)-RES1D%X(1))/(RES1D%N-1)
      PAR=0.0
      ePAR=0.0
      if (dbg) DBG_NY=0

! fit separately peaks in all bunches
      do ib=1,NB
        NP=BS(ib)%NP
        if (NP>0) then
        KMIN=MIN(BS(ib)%REF(1),BS(ib)%REF(NP))
        KMAX=MIN(BS(ib)%REF(1),BS(ib)%REF(NP))
      ! collect simulated data into a histogram with appropriate sampling
        RES1D%X=0.0
        RES1D%Y=0.0
        RES1D%dY=0.0
        RES1D%PKMIN=BS(ib)%REF(1)
        RES1D%PKMAX=BS(ib)%REF(NP)
        NX=(NP+1)*10
        NX=max(24,NX)
        NX=min(2048,NX)
        call FILLBEAM1D(IVAR,show_detector,1.D0*BS(ib)%XMIN,1.D0*BS(ib)%XMAX,NX,RES1D)
        XSTEP=(RES1D%X(RES1D%N)-RES1D%X(1))/(RES1D%N-1)
      ! calculate initial peaks parameters
        call MPEAK_GETPARAM(BS(ib),RES1D%PKINT,MPAR)
        do I=1,NP
          K=2*(I-1)+1
        ! MPEAK_GETPARAM returns widths in steps => re-scale to physical units
          MPAR(K+1)=MPAR(K+1)*XSTEP/XSTEP0
          MPAR(K+2)=MPAR(K+2)*XSTEP0
        ! limits for intensity
          PMI(K+1)=MPAR(K+1)*1e-6
          PMA(K+1)=MPAR(K+1)*1e6
        ! limits for width
          PMI(K+2)=MPAR(K+2)*0.1
          PMA(K+2)=MPAR(K+2)*10.0
          CTR(i)=RES1D%PKPOS(BS(ib)%REF(i))
        enddo
      ! use non-zero error bars outside peak range
        DY0=0.0
        do i=1,RES1D%N
          if ((DY0==0.0).and.(RES1D%DY(i)>0.0)) DY0=RES1D%DY(i)
        enddo
        do i=1,RES1D%N
          if (RES1D%DY(i)==0.0) RES1D%DY(i)=DY0
        enddo

        ! *** DBG ***
        if (dbg) then
        DBG_NY=DBG_NY+1
        DBG_NX(DBG_NY)=RES1D%N
        K=3*(DBG_NY-1)
        do i=1,RES1D%N
          DBG_XYF(i,K+1)=RES1D%X(i)
          DBG_XYF(i,K+2)=RES1D%Y(i)
          DBG_XYF(i,K+3)=MGAUSS1(MPAR,CTR,NP*2+1,RES1D%X(i))
        enddo
        endif
        ! *** DBG ***

        ! limits for x-scale factor
        PMI(1)=0.5
        PMA(1)=1.5
        ! initial increment and tolerance
        DPAR=0.0
        TOL=0.005
        !write(*,1) 'FITMGAUSS  :',ib,BS(ib)%REF(1),BS(ib)%REF(NP),MPAR(2:3)
        !write(*,1) '   X-range :',RES1D%X(1),RES1D%X(RES1D%N)
        IERR=LMFIT(MGAUSS1,RES1D%X,RES1D%Y,RES1D%DY,RES1D%N,CTR,MPAR,eMPAR,PMI,PMA,DPAR,NP*2+1,TOL,ISIL)
        if (IERR<-2) then
          CHISQR=1.D10
          if (isil.eq.3) write(*,1)   'bunch fit failed: ',ib,np
        else
          CHISQR = GETCHI2(MGAUSS1,RES1D%X,RES1D%Y,RES1D%DY,RES1D%N,CTR,MPAR,NP*2+1)
          do I=1,NP
            K=2*(I-1)+1
            M=3*(NPK+I-1)
          ! we need to rescale peak intensity according to the original bin width
            PAR(M+1)=MPAR(K+1)*XSTEP0/XSTEP
            PAR(M+3)=MPAR(K+2)
            ePAR(M+1)=eMPAR(K+1)*XSTEP0/XSTEP
            ePAR(M+3)=eMPAR(K+2)
            PAR(M+2)=CTR(i)*MPAR(1)
            ePAR(M+2)=CTR(i)*eMPAR(1)
            if (isil.eq.3) write(*,1) 'peaks : ',i,BS(ib)%REF(i),MPAR(K+1),MPAR(K+2),CTR(i)*MPAR(1)
          enddo
          NPK=NPK+NP
          CHISQR1=max(CHISQR1,CHISQR)
          if (isil.eq.3) write(*,1)   'bunch fit  : ',ib,np,CHISQR
        endif
        endif
      endDO
      if (isil.eq.3) write(*,*) 'max. fitted CHI2: ',CHISQR1

      ! *** DBG ***
      if (dbg) then
        Open(Unit=22,File='dbg_FITMGAUSS.txt',Status='Unknown')
        MAXROW=0
        do i=1,DBG_NY
          MAXROW=MAX(MAXROW,DBG_NX(i))
        enddo
2       format(G13.5,a1,G13.5,a1,G13.5,a1,$)
3       format(a1,a1,a1,$)
        do i=1,MAXROW
          do j=1,DBG_NY
            K=3*(j-1)
            if (DBG_NX(j)>=i) then
              write(22,2) DBG_XYF(i,k+1),char(9),DBG_XYF(i,k+2),char(9),DBG_XYF(i,k+3),char(9)
            else
              write(22,3) char(9),char(9),char(9)
            endif
          enddo
          write(22,*)
        enddo
        CLOSE(22)
      endif
      ! *** DBG ***


    ! restore RES1D
      RES1D=RES0
      DO I=1,N
        Y(I)=MGAUSS(PAR,NPK*3,X(I))
      enddo
      FITMGAUSS=CHISQR1
      END function FITMGAUSS

!-------------------------------------------------------------------
      REAL(KIND(1.D0)) function FITMGAUSS2(IVAR,PAR,ePAR,X,Y,N,NPK)
! Fit Gaussian curves to multiple peaks in RES1D%X,Y
! PAR(3*NP) ... amplitude, center and fwhm of the gaussian
!            input: initial guess, output: fitted values
! X(N),Y(N),N ... fitted curve
!-------------------------------------------------------------------
      USE OPTIMIZATION
      INTEGER,intent(in) :: IVAR,N
      INTEGER,intent(out) :: NPK
      REAL :: PAR(:),ePAR(:)
      REAL :: X(:),Y(:)
      INTEGER :: I,J,K,M,ISIL,IERR,ib,KMIN,KMAX,NX,NP
      REAL :: DPAR(MAXPKS+2),PMI(MAXPKS+2),PMA(MAXPKS+2),TOL,DY0
      REAL :: MPAR(MAXPKS+2),eMPAR(MAXPKS+2),CTR(MAXPKS),XSTEP,XSTEP0
      REAL(kind(1.D0)) :: CHISQR,CHISQR1

      TYPE(TRES1D) :: RES0
      TYPE(TPKBUNCH) :: BS(MAXPKS)
      integer :: NB

      logical :: dbg=.false.
      REAL :: DBG_XYF(SPCM,3*MAXPKS)  ! values of dhkl, INT, FIT for all peaks
      integer :: DBG_NX(MAXPKS),DBG_NY,MAXROW
      ISIL=4
      if (REPOPT%MUTE) ISIL=4
      if (REPOPT%PROG) ISIL=4

1     format(a,1x,6(G13.5,2x))
    ! estimate peak parameters
      call MPEAK_GETBUNCHES(RES1D,BS,MAXPKS,NB)
      !write(*,1) 'FITMGAUSS  bunches, peaks:',NB,RES1D%PKMIN,RES1D%PKMAX
    ! store temporarily RES1D
      RES0=RES1D
      CHISQR1=0.0
      NPK=0
      RES1D%APPEND=.false.
      XSTEP0=(RES1D%X(RES1D%N)-RES1D%X(1))/(RES1D%N-1)
      PAR=0.0
      ePAR=0.0
      if (dbg) DBG_NY=0
! fit separately peaks in all bunches
      do ib=1,NB
        NP=BS(ib)%NP
        if (NP>0) then
        KMIN=MIN(BS(ib)%REF(1),BS(ib)%REF(NP))
        KMAX=MIN(BS(ib)%REF(1),BS(ib)%REF(NP))
      ! collect simulated data into a histogram with appropriate sampling
        RES1D%X=0.0
        RES1D%Y=0.0
        RES1D%dY=0.0
        RES1D%PKMIN=BS(ib)%REF(1)
        RES1D%PKMAX=BS(ib)%REF(NP)
        NX=(NP+1)*10
        NX=max(24,NX)
        NX=min(2048,NX)
        call FILLBEAM1D(IVAR,show_detector,1.D0*BS(ib)%XMIN,1.D0*BS(ib)%XMAX,NX,RES1D)
        XSTEP=(RES1D%X(RES1D%N)-RES1D%X(1))/(RES1D%N-1)
      ! calculate initial peaks parameters
        call MPEAK_GETPARAM2(BS(ib),RES1D%PKINT,MPAR)
        if (isil.eq.3) write(*,1) 'bunch : ',ib,NP,MPAR(1),MPAR(2)*XSTEP0
        MPAR(2)=MPAR(2)*XSTEP0
        PMI(2)=MPAR(2)*0.1
        PMA(2)=MPAR(2)*10.0
        do I=1,NP
        !  K=I+1
        ! MPEAK_GETPARAM returns widths in steps => re-scale to physical units
          MPAR(I+2)=MPAR(I+2)*XSTEP/XSTEP0
        ! limits for intensity
          PMI(I+2)=MPAR(I+2)*1e-6
          PMA(I+2)=MPAR(I+2)*1e6
        ! limits for width
          CTR(i)=RES1D%PKPOS(BS(ib)%REF(i))
          if (isil.eq.3) write(*,1) 'peaks : ',i,BS(ib)%REF(i),MPAR(I+2),CTR(i)
        enddo
      ! use non-zero error bars outside peak range
        DY0=0.0
        do i=1,RES1D%N
          if ((DY0==0.0).and.(RES1D%DY(i)>0.0)) DY0=RES1D%DY(i)
        enddo
        do i=1,RES1D%N
          if (RES1D%DY(i)==0.0) RES1D%DY(i)=DY0
        enddo

        ! *** DBG ***
        if (dbg) then
        DBG_NY=DBG_NY+1
        DBG_NX(DBG_NY)=RES1D%N
        K=3*(DBG_NY-1)
        do i=1,RES1D%N
          DBG_XYF(i,K+1)=RES1D%X(i)
          DBG_XYF(i,K+2)=RES1D%Y(i)
          DBG_XYF(i,K+3)=MGAUSS2(MPAR,CTR,NP+2,RES1D%X(i))
          !if (ib==2) then
          !  write(*,1) 'MGAUSS2 ',RES1D%X(i),MPAR(1),MPAR(2)
          !  do j=1,NP
          !    write(*,1) 'MGAUSS2 ',j,CTR(j),MPAR(j+2)
          !  enddo
          !endif
        enddo
        endif
        ! *** DBG ***

        ! limits for x-scale factor
        PMI(1)=0.5
        PMA(1)=1.5
        ! initial increment and tolerance
        DPAR=0.0
        TOL=0.005
        !write(*,1) 'FITMGAUSS  :',ib,BS(ib)%REF(1),BS(ib)%REF(NP),MPAR(2:3)
        !write(*,1) '   X-range :',RES1D%X(1),RES1D%X(RES1D%N)
        IERR=LMFIT(MGAUSS2,RES1D%X,RES1D%Y,RES1D%DY,RES1D%N,CTR,MPAR,eMPAR,PMI,PMA,DPAR,NP+2,TOL,ISIL)
        if (IERR<-2) then
          CHISQR=1.D10
          if (isil.eq.3) write(*,1)   'bunch fit failed: ',ib,np
        else
          CHISQR = GETCHI2(MGAUSS2,RES1D%X,RES1D%Y,RES1D%DY,RES1D%N,CTR,MPAR,NP+2)
          if (isil.eq.3) write(*,1) 'bunch  fit: ',ib,NP,MPAR(1),MPAR(2),CHISQR
          do I=1,NP
            M=3*(NPK+I-1)
          ! we need to rescale peak intensity according to the original bin width
            PAR(M+1)=MPAR(I+2)*XSTEP0/XSTEP
            PAR(M+3)=MPAR(2)
            ePAR(M+1)=eMPAR(I+2)*XSTEP0/XSTEP
            ePAR(M+3)=eMPAR(2)
            PAR(M+2)=CTR(i)*MPAR(1)
            ePAR(M+2)=CTR(i)*eMPAR(1)
            !if (isil.eq.3) write(*,1) 'peaks : ',i,BS(ib)%REF(i),MPAR(K+1),MPAR(2),CTR(i)*MPAR(1)
            if (isil.eq.3) write(*,1) 'peak fit : ',i,BS(ib)%REF(i),MPAR(I+2),CTR(i)
          enddo
          NPK=NPK+NP
          CHISQR1=max(CHISQR1,CHISQR)
          !if (isil.eq.3) write(*,1)   'bunch fit  : ',ib,np,CHISQR
        endif
        endif
      endDO
      if (isil.eq.3) write(*,*) 'max. fitted CHI2: ',CHISQR1

      ! *** DBG ***
      if (dbg) then
        Open(Unit=22,File='dbg_FITMGAUSS.txt',Status='Unknown')
        MAXROW=0
        do i=1,DBG_NY
          MAXROW=MAX(MAXROW,DBG_NX(i))
        enddo
2       format(G13.5,a1,G13.5,a1,G13.5,a1,$)
3       format(a1,a1,a1,$)
        do j=1,DBG_NY
          K=3*(j-1)
          write(*,*) j,DBG_NX(j),DBG_XYF(1,k+1)
        enddo
        do i=1,MAXROW
          do j=1,DBG_NY
            K=3*(j-1)
            if (DBG_NX(j)>=i) then
              write(22,2) DBG_XYF(i,k+1),char(9),DBG_XYF(i,k+2),char(9),DBG_XYF(i,k+3),char(9)
            else
              write(22,3) char(9),char(9),char(9)
            endif
          enddo
          write(22,*)
        enddo
        CLOSE(22)
      endif
      ! *** DBG ***


    ! restore RES1D
      RES1D=RES0
      DO I=1,N
        Y(I)=MGAUSS(PAR,NPK*3,X(I))
      enddo
      FITMGAUSS2=CHISQR1
      END function FITMGAUSS2

!-------------------------------------------------------------------
      REAL FUNCTION MGAUSS1(PAR,CTR,NP,X)
! Get value of multiple Gauss function at X
! NP  .. number of parameters
! PAR .. array with dimension NP
! CTR .. nominal peak positions (not fitted)
! NOTE: NP must be odd integer, np=2*npk+1, where nkp is the number of peaks
! no check is made in this procedure !!
!-------------------------------------------------------------------
      REAL,intent(inout) :: PAR(:)
      REAL,intent(in) :: CTR(:),X
      integer,intent(in) :: NP
      REAL(KIND(1.D0)) :: Z,Z1
      INTEGER :: K,M,NPK
      Z=0.D0
      NPK=(NP-1)/2
      DO M=1,NPK
        K=2*(M-1)+1
        PAR(K+2)=ABS(PAR(K+2))
        IF (ABS(PAR(K+2)).LT.1E-10) PAR(K+2)=1E-10
        Z1=((X-CTR(M)*PAR(1))/PAR(K+2)*R8LN2)**2
        if (Z1<20.D0) Z=Z+PAR(K+1)*EXP(-0.5*Z1)
      enddo
      MGAUSS1=Z
      END FUNCTION MGAUSS1


!-------------------------------------------------------------------
      REAL FUNCTION MGAUSS2(PAR,CTR,NP,X)
! Like MGAUSS1, but with only one parameter for peak width
!-------------------------------------------------------------------
      REAL,intent(inout) :: PAR(:)
      REAL,intent(in) :: CTR(:),X
      integer,intent(in) :: NP
      REAL(KIND(1.D0)) :: Z,Z1
      INTEGER :: M,NPK
      Z=0.D0
      NPK=NP-2
      DO M=1,NPK
        PAR(2)=ABS(PAR(2))
        IF (PAR(2).LT.1E-10) PAR(2)=1E-10
        Z1=((X-CTR(M)*PAR(1))/PAR(2)*R8LN2)**2
        if (Z1<20.D0) Z=Z+PAR(M+2)*EXP(-0.5*Z1)
        ! if (X>1.64) write(*,*) 'MGAUSS2 ',X,M,PAR(1),PAR(2),CTR(M)
      enddo
      MGAUSS2=Z
      END FUNCTION MGAUSS2




      END MODULE BEAM1D
