!//////////////////////////////////////////////////////////////////////
!////  $Id: grfplot.f,v 1.11 2019/06/22 20:48:05 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.11 $
!////     $Date: 2019/06/22 20:48:05 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Module with graphics subroutines for plotting basic graphs in SIMRES
!////  inherited from GRFDAT.mod RESTRAX version, Revision: 1.23
!////  Requires: PGPLOT
!//////////////////////////////////////////////////////////////////////
      MODULE GRFPLOT
      use GRFDATA
      USE COMMANDS
      implicit NONE

      contains


!-------------------------------------------------
      logical function INITGRF(TOPRINT)
!  Initialization of the PGPlot graphics device.
!  TOPRINT=0 ... output to the current device
!  TOPRINT=1 ... output to the PostScript file "restrax.ps"
!-------------------------------------------------
      integer,intent(in) :: TOPRINT
      CHARACTER*60 DST,AUX
      INTEGER :: IDEVSTR,PGOPEN,IFIRST
      SAVE IFIRST

      call getDeviceStr(GRFSEL,DST)
      IF(TOPRINT.EQ.1) DST='"restrax.ps"/VPS'

      IDEVSTR=PGOPEN(TRIM(DST))
      IF (IDEVSTR.LE.0) THEN
          write(*,*) 'PGOPEN not successful, graphics disabled'
          call setCmdParam('PLOT','ENABLED',0.D0)
          INITGRF=.false.
      else
        INITGRF=.true.
        call getDeviceStr(GRFSEL,DST)
      END IF

! set /XSERV as portrait at the beginning
      IF (IFIRST.EQ.0) THEN
        AUX=DST
        CALL MKUPCASE(AUX)
        IF (INDEX(TRIM(AUX),'/XSERV').GT.0) THEN
          CALL PGPAP(4.0,1.4)
          IFIRST=1
        ENDIF
      ENDIF
! open the page
      CALL PGPAGE
      END function INITGRF

!-------------------------------------------------
      SUBROUTINE CLOSEGRF(TOPRINT)
!  Initialization of the PGPlot graphics device.
!-------------------------------------------------
      integer,intent(in) :: TOPRINT
      CALL PGCLOS
      IF (TOPRINT.GT.0) THEN
  !       CALL PRINTFILE('restrax.ps')
      ENDIF
      end SUBROUTINE CLOSEGRF

!----------------------------------------------------
      SUBROUTINE LISTGRFDEV
!  Fill the array GRFDEV with available PGPLOT devices
!-----------------------------------------------------
      INTEGER*4 I,TLEN,DLEN,INTER,IS,IL
      CHARACTER*8 SELTYPE,DTYPE
      CHARACTER*80 DESCR
      SELTYPE=' '
      if (GRFSEL.gt.0) SELTYPE=GRFDEV(GRFSEL)%DTYPE
      GRFSEL=0
      CALL PGQNDT(GRFNUM)
      GRFNUM=MIN(GRFNUM,MAXGDEV)
      DO I=1,GRFNUM
        CALL PGQDT(I, DTYPE, TLEN, DESCR, DLEN, INTER)
        DLEN=MIN(DLEN,80)
        TLEN=MIN(TLEN,8)
        CALL BOUNDS(DTYPE(1:TLEN),IS,IL)
        GRFDEV(I)%IORD=I
        GRFDEV(I)%INTER=INTER
        call MKUPCASE(DTYPE)
        IF (trim(DTYPE).EQ.trim(SELTYPE)) GRFSEL=I
        GRFDEV(I)%DTYPE=DTYPE(IS:IS+IL-1)
!// remove brackets from DESCR
        CALL BOUNDS(DESCR(1:DLEN),IS,IL)
        IF (DESCR(IS:IS).EQ.'(') THEN
          IS=IS+1
          IL=IL-1
        endif
        IF (DESCR(IS+IL-1:IS+IL-1).EQ.')') IL=IL-1
        GRFDEV(I)%DESCR=DESCR(IS:IS+IL-1)
c! write(*,*) 'LISTGRFDEV ',trim(GRFDEV(I)%DTYPE)
      ENDDO
      if (GRFSEL.le.0) GRFSEL=1
      END SUBROUTINE LISTGRFDEV

!-------------------------------------------------------------
      SUBROUTINE EraseViewPort(iport)
!  clear viewport
!-------------------------------------------------------------
      IMPLICIT NONE
      integer, intent(in) :: iport
      real :: D(4)
      D=GRFPORT(iport)%DC
      CALL PGSVP(D(1),D(2),D(3),D(4))
      CALL PGSWIN(0.,1.,0.,1.)
      CALL PGSCI(0)
      CALL PGRECT(0.,1.,0.,1.)
      CALL PGSCI(1)
      END SUBROUTINE EraseViewPort

!--------------------------------------------------------------------
      subroutine AllocPort(iport)
! Allocates graph space by calling PGPLOT commands PGSVP and PGSWIN
!--------------------------------------------------------------------
      IMPLICIT NONE
      integer,intent(in) :: iport
      real :: DC(4),WC(4)

      DC=GRFPORT(iport)%DC
      WC=GRFPORT(iport)%WC
      CALL PGSVP(DC(1),DC(2),DC(3),DC(4))
      if (GRFPORT(iport)%LOGSCALE(2).ne.0) then
        if ((WC(3).GT.0).and.(WC(4).gt.0)) then
          WC(3)=log10(abs(WC(3)))
          WC(4)=log10(abs(WC(4)))
        else
          write(*,*) 'AllocPort: wrong log-scale ',WC(3),WC(4)
        endif
      endif
      CALL PGSWIN(WC(1),WC(2),WC(3),WC(4))
      CALL PGSCI(GRFPORT(iport)%COLOR)
      CALL PGSCH(GRFPORT(iport)%CHSIZE)
      call PGSLW(GRFPORT(iport)%LWIDTH)
      CALL PGSLS(1)
      end subroutine AllocPort

!--------------------------------------------------------------------
      SUBROUTINE PlotGrfData(iport,imin,imax)
!  Plot a curve on GRFPORT(iport).
!   iport      ... index of GRFPORT
!   imin,imax  ... index range of graph data to plot
! CALLED BY: VIEWSCAN
!--------------------------------------------------------------------
      INTEGER,intent(in) :: iport,imin,imax
      REAL :: XX(MAXPNTS),YY(MAXPNTS),Y1(MAXPNTS),Y2(MAXPNTS),REF(2)

      INTEGER :: igrf,i
      logical :: LGY

! nothing to plot ?
      IF (GRFPORT(iport)%WC(1).GE.GRFPORT(iport)%WC(2)) RETURN
      IF (GRFPORT(iport)%WC(3).GE.GRFPORT(iport)%WC(4)) RETURN
      LGY=(GRFPORT(iport)%LOGSCALE(2).NE.0.0)
      if (LGY.and.(GRFPORT(iport)%WC(3).le.0)) then
        call MSG_ERROR('PlotGrfData','Log-scale requires lower limit > 0',0,1)
        return
      endif
! don't forget to allocate space on the graph
      call AllocPort(iport)
      REF(1)=GRFPORT(iport)%REF(GRFPORT(iport)%IX)
      REF(2)=GRFPORT(iport)%REF(GRFPORT(iport)%IY)
! plot curves
      DO igrf=imin,imax
      if (GRF1D(igrf)%NP.GT.0) then
        !write(*,*) 'PlotGrfData ',igrf,GRF1D(igrf)%NPEAKS,GRF1D(igrf)%PKPAR(1:3)
        !write(*,*) 'PlotGrfData ',igrf,GRF1D(igrf)%LINE,GRF1D(igrf)%POINT

        if (LGY) then
          do i=1,GRF1D(igrf)%NP
            YY(i)=log10(MAX(GRF1D(igrf)%Y(i)-REF(2),GRFPORT(iport)%WC(3)/10.0))
            Y1(i)=log10(MAX(GRF1D(igrf)%Y(i)-GRF1D(igrf)%DY(i)-REF(2),GRFPORT(iport)%WC(3)/10.0))
            Y2(i)=log10(MAX(GRF1D(igrf)%Y(i)+GRF1D(igrf)%DY(i)-REF(2),GRFPORT(iport)%WC(3)/10.0))
          enddo
        else
          do i=1,GRF1D(igrf)%NP
            YY(i)=GRF1D(igrf)%Y(i)-REF(2)
            Y1(i)=YY(i)-GRF1D(igrf)%DY(i)
            Y2(i)=YY(i)+GRF1D(igrf)%DY(i)
          enddo
        endif
        do i=1,GRF1D(igrf)%NP
          XX(i)=GRF1D(igrf)%X(i)-REF(1)
        enddo
        CALL PGSCI(GRF1D(igrf)%COLOR)
        IF (GRF1D(igrf)%ER.NE.0) then
          CALL PGSLS(GRF1D(igrf)%LINE)
          CALL PGERRY(GRF1D(igrf)%NP,XX(1),Y1(1),Y2(1),0.0)
        endif
        IF (GRF1D(igrf)%LINE.NE.0) then
          CALL PGSLS(GRF1D(igrf)%LINE)
          CALL PGLINE(GRF1D(igrf)%NP,XX(1),YY(1))
        endif
        IF (GRF1D(igrf)%POINT.NE.0) then
        ! write(*,*) 'PlotGrfData ',igrf,GRF1D(igrf)%POINT,GRF1D(igrf)%PSIZE
          call PGSLS(1)
          CALL PGSCH(GRF1D(igrf)%PSIZE)
          CALL PGPOINT(GRF1D(igrf)%NP,XX(1),YY(1),GRF1D(igrf)%POINT)
        endif
      endif
      enddo
      END SUBROUTINE PlotGrfData


!--------------------------------------------------------------------
      SUBROUTINE PlotGrfData2D(iport,imin,imax)
! Plot gray-scale maps for all 2D graphs from the specified range.
! on given viewport
!--------------------------------------------------------------------
      INTEGER,intent(in) :: iport,imin,imax
      REAL :: logz
      INTEGER :: igrf,nx,ny,minz
1     format('PlotGrfData2D: igrf=',I2,' port=',I2)
2     format('PlotGrfData2D: range X=',G10.4,' .. ',G10.4,' range Y=',G10.4,' .. ',G10.4)

! don't forget to allocate space on the graph
      call AllocPort(iport)

! scan through the specified range of maps/viewports
      DO igrf=MAX(1,imin),MIN(MAX2D,imax)

        nx=GRF2D(igrf)%NX
        ny=GRF2D(igrf)%NY
        logz=GRFPORT(iport)%LOGSCALE(3)
        minz=GRFPORT(iport)%MINZ
        if (nx*ny.GT.0) then
! nothing to plot ?
          IF (GRFPORT(iport)%WC(1).GE.GRFPORT(iport)%WC(2)) RETURN
          IF (GRFPORT(iport)%WC(3).GE.GRFPORT(iport)%WC(4)) RETURN
! plot curves
          !write(*,1) igrf,iport
          !write(*,2) GRFPORT(iport)%WC(1:4)
          CALL PLOT2D(iport,GRF2D(igrf)%Z(1,1),nx,ny,MAXPIX,MAXPIX,logz,minz)
        endif
      enddo
      END SUBROUTINE PlotGrfData2D

!--------------------------------------------------------------------
      SUBROUTINE PlotGrfEll(iport,imin,imax)
!  Plot projection and section ellipses on GRFPORT(iport).
!   iport      ... index of GRFPORT
!   imin,imax  ... index range of graph data to plot
!--------------------------------------------------------------------
      INTEGER,intent(in) :: iport,imin,imax
      real, parameter :: C2LN2=1.38629436112 ! 2*ln(2)
      INTEGER :: igrf,I,N,NP,ip,iline,IX,IY
      REAL :: X(MAXPNTS),Y(MAXPNTS),B(2,2),U(2),PROJ(2,2),SECT(2,2)
      REAL :: WX,WY,WEX,WEY,D,STEP,DET,Z
      integer :: IDX(2)
1     format('Plot ellipse No ',I2,' at ',2(G10.4,1x),' port ',I2)
2     format(a,': ',4(G10.4,1x))
! nothing to plot ?
      IF (GRFPORT(iport)%WC(1).GE.GRFPORT(iport)%WC(2)) RETURN
      IF (GRFPORT(iport)%WC(3).GE.GRFPORT(iport)%WC(4)) RETURN

! don't forget to allocate space on the graph
      call AllocPort(iport)
  ! with of the plot area
      WX=GRFPORT(iport)%WC(2)-GRFPORT(iport)%WC(1)
      WY=GRFPORT(iport)%WC(4)-GRFPORT(iport)%WC(3)
  ! local copy of projection indices
      IX=GRFPORT(iport)%IX
      IY=GRFPORT(iport)%IY
  ! projection indices
      IDX=(/IX,IY/)
! plot projection and section ellipses
      DO igrf=imin,imax
      ! ellipsoid center projected on the viewport
        U(1)=GRFELL(igrf)%RZERO(IX)-GRFPORT(iport)%REF(IX)
        U(2)=GRFELL(igrf)%RZERO(IY)-GRFPORT(iport)%REF(IY)
        write(*,*) IDX
        write(*,1) igrf,U,iport
        call EllipsoidView4(GRFELL(igrf)%RMAT(1,1),IX,IY,PROJ,SECT)
        do ip=1,2
          if (ip.eq.1) then
      ! projection
            B=PROJ
            iline=GRFELL(igrf)%PLINE
          else
      ! section
            B=SECT
            iline=GRFELL(igrf)%SLINE
          endif

      ! plot ellipse defined by the B(2,2) matrix
          DET=B(2,2)*B(1,1)-B(1,2)**2
          write(*,*) 'part ',ip,det
          write(*,2) 'B',B(1:2,1)
          write(*,2) 'B',B(1:2,2)
          if ((DET.GT.0.0).and.(iline.gt.0)) then
      ! get ellipse widths
            WEX=2.0*SQRT(C2LN2*B(2,2)/DET)
            WEY=2.0*SQRT(C2LN2*B(1,1)/DET)
      ! get optimum number of points for the curve
      ! full viewport = MAXPNTS points
            NP=INT(MAXPNTS*MAX(WEX/WX,WEY/WY,0.99))
            N=(NP+1)/2
            STEP=WEX/(N-1)
      ! fill X,Y with ellipse data
            DO I=1,N
              D=I-1
              Z=-WEX/2.D0+D*STEP
              X(I)=Z+U(1)
              DET=(B(1,2)*Z)**2-B(2,2)*(B(1,1)*Z**2-C2LN2)
              Y(I)=(-B(1,2)*Z+SQRT(ABS(DET)))/B(2,2)+U(2)
            enddo
            DO I=N+1,2*N-1
              X(I)=-X(I-N+1)+2.*U(1)
              Y(I)=-Y(I-N+1)+2.*U(2)
            enddo
            if ((2*N-1).LT.NP) then
              X(NP)=X(1)
              Y(NP)=Y(1)
            endif
            CALL PGSCI(GRFELL(igrf)%COLOR)
            CALL PGSLS(iline)
            CALL PGLINE(NP,X(1),Y(1))
          ENDIF
        endDO
      enddo
      END SUBROUTINE PlotGrfEll


!--------------------------------------------------------------------
      SUBROUTINE PlotGrfTxt(iport,imin,imax)
!  Plot text assigned to GRFPORT(iport)
!   iport      ... index of GRFPORT
!   imin,imax  ... index range of text data
!--------------------------------------------------------------------
      INTEGER,intent(in) :: iport,imin,imax
      integer :: igrf,i,is,il,nc
      real :: di,sh,an
! don't forget to allocate space on the graph
      call AllocPort(iport)
      DO igrf=imin,imax
        if (GRFTXT(igrf)%TEXT.NE.' ') then
          di=GRFTXT(igrf)%DIST
          sh=GRFTXT(igrf)%SHIFT
          an=GRFTXT(igrf)%angle
          !write(*,*) 'PlotGrfTxt ',igrf,trim(GRFTXT(igrf)%TEXT)
          !write(*,*) 'PlotGrfTxt ',igrf,di,sh,an,GRFTXT(igrf)%SIDE

          call PGSLW(GRFTXT(igrf)%LWIDTH)
          call PGSCH(GRFTXT(igrf)%CHSIZE)
          call PGSCI(GRFTXT(igrf)%COLOR)
          call COUNTPAR(trim(GRFTXT(igrf)%TEXT),'|',NC)
          do i=1,NC
            CALL FINDSTRPAR(GRFTXT(igrf)%TEXT,'|',i,is,il)
            CALL PGMTEXT(GRFTXT(igrf)%SIDE,di+(i-1)*1.0,sh,an,trim(GRFTXT(igrf)%TEXT(is:is+il-1)))
          enddo
        endif
      enddo
      end SUBROUTINE PlotGrfTxt

!--------------------------------------------------------------------
      SUBROUTINE PlotAllGraphs
! Plot all graphs defined in the module
!--------------------------------------------------------------------
      INTEGER :: igrf,iport,ig,ng,np
1     format('Plot port=',I3,' graph=',I3,' type=',I3)
2     format(a8,': ',6(G10.4,1x))
    ! scan through viewports
        np=0
        do iport=1,MAXPORTS
          ng=GRFPORT(iport)%NGRF
          if (ng.le.0) goto 100
      !    write(*,2) 'GRFPORT',iport,GRFPORT(iport)%REF
        ! set world scale to the viewport if not yet done and clear it
          call SetViewPortScale(iport)
          call EraseViewPort(iport)
          do igrf=1,ng
            ig=GRFPORT(iport)%IGRF(igrf)
c            write(*,1) iport,igrf,GRFPORT(iport)%TGRF(igrf)
            select case(GRFPORT(iport)%TGRF(igrf))
            case(grtyp_1D)
              call PlotGrfData(iport,ig,ig)
            case(grtyp_2D)
              call PlotGrfData2D(iport,ig,ig)
            case(grtyp_ell)
              call PlotGrfEll(iport,ig,ig)
            case(grtyp_txt)
              call PlotGrfTxt(iport,ig,ig)
            end select
          enddo
        ! plot frame at the end
          CALL PLOTNEWFRAME(iport)
          call PLOTCOMMENT(iport)
          np=np+1 ! count number of ports
100       continue
        enddo
        CALL PGIDEN
    ! reset color, line and character style to defaults
        CALL PGSCI(1)   ! color
        CALL PGSLS(1)   ! line style
        call PGSLW(1)   ! line width
        call PGSCH(1.0) ! character size
      END SUBROUTINE PlotAllGraphs

!---------------------------------------------------------------
      SUBROUTINE PLOTCOMMENT(iport)
!  plots frame, axes and titles of a viewport from GRFDATA module
!---------------------------------------------------------------
      integer,intent(in) :: iport
      integer :: NC,i,is,il
      call COUNTPAR(trim(GRFPORT(iport)%comment),'|',NC)
      if (NC.gt.0) then
        call AllocPort(iport)
        call PGSCH(0.8*GRFPORT(iport)%CHSIZE)
        do i=1,NC
          CALL FINDSTRPAR(GRFPORT(iport)%comment,'|',i,is,il)
          CALL PGMTEXT('B',6.0+(i-1)*1.0,-0.05,0.0,trim(GRFPORT(iport)%comment(is:is+il-1)))
        enddo
      endif
      end SUBROUTINE PLOTCOMMENT

!---------------------------------------------------------------
      SUBROUTINE PLOTNEWFRAME(iport)
!  plots frame, axes and titles of a viewport from GRFDATA module
!---------------------------------------------------------------
      integer,intent(in) :: iport
      character*64 STYLEX,STYLEY
  ! no frame for style<0
      IF (GRFPORT(iport)%STYLE.LT.0) RETURN
      call AllocPort(iport)
      STYLEX='BCNST'
      STYLEY='BCNST'
      if (GRFPORT(iport)%LOGSCALE(2).NE.0) then
        STYLEY='BCLNST'
      endif
      IF(GRFPORT(iport)%STYLE.NE.0) THEN
        STYLEX='A'//trim(STYLEX)
        STYLEY='A'//trim(STYLEY)
      endif
      CALL PGBOX(trim(STYLEX),0.0,0,trim(STYLEY),0.0,0)
      CALL PGLAB(GRFPORT(iport)%XTIT,GRFPORT(iport)%YTIT,GRFPORT(iport)%HEAD)
      END SUBROUTINE PLOTNEWFRAME

!-------------------------------------------------------------
      SUBROUTINE PLOT2D(iport,A,NX,NY,NDX,NDY,LOGZ,MINZ)
! plots a gray-scale map of the array A to the viewport PORT
! assumes positive (>=0) values in A matrix, negative ones are ignored
! LOGZ = true to plot on log scale
! MINZ = how to get minimum z:
!         minimum (0), zero (1), minimum positive (2)
!-------------------------------------------------------------
      integer,intent(in) :: iport,NX,NY,NDX,NDY,MINZ
      REAL,intent(in) :: A(NDX,NDY)
      REAL,intent(in) :: LOGZ ! base for logarithmic scale
      INTEGER :: I,J
      REAL :: CX,DX,CY,DY,TR(6),AMIN,AMINP,AMAX,WC(4),DC(4)
      REAL :: AUX(MIMAX,MIMAX),L0,RANGE
1     format('Z-scale: ',G12.4,', ',G12.4,', ',G12.4)

      DC=GRFPORT(iport)%DC
      WC=GRFPORT(iport)%WC
      DX=(WC(2)-WC(1))/NX
      CX=(WC(2)+WC(1))/NX

      DY=(WC(4)-WC(3))/NY
      CY=(WC(4)+WC(3))/NY

      TR(1)=-0.5*DX+WC(1)
      TR(2)=DX
      TR(3)=0.
      TR(4)=-0.5*DY+WC(3)
      TR(5)=0.
      TR(6)=DY
      AMIN=1.D20
      AMINP=1.D20
      AMAX=-1.D20
      DO I=1,NX
      DO J=1,NY
        if (A(I,J).GT.0) then
          IF (A(I,J).GT.AMAX) AMAX=A(I,J)
          if (A(I,J).LT.AMIN) AMIN=A(I,J)
          if ((A(I,J).LT.AMINP).and.(A(I,J).GT.0)) AMINP=A(I,J)
        endif
      enddo
      enddo
      if (MINZ==1) then
        AMIN=0.0
      else if (MINZ==2) then
        AMIN=AMINP
      endif

      IF (LOGZ.LE.0.0) THEN
        CALL PGSVP(DC(1),DC(2),DC(3),DC(4))
        CALL PGSWIN(WC(1),WC(2),WC(3),WC(4))
        CALL PGGRAY(A,NDX,NDY,1,NX,1,NY,AMAX,AMIN,TR)
      ELSE IF ((AMAX.GT.0).AND.(AMIN.LT.AMAX)) THEN
        L0=LOG(10.)
        IF (AMINP.LE.0) AMINP=1.E-20
        RANGE=LOG(AMAX/AMINP)/L0
        IF (RANGE.GT.16.0) THEN
          AMIN=AMAX*exp(-16.0*L0)
          RANGE=16.0
        ENDIF
        IF (RANGE.LT.1.0) THEN
          AMIN=AMAX*exp(-1.0*L0)
          RANGE=1.0
        ENDIF
        DO I=1,NX
        DO J=1,NY
          IF (A(I,J).GT.0) THEN
            AUX(I,J)=LOG(A(I,J)/AMIN)/L0+RANGE/16
          ELSE
            AUX(I,J)=0.D0
          ENDIF
        ENDDO
        ENDDO
        CALL PGSVP(DC(1),DC(2),DC(3),DC(4))
        CALL PGSWIN(WC(1),WC(2),WC(3),WC(4))
        CALL PGGRAY(AUX,MIMAX,MIMAX,1,NX,1,NY,RANGE,0.,TR)
      ENDIF
      ! write(*,1) AMIN,AMINP,AMAX
      END SUBROUTINE PLOT2D

      end MODULE GRFPLOT