!//////////////////////////////////////////////////////////////////////
!////  $Id: grfdata.f,v 1.15 2019/08/16 17:16:26 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.15 $
!////     $Date: 2019/08/16 17:16:26 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Module with graphics data
!////  inherited from RESTRAX version, Revision: 1.23
!////
!//////////////////////////////////////////////////////////////////////
      MODULE GRFDATA
      use XMLINFO
      use SIMATH
      implicit none

      integer,parameter :: MIMAX=256    ! Max. & default dimension of image arrays

c array for gauge volume map
      REAL :: SVOL(MIMAX,MIMAX),ARRIMA(MIMAX,MIMAX)

! type describing PGPLOT output device
      integer,parameter :: MAXGDEV=32
      TYPE TGRFDEV; SEQUENCE
        integer :: IORD
        integer :: INTER
        character*8 :: DTYPE
        character*64 :: DESCR
      end type TGRFDEV
      TYPE(TGRFDEV) :: GRFDEV(MAXGDEV)
      character*128 :: GRFFILE='restrax.ps' ! filename for a file device
      integer :: GRFNUM=0 ! number of avaliable devices
      integer :: GRFSEL=0 ! currently selected device

!       integer, PARAMETER :: TOPRINT=0
!       type  VIEWSET ; sequence
!          INTEGER*4 IX,IY
!          REAL*4 DX1,DX2,DY1,DY2    ! Viewport in device coordinates
!          REAL*4 WX1,WX2,WY1,WY2    ! World coordinates
!          CHARACTER*32 XTIT,YTIT,HEAD
!       END type  VIEWSET
!       type(VIEWSET) :: PORT

! type for label information
      type GRFLABEL; SEQUENCE
        character*128 :: TEXT
        real :: DIST,SHIFT,CHSIZE,ANGLE
        integer :: COLOR,LWIDTH
        character*1 :: DUM_ALIGN
        CHARACTER*1 :: SIDE
      END TYPE GRFLABEL

! Viewport definition
      INTEGER,parameter :: VPORT_MAXGRF=255 ! max. number of graphs allocated to the viewport
      type VIEWPORT; sequence
        real(KIND(1.D0)) :: RL2PORT(4,4)  ! transformation matrix rec. lattice -> viewport coordinates
        real(KIND(1.D0)) :: PORT2RL(4,4)  ! inverse to RL2PORT
        REAL :: REF(4)            ! reference point of the viewport coordinate system
        INTEGER :: IX,IY          ! indices to QE(4) vector (for resol. function projections)
        INTEGER :: IGRF(VPORT_MAXGRF)  ! indices of the graphs allocated to the viewport
        INTEGER :: TGRF(VPORT_MAXGRF)  ! types of the graphs allocated to the viewport
        INTEGER :: NGRF           ! number of actually allocated graphs
        INTEGER :: MINZ           ! how to determine minimum z: (0) minimum, (1) zero, (2) minimum positive
        REAL :: DC(4)             ! Device coordinates
        REAL :: WC(4)             ! World coordinates
        real :: LOGSCALE(3)       ! base for log-scale (0=linear)
        REAL :: CHSIZE            ! character size
        integer :: STYLE          ! frame style: default(0),cross(1)
        integer :: LWIDTH         ! line width
        integer :: COLOR          ! color index
        CHARACTER*64 :: XTIT,YTIT,HEAD
        character*128 :: COMMENT
        character*128 :: LEGEND
      END type  VIEWPORT

! defined values of VIEWPORT%TYP
      INTEGER,parameter :: grtyp_none=0
      INTEGER,parameter :: grtyp_1D=1
      INTEGER,parameter :: grtyp_2D=2
      INTEGER,parameter :: grtyp_ell=3
      INTEGER,parameter :: grtyp_txt=4

! Define space for graph data:
      integer, PARAMETER :: MAXPORTS=16   ! number of viewports
      integer, PARAMETER :: MAXPIX=512    ! pixel resolution
      integer, PARAMETER :: MAX2D=4       ! number of 2D plots
      integer, PARAMETER :: MAX1D=32      ! number of curves
      integer, PARAMETER :: MAXELL=256    ! number of ellipsoids
      integer, PARAMETER :: MAXTXT=256    ! number of text rows
      integer, PARAMETER :: MAXPNTS=4096  ! number of points for 1 curve
      integer, PARAMETER :: MAXPEAKS=64   ! number of points for 1 curve

! viewports for standard graphs in RESTRAX:
      type(VIEWPORT) :: GRFPORT(MAXPORTS)

! graph types
  ! ellipsoid projection and section
      TYPE TPLOTELL;SEQUENCE
        REAL :: RMAT(4,4)
        REAL :: RZERO(4)
        INTEGER :: COLOR
        INTEGER :: PLINE,SLINE
      END TYPE TPLOTELL

  ! 1D plot
      TYPE TPLOT1D;SEQUENCE
        REAL :: X(MAXPNTS),Y(MAXPNTS),DY(MAXPNTS)
        REAL :: PKPAR(3*MAXPEAKS)
        REAL :: PKERR(3*MAXPEAKS)
        integer :: NPEAKS
        INTEGER :: COLOR,POINT,LINE
        real :: PSIZE
        integer :: NP
        INTEGER :: ER
        CHARACTER*20 :: X_CAP,Y_CAP
      END TYPE TPLOT1D
  ! 2D map
      TYPE TPLOT2D;SEQUENCE
        REAL :: Z(MAXPIX,MAXPIX)
        integer :: NX,NY
        REAL :: SCALE(4)
      END TYPE TPLOT2D

! define variables for graph data:
  ! 1D plot
      integer :: NGRF1D
      TYPE(TPLOT1D) :: GRF1D(MAX1D)
  ! 2D maps
      integer :: NGRF2D
      TYPE(TPLOT2D) :: GRF2D(MAX2D)
  ! ELLIPSES
      integer :: NGRFELL
      TYPE(TPLOTELL) :: GRFELL(MAXELL)
  ! text
      integer :: NGRFTXT
      TYPE(GRFLABEL) :: GRFTXT(MAXTXT)

      CONTAINS


!----------------------------------------------------------------------
      SUBROUTINE getDeviceStr(IDEV,DSTR)
! Clear 1D graph data arrays
!----------------------------------------------------------------------
      integer,intent(in) :: IDEV
      character(*),intent(out) :: DSTR
        if ((IDEV.gt.0).and.(IDEV.le.GRFNUM)) then
          if (GRFDEV(IDEV)%INTER.gt.0) then
            DSTR=trim(GRFDEV(IDEV)%DTYPE)
          else
            DSTR=trim(GRFFILE)//trim(GRFDEV(IDEV)%DTYPE)
          endif
        else
          DSTR='/NULL'
        endif
      end SUBROUTINE getDeviceStr

!----------------------------------------------------------------------
      SUBROUTINE ClrAllGraphs
! Clear 1D graph data arrays
!----------------------------------------------------------------------
      INTEGER :: ig
      do ig=1,MAX1D
        call ClrData(ig)
      enddo
      do ig=1,MAX2D
        call ClrData2D(ig)
      enddo
      do ig=1,MAXELL
        call ClrEll(ig)
      enddo
      do ig=1,MAXTXT
        call ClrTxt(ig)
      enddo
      do ig=1,MAXPORTS
        call ClrViewPort(ig)
      enddo
      NGRF1D=0
      NGRF2D=0
      NGRFELL=0
      NGRFTXT=0

      end SUBROUTINE ClrAllGraphs

!----------------------------------------------------------------------
      SUBROUTINE ClrData2D(ig)
! Clear 1D graph data arrays
!----------------------------------------------------------------------
      INTEGER, intent(in) :: ig
      GRF2D(ig)%NX=0
      GRF2D(ig)%NY=0
      GRF2D(ig)%Z(1:MAXPIX,1:MAXPIX)=0.0
      end SUBROUTINE ClrData2D

!----------------------------------------------------------------------
      SUBROUTINE ClrData(ig)
! Clear 1D graph data arrays
!----------------------------------------------------------------------
      INTEGER, intent(in) :: ig
      GRF1D(ig)%NP=0
      GRF1D(ig)%X(1:MAXPNTS)=0.0
      GRF1D(ig)%Y(1:MAXPNTS)=0.0
      GRF1D(ig)%PSIZE=1.0
      GRF1D(ig)%COLOR=1
      GRF1D(ig)%POINT=1
      GRF1D(ig)%LINE=1
      GRF1D(ig)%X_CAP=' '
      GRF1D(ig)%Y_CAP=' '
      GRF1D(ig)%NPEAKS=0
      GRF1D(ig)%PKPAR=0.0
      end SUBROUTINE ClrData

!----------------------------------------------------------------------
      SUBROUTINE ClrEll(ig)
! Clear 1D graph data arrays
!----------------------------------------------------------------------
      INTEGER, intent(in) :: ig
      GRFELL(ig)%RMAT=0.0
      GRFELL(ig)%RZERO=0.0
      GRFELL(ig)%COLOR=1
      GRFELL(ig)%PLINE=1
      GRFELL(ig)%SLINE=2
      end SUBROUTINE ClrEll

!----------------------------------------------------------------------
      SUBROUTINE ClrTxt(ig)
! Clear text data arrays
!----------------------------------------------------------------------
      INTEGER, intent(in) :: ig
        GRFTXT(ig)%TEXT=' '
        GRFTXT(ig)%DIST=-1.5
        GRFTXT(ig)%SHIFT=0.05
        GRFTXT(ig)%CHSIZE=1.0
        GRFTXT(ig)%ANGLE=0.0
        GRFTXT(ig)%COLOR=1
        GRFTXT(ig)%LWIDTH=1
        GRFTXT(ig)%SIDE='T'
      end SUBROUTINE ClrTxt

!----------------------------------------------------------------------
      SUBROUTINE ClrViewPort(iport)
! Clear viewport for future plots
!----------------------------------------------------------------------
      INTEGER, intent(in) :: iport
      GRFPORT(iport)%TGRF=grtyp_none
      GRFPORT(iport)%IGRF=0
      GRFPORT(iport)%NGRF=0
      GRFPORT(iport)%DC=0.0
      GRFPORT(iport)%WC=0.0
      GRFPORT(iport)%REF=0.0
      GRFPORT(iport)%XTIT=' '
      GRFPORT(iport)%YTIT=' '
      GRFPORT(iport)%HEAD=' '
      GRFPORT(iport)%COMMENT=' '
      GRFPORT(iport)%LEGEND=' '
      GRFPORT(iport)%IX=1
      GRFPORT(iport)%IY=2
      GRFPORT(iport)%MINZ=0
      GRFPORT(iport)%LOGSCALE=0.0
      GRFPORT(iport)%STYLE=0
      GRFPORT(iport)%COLOR=1
      GRFPORT(iport)%LWIDTH=2
      GRFPORT(iport)%CHSIZE=1.2
      call UNIMAT(GRFPORT(iport)%RL2PORT,4,4)
      call UNIMAT(GRFPORT(iport)%PORT2RL,4,4)

      end SUBROUTINE ClrViewPort

!----------------------------------------------------------------------
      SUBROUTINE CollectGraphs(iport,itype,maxg,ig,ng)
! Collect indices of all graphs of given type in given viewport
!----------------------------------------------------------------------
      integer,intent(in) :: iport,itype,maxg
      integer,intent(out) :: ig(maxg),ng
      integer :: ngraphs,igrf
! get total number of graphs assigned to the viewport
      ngraphs=GRFPORT(iport)%NGRF
      ng=0
      do igrf=1,ngraphs
        if (GRFPORT(iport)%TGRF(igrf).eq.itype) then
          ng=ng+1
          ig(ng)=GRFPORT(iport)%IGRF(igrf)
        endif
      endDO
      end SUBROUTINE CollectGraphs

!----------------------------------------------------------------------
      INTEGER FUNCTION getNextGraph(itype)
! Return the index of next non-allocated graph of given type
!----------------------------------------------------------------------
      INTEGER,intent(in) :: itype
      INTEGER :: ig
      ig=0
      select case(itype)
      case(grtyp_1D)
        ig=min(MAX1D,NGRF1D+1)
      case(grtyp_2D)
        ig=min(MAX2D,NGRF2D+1)
      case(grtyp_ell)
        ig=min(MAXELL,NGRFELL+1)
      case(grtyp_txt)
        ig=min(MAXTXT,NGRFTXT+1)
      end select
      getNextGraph=ig
      end FUNCTION getNextGraph

!--------------------------------------------------------------------
      SUBROUTINE AssignGraph(igrf,iport,itype)
! Assigns a graph to viewport.
! Also updates the maximum number of allocated graphs of given type
! igrf ... index of the graph in corresponding array (GRF1D, GRF2D or GRFELL)
! iport ... index of the viewport to be allocate to
! itype ... graph type: 1D (1) 2D (2) or ellipse (3)
!--------------------------------------------------------------------
      INTEGER,intent(in) :: igrf,iport,itype
      INTEGER :: ng
      logical :: canAssign
3     format('Assign graph ',I2,' to port ',I2,' at position ',I2,' typ ',I2)
1     format(a8,': ',6(G10.4,1x))
        canAssign=.FALSE.
        if ((iport.gt.0).and.(iport.le.MAXPORTS)) then
          if (GRFPORT(iport)%NGRF.LT.VPORT_MAXGRF) then
            select case(itype)
            case(grtyp_1D)
              canAssign=(igrf.le.MAX1D)
            case(grtyp_2D)
              canAssign=(igrf.le.MAX2D)
            case(grtyp_ell)
              canAssign=(igrf.le.MAXELL)
            case(grtyp_txt)
              canAssign=(igrf.le.MAXTXT)
            end select
          endif
        endif
        if (canAssign) then
          ng=GRFPORT(iport)%NGRF+1
          GRFPORT(iport)%NGRF=ng
          GRFPORT(iport)%IGRF(ng)=igrf
          GRFPORT(iport)%TGRF(ng)=itype
          select case(itype)
          case(grtyp_1D)
            NGRF1D=max(NGRF1D,igrf)
          case(grtyp_2D)
            NGRF2D=max(NGRF2D,igrf)
          case(grtyp_ell)
            NGRFELL=max(NGRFELL,igrf)
      !      write(*,3) igrf,iport,ng,GRFPORT(iport)%TGRF(ng)
      !      write(*,1) 'center',GRFELL(igrf)%RZERO
          case(grtyp_txt)
            NGRFTXT=max(NGRFTXT,igrf)
          end select
        endif
      END SUBROUTINE AssignGraph


!----------------------------------------------------------------------
      SUBROUTINE WriteGrfData2D(iunit,iport)
! Writes all 2D graph data assigned to GRFPORT(iport) to given unit
! writes warning message if iwarn>0
!----------------------------------------------------------------------
      INTEGER, intent(in) :: iunit,iport
      integer :: i,j,k,igrf,jmax
      integer :: ig(MAX2D),ng
      character*256 :: MSG
      character*32 :: CNUM(4)
1     FORMAT(a)
10    FORMAT(G10.4,' ',$)
11    format('dimension (col,row): (',a,',',a,')')
13    format('scale (col,row): ((',a,',',a,'),(',a,',',a,'))')

      call CollectGraphs(iport,grtyp_2D,MAX2D,ig,ng)
      if (ng.le.0) RETURN

!c// get max. number of points
      jmax=0
      do i=1,ng
        jmax=MAX(jmax,GRF2D(ig(i))%NX*GRF2D(ig(i))%NY)
      enddo
      if (jmax.le.0) RETURN

!c// write all 2D graph data sequentially
      do k=1,ng
        igrf=ig(k)
        call int2str(GRF2D(igrf)%NX,CNUM(1))
        call int2str(GRF2D(igrf)%NY,CNUM(2))
  ! write header info
        write(iunit,11,err=998) trim(CNUM(1)),trim(CNUM(2))
        do i=1,4
          call float2str(1.D0*GRFPORT(iport)%WC(i),CNUM(i))
        enddo
        WRITE(iunit,13,err=998) (trim(CNUM(i)),i=1,4)
        write(iunit,1) 'DATA_2D:'
  ! write data table
        do j=1,GRF2D(igrf)%NY
          do i=1,GRF2D(igrf)%NX
            write(iunit,10,err=998) GRF2D(igrf)%Z(i,j)
          enddo
          write(iunit,*)
        enddo
      enddo
      return

998   write(MSG,*) 'Cannot write graph data to unit ',iunit
      call MSG_ERROR('WriteGrfData',trim(MSG),1,1)
      END SUBROUTINE WriteGrfData2D

!----------------------------------------------------------------------
      SUBROUTINE WriteGrfData(iunit,iport)
! Writes all 1D graph data assigned to GRFPORT(iport) to given unit
! writes warning message if iwarn>0
!----------------------------------------------------------------------
      INTEGER, intent(in) :: iunit,iport
      integer :: i,j,k,jmax,igrf
      integer :: ig(MAX1D),ng
      character*256 :: MSG
      character*32 :: CNUM,XLAB,YLAB

6     FORMAT('"',a,'"  ',$)
10    FORMAT(2(G13.5,2x),$)
11    FORMAT(G13.5,2x,$)
12    FORMAT(I3,2x,$)
13    FORMAT(G13.5,2x,G13.5,2x,$)

      call CollectGraphs(iport,grtyp_1D,MAX1D,ig,ng)
      if (ng.le.0) RETURN

!c// get max. number of points
      jmax=0
      do i=1,ng
        jmax=MAX(jmax,GRF1D(ig(i))%NP)
      enddo
      if (jmax.le.0) RETURN

! write peak parameters
      do k=1,ng
        igrf=ig(k)
        !write(*,*) 'GRFDATA:  ',igrf,GRF1D(igrf)%NPEAKS,GRF1D(igrf)%PKPAR(1:3)
        if (GRF1D(igrf)%NPEAKS>0) then
          write(iunit,*) 'PEAKPARAM:'
          write(iunit,*) 'INTENSITY  eINT  WIDTH  eWIDTH  POSITION  ePOS'
          do i=1,GRF1D(igrf)%NPEAKS
            j=3*(i-1)
            write(iunit,12,err=998) i
            write(iunit,13,err=998) GRF1D(igrf)%PKPAR(j+1),GRF1D(igrf)%PKERR(j+1)
            write(iunit,13,err=998) GRF1D(igrf)%PKPAR(j+3),GRF1D(igrf)%PKERR(j+3)
            write(iunit,13,err=998) GRF1D(igrf)%PKPAR(j+2),GRF1D(igrf)%PKERR(j+2)
            write(iunit,*,err=998)
          enddo
        endif
      endDO
!c// write all 1D graph data in one table
  ! write header line
      write(iunit,*) 'DATA_1D:'
      do k=1,ng
        igrf=ig(k)
        CNUM=' '
        if (ng.gt.1) call INT2STR(k,CNUM)
        XLAB='X'//trim(CNUM)
        YLAB='Y'//trim(CNUM)
        write(iunit,6,err=998) trim(XLAB)
        write(iunit,6,err=998) trim(YLAB)
        if (GRF1D(igrf)%ER.gt.0) then
          YLAB='dY'//trim(CNUM)
          write(iunit,6,err=998) trim(YLAB)
        endif
      enddo
      write(iunit,*)
  ! write data table
      do j=1,jmax
        do k=1,ng
          igrf=ig(k)
          if (GRF1D(igrf)%NP.GE.j) then
            write(iunit,10,err=998) GRF1D(igrf)%X(j),GRF1D(igrf)%Y(j)
            if (GRF1D(igrf)%ER.gt.0) write(iunit,11,err=998) GRF1D(igrf)%dY(j)
          endif
        enddo
        write(iunit,*)
      enddo
      return

998   write(MSG,*) 'Cannot write graph data to unit ',iunit
      call MSG_ERROR('WriteGrfData',trim(MSG),1,1)
      END SUBROUTINE WriteGrfData

!----------------------------------------------------------------------
      SUBROUTINE WriteGrfEll(iunit,iport)
! Writes all ellipses assigned to GRFPORT(iport) to given unit
! writes warning message if iwarn>0
!----------------------------------------------------------------------
      INTEGER, intent(in) :: iunit,iport
      integer :: k,igrf
      integer :: ig(MAXELL),ng,IX,IY
      REAL :: PROJ(2,2),SECT(2,2),ZERO(2)
      character*256 :: MSG
      character*32 :: CNUM

1     FORMAT('ELLIPSE  X0  Y0  proj_XX  proj_YX  proj_XY  proj_YY sec_XX  sec_YX  sec_XY  sec_YY')
10    FORMAT(a,' ',10(G13.5,2x))

      call CollectGraphs(iport,grtyp_ell,MAXELL,ig,ng)
      if (ng.le.0) RETURN
      IX=GRFPORT(iport)%IX
      IY=GRFPORT(iport)%IY
!c// write all graph data in one table
  ! write header line
      write(iunit,*) 'DATA_ELL: (center position and matrix elements)'
      write(iunit,1)
  ! write data table
      do k=1,ng
        igrf=ig(k)
        call EllipsoidView4(GRFELL(igrf)%RMAT(1,1),IX,IY,PROJ,SECT)
        call INT2STR(k,CNUM)
        ZERO(1)=GRFELL(igrf)%RZERO(IX)
        ZERO(2)=GRFELL(igrf)%RZERO(IY)
        write(iunit,10,err=998) trim(CNUM),ZERO,PROJ,SECT
      enddo
      return

998   write(MSG,*) 'Cannot write graph data to unit ',iunit
      call MSG_ERROR('WriteGrfEll',trim(MSG),1,1)
      END SUBROUTINE WriteGrfEll


!----------------------------------------------------------------------
      SUBROUTINE WriteGrfTxt(iunit,iport)
! Writes all text data assigned to GRFPORT(iport) to given unit
!----------------------------------------------------------------------
      INTEGER, intent(in) :: iunit,iport
      integer :: igrf,k
      integer :: ig(MAXTXT),ng
      character*32 :: CNUM
6     FORMAT('LABEL(',a,') : "',a,'"')
      call CollectGraphs(iport,grtyp_txt,MAXTXT,ig,ng)
      do k=1,ng
        igrf=ig(k)
        if (len_trim(GRFTXT(igrf)%text).gt.0) then
          call INT2STR(k,CNUM)
          write(iunit,6) trim(CNUM),trim(GRFTXT(igrf)%text)
        endif
      enddo
      end SUBROUTINE WriteGrfTxt

!----------------------------------------------------------------------
      SUBROUTINE WritePortData(iunit,iport)
! Writes all data assigned to GRFPORT(iport) to given unit
! Includes header info for the viewport
!----------------------------------------------------------------------
      INTEGER, intent(in) :: iunit,iport
      integer :: i
      character*256 :: MSG
      character*32 CNUM(4)
      REAL :: REF(2)

1     FORMAT('TITLE:    "',a,'"')
2     FORMAT('COMMENT:  "',a,'"')
3     FORMAT('X-TITLE:  "',a,'"')
4     FORMAT('Y-TITLE:  "',a,'"')
31    FORMAT('X-RANGE:  ',a,' ',a)
41    FORMAT('Y-RANGE:  ',a,' ',a)
5     FORMAT('LEGEND :  "',a,'"')
7     FORMAT('LOGSCALE(x,y,z): ',a,' ',a,' ',a)

      REF(1)=GRFPORT(iport)%REF(GRFPORT(iport)%IX)
      REF(2)=GRFPORT(iport)%REF(GRFPORT(iport)%IY)
      do i=1,2
        call float2str(1.D0*(GRFPORT(iport)%WC(i)+REF(1)),CNUM(i))
        call float2str(1.D0*(GRFPORT(iport)%WC(i+2)+REF(2)),CNUM(i+2))
      enddo
      write(iunit,1,err=998) trim(GRFPORT(iport)%HEAD)
      WRITE(iunit,2,err=998) trim(GRFPORT(iport)%COMMENT)
      write(iunit,3,err=998) trim(GRFPORT(iport)%XTIT)
      write(iunit,4,err=998) trim(GRFPORT(iport)%YTIT)
      write(iunit,31,err=998) trim(CNUM(1)),trim(CNUM(2))
      write(iunit,41,err=998) trim(CNUM(3)),trim(CNUM(4))
      do i=1,3
        call float2str(GRFPORT(iport)%LOGSCALE,CNUM(i))
      enddo
      write(iunit,7,err=998) (trim(CNUM(i)),i=1,3)
      write(iunit,5,err=998) trim(GRFPORT(iport)%LEGEND)
      call WriteGrfTxt(iunit,iport)
      call WriteGrfData(iunit,iport)
      call WriteGrfData2D(iunit,iport)
      call WriteGrfEll(iunit,iport)
      return

998   write(MSG,*) 'Cannot write graph data to unit ',iunit
      call MSG_ERROR('WritePortData',trim(MSG),1,1)
      END SUBROUTINE WritePortData


!--------------------------------------------------------------------
      SUBROUTINE SetYscale(ip)
! set automatic y-scale for GRFPORT(ip)
!--------------------------------------------------------------------
      INTEGER,intent(in) :: ip
      real :: YMIN,YMAX,rang,centre
      integer :: i,ig,ng,IY
      real, parameter :: C2LN2=1.38629436112 ! 2*ln(2)
      real :: AUX(4,4)

      GRFPORT(ip)%WC(3)=1.E30
      GRFPORT(ip)%WC(4)=-1.E30
      IY=GRFPORT(ip)%IY
      do ng=1,GRFPORT(ip)%NGRF
        ig=GRFPORT(ip)%IGRF(ng)
        YMAX=-1.E30
        YMIN=1.E30
        select case (GRFPORT(ip)%TGRF(ng))
        case(grtyp_1D)
          DO I=1,GRF1D(ig)%NP
            YMAX=MAX(YMAX,1.0*GRF1D(ig)%Y(i))
            YMIN=MIN(YMIN,1.0*GRF1D(ig)%Y(i))
          ENDDO
          if (YMIN.gt.0) then ! stick to zero
            if (YMIN.lt.(YMAX-YMIN)*0.1) YMIN=0.0
          endif
        case(grtyp_2D)
          YMIN=GRF2D(ig)%SCALE(1)-GRFPORT(ip)%REF(IY)
          YMAX=GRF2D(ig)%SCALE(2)-GRFPORT(ip)%REF(IY)
        case(grtyp_ell)
          I=INVERT(4,GRFELL(ig)%RMAT,4,AUX,4)
          if (I.EQ.0) then
            rang=SQRT(C2LN2*AUX(IY,IY))
            CENTRE=GRFELL(ig)%RZERO(IY)-GRFPORT(ip)%REF(IY)
            YMIN=CENTRE-rang
            YMAX=CENTRE+rang
          endif
        end select
        if (GRFPORT(ip)%LOGSCALE(2).EQ.0) then
          GRFPORT(ip)%WC(3)=MIN(GRFPORT(ip)%WC(3),YMIN)
          GRFPORT(ip)%WC(4)=MAX(GRFPORT(ip)%WC(4),YMAX+(YMAX-YMIN)*0.1)
        else
          if (YMAX.GT.0) then
            i=INT(LOG10(YMAX*1.2))
            if (YMAX.GT.1) i=i+1
            GRFPORT(ip)%WC(4)=max(GRFPORT(ip)%WC(4),10.0**i)
          else
            GRFPORT(ip)%WC(4)=max(GRFPORT(ip)%WC(4),1.0)
          endif
          if (YMIN.GT.0) then
            i=INT(LOG10(YMIN))
            if (YMIN.LT.1) i=i-1
            GRFPORT(ip)%WC(3)=min(GRFPORT(ip)%WC(3),10.0**i)
          else
            GRFPORT(ip)%WC(3)=min(GRFPORT(ip)%WC(3),ABS(YMAX)/1.E5)
          endif
        endif
      enddo

      end SUBROUTINE SetYscale

!--------------------------------------------------------------------
      SUBROUTINE SetViewPortScale(iport)
! Set physical scale of the viewport automatically if not yet defined
!--------------------------------------------------------------------
      INTEGER,intent(in) :: iport
      INTEGER :: ng,igrf
      logical :: needx,needy
      real :: W1,W2,V1,V2,centre,rang
      INTEGER :: i,j1,j2,jr,IX,IY
      real, parameter :: C2LN2=1.38629436112 ! 2*ln(2)
      real :: AUX(4,4)

      needx=(GRFPORT(iport)%WC(1).GE.GRFPORT(iport)%WC(2))
      needy=(GRFPORT(iport)%WC(3).GE.GRFPORT(iport)%WC(4))
      IX=GRFPORT(iport)%IX
      IY=GRFPORT(iport)%IY
!c set horizontal scale
      if (needx) then
        GRFPORT(iport)%REF(IX)=0.0
        write(*,*) 'Set VPORT scale: ',iport,' title=',trim(GRFPORT(iport)%YTIT)
        W1=+1E+20
        W2=-1E+20
        do ng=1,GRFPORT(iport)%NGRF
          igrf=GRFPORT(iport)%IGRF(ng)
          select case (GRFPORT(iport)%TGRF(ng))
    ! 1D data
          case(grtyp_1D)
!       write(*,*) '1D data hor. ',igrf,' np=',GRF1D(igrf)%NP
            if (GRF1D(igrf)%NP.GT.0) then
              CENTRE=(GRF1D(igrf)%X(1)+GRF1D(igrf)%X(GRF1D(igrf)%NP))/2.-GRFPORT(iport)%REF(IX)
      ! truncate points with zero values
              i=1
              do while ((GRF1D(igrf)%Y(i).eq.0).and.(i.lt.GRF1D(igrf)%NP))
                i=i+1
              enddo
              if(i.gt.1) i=i-1
              j1=i
              i=GRF1D(igrf)%NP
              do while ((GRF1D(igrf)%Y(i).eq.0).and.(i.gt.1))
                i=i-1
              enddo
              if(i.lt.GRF1D(igrf)%NP) i=i+1
              j2=i
      ! rang is the number of non-zero points + 5
              jr=j2-j1+5
      ! rang is at least 20 points
              jr=MAX(jr,20)
              j1=MAX(j1-jr/2,1)
              j2=MIN(j2+jr/2,GRF1D(igrf)%NP)
              rang=ABS(GRF1D(igrf)%X(j2)-GRF1D(igrf)%X(j1))
              W1=MIN(CENTRE-rang/2.,W1)
              W2=MAX(CENTRE+rang/2.,W2)
!         write(*,*) '1D data  hor. W1=',W1,' W2=',W2
            endif
    ! 2D data
          case(grtyp_2D)
!       write(*,*) '2D data  hor. ',igrf,' nx=',GRF2D(igrf)%NX
            if (GRF2D(igrf)%NX.GT.0) then
              W1=MIN(GRF2D(igrf)%SCALE(1)-GRFPORT(iport)%REF(IX),W1)
              W2=MAX(GRF2D(igrf)%SCALE(2)-GRFPORT(iport)%REF(IX),W2)
!       write(*,*) '2D data  hor. W1=',W1,' W2=',W2
            endif
    ! ellipse
          case(grtyp_ell)
            I=INVERT(4,GRFELL(igrf)%RMAT,4,AUX,4)
            if (I.EQ.0) then
              rang=SQRT(C2LN2*AUX(IX,IX))
              CENTRE=GRFELL(igrf)%RZERO(IX)-GRFPORT(iport)%REF(IX)
              W1=MIN(CENTRE-rang,W1)
              W2=MAX(CENTRE+rang,W2)
            endif
          end select
        enddo
        if (W1.LT.W2) then
          GRFPORT(iport)%WC(1)=W1
          GRFPORT(iport)%WC(2)=W2
!           write(*,*) 'New horizontal scale: ',W1,W2
        endif
      endif
!c set vertical scale
      if (needy) then
        W1=+1E+20
        W2=-1E+20
        do ng=1,GRFPORT(iport)%NGRF
          igrf=GRFPORT(iport)%IGRF(ng)
    ! 1D data
          select case (GRFPORT(iport)%TGRF(ng))
          case(grtyp_1D)
!       write(*,*) '1D data ver. ',igrf,' np=',GRF1D(igrf)%NP
            V1=+1E+20
            V2=-1E+20
            if (GRF1D(igrf)%NP.GT.0) then
              DO I=1,GRF1D(igrf)%NP
                V2=MAX(V2,GRF1D(igrf)%Y(i)-GRFPORT(iport)%REF(IY))
                V1=MIN(V1,GRF1D(igrf)%Y(i)-GRFPORT(iport)%REF(IY))
              END DO
              rang=ABS(V2-V1)
              IF (V1.GE.0) then
                if (V1.LT.0.2*rang) then
                  V1=0.0
                else
                  V1=V1-0.05*rang
                endif
              else
                V1=V1-0.1*rang
              endif
              V2=V2+0.15*rang
!       write(*,*) '1D data  ver. V1=',V1,' V2=',V2
            endif
            W1=MIN(W1,V1)
            W2=MAX(W2,V2)
    ! 2D data
          case(grtyp_2D)
            if (GRF2D(igrf)%NY.GT.0) then
              W1=MIN(GRF2D(igrf)%SCALE(3)-GRFPORT(iport)%REF(IY),W1)
              W2=MAX(GRF2D(igrf)%SCALE(4)-GRFPORT(iport)%REF(IY),W2)
!       write(*,*) '2D data  ver. W1=',W1,' W2=',W2
            endif
    ! ellipse
          case(grtyp_ell)
            I=INVERT(4,GRFELL(igrf)%RMAT,4,AUX,4)
            if (I.EQ.0) then
              rang=SQRT(C2LN2*AUX(IY,IY))
              CENTRE=GRFELL(igrf)%RZERO(IY)-GRFPORT(iport)%REF(IY)
              W1=MIN(CENTRE-rang,W1)
              W2=MAX(CENTRE+rang,W2)
            endif
          end select
        enddo
        if (W1.LT.W2) then
          GRFPORT(iport)%WC(3)=W1
          GRFPORT(iport)%WC(4)=W2
!           write(*,*) 'New vertical scale: ',W1,W2
        endif
      endif
      end subroutine SetViewPortScale

      END MODULE GRFDATA
