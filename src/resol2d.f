!////  $Id: resol2d.f,v 1.9 2019/08/15 15:02:06 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2009, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.9 $
!////     $Date: 2019/08/15 15:02:06 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Module for 2D resolution functions
!////
!//////////////////////////////////////////////////////////////////////
      MODULE RESOL2D
      use RESULTS
      use NSTORE
      USE TRACINGDATA
      use RESOL1D
      use DATASETS
      implicit none
      SAVE
      PRIVATE

      public GETQRANGE2D,DORESOL2D

      contains

!---------------------------------------------------
      subroutine GETQRANGE2D(IX,IY,X0,DX,Y0,DY)
! get range of a resolution variable given by index
!---------------------------------------------------
      integer,intent(in) :: IX,IY
      real(kind(1.D0)),intent(out) :: X0,DX,Y0,DY
      real(kind(1.D0)) :: xmin,xmax,xval,ymin,ymax,yval
      integer :: i,NCNT2,NS2,ND2
      type(NEUTRON) :: NEU1,NEU2
      integer,parameter :: MC=10
      type(PDATASET) :: DT
      DT=DSET_DATA(DSET_ISEL)
      NS2=CSUBSET(2)
      ND2=DSET_ISEL
      NCNT2=NSTORE_MAXN(DT%P%REG2) ! get number of events NCNT
      xmin=1.D20
      xmax=-1.D20
      ymin=1.D20
      ymax=-1.D20
      X0=0.D0
      Y0=0.D0
      DX=10.D0
      DY=10.D0
      IF (NCNT2.GT.MC) THEN
        DO I=1,NCNT2
          !call NSTORE_REFERENCE(I,NS2,ND2,IREF)
          !call NSTORE_GET_KF(I,NS2,ND2,NEU2)
          !call NSTORE_GET_KI2(IREF,NEU1)

          call NSTORE_GETEVENT(DT%P%REG2,i,NEU2)
          call NSTORE_GETEVENT(NEU2%REF(1),NEU2%REF(2),NEU1)
          XVAL=GETQVAR(IX,NEU1,NEU2)
          xmin=min(xmin,XVAL)
          xmax=max(xmax,XVAL)
          YVAL=GETQVAR(IY,NEU1,NEU2)
          ymin=min(ymin,YVAL)
          ymax=max(ymax,YVAL)
        enddo
        X0=(xmin+xmax)/2.D0
        DX=(xmax-xmin)
        Y0=(ymin+ymax)/2.D0
        DY=(ymax-ymin)
      endif
      end subroutine GETQRANGE2D

C----------------------------------------------------
      SUBROUTINE DORESOL2D(APP,IX,IY,X0,Y0,DX,DY,NX,NY,RES)
C get beam profile integrated along vertical coordinate
C----------------------------------------------------
      logical,intent(in) :: APP ! append
      integer,intent(in) :: IX,IY,NX,NY
      REAL(KIND(1.D0)),intent(in) :: X0,Y0,DX,DY
      TYPE(TRES2D) :: RES
      REAL(KIND(1.D0)) ::  xmax,xmin,ddx,xval,ymax,ymin,ddy,yval,CNORM
      integer :: I,JX,JY,NCNT2,NCNT1,NNX,NNY !,NS2,ND2
      type(NEUTRON) :: NEU1,NEU2
      type(PDATASET) :: DT
      DT=DSET_DATA(DSET_ISEL)
1     format(a12,': ',6(G10.4,1x),I6)


      NCNT1=NSTORE_MAXN(DT%P%REG1)
      NCNT2=NSTORE_MAXN(DT%P%REG2)
      CNORM=NSTORE_GETNORM(DT%P%REG2)
      if ((NCNT2<=0).or.(NCNT1<=0).or.(CNORM<=0.D0)) then
        call MSG_WARN('Can''t calculate resolution function, no neutrons stored in the scattering channel.',1)
        return
      endif

      !NS2=CSUBSET(2)
      !ND2=DSET_ISEL
      NNX=MIN(NX,MAX2DIM)
      NNY=MIN(NY,MAX2DIM)
  !  write(*,*) 'DORESOL2D NCNT=',NCNT2,IX,IY
      IF (NCNT2.GT.10) THEN
        xmin=X0-DX/2
        xmax=X0+DX/2
        ymin=Y0-DY/2
        ymax=Y0+DY/2
  !    write(*,1) 'DORESOL2D limits =',xmin,xmax,ymin,ymax
        ddx=(xmax-xmin)/NNX
        ddy=(ymax-ymin)/NNY
        if (.not.APP) RES%Z=0.D0
        DO I=1,NCNT2
          !call NSTORE_REFERENCE(I,NS2,ND2,IREF)
          !call NSTORE_GET_KF(I,NS2,ND2,NEU2)
          !call NSTORE_GET_KI2(IREF,NEU1)
          call NSTORE_GETEVENT(DT%P%REG2,i,NEU2)
          call NSTORE_GETEVENT(NEU2%REF(1),NEU2%REF(2),NEU1)

  !      write(*,1) 'DORESOL2D IREF=',IREF
  !      write(*,1) 'DORESOL2D P=',NEU2%P
  !      write(*,1) 'DORESOL2D NEU1=',NEU1%R,NEU1%K
  !      write(*,1) 'DORESOL2D NEU2=',NEU2%R,NEU2%K
          XVAL=GETQVAR(IX,NEU1,NEU2)
          YVAL=GETQVAR(IY,NEU1,NEU2)
  !      write(*,1) 'DORESOL2D X,Y=',XVAL,YVAL

          JX=INT((XVAL-xmin)/ddx)+1
          JY=INT((YVAL-ymin)/ddy)+1
  !      write(*,*) 'DORESOL2D IX,IY=',JX,JY
  !      READ(*,*)
          IF ((JX.GE.1).AND.(JX.LE.NNX)) THEN
          IF ((JY.GE.1).AND.(JY.LE.NNY)) THEN
            RES%Z(JX,JY)=RES%Z(JX,JY)+CNORM*NEU2%P
          ENDIF
          ENDIF
        enddo
        RES%NX=NNX
        RES%NY=NNY
      else
        RES%NX=0
        RES%NY=0
      ENDIF
      END SUBROUTINE DORESOL2D

      END MODULE RESOL2D
