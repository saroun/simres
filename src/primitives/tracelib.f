!//////////////////////////////////////////////////////////////////////
!////  $Id: tracelib.f,v 1.24 2019/08/15 17:24:08 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2011, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.24 $
!////     $Date: 2019/08/15 17:24:08 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Library of basic ray-tracing subroutines
!////
!////////////////////////////////////////////////////////////////////////
      MODULE TRACELIB
      use CONSTANTS
      implicit none

      REAL(kind(1.D0)), parameter :: TLIB_INF=1.D30
      REAL(kind(1.D0)), parameter :: TLIB_EPS=1.D-8

      logical,private :: dbg=.true.

      contains

!-----------------------------------------------------------------------
      REAL(kind(1.D0)) function MIN_POSITIVE(A,NA)
! return minimum positive value from the array A
! return TLIB_INF if there is no positive value
!-----------------------------------------------------------------------
      real(kind(1.D0)),intent(in) :: A(:)
      integer, intent(in) :: NA
      real(kind(1.D0)) :: t
      integer :: i
      t=TLIB_INF
      do i=1,NA
        if (A(i)>0.D0) t=min(t,A(i))
      enddo
      MIN_POSITIVE=t
      end function MIN_POSITIVE

!-----------------------------------------------------------------------
      REAL(kind(1.D0)) function MAX_NEGATIVE(A,NA)
! return maximum negative value from the array A
! return -TLIB_INF if there is no positive value
!-----------------------------------------------------------------------
      real(kind(1.D0)),intent(in) :: A(:)
      integer, intent(in) :: NA
      real(kind(1.D0)) :: t
      integer :: i
      t=-TLIB_INF
      do i=1,NA
        if (A(i)<0.D0) t=max(t,A(i))
      enddo
      MAX_NEGATIVE=t
      end function MAX_NEGATIVE

!-------------------------------------------------------------------
      SUBROUTINE TLIB_QUAD(EPS,A,B,C,X1,X2)
! Solve quadratic equation A*X^2 + B*X + C = 0
! Sort so that X1<=X2
! Require solution |X|>EPS (avoid zero transport times)
! If there is no solution, return +- TLIB_INF
!-------------------------------------------------------------------
      REAL(kind(1.D0)),intent(in) :: EPS,A,B,C
      REAL(kind(1.D0)),intent(out) :: X1,X2
      REAL(kind(1.D0)) :: DET,X
      integer :: i
! get number of solutions
      i=2
    ! linear equation
      IF (ABS(A).LE.TLIB_EPS**2) THEN
        i=3
        IF (ABS(B).LE.TLIB_EPS) i=0 ! no cross-section
    ! quadratic equation
      else
        DET=B**2-4*A*C
        IF (DET.EQ.0.D0) then
          i=1
        else IF (DET.LT.0.D0) then
          i=0
        endif
      endif
      select case(i)
    ! no sloution
      case(0)
        X1=-TLIB_INF
        X2=TLIB_INF
    ! single solution
      case(1,3)
        if (i.eq.1) then
          X=-B/2.D0/A
        else
          X=-C/B
        endif
        if (X.gt.0.D0) then
          X2=X
          X1=-TLIB_INF
        else
          X1=X
          X2=TLIB_INF
        endif
    ! 2 solutions
      case default
        DET=SQRT(DET)
        X1=(-B+DET)/2.D0/A
        X2=(-B-DET)/2.D0/A
        if (X2.LT.X1) then
          X=X1
          X1=X2
          X2=X
        endif
      end select
      if (abs(X1).lt.EPS) X1=-TLIB_INF
      if (abs(X2).lt.EPS) X2=TLIB_INF
      END SUBROUTINE TLIB_QUAD



!----------------------------------------------------------------------------
      subroutine TLIB_CROSS(kx,kz,x,z,d,gx,t1,t2)
! Calculate time of flight to a cross-point with a parabolic surface
! in the form x=d(1)*z**2 + d(2)*z + d(3)
! The results are ordered so that t1<=t2
! If there is no solution, return +- TLIB_INF
! x,z    .. starting point
! kx,kz  .. velocity components
! gx     .. accelleration in x-direction
!----------------------------------------------------------------------------
      REAL(kind(1.D0)),intent(in) :: kx,kz,x,z,d(3),gx
      REAL(kind(1.D0)),intent(out) :: t1,t2
      REAL(kind(1.D0)) :: aa,bb,cc
      aa=d(1)*kz**2 - 0.5D0*gx
      bb=(2.D0*d(1)*z + d(2))*kz - kx
      cc=d(1)*z**2 + d(2)*z - x + d(3)
      CALL TLIB_QUAD(TLIB_EPS,aa,bb,cc,t1,t2)
      END subroutine TLIB_CROSS

!----------------------------------------------------------------------------
      subroutine TLIB_DBG(arg)
!----------------------------------------------------------------------------
      logical,intent(in)  :: arg
        dbg=arg
      end subroutine TLIB_DBG

!----------------------------------------------------------------------------
      subroutine TLIB_CROSS3(kx0,kz0,x0,z0,tmin,tmax,d,gx,tol,tc)
! Calculate time of flight to a cross-point with a curved surface, in constant acceleration field
! Surface profile function is a polynom, x=d(1)*z**3 + d(2)*z**2 + d(3)*z + d(4)
! Use iteration method with parabolic approximation.
! Touch condition is limited by zmin+kz*TOL<z<zmax
! If there is no solution, return -TLIB_INF
! x0,z0    .. starting point
! kx0,kz0  .. k-vector components
! tmin,tmax .. time limits
! d(4)   .. wall shape parameters (cubic form)
! gx     .. accelleration in x-direction
! tol    .. touch tolerance [mm]
! tc time-of-flight to the cross point, or TLIB_INF if there is no touch
!----------------------------------------------------------------------------
      REAL(kind(1.D0)),intent(in) :: kx0,kz0,x0,z0,tmin,tmax,d(4),gx,tol
      REAL(kind(1.D0)),intent(out) :: tc
      REAL(kind(1.D0)) :: kx,x,z,dt,dif,xw,tsum,ttol,t1,t2
      logical :: touch
      integer :: iter
      integer,parameter :: imax=10
1     format(a,8(1x,G12.5))

      if (dbg) then
        write(*,1) 'TLIB_CROSS3 ',kx0,z0,x0,d(1:4)
      endif

      tc=-TLIB_INF
      x=x0
      z=z0
      kx=kx0
      ! wall position
      xw=d(1)*z**3+d(2)*z**2+d(3)*z+d(4)
      ! x-distance of neutron from the wall position
      dif=x0-xw
      !sg=sign(1.D0,dif)
      !sgt=1.D0
      dt=-TLIB_INF
      touch=.false.
      tsum=0.D0
      iter=0
      ttol=tol/abs(kz0)
      do while((iter<imax).and.(.not.touch))
        iter=iter+1
        ! 1st approximation - assume quadratic model
        call TLIB_CROSS(kx,kz0,x,z,d(2:4),gx,t1,t2)
        if (dbg) then
          write(*,1) 'TLIB_CROSS3 iter=',iter,xw,z,x,t1,t2
        endif
        ! exclude solutions outside limits
        if ((t1+tsum<ttol+tmin).or.(t1+tsum>tmax)) t1=TLIB_INF
        if ((t2+tsum<ttol+tmin).or.(t2+tsum>tmax)) t2=TLIB_INF
        if (abs(t1)<abs(t2)) then
          dt=t1
        else
          dt=t2
        endif
      ! no touch - exit
        if (dt+ttol>=TLIB_INF) return
      ! shift by dt
        tsum=tsum+dt
        z=z0+kz0*tsum
        x=x0+kx0*tsum+0.5D0*gx*tsum**2
        kx=kx0+gx*tsum
      ! wall distance at the new z-value
        xw=d(1)*z**3+d(2)*z**2+d(3)*z+d(4)
        dif=x-xw
        if (dbg) then
          write(*,1) 'TLIB_CROSS3 dif=',xw,z,x,dif,dt
        endif
        touch=(abs(dif)<ttol)
      enddo
      if (touch) tc=tsum
      END subroutine TLIB_CROSS3


!-------------------------------------------------------------------
      SUBROUTINE TLIB_ROT_POS(R0,K0,FRQ,PHI,IR,RC,IC,T,R,K)
! Get position in rotating frame. Rotation axis // ir-th axis, displaced at -RC along ic-th axis.
! No input validity check !!
! R0(3) ... position at time=T
! K0(3) ... propagation vector at time=T
! RC    ... rotation center distance
! IC    ... direction of rotation center (1..3)
! FRQ   ... frequency [MHz]
! PHI   ... phase [2PI], i.e. as a fraction of the period
! IR    ... rotation axis index (1..3)
! T     ... time of flight [in units of m/hbar, as always in SIMRES)
! R ... resulting position at time=T
! K ... resulting propagation vector at time=T
! assume phase=PHI at T=0
!-------------------------------------------------------------------
      REAL(kind(1.D0)),intent(in) :: R0(3),K0(3),RC,FRQ,PHI,T
      integer,intent(in) :: IR,IC
      REAL(kind(1.D0)),intent(out) :: R(3),K(3)
      integer :: IAX(6) = (/2,3,3,1,1,2/) ! axes permutations
      REAL(kind(1.D0)) :: R1(3),SS,CC,PHASE,PH
      integer :: IX,IY
      PH=FRQ*T/HOVM+PHI
      PHASE=TWOPI*PH
      CC=COS(PHASE)
      SS=SIN(PHASE)
      IX=IAX((IR-1)*2+1)
      IY=IAX((IR-1)*2+2)
      R1=R0
      if (RC.NE.0.D0) R1(IC)=R1(IC)+RC
      R(IX)=R1(IX)*CC + R1(IY)*SS
      R(IY)=R1(IY)*CC - R1(IX)*SS
      R(IR)=R1(IR)
      K(IX)=K0(IX)*CC + K0(IY)*SS
      K(IY)=K0(IY)*CC - K0(IX)*SS
      K(IR)=K0(IR)
      if (RC.NE.0.D0) R(IC)=R(IC)-RC
      end SUBROUTINE TLIB_ROT_POS



!-------------------------------------------------------------------
      SUBROUTINE TLIB_ROT_POSMV(R0,K0,FRQ,PHI,IR,RC,IC,T,R,K)
! R0(3) ... position at time=0
! K0(3) ... propagation vector at time=0
! Move by time T and convert to rotating frame
!-------------------------------------------------------------------
      REAL(kind(1.D0)),intent(in) :: R0(3),K0(3),RC,FRQ,PHI,T
      integer,intent(in) :: IR,IC
      REAL(kind(1.D0)),intent(out) :: R(3),K(3)
      REAL(kind(1.D0)) :: R1(3)
      integer :: i
      do i=1,3
        R1(i)=R0(i)+K0(i)*T
      enddo
      call TLIB_ROT_POS(R1,K0,FRQ,PHI,IR,RC,IC,T,R,K)
      end SUBROUTINE TLIB_ROT_POSMV


!-------------------------------------------------------------------
      REAL(KIND(1.D0)) FUNCTION TLIB_BLADE_DIST(R,IR,rho,alpha,delta)
! Get distance from the point R to a curved surface defined as
! y = 0.5*rho*x**2 + alpha*x + delta, when IR=3
! the distance is measured along y.
! Permutation of axes is possible by setting different IR value.
! e.g. for a Fermi chopper, we would use IR=2, so that the above equation changes to
! x = 0.5*rho*z**2 + alpha*z + delta, and the distance is measured along x
! Arguments:
! R     ... current neutron position
! IR    ... rotation axis index (1..3)
! rho,alpha,delta ... blade profile parameters
!-------------------------------------------------------------------
      REAL(kind(1.D0)),intent(in) :: R(3),rho,alpha,delta
      integer,intent(in) :: IR
      integer :: IAX(6) = (/2,3,3,1,1,2/) ! axes permutations
      integer :: IX,IY
      IX=IAX((IR-1)*2+1)
      IY=IAX((IR-1)*2+2)
      TLIB_BLADE_DIST=R(IY)-(0.5D0*rho*R(IX)**2+alpha*R(IX)+delta)
      end FUNCTION TLIB_BLADE_DIST


!-------------------------------------------------------------------
      SUBROUTINE TLIB_ROT_CROSS_APX(what,R0,K0,FRQ,IR,rho,alpha,delta,T)
! Get time of cross-section with a curved plane defined as
! y = 0.5*rho*x**2 + alpha*x + delta
! rotation axis is z, when IR=3
! Permutation of axes is possible by setting different IR value.
! e.g. for a Fermi chopper, we would use IR=2, so that the above equation changes to
! x = 0.5*rho*z**2 + alpha*z + delta, while beam passes along z
! what=-1  ... return largest negative solution, if any
! what=+1  ... return the smallest positive soltion, if any
! what=0   ... return the solution with the smallest absolute value
! APPROXIMATION: use quadratic expansion for T
! No input validity check !!
! R0(3) ... position at time=0
! K0(3) ... propagation vector at time=0
! FRQ   ... frequency [MHz]
! IR    ... rotation axis index (1..3)
! rho,alpha,delta ... blade profile parameters
! T     ... time of flight to the cross point in [m/hbar]
!-------------------------------------------------------------------
      integer,intent(in) :: what,IR
      REAL(kind(1.D0)),intent(in) :: R0(3),K0(3),FRQ,rho,alpha,delta
      REAL(kind(1.D0)),intent(out) :: T
      integer :: IAX(6) = (/2,3,3,1,1,2/) ! axes permutations
      REAL(kind(1.D0)) :: P0,P1,P2,T1,T2,OM
      integer :: IX,IY
! note: time in [usec] = T/HOVM => we have to re-scale frequency
      OM=FRQ*TWOPI/HOVM
      IX=IAX((IR-1)*2+1)
      IY=IAX((IR-1)*2+2)
      P0=delta+alpha*R0(IX)-R0(IY)
      P1=alpha*(K0(IX)+OM*R0(IY))-K0(IY)+OM*R0(IX)
      P2=(alpha*R0(IX)-R0(IY))*OM**2 + (K0(IX)+alpha*K0(IY))*OM
      if (rho.ne.0.D0) then
        P0=P0+0.5D0*rho*R0(IX)**2
        P1=P1+R0(IX)*rho*(K0(IX)+OM*R0(IY))
        P2=P2+rho*OM**2*(R0(IX)**2+0.5D0*R0(IY)**2)
        P2=P2+rho*(R0(IY)*K0(IX)+R0(IX)*K0(IY) + 0.5D0*K0(IX)**2)
      endif
      select case(what)
      case(-1)
        CALL TLIB_QUAD(TLIB_EPS,P2,P1,P0,t1,t2)
        if (t2.le.0.D0) then
          T=max(t1,t2)
        else if (t1.le.0.D0) then
          T=t1
        else
          T=-TLIB_INF
        endif
      case(1)
        CALL TLIB_QUAD(TLIB_EPS,P2,P1,P0,t1,t2)
        if (t1.gt.0.D0) then
          T=min(t1,t2)
        else if (t2.gt.0.D0) then
          T=t2
        else
          T=TLIB_INF
        endif
      case default
        CALL TLIB_QUAD(0.D0,P2,P1,P0,t1,t2)
        if (abs(t1).lt.abs(t2)) then
          T=T1
        else
          T=T2
        endif
      end select
      end SUBROUTINE TLIB_ROT_CROSS_APX

!-------------------------------------------------------------------
      SUBROUTINE TLIB_ROT_CROSS(what,R0,K0,FRQ,IR,rho,alpha,delta,T,iter)
! Get time of cross-section with a curved plane defined as
! y = 0.5*rho*x**2 + alpha*x + delta
! rotation axis is z, when IR=3
! Permutation of axes is possible by setting different IR value.
! e.g. for a Fermi chopper, we would use IR=2, so that the above equation changes to
! x = 0.5*rho*z**2 + alpha*z + delta, while beam passes along z
! ITERATIVE method
! No input validity check !!
! what=-1  ... return largest negative solution, if any
! what=+1  ... return the smallest positive soltion, if any
! R0(3) ... position at time=0
! K0(3) ... propagation vector at time=0
! FRQ   ... frequency [Hz]
! IR    ... rotation axis index (1..3)
! rho,alpha,delta ... blade profile parameters
! T ... time of flight to the cross point
!-------------------------------------------------------------------
      integer,intent(in) :: what,IR
      REAL(kind(1.D0)),intent(in) :: R0(3),K0(3),FRQ,rho,alpha,delta
      REAL(kind(1.D0)),intent(out) :: T
      integer,intent(out) :: iter
      integer :: IAX(6) = (/2,3,3,1,1,2/) ! axes permutations
      REAL(kind(1.D0)) :: t1,t2,EPS,R(3),K(3),R1(3),K1(3),dd
      integer,parameter :: MAXITER=10
      logical :: loop
    ! precission ~ 10^-6 mm
      EPS=1.D-6
      iter=0
      T=0.D0
      call TLIB_ROT_CROSS_APX(what,R0,K0,FRQ,IR,rho,alpha,delta,t1)
      if (abs(t1).lt.TLIB_INF) then
        call TLIB_ROT_POSMV(R0,K0,FRQ,0.D0,IR,0.D0,3,t1,R,K)
        dd=TLIB_BLADE_DIST(R,IR,rho,alpha,delta)
        iter=1
        loop=(abs(dd).gt.EPS)
        do while(loop)
          iter=iter+1
          call TLIB_ROT_CROSS_APX(0,R,K,FRQ,IR,rho,alpha,delta,t2)
          loop=(abs(t2).lt.TLIB_INF)
          if (loop) then
            t1=t1+t2
            call TLIB_ROT_POSMV(R,K,FRQ,0.D0,IR,0.D0,3,t2,R1,K1)
            dd=TLIB_BLADE_DIST(R1,IR,rho,alpha,delta)
          endif
          loop=(loop.and.(abs(dd).gt.EPS))
          loop=(loop.and.(iter.lt.MAXITER))
        enddo
      endif
      T=t1
      end SUBROUTINE TLIB_ROT_CROSS


!-------------------------------------------
      LOGICAL FUNCTION TLIB_INFRAME(R,FSIZE,FSHAPE)
! Returns .TRUE. if R(3) is inside the OBJ
!-------------------------------------------
      use CLASSES
      real(kind(1.D0)),intent(in) :: R(3),FSIZE(3)
      integer,intent(in) :: FSHAPE
      LOGICAL :: LOG1
      select case (FSHAPE)
      case(FRAME_SHAPE_BOX)
        LOG1=(
     1     (ABS(R(1)).LT.FSIZE(1)/2.D0).AND.
     2     (ABS(R(2)).LT.FSIZE(2)/2.D0).AND.
     3     (ABS(R(3)).LT.FSIZE(3)/2.D0))
      case(FRAME_SHAPE_DISC)
        LOG1=(
     1     (((R(1)*2.D0/FSIZE(1))**2+
     2     (R(2)*2./FSIZE(2))**2).LT.1.D0).AND.
     3     (ABS(R(3)).LT.FSIZE(3)/2.D0))
      case(FRAME_SHAPE_CYLINDER)
        LOG1=(
     1      (((R(1)*2.D0/FSIZE(1))**2+
     2     (R(3)*2.D0/FSIZE(3))**2).LT.1.D0).AND.
     3     (ABS(R(2)).LT.FSIZE(2)/2.D0))
      case(FRAME_SHAPE_ELLIPSOID)
        LOG1=(
     1     ((R(1)*2.D0/FSIZE(1))**2+
     2      (R(2)*2.D0/FSIZE(2))**2+
     3      (R(3)*2.D0/FSIZE(3))**2).LT.1.D0)
      case DEFAULT
        LOG1=.FALSE.
      END select
      TLIB_INFRAME=LOG1
      END FUNCTION TLIB_INFRAME


!-------------------------------------------------
      LOGICAL FUNCTION TLIB_INSIDE(R,FSIZE,FSHAPE)
! Returns .TRUE.  if R(2) is inside the frame
! defined by the FSIZE and FSHAPE parameters
!-------------------------------------------------
      use CLASSES
      real(kind(1.D0)),intent(in) :: R(2),FSIZE(2)
      integer,intent(in) :: FSHAPE
      LOGICAL :: LOG1
      select case (FSHAPE)
      case(FRAME_SHAPE_DISC,FRAME_SHAPE_ELLIPSOID)
         LOG1=(((R(1)*2./FSIZE(1))**2+(R(2)*2./FSIZE(2))**2).LT.1)
      case DEFAULT
         LOG1=((ABS(R(1)).LT.FSIZE(1)/2.).AND.(ABS(R(2)).LT.FSIZE(2)/2.))
      END select
      TLIB_INSIDE=LOG1
      END FUNCTION TLIB_INSIDE

!-------------------------------------------------------
      LOGICAL FUNCTION TLIB_PASS_FRAME(R,K,FSIZE,FSHAPE)
! Returns .TRUE. if neutron passes through both entry
! and exit windows (normal to z-axis) of the frame defined
! by the FSIZE and FSHAPE parameters.
!-------------------------------------------------------
      real(kind(1.D0)),intent(in) :: R(3),K(3),FSIZE(3)
      integer,intent(in) :: FSHAPE
      LOGICAL :: LOG1
      real(kind(1.D0)) :: t1,t2,X(2)
      integer :: i
      call BORDER_WALL(R(3),K(3),FSIZE(3),t1,t2)
      LOG1=(t1.lt.t2)
      if (LOG1) then
        do i=1,2
          X(i)=R(i)+K(i)*t1
        enddo
        LOG1=TLIB_INSIDE(X,FSIZE(1:2),FSHAPE)
        if (LOG1) then
          do i=1,2
            X(i)=R(i)+K(i)*t2
          enddo
          LOG1=TLIB_INSIDE(X,FSIZE(1:2),FSHAPE)
        endif
      endif
      TLIB_PASS_FRAME=LOG1
      END FUNCTION TLIB_PASS_FRAME

!========================================================
! BORDER checking procedures
!========================================================


!---------------------------------------------------------
      LOGICAL FUNCTION TLIB_BORDER(shp,sz,r,k,tin,tout)
! test neutron cross-section with a frame of given shape and size
! return cross-section times, if any
! assume r,k in local coordinates
! Time is in units [m/h_bar] everywhere   i.e. length=time*K
!---------------------------------------------------------
      use CLASSES
      integer, intent(in) :: shp
      REAL(kind(1.D0)),intent(in)  :: sz(3), r(3) , k(3)
      REAL(kind(1.D0)),intent(out)  :: tin,tout
      SELECT CASE(shp)
! rectangular box
      CASE(FRAME_SHAPE_BOX)
        call BORDER_BOX(R,k,sz,tin,tout)
! cylinder
      CASE(FRAME_SHAPE_CYLINDER)
        call BORDER_CYL(R,k,sz,tin,tout)
! disc
      CASE(FRAME_SHAPE_DISC)
        call BORDER_DISC(R,k,sz,tin,tout)
! ellipsoid
      CASE(FRAME_SHAPE_ELLIPSOID)
        call BORDER_ELLIPS(R,k,sz,tin,tout)
      END SELECT
      TLIB_BORDER=((tin.lt.tout).and.(tin.gt.-TLIB_INF).and.(tout.lt.TLIB_INF))
      END FUNCTION TLIB_BORDER


!-----------------------------------------------------------------------
      SUBROUTINE BORDER_PLANE(x,kx,dx,a,t)
! get time of flight to the intersection with a plane with given normal angle w.r.t. x-axis
! INPUT:
!   x      ... position of neutron
!   kx     ... direction of neutron
!   dx     ... position of the plane
!   a      ... tangents of the plane
! OUTPUT
!   t      ... time to cross-sections
!------------------------------------------------------------------------
      REAL(kind(1.D0)),intent(in) :: x,kx,dx,a
      REAL(kind(1.D0)),intent(out) :: t
      REAL(kind(1.D0)) :: kxx
      kxx=kx-a
      if(abs(kxx).gt.1.D-20) then
        t=-(x-dx)/kxx
  ! parallel direction, t=+INF
      else
        t=sign(1.D0,dx)*TLIB_INF
      endif
      end SUBROUTINE BORDER_PLANE

!-----------------------------------------------------------------------
      SUBROUTINE BORDER_WALL(x,kx,dist,T1,T2)
! get time of flight between two planparallel walls
! INPUT:
!   x      ... position
!   kx     ... direction
!   dist   ... walls distance
! OUTPUT
!   T1,T2  ... time of cross-sections with the two walls in units [x/kx]
!------------------------------------------------------------------------
      REAL(kind(1.D0)),intent(in) :: x,kx,dist
      REAL(kind(1.D0)),intent(out) :: t1,t2
      REAL(kind(1.D0)) :: z
      if(abs(kx).gt.1.D-20) then
        t1=(-dist/2-x)/kx
        t2=(dist/2-x)/kx
        if(t2.lt.t1) then  ! only if dist<0, unexpected case
            z=t1
            t1=t2
            t2=z
        endif
  ! parallel direction, t=+-INF
      else
        t1=-1.D31
        t2=1.D31
      endif
      end SUBROUTINE BORDER_WALL


!-----------------------------------------------------------------------
      SUBROUTINE BORDER_WALL_QUAD(r,k,dist,curv,tin,tout,db)
! get time of flight between two planparallel quadratic surfaces
! radius vectors // -z
! surfaces are separated by +-dist/2 along z
! INPUT:
!   r         ... neutron position
!   k         ... wave vector
!   dist      ... surface separation
!   curv ... horizontal and vertical curvatures
! OUTPUT:
!   tin,tout  ... entry and exit times in units [sec*h_bar/m]
! tout<tin if there is no cross-section
!------------------------------------------------------------------------
      REAL(kind(1.D0)),intent(in) :: r(3),k(3),dist,curv(2)
      REAL(kind(1.D0)),intent(out) :: tin,tout
      REAL(kind(1.D0)) :: t1,t2,t3,t4
      REAL(kind(1.D0)) :: R1(3),dc
      logical :: db
1     format(a,' : ',6(1x,G12.5))
      if (db) write(*,1) 'BORDER_WALL_QUAD r,k',r,k
      if (db) write(*,1) '   curv',curv
  ! top surface
      R1=r+(/0.D0,0.D0,-0.5D0*dist/)
      call SURFACE_QUAD(R1,k,curv,t1,t2)
  ! botom surface
      R1=r+(/0.D0,0.D0,+0.5D0*dist/)
      call SURFACE_QUAD(R1,k,curv,t3,t4)
  ! distance from the wafer zero surface
      dc=r(3)-0.5*(curv(1)*r(1)**2+curv(2)*r(2)**2)
  ! above the wafer
      if (db) write(*,1) '   t1..4, dc',t1,t2,t3,t4,dc
      if (dc.gt.(0.5D0*dist)) then
        ! entry is the minimum positive time through the top surface
        ! exit is the minimum positive time through the bottom surface
        tin=MIN_POSITIVE((/t1,t2/),2)
        tout=MIN_POSITIVE((/t3,t4/),2)
        if (db) write(*,1) '  above',tin,tout
  ! below the wafer
      else if (dc.lt.(-0.5D0*dist)) then
      ! entry is the minimum positive time through the bottom surface
      ! exit is the minimum positive time through the top surface
        tin=MIN_POSITIVE((/t3,t4/),2)
        tout=MIN_POSITIVE((/t1,t2/),2)
        if (db) write(*,1) '  below',tin,tout
  ! inside the wafer
      else
      ! entry is the maximum negative time through any surface
      ! exit is the minimum positive time through any surface
        tin=MAX_NEGATIVE((/t1,t2,t3,t4/),4)
        tout=MIN_POSITIVE((/t1,t2,t3,t4/),4)
        if (db) write(*,1) '  inside',tin,tout
      endif
      end SUBROUTINE BORDER_WALL_QUAD

!-----------------------------------------------------------------------
      SUBROUTINE BORDER_BOX(r,k,csize,tin,tout)
! return neutron entry and exit times for a cuboid
! INPUT:
!   r         ... neutron position
!   k         ... wave vector
!   csize     ... box size
! OUTPUT:
!   tin,tout  ... entry and exit times in units [h_bar/m]
! tout<tin if there is no cross-section
!------------------------------------------------------------------------
      REAL(kind(1.D0)),intent(in) :: r(3),k(3),csize(3)
      REAL(kind(1.D0)),intent(out) :: tin,tout
      REAL(kind(1.D0)) :: t1,t2,t3,t4,t5,t6
      call BORDER_WALL(r(1),k(1),csize(1),t1,t2)
      call BORDER_WALL(r(2),k(2),csize(2),t3,t4)
      call BORDER_WALL(r(3),k(3),csize(3),t5,t6)
      tin=max(t1,t3,t5)
      tout=min(t2,t4,t6)
      if ((tin.ge.tout).or.(abs(tout-tin).gt.TLIB_INF)) then
         tin=1D31
         tout=-1D31
      endif
      end SUBROUTINE BORDER_BOX


!-----------------------------------------------------------------------
      SUBROUTINE BORDER_WAFER(r,k,curv,csize,tin,tout)
! return neutron entry and exit times for a bent wafer
! INPUT:
!   r         ... neutron position
!   k         ... wave vector
!   csize     ... box size
!   rhoX, rhoY ... horizontal and vertical curvatures
! OUTPUT:
!   tin,tout  ... entry and exit times in units [h_bar/m]
! tout<tin if there is no cross-section
!------------------------------------------------------------------------
      REAL(kind(1.D0)),intent(in) :: r(3),k(3),csize(3),curv(2)
      REAL(kind(1.D0)),intent(out) :: tin,tout
      REAL(kind(1.D0)) :: t1,t2,t3,t4,t5,t6
      call BORDER_WALL(r(1),k(1),csize(1),t1,t2)
      call BORDER_WALL(r(2),k(2),csize(2),t3,t4)
      call BORDER_WALL_QUAD(r,k,csize(3),curv,t5,t6,.false.)
      tin=max(t1,t3,t5)
      tout=min(t2,t4,t6)
      if ((tin.ge.tout).or.(abs(tout-tin).gt.TLIB_INF)) then
         tin=1D31
         tout=-1D31
      endif
      end SUBROUTINE BORDER_WAFER

!-----------------------------------------------------------------------
      SUBROUTINE BORDER_WAFER_DBG(r,k,curv,csize,tin,tout)
! return neutron entry and exit times for a bent wafer
! INPUT:
!   r         ... neutron position
!   k         ... wave vector
!   csize     ... box size
!   rhoX, rhoY ... horizontal and vertical curvatures
! OUTPUT:
!   tin,tout  ... entry and exit times in units [h_bar/m]
! tout<tin if there is no cross-section
!------------------------------------------------------------------------
      REAL(kind(1.D0)),intent(in) :: r(3),k(3),csize(3),curv(2)
      REAL(kind(1.D0)),intent(out) :: tin,tout
      REAL(kind(1.D0)) :: t1,t2,t3,t4,t5,t6
1     format(a,' : ',6(2x,G12.5))
      call BORDER_WALL(r(1),k(1),csize(1),t1,t2)
      call BORDER_WALL(r(2),k(2),csize(2),t3,t4)
      call BORDER_WALL_QUAD(r,k,csize(3),curv,t5,t6,.true.)
      write(*,1) 'BORDER_WAFER',t1,t2,t3,t4,t5,t6
      tin=max(t1,t3,t5)
      tout=min(t2,t4,t6)
      write(*,1) 'BORDER_WAFER',tin,tout
      if ((tin.ge.tout).or.(abs(tout-tin).gt.TLIB_INF)) then
         tin=1D31
         tout=-1D31
      endif
      end SUBROUTINE BORDER_WAFER_DBG

!-----------------------------------------------------------------------
      SUBROUTINE BORDER_SANDWICH(r,k,astack,curv,csize,tin,tout,dbg)
! return neutron entry and exit times for a bent sandwich of wafers
! like BENT_WAFER, but includes also stacking angles
! INPUT:
!   r         ... neutron position
!   k         ... wave vector
!   csize     ... box size
!   astack    ... tg of stackig angles for horizontal and vertical directions
!   curv      ... horizontal and vertical curvatures
! OUTPUT:
!   tin,tout  ... entry and exit times in units [h_bar/m]
! tout<tin if there is no cross-section
!------------------------------------------------------------------------
      REAL(kind(1.D0)),intent(in) :: r(3),k(3),astack(2),curv(2),csize(3)
      REAL(kind(1.D0)),intent(out) :: tin,tout
      REAL(kind(1.D0)) :: t1,t2,t3,t4,t5,t6,z
      logical :: dbg
1     format(a,': ',7(1x,G12.5))

      ! call BORDER_WALL(r(1),k(1),csize(1),t1,t2)
      ! call BORDER_WALL(r(2),k(2),csize(2),t3,t4)

      if (dbg) write(*,1) 'BORDER_SANDWICH ',r,k,astack(1)

      call BORDER_PLANE(r(1),k(1),-0.5*csize(1)+astack(1)*r(3),astack(1),t1)
      call BORDER_PLANE(r(1),k(1),0.5*csize(1)+astack(1)*r(3),astack(1),t2)
      if(t2.lt.t1) then
            z=t1
            t1=t2
            t2=z
      endif
      call BORDER_PLANE(r(2),k(2),-0.5*csize(2)+astack(2)*r(3),astack(2),t3)
      call BORDER_PLANE(r(2),k(2),0.5*csize(2)+astack(2)*r(3),astack(2),t4)
      if(t4.lt.t3) then
            z=t3
            t3=t4
            t4=z
      endif
      if (dbg) write(*,1) '   t1 .. t4',t1,t2,t3,t4

      call BORDER_WALL_QUAD(r,k,csize(3),curv,t5,t6,.false.)
      if (dbg) write(*,1) '   t5,t6',t5,t6
      tin=max(t1,t3,t5)
      tout=min(t2,t4,t6)
      if ((tin.ge.tout).or.(abs(tout-tin).gt.TLIB_INF)) then
         tin=1D31
         tout=-1D31
      endif
      end SUBROUTINE BORDER_SANDWICH

!-----------------------------------------------------------------------
      SUBROUTINE BORDER_BLADE(r,k,w,h1,h2,length,tin,tout)
! return neutron entry and exit times for a blade
! INPUT:
!   r         ... neutron position
!   k         ... wave vector
!   w,h1,h2,length     ... thickness, entry and exit height, length
! OUTPUT:
!   tin,tout  ... entry and exit times in units [h_bar/m]
! T=1e30 if there is no cross-section
!------------------------------------------------------------------------
      REAL(kind(1.D0)),intent(in) :: r(3),k(3),w,h1,h2,length
      REAL(kind(1.D0)),intent(out) :: tin,tout
      REAL(kind(1.D0)) :: t1,t2,t3,t4,t5,t6,z,a,dx
      ! side surfaces
      call BORDER_WALL(r(1),k(1),w,t1,t2)
      ! top edge
      a=0.5D0*(h2-h1)/length
      dx=0.5D0*h1
      call BORDER_PLANE(r(2),k(2),dx,a,t3)
      ! bottom edge
      a=-a
      dx=-dx
      call BORDER_PLANE(r(2),k(2),dx,a,t4)
      ! front and rear edge
      if(t4.lt.t3) then  ! ensure t3<=t4
        z=t3
        t3=t4
        t4=z
      endif
      call BORDER_WALL(r(3)-0.5D0*length,k(3),length,t5,t6)
      tin=max(t1,t3,t5)
      tout=min(t2,t4,t6)
      if ((tin.ge.tout).or.(abs(tout-tin).gt.TLIB_INF)) then
         tin=1D31
         tout=-1D31
      endif
      end SUBROUTINE BORDER_BLADE

!-----------------------------------------------------------------------
      SUBROUTINE BORDER_CYL(r,k,csize,tin,tout)
! return neutron entry and exit times for a cylindrical container
! INPUT:
!   r         ... neutron position
!   k         ... wave vector
!   csize     ... cylinder size, height=csize(2)
! OUTPUT:
!   tin,tout  ... entry and exit times in units [sec*h_bar/m]
! T=-+1e31 if there is no cross-section
!------------------------------------------------------------------------
      REAL(kind(1.D0)),intent(in)  :: r(3),k(3),csize(3)
      REAL(kind(1.D0)),intent(out)  :: tin,tout
      REAL(kind(1.D0)) :: a,b,c,dd,t1,t2,t3,t4

      a=4*((k(1)/csize(1))**2+(k(3)/csize(3))**2)
      b=8*(k(1)*r(1)/csize(1)**2+k(3)*r(3)/csize(3)**2)
      c=4*((r(1)/csize(1))**2+(r(3)/csize(3))**2)-1.D0
      dd=b**2-4*a*c
      if ((dd.gt.0.D0).and.(a.gt.1.D-30)) then
        t1=(-b-sqrt(dd))/2/a
        t2=(-b+sqrt(dd))/2/a
        call BORDER_WALL(r(2),k(2),csize(2),t3,t4)
        tin=max(t1,t3)
        tout=min(t2,t4)
      else
        tin=1.D31
        tout=-1.D31
      endif
      end SUBROUTINE BORDER_CYL

!-----------------------------------------------------------------------
      SUBROUTINE BORDER_DISC(r,k,csize,tin,tout)
! return neutron entry and exit times for a DISC container
! similar as cylinder, but coordinates 2 and 3 are exchanged
! INPUT:
!   r         ... neutron position
!   k         ... wave vector
!   csize     ... cylinder size, height=csize(2)
! OUTPUT:
!   tin,tout  ... entry and exit times in units [sec*h_bar/m]
! T=-+1e31 if there is no cross-section
!------------------------------------------------------------------------
      REAL(kind(1.D0)),intent(in)  :: r(3),k(3),csize(3)
      REAL(kind(1.D0)),intent(out)  :: tin,tout
      REAL(kind(1.D0)) :: a,b,c,dd,t1,t2,t3,t4
      a=4*((k(1)/csize(1))**2+(k(2)/csize(2))**2)
      b=8*(k(1)*r(1)/csize(1)**2+k(2)*r(2)/csize(2)**2)
      c=4*((r(1)/csize(1))**2+(r(2)/csize(2))**2)-1.D0
      dd=b**2-4*a*c
      if ((dd.gt.0.D0).and.(a.gt.1.D-30)) then
        t1=(-b-sqrt(dd))/2/a
        t2=(-b+sqrt(dd))/2/a
        call BORDER_WALL(r(3),k(3),csize(3),t3,t4)
        tin=max(t1,t3)
        tout=min(t2,t4)
      else
        tin=1.D31
        tout=-1.D31
      endif
      end SUBROUTINE BORDER_DISC

!-----------------------------------------------------------------------
      SUBROUTINE BORDER_ELLIPS(r,k,csize,tin,tout)
! return neutron entry and exit times for an ellipsoid
! INPUT:
!   r         ... neutron position
!   k         ... wave vector
!   csize     ... ellipse size, height=csize(2)
! OUTPUT:
!   tin,tout  ... entry and exit times in units [sec*h_bar/m]
! T=-+1e31 if there is no cross-section
!------------------------------------------------------------------------
      REAL(kind(1.D0)),intent(in)  :: r(3),k(3),csize(3)
      REAL(kind(1.D0)),intent(out)  :: tin,tout
      REAL(kind(1.D0)) :: a,b,c,dd

      a=4*((k(1)/csize(1))**2+(k(2)/csize(2))**2+(k(3)/csize(3))**2)
      b=8*(k(1)*r(1)/csize(1)**2+k(2)*r(2)/csize(2)**2+k(3)*r(3)/csize(2)**2)
      c=4*((r(1)/csize(1))**2+(r(2)/csize(2))**2+(r(3)/csize(3))**2)-1.D0
      dd=b**2-4*a*c
      if ((dd.gt.0.D0).and.(a.gt.1.D-30)) then
        tin=(-b-sqrt(dd))/2/a
        tout=(-b+sqrt(dd))/2/a
      else
        tin=1.D31
        tout=-1.D31
      endif
      end SUBROUTINE BORDER_ELLIPS

!-----------------------------------------------------------------------
      SUBROUTINE SURFACE_SPH(r,k,rad,tin,tout)
! cross-section times through as sphere surface
! simplyfied version of BORDER_ELLIPS
! r is relative to the spehere center
! INPUT:
!   r         ... neutron position
!   k         ... wave vector
!   rad       ... sphere radius
! OUTPUT:
!   tin,tout  ... entry and exit times in units [sec*h_bar/m]
! T=-+1e31 if there is no cross-section
!------------------------------------------------------------------------
      REAL(kind(1.D0)),intent(in)  :: r(3),k(3),rad
      REAL(kind(1.D0)),intent(out)  :: tin,tout
      REAL(kind(1.D0)) :: a,b,c,dd
      a=k(1)**2+k(2)**2+k(3)**2
      b=2.D0*(k(1)*r(1)+k(2)*r(2)+k(3)*r(3))
      c=r(1)**2+r(2)**2+r(3)**2-rad**2
      dd=b**2-4*a*c
      if (dd.gt.0.D0) then
        tin=(-b-sqrt(dd))/(2.D0*a)
        tout=(-b+sqrt(dd))/(2.D0*a)
      else
        tin=1.D31
        tout=-1.D31
      endif
      end SUBROUTINE SURFACE_SPH

!-----------------------------------------------------------------------
      SUBROUTINE SURFACE_CYL(r,k,rad,tin,tout)
! like BORDER_CYL with circular base, but does not limit the top and bottom
!------------------------------------------------------------------------
      REAL(kind(1.D0)),intent(in)  :: r(3),k(3),rad
      REAL(kind(1.D0)),intent(out)  :: tin,tout
      REAL(kind(1.D0)) :: a,b,c,dd
      a=k(1)**2+k(3)**2
      b=2.D0*(k(1)*r(1)+k(3)*r(3))
      c=r(1)**2+r(3)**2-rad**2
      dd=b**2-4*a*c
      if (dd.gt.0.D0) then
        tin=(-b-sqrt(dd))/(2.D0*a)
        tout=(-b+sqrt(dd))/(2.D0*a)
      else
        tin=1.D31
        tout=-1.D31
      endif
      end SUBROUTINE SURFACE_CYL


!-----------------------------------------------------------------------
      SUBROUTINE SURFACE_QUAD(r,k,curv,tin,tout)
! cross times with a quadratic surface
! rhoH, rhoV are horizotal (//x} and vertical (//y} curvatures
!------------------------------------------------------------------------
      REAL(kind(1.D0)),intent(in)  :: r(3),k(3),curv(2)
      REAL(kind(1.D0)),intent(out)  :: tin,tout
      REAL(kind(1.D0)) :: a,b,c,dd,t1
      a=0.5*(curv(1)*k(1)**2+curv(2)*k(2)**2)
      b=curv(1)*k(1)*r(1)+curv(2)*k(2)*r(2)-k(3)
      c=0.5*(curv(1)*r(1)**2+curv(2)*r(2)**2)-r(3)
      dd=b**2-4*a*c
      if((a.eq.0.D0).and.(abs(b)>0.D0)) then
        tin=-c/b
        tout=tin
      else if (dd.gt.0.D0) then
        tin=(-b-sqrt(dd))/(2.D0*a)
        tout=(-b+sqrt(dd))/(2.D0*a)
        if (tout<tin) then
          t1=tout
          tout=tin
          tin=t1
        endif
      else
        tin=1.D31
        tout=-1.D31
      endif
      end SUBROUTINE SURFACE_QUAD

!-----------------------------------------------------------------------
      SUBROUTINE SURFACE_PLANE(r,k,norm,t)
! get cross time with a plane defined by a normal vector
! (the plane crosses the coord. origin)
! t=1.D31 if k is parallel to the plane
!------------------------------------------------------------------------
      REAL(kind(1.D0)),intent(in)  :: r(3),k(3),norm(3)
      REAL(kind(1.D0)),intent(out)  :: t
      REAL(kind(1.D0)) :: kn,rn
      kn=norm(1)*k(1)+norm(2)*k(2)+norm(3)*k(3)
      rn=norm(1)*r(1)+norm(2)*r(2)+norm(3)*r(3)
      if (abs(kn)<1.D-30) then
        t=-1D31*sign(1.D0,kn*rn)
      else
        t=-rn/kn
      endif
      end SUBROUTINE SURFACE_PLANE

!-----------------------------------------------------------------------
      SUBROUTINE SURFACE_XZ(r,k,A,B,is,t)
! get cross time with a planar surface normal to the (is,z) plane
! the surface is defined by z=A+B*is
! is is either x=1 or y=2
! r,k are neutron position and velocity
! t=-TLIB_INF if k is parallel to the plane
!------------------------------------------------------------------------
      REAL(kind(1.D0)),intent(in)  :: r(3),k(3),A,B
      integer,intent(in) :: is
      REAL(kind(1.D0)),intent(out)  :: t
      REAL(kind(1.D0)) :: dnm
      dnm=K(3)-B*K(is)
      if (abs(dnm)>1.D-10) then
        t=(A - R(3) + B*R(is))/dnm
      else
        t=-TLIB_INF
      endif
      end SUBROUTINE SURFACE_XZ

!-----------------------------------------------------------------------
      SUBROUTINE SURFACE_SECTOR(r,k,cosa,sina,tin,tout)
! get entry and exit times for a sector defined by 2 planes
! with cross-section = y-xis
! r=(x,z), k=(kx,kz)
!------------------------------------------------------------------------
      REAL(kind(1.D0)),intent(in)  :: r(2),k(2),cosa,sina
      REAL(kind(1.D0)),intent(out)  :: tin,tout
      REAL(kind(1.D0)) :: kn1,kn2,rn1,rn2,t1,t2
      REAL(kind(1.D0)),parameter :: EPS=1.D-30
      tin=TLIB_INF
      tout=-TLIB_INF
      kn1=cosa*k(1)-sina*k(2)
      kn2=-cosa*k(1)-sina*k(2)
      rn1=cosa*r(1)-sina*r(2)
      rn2=-cosa*r(1)-sina*r(2)
      ! avoid tangential directions
      if (kn1==0.D0) then
        kn1=EPS
      else if (abs(kn1)<EPS) then
        kn1=sign(EPS,kn1)
      endif
      if (kn2==0.D0) then
        kn2=EPS
      else if (abs(kn2)<EPS) then
        kn2=sign(EPS,kn2)
      endif
    ! there is only one entry and one exit
      if (kn1*kn2<0.D0) then
        if (kn2<0.D0) then
          t1=-rn2/kn2
          t2=-rn1/kn1
        else
          t1=-rn1/kn1
          t2=-rn2/kn2
        endif
      endif
    ! we have to distinguish two cases:
    ! sector angle < PI; once you escape, you can't re-enter
      if (cosa>0.D0) then
        ! there is only one entry and one exit
        if (kn1*kn2<0.D0) then
          if (t1<t2) then
            tin=t1
            tout=t2
          else
          ! no crossing, always out of the sector
            tin=TLIB_INF
            tout=-TLIB_INF
          endif
        ! 2 entries => take the latest one
        else if ((kn1<0.D0).and.(kn2<0.D0)) then
          tin=max(-rn1/kn1,-rn2/kn2)
          tout=TLIB_INF
        ! 2 exits => take the 1st one
        else
          tin=-TLIB_INF
          tout=min(-rn1/kn1,-rn2/kn2)
        endif
    ! sector angle > PI; once you enter, you can't escape
      else
        ! there is only one entry and one exit
        if (kn1*kn2<0.D0) then
          if (t1>t2) then
            tin=t1
            tout=t2
          else
          ! no crossing, always inside the sector
            tin=-TLIB_INF
            tout=TLIB_INF
          endif
        ! 2 entries => take the 1st one
        else if ((kn1<0.D0).and.(kn2<0.D0)) then
          tin=min(-rn1/kn1,-rn2/kn2)
          tout=TLIB_INF
        ! 2 exits => take the latest one
        else
          tin=-TLIB_INF
          tout=max(-rn1/kn1,-rn2/kn2)
        endif
      endif
      end SUBROUTINE SURFACE_SECTOR

      end module TRACELIB
