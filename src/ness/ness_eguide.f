C $Id: ness_eguide.f,v 1.45 2019/08/16 17:16:26 saroun Exp $
C//////////////////////////////////////////////////////////////////////
C////                                                              ////
C////  NEutron Scattering Simulation - v.2.0, (c) J.Saroun, 2009   ////
C////                                                              ////
C//////////////////////////////////////////////////////////////////////
C////
C////  Subroutines describing objects - ELLIPTIC GUIDE
C////
C////
C//////////////////////////////////////////////////////////////////////

C---------------------------------------------------------------------
      SUBROUTINE EGUIDE_INIT(OBJ)
C----------------------------------------------------------------------
      USE GUIDES
      use FRAMES_TRACE
      IMPLICIT NONE

      REAL*4 EGUIDE_A
      REAL*8 f,ZL,DIF,H1,H2,V1,V2,W
      TYPE(GUIDE) :: OBJ
      INTEGER*4 I

      CALL FRAME_INIT(OBJ%FRAME)

C limit for number of slits = 255
      IF (OBJ%NLH.GT.255) OBJ%NLH=255
      IF (OBJ%NLV.GT.255) OBJ%NLV=255

C HORIZONTAL
c elliptic profile is determined by the dimensions
      H1=OBJ%FRAME%SIZE(1)
      H2=OBJ%W2
      IF (H2.EQ.H1.OR.OBJ%ROH.EQ.0) THEN  !  no ellipsa, flat walls
         OBJ%CH=0.D0
      ELSE IF (H2.GT.H1) THEN
       DIF=SQRT(H2**2-H1**2)
       OBJ%CH=-OBJ%FRAME%SIZE(3)*H2/DIF
      ELSE
       DIF=SQRT(H1**2-H2**2)
       OBJ%CH=OBJ%FRAME%SIZE(3)*H1/DIF
      ENDIF

      f=OBJ%CH
      IF (f.EQ.0) THEN ! flat lamellas, AH & LH are angles & positions, respectively
        DO I=0,OBJ%NLH
           ZL=I*1.D0/OBJ%NLH - 0.5D0
         OBJ%AH(I)=ZL*(H2-H1)/OBJ%FRAME%SIZE(3)
         OBJ%LH(I)=ZL*H1
        ENDDO
      ELSE  ! elliptic lamellas, AH & LH are parameters & lengths, respectively
      w=max(H1,H2)
        DO I=0,OBJ%NLH
           OBJ%AH(I)=EGUIDE_A(I,OBJ%NLH,w)
           OBJ%LH(I)=OBJ%FRAME%SIZE(3)
        ENDDO
      ENDIF

C VERTICAL
c elliptic profile is determined by the dimensions
      V1=OBJ%FRAME%SIZE(2)
      V2=OBJ%H2
      IF (V2.EQ.V1.OR.OBJ%ROV.EQ.0) THEN  !  no parabola, flat walls
       OBJ%CV=0.D0
      ELSE IF (V2.GT.V1) THEN
       DIF=SQRT(V2**2-V1**2)
       OBJ%CV=-OBJ%FRAME%SIZE(3)*V2/DIF
      ELSE
       DIF=SQRT(V1**2-V2**2)
       OBJ%CV=OBJ%FRAME%SIZE(3)*V1/DIF
      ENDIF
      f=OBJ%CV
      IF (f.EQ.0) THEN ! flat lamellas, AH & LH are angles & positions, respectively
        DO I=0,OBJ%NLV
           ZL=I*1.E0/OBJ%NLV - 0.5D0
         OBJ%AV(I)=ZL*(V2-V1)/OBJ%FRAME%SIZE(3)
         OBJ%LV(I)=ZL*V1
        ENDDO
      ELSE  ! elliptic lamellas, AH & LH are parameters & lengths, respectively
      w=max(V1,V2)
        DO I=0,OBJ%NLV
           OBJ%AV(I)=EGUIDE_A(I,OBJ%NLV,w)
           OBJ%LV(I)=OBJ%FRAME%SIZE(3)
        ENDDO
      ENDIF

      END SUBROUTINE EGUIDE_INIT

C---------------------------------------------------------------------
      REAL*8 FUNCTION EGUIDE_LAM(OBJ,ID,IL,Z,IOR)
C function describing lamella profile
C ID   ... derivative
C IL   ... lamella index
C Z    ... distance from guide entry
C IOR  ... horizontal (1) or vertical (2)
C----------------------------------------------------------------------
      USE GUIDES
      IMPLICIT NONE
      TYPE(GUIDE) :: OBJ

      REAL*8 a,b,zz,aa,x0
      REAL*8 Z
      INTEGER*4 ID,IL,IOR
C//  a is the ellipse main axis // z (or slope for a planar lamellla)
C//  b is the ellipse smaller axis

      IF (IOR.EQ.2) THEN   ! vertical slit
        a=OBJ%CV
        b=OBJ%AV(IL)
        x0=OBJ%LV(IL)
      ELSE                 ! horizontal slit
        a=OBJ%CH
        b=OBJ%AH(IL)
        x0=OBJ%LH(IL)
      ENDIF

      IF (a.LT.0) THEN     ! focal point before the guide entry
        zz=z-OBJ%FRAME%SIZE(3)
      ELSE    ! focal point behind the guide
        zz=z
      ENDIF
      aa=a**2-zz**2

      IF(a.NE.0) THEN
        IF(aa.LE.0) THEN
          EGUIDE_LAM=0.D0
          RETURN
        ELSE
          aa=SQRT(aa)
        ENDIF
      ENDIF

c zero deriv.
      IF (ID.LE.0) THEN
        IF (a.EQ.0) THEN
          EGUIDE_LAM=x0+b*z
        ELSE
          EGUIDE_LAM=b/abs(a)*aa
        ENDIF
c 1st deriv.
      ELSE IF (ID.EQ.1) THEN
        IF (a.EQ.0) THEN
        EGUIDE_LAM=b
      ELSE
          EGUIDE_LAM=-b/abs(a)/aa*zz
      ENDIF
c 2nd deriv.
      ELSE IF (ID.EQ.2) THEN
        IF (a.EQ.0) THEN
        EGUIDE_LAM=0.D0
      ELSE
          EGUIDE_LAM=-b*abs(a)/aa/aa**2
      ENDIF
      ELSE
        EGUIDE_LAM=0.D0
      ENDIF

      END FUNCTION EGUIDE_LAM


C-------------------------------------------------------------------------------
      REAL*4 FUNCTION EGUIDE_A(il,nl,w)
C ellipsa parameter for il-th lamella = smaller axis
C il ... lamella index, counted from right, right side il=0
C w  ... small axis of the outer profile
C nl ... number of slits (number of lamellae + 1)
C sign(A) determines which side from the guide center: right/bottom(<0) or left/top(>0)
C--------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 w
      INTEGER*4 il,nl

      EGUIDE_A=w*(il*1.E0/nl-0.5E0)

      END FUNCTION EGUIDE_A


C----------------------------------------------------------------------------
      REAL*8 FUNCTION EGUIDE_CROSS(kx,kz,x,z,b,a,tmax,glen)
C return TOF to a cross-point with elliptic surface
C if there are 2 solutions, take the minimum one which is > 0
C return 10^30 if there is no solution 0 < result < tmax
C b      .. smaller ellipsa axis
C a      .. main ellipsa axis
C lmax   .. limit time (TOF to exit)
C x,z    .. starting point (z // guide axis)
C kx,kz  .. ray direction
C glen   .. guide length
C----------------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 EPS
      PARAMETER (EPS=1.D-7)
      REAL*8 kx,kz,x,z,b,a,tmax,zz,glen
      REAL*8 T,T1,T2,dtm,dnom,t0

      IF (a.LT.0.0) THEN     ! focal point before the guide entry
        zz=z-glen
      ELSE    ! focal point behind the guide
        zz=z
      ENDIF
      T=1.D30
      dnom=(kx*a)**2 + (kz*b)**2
      dtm=(a*b)**2*(dnom-(kx*zz-kz*x)**2)
      t0=x*kx*a**2 + zz*kz*b**2
      IF (dtm.GE.0.D0) THEN 
        T1=(-t0+SQRT(dtm))/dnom
        T2=(-t0-SQRT(dtm))/dnom
c        write(*,*) '    EGUIDE_CROS zz,a,b',zz,a,b
c        write(*,*) '    EGUIDE_CROS T1,T2',T1,T2,SQRT(dtm)
! choose the minimum positive value
        if ((T1.GT.EPS).and.(T2.GT.EPS)) then
          T=MIN(T1,T2)
        else if (T1.GT.EPS) then
          T=T1
        else if (T2.GT.EPS) then
          T=T2
        ENDIF
      ENDIF
      IF (T.GT.tmax) T=1.D30 ! no reflection behind the guide
      EGUIDE_CROSS=T

      END FUNCTION EGUIDE_CROSS


C---------------------------------------------------------------------
      SUBROUTINE EGUIDE_CON(OBJ,R,K,IH,IV,IC,DT,Q,NW,KG,P)
C find next contact with a slit side
C moves to the next contact point, or to the exit
C turns K, i.e. returns new K=K+Q !
C *** INPUT:
C IH,IV   ... slit indices
C R(3)    ... neutron coordinate
C K(3)    ... neutron k-vector
C *** RETURN:
C Q     .. reflection vector (magnitude)
C DT    .. time to reach the next contact point
C IC    .. which wall,  left(1), right(2), top(3), bottom (4)
C NW    .. normal to the reflecting surface (with waviness)
C KG    .. K+Q*NW
C IC=0  .. no contact, pass through
C----------------------------------------------------------------------
      USE GUIDES
      USE TRACINGDATA
      USE CONSTANTS
      use MIRROR_TABLE
      use SIMATH
      IMPLICIT NONE
      TYPE(GUIDE) :: OBJ
      REAL(kind(1.D0)),intent(in) :: R(3),K(3)
      REAL(kind(1.D0)),intent(out) :: Q,DT,NW(3),KG(3),P
      INTEGER :: IH,IV,IC
      REAL(kind(1.D0)) :: EGUIDE_CROSS,GUIDE_CROSS1,EGUIDE_LAM
      REAL(kind(1.D0)) :: kx,kz,x,z,a,f,tmax,glen,ZL
      REAL(kind(1.D0)) :: LZ(4),ANG,T0,GX,N(3)
      INTEGER :: I,IX
      REAL(kind(1.D0)) :: OSCD,OSCA
      integer,parameter :: ID(4)=(/1,1,2,2/)
      integer,parameter :: IDINV(4)=(/2,2,1,1/)
      COMMON /OSCBEND/ OSCD,OSCA
      LOGICAL*4 DBG
      common /dbg_eg/ dbg

1     FORMAT(a,1x,9(1x,G12.5))
c10    FORMAT(a11,1x,6(1x,G11.5),2(1x,I4),2(1x,G10.4))
c      write(*,10) 'CON START: ',R,K,IH,IV,OBJ%CH
      KG=K
      z=R(3)
      kz=K(3)
      gx=NEUT_GR(1)

! TOF to the guide end
      if (TRACING_UP) then
        T0=-z/kz
      else
        T0=(OBJ%FRAME%SIZE(3)-z)/kz
      endif
      glen=OBJ%FRAME%SIZE(3)
C HORIZONTAL RIGHT:
      f=OBJ%CH
      kx=K(1)
      tmax=T0
      IF (f.EQ.0) THEN
        a=OBJ%AH(IH)-OSCA
        x=R(1)-OBJ%DLH/2.D0-OBJ%LH(IH)-OSCD-R(3)*OSCA
        !LZ(2)=GUIDE_CROSS(kx,kz,x,z,a,f)
        LZ(2)=GUIDE_CROSS1(kx,kz,x,z,gx,a,f)
      ELSE
        a=OBJ%AH(IH)
        x=R(1)-OBJ%DLH/2.D0
        LZ(2)=EGUIDE_CROSS(kx,kz,x,z,a,f,tmax,glen)
        if (dbg)  write(*,1) 'RIGHT a,x,LZ ',a,x,LZ(2)
      ENDIF
C HORIZONTAL LEFT:
      IF (f.EQ.0) THEN
         a=OBJ%AH(IH+1)-OSCA
         x=R(1)+OBJ%DLH/2.D0-OBJ%LH(IH+1)-OSCD-R(3)*OSCA
         !LZ(1)=GUIDE_CROSS(kx,kz,x,z,a,f)
         LZ(1)=GUIDE_CROSS1(kx,kz,x,z,gx,a,f)
      ELSE
         a=OBJ%AH(IH+1)
         x=R(1)+OBJ%DLH/2.D0
         LZ(1)=EGUIDE_CROSS(kx,kz,x,z,a,f,tmax,glen)
        if (dbg)  write(*,1) 'LEFT  a,x,LZ ',a,x,LZ(1)
      ENDIF
C VERTICAL BOTTOM:
      f=OBJ%CV
      kx=K(2)
      gx=NEUT_GR(2)
      tmax=T0
      IF (f.EQ.0) THEN
         a=OBJ%AV(IV)
         x=R(2)-OBJ%DLV/2.D0-OBJ%LV(IV)
         !LZ(4)=GUIDE_CROSS(kx,kz,x,z,a,f)
         LZ(4)=GUIDE_CROSS1(kx,kz,x,z,gx,a,f)
      ELSE
         a=OBJ%AV(IV)
         x=R(2)-OBJ%DLV/2.D0
         LZ(4)=EGUIDE_CROSS(kx,kz,x,z,a,f,tmax,glen)
      ENDIF
C VERTICAL TOP:
      IF (f.EQ.0) THEN
         a=OBJ%AV(IV+1)
         x=R(2)+OBJ%DLV/2.D0-OBJ%LV(IV+1)
         !LZ(3)=GUIDE_CROSS(kx,kz,x,z,a,f)
         LZ(3)=GUIDE_CROSS1(kx,kz,x,z,gx,a,f)
      ELSE
         a=OBJ%AV(IV+1)
         x=R(2)+OBJ%DLV/2.D0
         LZ(3)=EGUIDE_CROSS(kx,kz,x,z,a,f,tmax,glen)
      ENDIF
      if (dbg)    write(*,1) 'LZ ',LZ,T0
c      write(*,10) 'kz*times: ',LZ(1),LZ(2)

! select the shortest time to cross a surface
      DT=MIN(LZ(1),LZ(2),LZ(3),LZ(4),T0)
      IF (DT.LE.0D0) goto 99 ! signal to stop


! and find which of the 4 surfaces it was
      IC=0
      DO I=1,4
        IF (DT.EQ.LZ(I)) then
          if (IC.eq.0) then
            IC=I
      ! both sides belong to the same ellipse
      ! move to the cross point and look which side is it
          else
            x=R(ID(IC))+K(ID(IC))*DT
            if (dbg) write(*,1) 'single ellipse: x,IC,I ',x,IC,I
            if (x.gt.0) then
              IC=MIN(IC,I)
            else
              IC=MAX(IC,I)
            endif
            if (dbg) write(*,1) 'the result is IC=',IC
          endif
        endif
      ENDDO
! no contact, passed through
      IF (IC.EQ.0) THEN
        DT=T0
        Q=0.D0
        return
      ENDIF

C get the surface normal vector
      N=0.D0
      ZL=0.D0
      select case(IC)
      case(2)
        ANG=EGUIDE_LAM(OBJ,1,IH,Z+DT*kz,1)
        ZL=IH*1.D0/OBJ%NLH - 0.5D0
        N(3)=-ANG/SQRT((1.D0+ANG**2))
        N(1)=SQRT(1.D0-N(3)**2)
      case(1)
        ANG=EGUIDE_LAM(OBJ,1,IH+1,Z+DT*kz,1)
        ZL=(IH+1)*1.D0/OBJ%NLH - 0.5D0
        N(3)=ANG/SQRT((1.D0+ANG**2))
        N(1)=-SQRT(1.D0-N(3)**2)
      case(4)
        ANG=EGUIDE_LAM(OBJ,1,IV,Z+DT*kz,2)
        ZL=IV*1.D0/OBJ%NLV - 0.5D0
        N(3)=-ANG/SQRT((1.D0+ANG**2))
        N(2)=SQRT(1.D0-N(3)**2)
      case(3)
        ANG=EGUIDE_LAM(OBJ,1,IV+1,Z+DT*kz,2)
        ZL=(IV+1)*1.D0/OBJ%NLV - 0.5D0
        N(3)=ANG/SQRT((1.D0+ANG**2))
        N(2)=-SQRT(1.D0-N(3)**2)
      END select

      IX=ID(IC)
! stop events reflected from concave side
      IF (OBJ%ONESIDE.eq.1) then
        IF (N(IX)*ZL.GT.0) goto 99
      endif
      call MIRROR_REFLECT(OBJ%WAV,K,N,IC,KG,NW,Q,P)
      if (P<=0.D0) goto 99
      Return
! signal to stop the event
99    Q=-10.D0
      DT=0.D0
      IC=0
      P=0.D0
      END SUBROUTINE EGUIDE_CON


