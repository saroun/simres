C $Id: ness_pguide.f,v 1.45 2019/08/16 17:16:26 saroun Exp $
C//////////////////////////////////////////////////////////////////////
C////                                                              ////
C////  NEutron Scattering Simulation - v.2.0, (c) J.Saroun, 2004   ////
C////                                                              ////
C//////////////////////////////////////////////////////////////////////
C////
C////  Subroutines describing objects - PARABOLIC GUIDE
C////
C////
C//////////////////////////////////////////////////////////////////////

C---------------------------------------------------------------------
      SUBROUTINE PGUIDE_INIT(OBJ)
C----------------------------------------------------------------------
      USE GUIDES
      USE FRAMES_TRACE
      IMPLICIT NONE

      REAL*4 PGUIDE_A
      REAL*4 XX
      REAL*8 f,ZL,DIF,A,H1,H2,V1,V2
      TYPE(GUIDE) :: OBJ
      INTEGER*4 I
1     format(a,' ',I5,' ',3(G12.5,2x))
2     format(a,7(1x,G11.5))
      CALL FRAME_INIT(OBJ%FRAME)

C limit for number of slits = 255
      IF (OBJ%NLH.GT.255) OBJ%NLH=255
      IF (OBJ%NLV.GT.255) OBJ%NLV=255

C HORIZONTAL
c parabolic profile is determined by the dimensions
      H1=OBJ%FRAME%SIZE(1)
      H2=OBJ%W2
      IF (H2.EQ.H1.OR.OBJ%ROH.EQ.0) THEN  !  no parabola, flat walls
        OBJ%CH=0.D0
      ELSE IF (H2.GT.H1) THEN
        DIF=H2**2-H1**2
    ! ERROR: should be
        A=DIF/8.D0/OBJ%FRAME%SIZE(3)
    !   A=DIF/4.D0/OBJ%FRAME%SIZE(3)
        OBJ%CH=A-H1**2*OBJ%FRAME%SIZE(3)/DIF
      ELSE IF (H2.LT.H1) THEN
        DIF=H1**2-H2**2
        A=DIF/8.D0/OBJ%FRAME%SIZE(3)
        OBJ%CH=H2**2*OBJ%FRAME%SIZE(3)/DIF-A
      ENDIF

      f=OBJ%CH
  !    write(*,2) 'PGUIDE_INIT ID='//trim(OBJ%FRAME%ID)//' f=',f
      IF (f.EQ.0) THEN ! flat lamellas, AH & LH are angles & positions, respectively
        DO I=0,OBJ%NLH
          ZL=I*1.D0/OBJ%NLH - 0.5D0
          OBJ%AH(I)=ZL*(OBJ%W2-OBJ%FRAME%SIZE(1))/OBJ%FRAME%SIZE(3)
          OBJ%LH(I)=ZL*OBJ%FRAME%SIZE(1)
  !        write(*,1) 'PGUIDE_INIT AH,LH=',I,OBJ%AH(I),OBJ%LH(I)
        ENDDO
      ELSE  ! parabolic lamellas, AH & LH are parameters & lengths, respectively
        IF (f.GT.0) f=f+OBJ%FRAME%SIZE(3)
        DO I=0,OBJ%NLH
          OBJ%AH(I)=PGUIDE_A(I,OBJ%NLH,OBJ%FRAME%SIZE(1),f)
        !  write(*,1) 'PGUIDE_INIT AH=',I,OBJ%AH(I)
        ENDDO
        IF ((OBJ%TYP.EQ.3).AND.(OBJ%NLH.GT.1)) THEN  ! set optimum lengths of lamellae
          DO I=0,OBJ%NLH
            XX=ABS(I-0.5*OBJ%NLH)
            IF (XX.EQ.0) THEN
              OBJ%LH(I)=OBJ%FRAME%SIZE(3)
            ELSE
              OBJ%LH(I)=MIN(2.0*ABS(OBJ%CH)/XX,OBJ%FRAME%SIZE(3))
            ENDIF
c      write(*,*) 'lam: ',I,'  ',OBJ%LH(I),'  ',XX*OBJ%FRAME%SIZE(1)/OBJ%NLH
          ENDDO
        ELSE
          DO I=0,OBJ%NLH
            OBJ%LH(I)=OBJ%FRAME%SIZE(3)
          ENDDO
        ENDIF
      ENDIF

C VERTICAL
c parabolic profile is determined by the dimensions
      V1=OBJ%FRAME%SIZE(2)
      V2=OBJ%H2
      IF (V2.EQ.V1.OR.OBJ%ROV.EQ.0) THEN  !  no parabola, flat walls
        OBJ%CV=0.D0
      ELSE IF (V2.GT.V1) THEN
        DIF=V2**2-V1**2
        A=DIF/8.D0/OBJ%FRAME%SIZE(3)
        OBJ%CV=A-V1**2*OBJ%FRAME%SIZE(3)/DIF
      ELSE IF (V2.LT.V1) THEN
        DIF=V1**2-V2**2
        A=DIF/8.D0/OBJ%FRAME%SIZE(3)
        OBJ%CV=V2**2*OBJ%FRAME%SIZE(3)/DIF-A
      ENDIF
      f=OBJ%CV
      IF (f.EQ.0) THEN ! flat lamellas, AH & LH are angles & positions, respectively
        DO I=0,OBJ%NLV
           ZL=I*1.D0/OBJ%NLV - 0.5D0
         OBJ%AV(I)=ZL*(OBJ%H2-OBJ%FRAME%SIZE(2))/OBJ%FRAME%SIZE(3)
         OBJ%LV(I)=ZL*OBJ%FRAME%SIZE(2)
        ENDDO
      ELSE  ! parabolic lamellas, AH & LH are parameters & lengths, respectively
      IF (f.GT.0) f=f+OBJ%FRAME%SIZE(3)
        DO I=0,OBJ%NLV
           OBJ%AV(I)=PGUIDE_A(I,OBJ%NLV,OBJ%FRAME%SIZE(2),f)
        ENDDO
        IF ((OBJ%TYP.EQ.3).AND.(OBJ%NLV.GT.1)) THEN  ! set optimum lengths of lamellae
          DO I=0,OBJ%NLV
            XX=ABS(I-0.5*OBJ%NLV)
            IF (XX.EQ.0) THEN
              OBJ%LV(I)=OBJ%FRAME%SIZE(3)
            ELSE
              OBJ%LV(I)=MIN(2.0*ABS(OBJ%CV)/XX,OBJ%FRAME%SIZE(3))
            ENDIF
          ENDDO
        ELSE
          DO I=0,OBJ%NLV
           OBJ%LV(I)=OBJ%FRAME%SIZE(3)
          ENDDO
        ENDIF
      ENDIF

      ! write(*,2) 'PGUIDE_INIT'
      ! write(*,2) 'HOR: ', OBJ%NLH, OBJ%CH, OBJ%AH(0:OBJ%NLH),OBJ%LH(0:OBJ%NLH)
      ! write(*,2) 'VER: ', OBJ%NLV, OBJ%CV, OBJ%AV(0:OBJ%NLV),OBJ%LV(0:OBJ%NLV)


      END SUBROUTINE PGUIDE_INIT

C---------------------------------------------------------------------
      REAL*8 FUNCTION PGUIDE_LAM(OBJ,ID,IL,Z,IORI)
C function describing lamella profile
C ID   ... derivative
C IL   ... lamella index
C Z    ... distance from guide entry
C IORI  ... horizontal (1) or vertical (2)
C----------------------------------------------------------------------
      USE GUIDES
      IMPLICIT NONE
      TYPE(GUIDE) :: OBJ

      REAL*8 a,f,zz,aa,x0
      REAL*8 Z
      INTEGER*4 ID,IL,IORI

      IF (IORI.EQ.2) THEN   ! vertical slit
        f=OBJ%CV
        a=OBJ%AV(IL)
        x0=OBJ%LV(IL)
      ELSE                 ! horizontal slit
        f=OBJ%CH
        a=OBJ%AH(IL)
        x0=OBJ%LH(IL)
      ENDIF

      IF (f.LT.0) THEN     ! focal point before the guide
        zz=z-f+abs(a)
      ELSE IF (f.GT.0) THEN   ! focal point behind the guide
        zz=-z+f+OBJ%FRAME%SIZE(3)+abs(a)
      ELSE
        zz=0
      ENDIF
      aa=SIGN(1.D0,a)*SQRT(abs(a))

      IF (zz.LE.0.D0.AND.f.NE.0) THEN
        PGUIDE_LAM=0.D0
      RETURN
      ENDIF

c zero deriv.
      IF (ID.LE.0) THEN
        IF (f.EQ.0) THEN
          PGUIDE_LAM=x0+a*z
        ELSE
          PGUIDE_LAM=2.D0*aa*SQRT(zz)
        ENDIF
c 1st deriv.
      ELSE IF (ID.EQ.1) THEN
        IF (f.EQ.0) THEN
          PGUIDE_LAM=a
        ELSE
          PGUIDE_LAM=-aa*SIGN(1.D0,f)/SQRT(zz)
        ENDIF
c 2nd deriv.
      ELSE IF (ID.EQ.2) THEN
        IF (f.EQ.0) THEN
          PGUIDE_LAM=0.D0
        ELSE
          PGUIDE_LAM=-aa/SQRT(zz)/zz/2.D0
        ENDIF
      ELSE
        PGUIDE_LAM=0.D0
      ENDIF

      END FUNCTION PGUIDE_LAM


C-------------------------------------------------------------------------------
      REAL*4 FUNCTION PGUIDE_A(il,nl,w,f)
C parabola parameter for il-th lamella
C il ... lamella index, counted from right, right side il=0
C f  ... focal distance
C w  ... width of the entry
C nl ... number of slits (number of lamellae + 1)
C sign(a) determines which side from the guide center: right/bottom(<0) or left/top(>0)
C--------------------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 f,w
      REAL*4 xx
      INTEGER*4 il,nl

      xx=w*(il*1.E0/nl-0.5D0)
      IF (abs(xx).LE.1.E-10) THEN
         PGUIDE_A=0.D0
      ELSE
         PGUIDE_A=SIGN(1.E0,xx)*(SQRT(f**2+xx**2)-abs(f))/2.E0
      ENDIF
c      Z=SIGN(1.E0,xx)*(SQRT(f**2+xx**2)-abs(f))/2.E0

      END FUNCTION PGUIDE_A


C----------------------------------------------------------------------------
      REAL*8 FUNCTION PGUIDE_CROSS(kx,kz,x,z,a,f,tmax)
C return TOF to a cross-point with parabolic surface
C if there are 2 solutions, take the minimum one which is > 0
C return 10^30 if there is no solution 0 < result < lmax
C a      .. abs(parabola parameter)
C f      .. focal distance from guide entry (origin of z coordinate)
C tmax   .. limit time (TOF to guide exit)
C x,z    .. starting point (z // guide axis)
C kx,kz  .. ray direction
C----------------------------------------------------------------------------
      use TRACINGDATA
      IMPLICIT NONE
      REAL*8 EPS
      PARAMETER (EPS=1.D-7)
      REAL*8 kx,kz,x,z,a,f,tmax
      REAL*8 T,T1,T2,tang,s,aa,dtm,ab,SC1,SC2

c10    FORMAT(a10,1x,8(1x,G13.7))
c      write(*,10) 'cross',kx,kz,x,z,a,f,tmax

      T=1.D30
      ab=abs(a)
! beam along z
      IF (ABS(kx).LT.1.D-6*abs(kz)) THEN
        IF (ab.GT.EPS) THEN
          T=(f+SIGN(1.D0,f)*(ab-x**2/ab/4.D0)-z)/kz
        ENDIF
      ELSE
        tang=kz/kx
        dtm=ab*(1+tang**2) + SIGN(1.D0,f)*(x*tang+f-z)
! no solution
        IF (dtm.GT.0.D0) THEN
! two solutions
          s=2*SIGN(1.D0,f)*SIGN(1.D0,kx)
          aa=SQRT(ab)
          T1=(s*aa*(-aa*kz/abs(kx)+SQRT(dtm))-X)/kx
          T2=(s*aa*(-aa*kz/abs(kx)-SQRT(dtm))-X)/kx
          SC1=a*(X+kx*T1)
          SC2=a*(X+kx*T2)
          if (SC1.LT.0.D0) T1=1.D30
          if (SC2.LT.0.D0) T2=1.D30
! choose the minimum positive value
          if ((T1.GT.EPS).and.(T2.GT.EPS)) then
            ! apply sign condition
            T=MIN(T1,T2)
          else if (T1.GT.EPS) then
            T=T1
          else if (T2.GT.EPS) then
            T=T2
          ENDIF
        ENDIF
      ENDIF
c      write(*,10) 'cross',T1,T2,T,tmax
c      pause
      IF (T.GT.tmax) T=1.D30
      PGUIDE_CROSS=T

      END FUNCTION PGUIDE_CROSS


C---------------------------------------------------------------------
      SUBROUTINE PGUIDE_CON(OBJ,R,K,IH,IV,IC,DT,Q,NW,KG,P)
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
C IC    .. which wall,  right(1), left(2), bottom (3) or top(4)
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
      INTEGER*4 IH,IV,IC,IX
      REAL(kind(1.D0)) :: PGUIDE_CROSS,GUIDE_CROSS1,PGUIDE_LAM
      REAL(kind(1.D0)) :: kx,kz,x,z,a,f,tmax,ZL
      REAL(kind(1.D0)) :: LZ(4),ANG,T0,GX,N(3)
      INTEGER*4 I
      REAL(kind(1.D0)) :: OSCD,OSCA
      COMMON /OSCBEND/ OSCD,OSCA
      LOGICAL :: DBG=.false.
      integer,parameter :: ID(4)=(/1,1,2,2/)
      integer,parameter :: IDINV(4)=(/2,2,1,1/)

    !  DBG=(OBJ%FRAME%COUNT.gt.99900)
    !  DBG=(DBG.and.(R(1).gt.20.D0))

1     FORMAT(a,1x,6(1x,G11.5))
10    FORMAT(a11,1x,6(1x,G11.5),2(1x,I4),2(1x,G10.4))
      ZL=0.D0
      KG=K
      z=R(3)
      kz=K(3)

! TOF to the guide end
      if (TRACING_UP) then
        T0=-z/kz
      else
        T0=(OBJ%FRAME%SIZE(3)-z)/kz
      endif
      if (DBG) write(*,10) 'CON START: ',R,K,IH,IV,OBJ%CH
      if (DBG) write(*,10) '   AH, T0: ',OBJ%AH(0),OBJ%AH(1),T0

C HORIZONTAL RIGHT:
      f=OBJ%CH
      kx=K(1)
      gx=NEUT_GR(1)
      tmax=T0
      IF (f.EQ.0) THEN
        a=OBJ%AH(IH)-OSCA
        x=R(1)-OBJ%DLH/2.D0-OBJ%LH(IH)-OSCD-R(3)*OSCA
        !LZ(2)=GUIDE_CROSS(kx,kz,x,z,a,f)
        LZ(2)=GUIDE_CROSS1(kx,kz,x,z,gx,a,f)
      ELSE
        a=OBJ%AH(IH)
        x=R(1)-OBJ%DLH/2.D0
        if (f.GT.0) f=f+OBJ%FRAME%SIZE(3)
        if ((f.LT.0).and.(.not.TRACING_UP)) tmax=(OBJ%LH(IH)-z)/kz
        LZ(2)=PGUIDE_CROSS(kx,kz,x,z,a,f,tmax)
        if (DBG) write(*,1) 'right: ',a,x,f,LZ(2)
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
         if ((f.LT.0).and.(.not.TRACING_UP)) tmax=(OBJ%LH(IH+1)-z)/kz
         LZ(1)=PGUIDE_CROSS(kx,kz,x,z,a,f,tmax)
         if (DBG) write(*,1) 'left: ',a,x,f,LZ(1)
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
         if (f.GT.0) f=f+OBJ%FRAME%SIZE(3)
         if ((f.LT.0).and.(.not.TRACING_UP)) tmax=(OBJ%LV(IV)-z)/kz
         LZ(4)=PGUIDE_CROSS(kx,kz,x,z,a,f,tmax)
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
         if ((f.LT.0).and.(.not.TRACING_UP)) tmax=(OBJ%LV(IV+1)-z)/kz
         LZ(3)=PGUIDE_CROSS(kx,kz,x,z,a,f,tmax)
      ENDIF

      if (DBG) write(*,1) 'times: ',LZ(1:4)

      DT=MIN(LZ(1),LZ(2),LZ(3),LZ(4),T0)

      IC=0
      DO I=1,4
        IF (DT.EQ.LZ(I)) IC=I
      ENDDO
      IF (IC.EQ.0) THEN ! no contact, passed through
        DT=T0
        Q=0.D0
        return
      ENDIF

C get the surface normal vector
      N=0.D0
      select case(IC)
      case(2)
        ANG=PGUIDE_LAM(OBJ,1,IH,Z+DT*kz,1)
        N(3)=-ANG/SQRT((1.D0+ANG**2))
        N(1)=SQRT(1.D0-N(3)**2)
      case(1)
        ANG=PGUIDE_LAM(OBJ,1,IH+1,Z+DT*kz,1)
        N(3)=ANG/SQRT((1.D0+ANG**2))
        N(1)=-SQRT(1.D0-N(3)**2)
      case(4)
        ANG=PGUIDE_LAM(OBJ,1,IV,Z+DT*kz,2)
        N(3)=-ANG/SQRT((1.D0+ANG**2))
        N(2)=SQRT(1.D0-N(3)**2)
      case(3)
        ANG=PGUIDE_LAM(OBJ,1,IV+1,Z+DT*kz,2)
        N(3)=ANG/SQRT((1.D0+ANG**2))
        N(2)=-SQRT(1.D0-N(3)**2)
      END select

      IX=ID(IC)
! stop events reflected from concave side
      IF ((OBJ%ONESIDE.eq.1).and.(IC.GT.0)) then
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

      END SUBROUTINE PGUIDE_CON


!----------------------------------------------------
      SUBROUTINE PGUIDE_MOVE_ENTRY(OBJ)
! Move neutron to the guide entry point (for TYPE=3, variable lamella length)
!----------------------------------------------------
      USE GUIDES
      USE TRACINGDATA
      TYPE(GUIDE) :: OBJ
      REAL(kind(1.D0)) :: DT
      REAL(kind(1.D0)), parameter :: EPS=1.D-10
      REAL(kind(1.D0)) :: lx,xx,fx,sx,sk,LD,Z,DETM,DTH,DTV
! Guide with varying lamellae lengths, traced from the wider side
      IF ((OBJ%TYP.EQ.3).AND.(.not.TRACING_UP)) then
        DTH=0.D0
        DTV=0.D0
        LD=OBJ%FRAME%SIZE(3)-NEUT%R(3) ! distance to the guide exit
! try horizontal
        IF ((OBJ%CH.GT.0).AND.(OBJ%NLH.GT.2)) THEN ! focus is behind exit
          xx=NEUT%R(1)
          sx=sign(1.D0,xx)
          sk=sign(1.D0,NEUT%K(1))
          fx=2.D0*OBJ%CH*OBJ%FRAME%SIZE(1)/OBJ%NLH
          lx=OBJ%FRAME%SIZE(3)
          IF(ABS(xx).GT.EPS) lx=MIN(fx/ABS(xx),lx) ! "optimum" lamella length
          IF (lx.LT.LD) THEN
            Z=LD/NEUT%K(3)-xx/NEUT%K(1)
            DETM=Z**2+4.D0*(LD*xx-fx*sx)/NEUT%K(1)/NEUT%K(3)
            DTH=(Z+sx*sk*SQRT(DETM))/2.D0   ! positive solution of quadratic equation
          ENDIF
        ENDIF
! try vertical
        IF ((OBJ%CV.GT.0).AND.(OBJ%NLV.GT.2)) THEN ! focus is behind exit
          xx=NEUT%R(2)
          sx=sign(1.D0,xx)
          sk=sign(1.D0,NEUT%K(2))
          fx=2.D0*OBJ%CV*OBJ%FRAME%SIZE(2)/OBJ%NLV
          lx=OBJ%FRAME%SIZE(3)
          IF(ABS(xx).GT.EPS) lx=MIN(fx/ABS(xx),lx) ! "optimum" lamella length
          IF (lx.LT.LD) THEN
            Z=LD/NEUT%K(3)-xx/NEUT%K(2)
            DETM=Z**2 + 4.D0*(LD*xx-fx*sx)/NEUT%K(2)/NEUT%K(3)
            DTV=(Z+sx*sk*SQRT(DETM))/2.D0   ! positive solution of quadratic equation
          ENDIF
        ENDIF
        DT=MIN(DTH,DTV)
! move to the real guide entry
        call TransportG(DT)
      endif
      END SUBROUTINE PGUIDE_MOVE_ENTRY

