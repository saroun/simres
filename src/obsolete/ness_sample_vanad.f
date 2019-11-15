
C      ------------------------------------------
      LOGICAL*4 FUNCTION VAN_GO(VAN,NEUI,NEUF,Q)
C      ------------------------------------------
      implicit NONE

        INCLUDE 'const.inc'
        INCLUDE 'structures.inc'
        INCLUDE 'randvars.inc'
        INCLUDE 'tracing.inc'

        REAL*8 SSCAT,SABS
c// scattering and absorption cross-section for k=2.664A-1
        PARAMETER (SSCAT=0.0362,SABS=0.0476)
c        PARAMETER (SSCAT=0.0362,SABS=0.)

      RECORD /SLIT/ VAN
      RECORD /NEUTRON/ NEUI,NEUF
      LOGICAL*4 LOG1,SAM_BORDER
      REAL*8 V(3),K(3),R(3),Q,X1(3)
      REAL*8 KF,PP,P0,T1,T2,DELTA,stb,ctb,s2tb,c2tb,sinn,coss,ksi
        REAL*4 DT
      REAL*4 RAN1
        INTEGER*4 I

      LOG1=.TRUE.
        PP=1
        CALL SLIT_PRE(VAN,NEUI.R,NEUI.K,V,K)

C/// move neutron to the centre of the SLIT:
      NEUF.T=NEUI.T-V(3)/HOVM/K(3)
      DO 10 I=1,2
10      R(I)=V(I)-V(3)/K(3)*K(I)
        R(3)=0.

      IF (SAM_BORDER(VAN,R,K,T1,T2)) THEN
C/// move neutron to the sample entry:
        NEUF.T=NEUI.T+T1/HOVM
        DO I=1,3
          R(I)=R(I)+K(I)*T1
        ENDDO
        DELTA=T2-T1
        kf=SQRT(K(1)**2+K(2)**2+K(3)**2)
           ksi=RAN1()
           P0=1.-exp(-SSCAT*kf*DELTA) ! scattering probability
         DT=-Log(1-ksi*P0)/SSCAT/kf
C/// move to the scattering point:
           NEUF.T=NEUF.T+DT/HOVM
         DO I=1,3
           R(I)=R(I)+K(I)*DT
           ENDDO

           stb=Q/2/kf
           if (abs(stb).ge.1) then
              NEUF.P=0
              VAN_GO=.FALSE.
              RETURN
           endif
           ctb=SQRT(1-stb**2)
           s2tb=sin(XRND(8)+2*ATAN(stb/ctb))
           c2tb=SQRT(1-s2tb**2)
           if (stb**2.gt.0.5) c2tb=-c2tb
           sinn=sin(XRND(7))
           coss=SQRT(1-sinn**2)
           IF(ABS(XRND(7)).GT.PI/2) COSS=-COSS
           X1(1)=kf*s2tb*coss
           X1(2)=kf*sinn
           x1(3)=kf*c2tb*coss
           CALL SLIT_POST(VAN,R,X1,NEUF.R,NEUF.K)

           IF (SAM_BORDER(VAN,R,X1,T1,T2)) THEN
              PP=EXP(-SABS*(DT+T2)*2.664)  ! absorption at k=kf,  (DT+T2)*kf = flight path
              NEUF.P=NEUI.P*P0*PP/4/PI  ! convert Sigma to dSigma/dOmega
            NEUF.S=NEUI.S
              LOG1=(NEUF.P.GT.0.D0)
            IF(LOG1) VAN.COUNT=VAN.COUNT+1
           ENDIF
      ELSE
          LOG1=.FALSE.
          NEUF.P=0
      END IF

      VAN_GO=LOG1
      RETURN
      END

C      ---------------------------------------------
      LOGICAL*4 FUNCTION VAN_TRANS(VAN,NEUI,NEUF)
C      ---------------------------------------------
      implicit NONE

        INCLUDE 'const.inc'
        INCLUDE 'structures.inc'

        REAL*8 SSCAT,SABS
        PARAMETER (SSCAT=0.0362,SABS=0.0476)
c        PARAMETER (SSCAT=0.0362,SABS=0.)

      RECORD /NEUTRON/ NEUI,NEUF
      RECORD /SLIT/ VAN
      LOGICAL*4 SAM_BORDER
      REAL*8 V(3),K(3),R(3),DT
      REAL*8 KF,T1,T2,EXP
        INTEGER*4 I

        CALL SLIT_PRE(VAN,NEUI.R,NEUI.K,V,K)

C/// move neutron to the centre of the SLIT:
        IF (K(3).GT.1D-5) THEN
          DT=-V(3)/K(3)
        ELSE
          DT=-V(1)/K(1)
        ENDIF
      NEUF.T=NEUI.T+DT/HOVM
      DO I=1,3
          R(I)=V(I)+DT*K(I)
        ENDDO

      IF (SAM_BORDER(VAN,R,K,T1,T2)) THEN
C/// transmission through the sample:
           kf=SQRT(K(1)**2+K(2)**2+K(3)**2)
           NEUF.P=NEUI.P*EXP(-(SABS+SSCAT)*kf*(T2-T1))
         NEUF.S=NEUI.S
      ELSE
          NEUF.P=NEUI.P
      END IF
        CALL SLIT_POST(VAN,R,K,NEUF.R,NEUF.K)

        VAN.COUNT=VAN.COUNT+1

      VAN_TRANS=.TRUE.

      RETURN
      END
