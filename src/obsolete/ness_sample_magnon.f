C----------------------------------------------------------
      REAL*8 FUNCTION SQE_AMAG(VKI,VKF)
C inelastic scattering cross-section
C spin waves in antiferromagnetics
C----------------------------------------------------------
      USE COMPONENTS
      USE TRACINGDATA
      USE CONSTANTS
      USE IO
      IMPLICIT NONE
      INCLUDE 'ness_common.inc'
      INCLUDE 'rescal.inc'

      REAL*8 bf
      PARAMETER(bf=11.609) !Bose factor conversion kT -> meV

      REAL*8 VKI(3),VKF(3)

      REAL*8 KF0,KI0,VQ(3),WQ(3),TAU(3),VE
      REAL*8 r1,nw,Eq,z,wqs
      INTEGER*4 I,J


      REAL*8 CS,EGAP,GAMMA,SIG,EPS,TEMP

      DATA CS,EGAP,GAMMA,SIG,EPS,TEMP /1.D-5,1.65,0.2,1000,0.6,1.6/
c      DATA ITMP/0/

      KF0=VKF(1)**2+VKF(2)**2+VKF(3)**2
      KI0=VKI(1)**2+VKI(2)**2+VKI(3)**2
c* get Q in lab. coord.
      DO I=1,3
        VQ(I)=VKF(I)-VKI(I)
      ENDDO
C* subtract nominal values from k vectors
      VQ(3)=VQ(3)-STP.KF*COMEGA+STP.KI
      VQ(1)=VQ(1)-STP.KF*SOMEGA
C*  transform to C&N coord.
      DO I=1,3
         WQ(I)=0
         DO J=1,3
            WQ(I)=WQ(I)+MLC(J,I)*VQ(J)
         ENDDO
      ENDDO
C*  transform to r.l. coord.
      DO I=1,3
         VQ(I)=0
         DO J=1,3
            VQ(I)=VQ(I)+MRC(I,J)*WQ(J)
         ENDDO
      ENDDO
C*  get energy transfer
      VE=HSQOV2M*(KI0-KF0)
c*  get TAU
      DO I=1,3
         TAU(I)=NINT(RES_DAT(i_QH+I-1))
      ENDDO
c*  get propagation vector
      DO I=1,3
         VQ(I)=VQ(I)+RES_DAT(i_QH+I-1)-TAU(I)
      ENDDO

C// bose factor
      z=exp(-ABS(VE)*bf/temp)
      nw=z/(1-z)
      if (VE.GT.0) nw=nw+1

C// excitation energy
      wqs=EGap**2*(1+SIG*(VQ(1)**2+VQ(2)**2+EPS*VQ(3)**2)) ! dispersion law

C// dynamic structure factor
      Eq=wqs+Gamma**2
      Z=(VE**2-Eq)**2+(VE*Gamma)**2
      r1=VE*Gamma/((VE**2-Eq)**2+(VE*Gamma)**2)
      Z=CS*nw*r1

      SQE_AMAG=Z
      END


C----------------------------------------------------------
      REAL*8 FUNCTION SQE_AMAG1(VKI,VKF)
C inelastic scattering cross-section
C dispersion with finite width and curvature along QH
C----------------------------------------------------------
      USE COMPONENTS
      USE TRACINGDATA
      USE CONSTANTS
      IMPLICIT NONE
      INCLUDE 'ness_common.inc'
      INCLUDE 'rescal.inc'

      REAL*8 bf
      PARAMETER(bf=11.609) !Bose factor conversion kT -> meV

      REAL*8 VKI(3),VKF(3)

      REAL*8 KF0,KI0,VQ(3),WQ(3),TAU(3),VE
      REAL*8 r1,nw,Eq,z,wqs
      INTEGER*4 I,J


      REAL*8 CS,EGAP,GAMMA,SIG,EPS,TEMP

      DATA CS,EGAP,GAMMA,SIG,EPS,TEMP /1.D-5,1.65,0.2,1000,0,1.6/
c      DATA ITMP/0/

      KF0=VKF(1)**2+VKF(2)**2+VKF(3)**2
      KI0=VKI(1)**2+VKI(2)**2+VKI(3)**2
c* get Q in lab. coord.
      DO I=1,3
        VQ(I)=VKF(I)-VKI(I)
      ENDDO
C* subtract nominal values from k vectors
      VQ(3)=VQ(3)-STP.KF*COMEGA+STP.KI
      VQ(1)=VQ(1)-STP.KF*SOMEGA
C*  transform to C&N coord.
      DO I=1,3
         WQ(I)=0
         DO J=1,3
            WQ(I)=WQ(I)+MLC(J,I)*VQ(J)
         ENDDO
      ENDDO
C*  transform to r.l. coord.
      DO I=1,3
         VQ(I)=0
         DO J=1,3
            VQ(I)=VQ(I)+MRC(I,J)*WQ(J)
         ENDDO
      ENDDO
C*  get energy transfer
      VE=HSQOV2M*(KI0-KF0)
c*  get TAU
      DO I=1,3
         TAU(I)=NINT(RES_DAT(i_QH+I-1))
      ENDDO
c*  get propagation vector
      DO I=1,3
         VQ(I)=VQ(I)+RES_DAT(i_QH+I-1)-TAU(I)
      ENDDO

C// bose factor
      z=exp(-ABS(VE)*bf/temp)
      nw=z/(1-z)
      if (VE.GT.0) nw=nw+1

C// excitation energy
      wqs=EGap**2*(1+SIG*(VQ(1)**2)) ! dispersion law

C// dynamic structure factor
      Eq=wqs+Gamma**2
      Z=(VE**2-Eq)**2+(VE*Gamma)**2
      r1=VE*Gamma/((VE**2-Eq)**2+(VE*Gamma)**2)
      Z=CS*nw*r1

      SQE_AMAG1=Z
      END
