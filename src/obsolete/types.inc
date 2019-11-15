C// Definition of TAS components - structures
C// $Id: types.inc,v 1.1 2007/08/03 15:46:24 saroun Exp $
C// shape ID's
      integer, PARAMETER :: SHP_ELLIPSOID = 0
      integer, PARAMETER :: SHP_CYLINDER = 1  ! cylinder, axis // y
      integer, PARAMETER :: SHP_DISC = 2      ! cylinder, axis // z
      integer, PARAMETER :: SHP_BOX = 3

C// shape identificator
      character*42 SHP_NAMES
      parameter(SHP_NAMES='ellipsoid:cylinder:disc:box')

C// class ID's
      integer, PARAMETER :: CLASS_SLIT = 1
      integer, PARAMETER :: CLASS_SOURCE = 2
      integer, PARAMETER :: CLASS_DETECTOR = 3
      integer, PARAMETER :: CLASS_SAMPLE = 4
      integer, PARAMETER :: CLASS_GUIDE = 5
      integer, PARAMETER :: CLASS_CRYSTAL = 6

C// class identificator (keep the same order as class ID !)
      character*42 CLASS_NAMES
      parameter(CLASS_NAMES='SLIT:SOURCE:DETECTOR:SAMPLE:GUIDE:CRYSTAL')

C SLIT:
      TYPE  SLIT; SEQUENCE
        CHARACTER*32  NAME
        REAL*8 AXX,AXY,AXZ,DIST          ! exit axis angles, distance
        REAL*8 SIZE(3),STA(3),GON(3)     ! dimensions, linear stage, gonio
        REAL*8 VEL(3)                    ! velocity (assumed constant)
        REAL*8 RLOC(3,3)           ! rotation from incident axis to local coordinates
        REAL*8 REXI(3,3)           ! rotation from incident to exit axis coordinates
        INTEGER*4 MLOC,MEXI        ! mask (=0 if RLOC or REXI are not required, else 1)
        INTEGER*4 COUNT                  ! event counter
        INTEGER*4  SHAPE     ! shape ID
        INTEGER*4 CLASS ! class ID
        INTEGER*4 DUM   ! dummy item for alignment
      END TYPE SLIT    ! size = 328

C DETECTOR
C TYP = 0        .. area, inf. thin, 100% efficiency, 0 mm resolution
C TYP = 1        .. array of ND tubes, efficiency=ALPHA
C TYP = 2        .. PSD with spatial resolution = RESX,RESY
      TYPE  DETECTOR; SEQUENCE
          TYPE(SLIT) :: FRAME
          INTEGER*4 TYP ! detector type
          REAL*8 RESX,RESY ! spatial resolution [mm], gaussian
          REAL*8 ALPHA ! efficiency coeff. [1/A/cm]
          REAL*8 SPACE ! distance between adjacent tubes
          INTEGER*4 ND ! number of tubes in an array
      END TYPE DETECTOR    ! size = 360

C BENDER:
C TYP < 0        .. ignore
C TYP = 0        .. normal course or soller collimator
C TYP = 1        .. guide or bender (reflecting walls)
C TYP = 2        .. parabolic guide
C TYP = 3        .. parabolic guide with optimum lamella lengths
C TYP = 4        .. elliptic guide, parallel lamellae at the wider end
      TYPE  BENDER; SEQUENCE
           TYPE(SLIT) :: FRAME
           REAL*8 W2,H2           ! width & height at the exit
           REAL*8 GHLU,GHLD,GHRU,GHRD,GVT,GVB  ! horizontal & vertical (H&V) critical angles
           REAL*8 RHLU,RHLD,RHRU,RHRD,RVT,RVB  ! horizontal& vertical reflectivities
           REAL*8 DLH,DLV         ! thickness of H&V lamellae
           REAL*8 ROH,ROV         ! curvatures
           real*8 MU              ! absorption coeff. (1/mm) for lamellae
           INTEGER*4 NLH,NLV      ! number of slits (H&V), max.=127
           INTEGER*4 TYP          ! normal(0), guide (1), parabolic (2,3), elliptic (4)
           INTEGER*4 NHLU,NHLD,NHRU,NHRD,NVT,NVB ! indexes pointing to the field containg reflectivity data
           INTEGER*4 OSCILATE     ! oscillating collimator, if > 0
           INTEGER*4 ONESIDE      ! reflection only from inner sides of lamellae, if >0
           INTEGER*4 TRFRONT      ! transparent front surface of lamellae, if >0
c calculated fields
           REAL*8 CH,CV           ! calculated curvatures, parabolic, elliptic, etc guides
           REAL*4 AH(0:127),AV(0:127) ! parameters of lamellae (for parabolic profile)
           REAL*4 LH(0:127),LV(0:127) ! lengths of lamellae (for parabolic profile)
      END TYPE BENDER   ! size = 2568

C CRYSTAL
      TYPE CRYSTAL; SEQUENCE
           TYPE(SLIT) :: FRAME
           character*8 REFNAME
           REAL*8 DHKL,CHI,POI,ANIZ ! dhkl,cutting angle, Poisson const., mos. anisotropy
           REAL*8 HMOS,VMOS,QML,VOL,MI,DW ! mosaicity (H,V), M.L. reflectivity, cell volume, absorption, Debye-Waller
           REAL*8 MAG !  magnetization
           REAL*8 RH,RV,RB,DH,DV,DB ! curvatures, segemnts spacing (H,V,z)
           REAL*8 sigmab, sigmaa, A, thetaD,C2 ! Freund parameters for absorption (see crystals.lib)
           REAL*8 DLam,Dext,Ext1,DELTA         ! extinction parameters: lamellae thickness, ext. length, primary ext., Darwin box
           REAL*8 THB,LAMBDA,STMCH,CTMCH,QHKL,REF,DETA,TANVOL ! Bragg angle, wavelength, sin(t), cos(t), kin. refl., ...
           REAL*8 DG_DR(3,3),G(3),GAMA(3),GTOT,DGR,DGA ! deformation gradient, diff. vector, d-gradient, gradient orientation
           LOGICAL*4 MAPG(3)
           INTEGER*4 DNRND ! base index to random numbers array
           INTEGER*4 nh,nv,nb,nref   ! number of segments (H,V,z), dummy
      END TYPE CRYSTAL            ! size = 776

C SAMPLE
C TYP = 0        .. generic, inelastic resolution
C TYP = 1        .. generic, elastic resolution
C TYP = 2        .. Vanad
C TYP = 3        ..
      TYPE  SAMPLE; SEQUENCE
          TYPE(SLIT) :: FRAME
          INTEGER*4 TYP ! sample type
          REAL*8 SIGSC  ! total coherent scattering cross-section [1/cm]
          REAL*8 SIGA   ! total absorption cross-section [1/cm]
      END TYPE SAMPLE   !

C SETUP - scattering triangle etc..
      TYPE SETUP; SEQUENCE
         INTEGER*4 NFX,SM,SS,SA
         REAL*8 KFIX,E,Q,KI,KF,EI0,EF0
      END TYPE SETUP
