      INTEGER,PARAMETER :: MSEG=128
      integer,PARAMETER :: ctyp_filter=0
      integer,PARAMETER :: ctyp_mosaic=1
      integer,PARAMETER :: ctyp_bent=2
      integer,PARAMETER :: ctyp_gradient=3
      integer,PARAMETER :: cmos_gauss=0
      integer,PARAMETER :: cmos_lorentz=1
      integer,PARAMETER :: cmos_voigt=2
      integer,PARAMETER :: cmos_box=3
      INTEGER :: NSEG,MDIST
      REAL(KIND(1.D0)) :: TSEG(0:MSEG),TSEG0(0:MSEG),ALPHA(0:MSEG),GRAD(0:MSEG)
      REAL(KIND(1.D0)) :: PSEG(0:MSEG),PSEG0(0:MSEG),PATH,RDIST,QKIN
      COMMON /CRBORDERS/ TSEG,TSEG0,ALPHA,GRAD,PSEG,PSEG0,PATH,RDIST,QKIN,NSEG,MDIST

