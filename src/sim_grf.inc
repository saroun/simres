      STRUCTURE /VIEWSET/
         INTEGER*4 IX,IY
         REAL*4 DX1,DX2,DY1,DY2    ! Viewport in device coordinates
         REAL*4 WX1,WX2,WY1,WY2    ! World coordinates
         CHARACTER*32 XTIT,YTIT,HEAD
      END STRUCTURE
      CHARACTER*72 GRFSTR(10)
      INTEGER*4 TOPRINT,GRFNARG
      REAL*8 GRFARG(40)
      COMMON /GRFSTAT/ GRFSTR,GRFARG,GRFNARG,TOPRINT

C// device settings
      CHARACTER*128 DEVSTR ! PGPLOT output device, e.g. "/xserv"
      INTEGER*4 GRFSAVE ! level of ASCII output
      INTEGER*4 DEVID  ! device ID (result of PGQINF)
      COMMON /GRAPH/ DEVSTR,GRFSAVE,DEVID  !//,EXTPLOT
