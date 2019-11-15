C///  Variables for the universal command line interpreter LINP
      INTEGER*4 linp_dim
      PARAMETER(linp_dim=256)
      INTEGER*4 linp_nc,linp_np,linp_in,linp_out,linp_err,linp_eof
      CHARACTER*5 linp_c(linp_dim)
      CHARACTER*60 linp_h(linp_dim)
      CHARACTER*20 linp_p
      COMMON /linp/ linp_nc,linp_np,linp_in,linp_out,linp_err,linp_eof,
     & linp_p,linp_c,linp_h


