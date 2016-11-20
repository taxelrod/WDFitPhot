PRO MOD_CUBIC,T_EFF,Log_g,TGRID,GGRID,WAVE,FLUX,FOUT
;************************************************************************
; MOD_CUBIC cubically interpolates a synthetic spectrum corresponding   *
; to a given T_eff at fixed Log_g from the input grid 'FLUX'.           *
; 11/04/99                                                              *
;************************************************************************
;========================================================================
; DEFINITIONS
;========================================================================
;T_eff .................. INPUT Effective Temperature  (scalar)
;Log_g .................. INPUT Log Surface Gravity    (scalar)
;TGRID .................. INPUT Temperature Vector     (1-D)
;GGRID .................. INPUT Gravity Vector         (1-D)
;WAVE.................... INPUT Wavelength Vector      (1-D)
;FLUX ................... INPUT Flux Grid              (2-D)
;FOUT ................... OUTPUT Flux Grid Vector      (1-D)
;========================================================================
; SUBROUTINES 
;========================================================================
; SPLINE ................ Spline Interpolation
;========================================================================

;========================================================================
; PART 1  Convert to Logs and Determine Temperature & Gravity Indices
;========================================================================

IMAX = N_ELEMENTS(TGRID) 
JMAX = N_ELEMENTS(GGRID) 
KMAX = N_ELEMENTS(WAVE)

XINDEX = FINDGEN(IMAX)         ;  vector of temperature indices
YINDEX = FINDGEN(JMAX)         ;  vector of gravity  indices
WINDEX = FINDGEN(KMAX)         ;  vector of wavelength indices

X    = ALOG10(TGRID)  
Y    = GGRID
Z    = ALOG10(FLUX)

XT   = ALOG10(T_eff)
YG   = Log_g

TVAL = SPLINE(X,XINDEX,XT)     ; Interpolated index for T_eff
GVAL = SPLINE(Y,YINDEX,YG)     ; Interpolated index for Log_g

;========================================================================
; PART 2  Interpolate a Model at T_eff and Log_g
;========================================================================

LHLAM = INTERPOLATE(Z,TVAL,GVAL,WINDEX,/GRID,CUBIC=-0.5)
LHLAM = REFORM(LHLAM,/OVERWRITE)
FOUT  = 10^(LHLAM)

RETURN
END



