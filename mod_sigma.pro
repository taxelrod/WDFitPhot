PRO MOD_SIGMA,VGRID,T_eff,Log_g,V_mag,SIGMA
;*************************************************************************
; MOD_SIGMA returns the interpolated stellar solid angle ( 4PI*{R/D}^2 ) *
; for a given T_eff, Log_g and V_mag.                                    *
; includes new Holbeg & Bergeron Flux normalization at lam=5423A.        *                                    *
; 09/18/07                                                               *
;*************************************************************************

;=========================================================================
; DEFINITIONS
;=========================================================================
;VGRID ..................... INPUT Grid File Name
;T_eff ..................... INPUT Effective Temperature
;Log_g ..................... INPUT Log Surface Gravity
;V_mag ..................... INPUT V Magnitude
;SIGMA ..................... OUTPUT Stellar Solid Angel...4PI*{R/D}^2
;========================================================================
; SUBROUTINES
;========================================================================
; SPLINE.................... Spline Interpolation
;========================================================================

path  = '/home/tsa/Dropbox/WD/GauthamCode/DATA/'
ngrid = path + VGRID

RESTORE,ngrid                       ; HEADLAB,PAR,TGRID,GGRID,FN

;========================================================================
; Part 1 Determine normalizing wavelength find observed flux
;========================================================================

WNORM = FLOAT( STRMID (VGRID,1,4) )
print,WNORM
wnor  = [5000, 5423, 5490, 5500, 5510, 5556 ]  ; Standard Wavelengths
fnor  = [4.834,3.804,3.674,3.655,3.633,3.547]  ; Absolute Fluxes
fnor  = fnor*1.0E-9

f0    = fnor(where(WNORM EQ wnor))  ; select normalizing flux
fobs  = f0*10^(-0.4*V_mag)          ; compute observed flux

;========================================================================
; Part 2 Interpolate Grid of Eddington Fluxes
;========================================================================

IMAX     = N_ELEMENTS(TGRID)  &  JMAX = N_ELEMENTS(GGRID)
XINDEX   = FINDGEN(IMAX)            ; vector of temperature indices
YINDEX   = FINDGEN(JMAX)            ; vector of gravity indices

X  = ALOG10(TGRID)
Y  = GGRID
Z  = ALOG10(FN)
XT = ALOG10(T_eff)
YG = Log_g

TVAL = SPLINE(X,XINDEX,XT)          ; interpolated index for T_eff
GVAL = SPLINE(Y,YINDEX,YG)          ; interpolated index for Log_g

LHLAM = INTERPOLATE(Z,TVAL,GVAL,/GRID,CUBIC=-0.5)

D = fobs/(10^LHLAM)
SIGMA = D(0)

RETURN
END










