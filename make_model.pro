PRO MAKE_MODEL,GRID,Teff,Logg,Vmag,FWHM,SAMP,WINDOW,W,F
;******************************************************************************
; MAKE_MODEL creates a calibrated model derived from GRID, defined by Vmag,   *
; Teff, Logg, at resolution FWHM, with sample size = SAMP and over window =   *
; Window.                                                                     *
; JBH  9/18/07                                                                *
;******************************************************************************

TGRID = GRID.t
GGRID = GRID.g
WAVE  = GRID.wl
FLUX  = GRID.fl
MOD_CUBIC,Teff,Logg,TGRID,GGRID,WAVE,FLUX,FOUT

GSMOOTH,WAVE,FOUT,FS,FWHM

MOD_SIGMA,'V5423.GRIDN',Teff,Logg,Vmag,SIGMA

NS = FIX( (WINDOW(1)-WINDOW(0))/SAMP )
WS = WINDOW(0) + SAMP*FINDGEN(NS+1)
LINTERP,WAVE,FS,WS,FX
TRIM,WINDOW(0),WINDOW(1),WS,FX

W = WS
F = FX*SIGMA

RETURN
END
