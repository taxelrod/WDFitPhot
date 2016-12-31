PRO GetModelFlux, T, logg, wl, flux

; Determine the Hubeny tgm file to use, based on T

modelFile = ''

IF (T ge 10000. and T lt 20000) THEN modelFile = 'RunData/da10-19.tgm.fits' ELSE $
IF (T ge 20000. and T lt 40000) THEN modelFile = 'RunData/da20-39.tgm.fits' ELSE $
IF (T ge 40000. and T lt 100000) THEN modelFile = 'RunData/da40-100.tgm.fits' ELSE BEGIN
   PRINT, 'Temp out of range: ', T
   RETURN
ENDELSE

wlFlux = uly_tgm_extr([T, logg, 0], MODEL_FILE=modelFile)

approxNorm = 4.980E-9 * 10^(-0.4*16) / max(wlFlux[1]) ; They're about mag 16
wl = wlFlux[0]
flux = wlFlux[1]*approxNorm

END
