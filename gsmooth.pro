PRO GSMOOTH,WIN,FIN,FOUT,FWHM
;**************************************************************************
; GSMOOTH performs a gaussian smooth of Full-Width-Half-Max = FWHM on an  *
; input vector FIN and returns FOUT of the same size.  The convolution    *
; assumes equally spaced data and centers the convolution kernel over     * 
; each data point. The convolution kernel is +/- 5 sigma in extent.       *
; 11/04/99                                                                *
;**************************************************************************

;==========================================================================
; DEFINATIONS
;==========================================================================
; WIN .................  INPUT Wavelength
; FIN .................  INPUT Vector
; FOUT ................  OUTPUT Vector
; FWHM ................  Full-Width-Half-Max for Gaussian Kernel
;==========================================================================

;==========================================================================
; PART 1 Determine input sample size and vector size
;==========================================================================

N     = N_ELEMENTS(WIN)
FSAMP = (WIN(N-1) - WIN(0))/(N-1)

;==========================================================================
; PART 2 Define Gaussian Smoothing Kernel
;==========================================================================

GSIG = 0.5*FWHM/SQRT(ALOG(2.0))         ; gaussian sigma
M    = FIX(5.0*GSIG/FSAMP)              ; gaussian kernel half-size
S    = (FINDGEN(2*M+1) - M)*FSAMP       ; gaussian sampling vector
arg  = (S/GSIG)*(S/GSIG)
z    = WHERE(arg GT 20)                 ; correct potential underflow cond.
IF(MIN(z) NE -1) THEN arg(z) = 20.0
GW   = EXP(-arg)/(SQRT(!PI)*GSIG)       ; gaussian weight kernel
GTOT = TOTAL(GW)

;==========================================================================
; PART 3 Perform Convolution
;==========================================================================

FOUT = CONVOL(FIN,GW,/EDGE_TRUNCATE)/GTOT

RETURN
END

