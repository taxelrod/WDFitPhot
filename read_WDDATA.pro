PRO read_WDDATA,filename,star,T,G,EBV,HSTObs, HSTObsUnc, Tunc, Gunc, bandDict
;**************************************************************************
;* read_WDDATA reads the basic WDDATA file and selects data for a given
;* star.
;* 10/03/13
;**************************************************************************
CLOSE,1
OPENR,1,filename

N = FILE_LINES(filename) - 1               ; number of data lines

hline =''
READF,1,hline         ; header line
; get the bands in use and their order
keyw = STRSPLIT(hline,/EXTRACT)
nHSTbands = (N_ELEMENTS(keyw)-6)/2

HSTObs = FLTARR(nHSTbands)
HSTObsUnc = FLTARR(nHSTbands)

bandDict = LIST()
FOR I = 0, nHSTbands-1 DO BEGIN
   band = keyw(I+4)
   bandDict.add,band
ENDFOR

IF star EQ '' THEN BEGIN
   return
ENDIF

WD = STRARR(N)
Teff = FLTARR(N)
lgg  = FLTARR(N)
EXT  = FLTARR(N)
HST  = FLTARR(N,nHSTbands)
HSTUnc  = FLTARR(N,nHSTbands)
sigmaT = FLTARR(N)
sigmaG = FLTARR(N)

FOR I  = 0,N-1  DO BEGIN
READF,1,hline
WRD      = STRSPLIT(hline,/EXTRACT)
WD(I)    = WRD(0)
Teff(I)  = FLOAT(WRD(1))
lgg(I)   = FLOAT(WRD(2)) 
EXT(I)   = FLOAT(WRD(3))

J=4
FOR L=0,nHSTbands-1 DO BEGIN
HST(I,L) = FLOAT(WRD(J))
J=J+1
ENDFOR

FOR L=0,nHSTbands-1 DO BEGIN
HSTUnc(I,L) = FLOAT(WRD(J))
J=J+1
ENDFOR

sigmaT(I) = FLOAT(WRD(J))
sigmaG(I) = FLOAT(WRD(J+1))
ENDFOR

Z = WHERE(star EQ WD)
IF Z EQ -1 THEN BEGIN
   PRINT,star, '  not found'
   STOP
ENDIF

;*************  Convert Output into Scalars and Velocity Vector *****

vec1 = Teff(Z)  & T    = vec1(0)
vec2 = lgg(Z)   & G    = vec2(0)
vec4 = EXT(Z)   & EBV  = vec4(0)
HSTObs = HST(Z,*)
HSTObsUnc = HSTUnc(Z,*)
vec6 = sigmaT(Z)  & Tunc    = vec6(0)
vec7 = sigmaG(Z)  & Gunc    = vec7(0)

CLOSE,1
RETURN
END
