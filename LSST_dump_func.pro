PRO LSST_dump_func,T,G,gobs,EBV,Nor,V,Wmod,Fmod
;****************************************************************************
;* LSST_dump_func is based on LSST_dump.  Instead of writing the flux
;* to a file, returns it as a vector.  Code to read input data from a
;* star has been stripped out for efficiency.  These must be supplied
;* as args.
;* Only the model flux is returned. The reddened flux is ignored
;* 05/26/15
;***************************************************************************

COMMON GRIDS, GRIDNUV, GRIDBALMER, GRIDIR

lam = 4670    &   f0 = 4.980E-9    &   DM = -0.0023

 
;********************** FILES  *********************************
 
path1    = 'RunData/'

BALMER   = path1 + 'BALMER.GRIDN'
IR       = path1 + 'IR.GRIDN'
NUV      = path1 + 'NUV.GRIDN'

;***************  Set UP Make Model *****************************

Teff =  T                     ;  Teff 
Logg =  G                     ;  Gravity
Vmag =  V                     ;  V magnitude
FWHM =  4.0                   ;  R = 1200
SAMP =  1.0                   ;  Sampling
WINDOW = [1350,3500]          ;  Wavelength Window

IF ~ISA(GRIDNUV) THEN BEGIN
   RESTORE, NUV
   GRIDNUV={fl:FLUX,g:GGRID,hd:HEADLAB,p:PAR,t:TGRID,wl:WAVE}
ENDIF

IF ~ISA(GRIDBALMER) THEN BEGIN
   RESTORE, BALMER
   GRIDBALMER={fl:FLUX,g:GGRID,hd:HEADLAB,p:PAR,t:TGRID,wl:WAVE}
ENDIF

IF ~ISA(GRIDIR) THEN BEGIN
   RESTORE, IR
   GRIDIR={fl:FLUX,g:GGRID,hd:HEADLAB,p:PAR,t:TGRID,wl:WAVE}
ENDIF

MAKE_MODEL,GRIDNUV,Teff,Logg,gobs,FWHM,SAMP,WINDOW,WN,FN

WINDOW = [3501,7000]          ;  Wavelength Window
MAKE_MODEL,GRIDBALMER,Teff,Logg,gobs,FWHM,SAMP,WINDOW,WB,FB

WINDOW = [7001,25000]         ;  Wavelength Window
MAKE_MODEL,GRIDIR,Teff,Logg,gobs,FWHM,SAMP,WINDOW,WI,FI

wmod = [WN,WB,WI]
fmod = [FN,FB,FI]

EXTN = fitzpatrick(wmod,3.1)*EBV
frd  = 10.0^(0.4*EXTN)
fred = fmod/frd               ; reddened model
zg    = WHERE(lam EQ wmod)
arg   = fred(zg)
fmodg = arg(0)
gsyn  = gobs - DM             ; g on synthetic scale

fsyn  = f0*10^(-0.4*gsyn)

nor   = fsyn/fmodg
help,fsyn,fmodg
fred  = fred*nor
fmod  = fmod*nor

RETURN
END



