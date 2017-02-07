;
; This version processes all stars simultaneously
;
FUNCTION ExtinctFuncAll, paramVec
;
; function to use with CURVEFIT
;
; paramVec - (zp_i, R, Ebv)
; vals - synth_mag_{ni} + (zp_i + X_i(R, Ebv_n))
; 
  COMMON TopInfo,  nStars, nBp, nBpHST, bpData, bandList, bandDict, zp, Tstars, TMinStars, TMaxStars, Gstars, GMinStars, GMaxStars, zpStar, EBV, modelWl, modelFluxes, sampleHST

  COMMON ZpMinimizeInfo, fitResult, chisqRes, scaleFactor, idRS, fluxRS, obsMagRS

  COMMON ScaleInfo, scaleFactor2, scaleFactor2Inv, scaleFactor3, scaleFactor3Inv, scaleFactor275, scaleFactor275Inv

  iBp275w = bandDict['F275W']

  zpBand = DBLARR(nBpHST)
  zpBand(0:(nBpHST-2)) = scaleFactor*paramVec(0:(nBpHST-2))
  zpBand(nBpHST-1) = -TOTAL(zpBand(0:(nBpHST-2))) ; zpBand sums to zero
  R = scalefactor*paramVec(nBpHST-1)
  EBV = scaleFactor*paramVec(nBpHST:(nBpHST+nStars-1))
  zpStar = scaleFactor*paramVec((nBpHST+nStars):(nBpHST+2*nStars-1))
  Tstars = scaleFactor2*paramVec((nBpHST+2*nStars):(nBpHST+3*nStars-1))
  Gstars = scaleFactor3*paramVec((nBpHST+3*nStars):(nBpHST+4*nStars-1))
  bp275wShift = scalefactor275*paramVec(nBpHST+4*nStars)

;  PRINTF, 2, FORMAT='("bp275wShift: ",e15.4)', bp275wShift

  fitResult = DBLARR(nBpHST*nStars)
  synMag = DBLARR(nBpHST*nStars)
  modelFluxes = LIST()

  FOR n = 1, nStars DO BEGIN
     LSST_dump_func, Tstars[n-1], Gstars[n-1], 0, EBV[n-1], 1.0, 0, modelWL, flux
     modelFluxes.add, flux
     extincMag = ext_odonnell(modelWL,R)*EBV[n-1]*R
     extinc = 10^(-0.4*extincMag)

     FOR i = 1, nBpHST DO BEGIN
        idx = (n-1)*nBpHST+i-1
        bp = bpData[i + nBp - 1]
        IF i-1 eq ibp275w THEN BEGIN
           ishift = FIX(bp275wShift)
           frcshift = bp275wShift - ishift
           synMag1 = synphot2(modelWl, modelFluxes[n-1]*extinc, bp.wavelength, SHIFT(bp.throughput, ishift), 0)
           synMag2 = synphot2(modelWl, modelFluxes[n-1]*extinc, bp.wavelength, SHIFT(bp.throughput, ishift+1), 0)
           synMag(idx) = synMag1 + (synMag2 - synMag1)*frcShift
	   printf, 2, 'F275 slope:', synMag2 - synMag1
        ENDIF ELSE BEGIN
           synMag(idx) = synphot2(modelWl, modelFluxes[n-1]*extinc, bp.wavelength, bp.throughput, 0)
        ENDELSE
        fitResult(idx) = zpBand[i-1] + zpStar[n-1] + synMag(idx)

     ENDFOR
  ENDFOR

; Stddev of RS mags as proxy for reddening
;
  IF idRS NE !NULL THEN BEGIN
     fluxRS = modelFluxes[idRS-1]

     bp275w = bpData[nBp+bandDict['F275W']] 
     synMagRSf275w = synphot2(modelWl, fluxRS, bp275w.wavelength, SHIFT(bp275w.throughput, bp275wShift), 0)
     obsMagRSf275w = obsMagRS[bandDict['F275W']]

     bp336w = bpData[nBp+bandDict['F336W']]
;  synMagRSf336w = synphot2(modelWl, fluxRS, bp336w.wavelength, bp336w.throughput, zpBand[bandDict['F336W']])
     synMagRSf336w = synphot2(modelWl, fluxRS, bp336w.wavelength, bp336w.throughput, 0)
     obsMagRSf336w = obsMagRS[bandDict['F336W']]

     bp475w = bpData[nBp+bandDict['F475W']]
;  synMagRSf475w = synphot2(modelWl, fluxRS, bp475w.wavelength, bp475w.throughput, zpBand[bandDict['F475W']])
     synMagRSf475w = synphot2(modelWl, fluxRS, bp475w.wavelength, bp475w.throughput, 0)
     obsMagRSf475w = obsMagRS[bandDict['F475W']]

     bp625w = bpData[nBp+bandDict['F625W']]
;  synMagRSf625w = synphot2(modelWl, fluxRS, bp625w.wavelength, bp625w.throughput, zpBand[bandDict['F625W']])
     synMagRSf625w = synphot2(modelWl, fluxRS, bp625w.wavelength, bp625w.throughput, 0)
     obsMagRSf625w = obsMagRS[bandDict['F625W']]

     bp775w = bpData[nBp+bandDict['F775W']] 
     synMagRSf775w = synphot2(modelWl, fluxRS, bp775w.wavelength, bp775w.throughput, 0)
     obsMagRSf775w = obsMagRS[bandDict['F775W']]

     bp160w = bpData[nBp+bandDict['F160W']] 
     synMagRSf160w = synphot2(modelWl, fluxRS, bp160w.wavelength, bp160w.throughput, 0)
     obsMagRSf160w = obsMagRS[bandDict['F160W']]

     RSdelta = DBLARR(6)
     RSdelta[0] = synMagRSf275w - obsMagRSf275w + zpBand[bandDict['F275W']]
     RSdelta[1] = synMagRSf336w - obsMagRSf336w + zpBand[bandDict['F336W']]
     RSdelta[2] = synMagRSf475w - obsMagRSf475w + zpBand[bandDict['F475W']]
     RSdelta[3] = synMagRSf625w - obsMagRSf625w + zpBand[bandDict['F625W']]
     RSdelta[4] = synMagRSf775w - obsMagRSf775w + zpBand[bandDict['F775W']]
     RSdelta[5] = synMagRSf160w - obsMagRSf160w + zpBand[bandDict['F160W']]

     RSred = TOTAL((RSdelta - MEAN(RSdelta))^2) / 0.01^2
  ENDIF ELSE BEGIN
     RSred = 0
  ENDELSE

  chisqRes = TOTAL((sampleHST - fitResult)^2) + RSred
  PRINT, paramVec
  PRINT, 'chisq:', chisqRes, RSred
  RETURN, chisqRes
END


FUNCTION AllOpt, FIXEDRV=fixedRv, FIXEDZP=zplist, FIXED275=shift275

;
;****************************************************************************
;* samples from multiD gaussian distribution and calls LSST_dump_func
;
;***************************************************************************
  COMMON TopInfo,  nStars, nBp, nBpHST, bpData, bandList, bandDict, zp, Tstars, TMinStars, TMaxStars, Gstars, GMinStars, GMaxStars, zpStar, EBV, modelWl, modelFluxes, sampleHST
  COMMON ZpMinimizeInfo, fitResult, chisqRes, scaleFactor, idRS, fluxRS, obsMagRS
  COMMON ScaleInfo, scaleFactor2, scaleFactor2Inv, scaleFactor3, scaleFactor3Inv, scaleFactor275, scaleFactor275Inv

  
  nParams = nStars+nBpHST+nStars+nStars+nStars+1
  params = DBLARR(nParams)      ; Av for each star, zp offset for each band (sum to zero), R, zp for each star, temp for each star, logg for each star
  paramsUB = REPLICATE(1.0e30, nParams)
  paramsLB = REPLICATE(-1.0e30, nParams)

  scaleFactor = 1.0e3           ; to fool CONSTRAINED_MIN into doing the right stepsize
  scaleFactorInv = 1./scaleFactor

  scaleFactor2 = 1.0e6           
  scaleFactor2Inv = 1./scaleFactor2

  scaleFactor3 = 1.0e2
  scaleFactor3Inv = 1./scaleFactor3

  scaleFactor275 = 1.e7
  scaleFactor275Inv = 1./scaleFactor275
;
; Fit for zeropoint, ebv and (maybe) R
;
     params(0:(nBpHST-2)) = scaleFactorInv*zp(0:(nBpHST-2)) ; delta zp for each band
     params(nBpHST:(nBpHST+nStars-1)) = scaleFactorInv*EBV ; scaleFactor*EBV for each star
     params((nBpHST+nStars):(nBpHST+2*nStars-1)) = 0
     params((nBpHST+nStars):(nBpHST+2*nStars-1)) = scaleFactorInv*zpStar
     params((nBpHST+2*nStars):(nBpHST+3*nStars-1)) = scaleFactor2Inv*Tstars  ; temp for each star
     params((nBpHST+3*nStars):(nBpHST+4*nStars-1)) = scaleFactor3Inv*Gstars  ; logg for each star
     params(nBpHST+4*nStars) = 0  ; shift for F275W
  IF KEYWORD_SET(zplist) THEN BEGIN
     paramsLB(0:(nBpHST-2)) = scaleFactorInv*zp(0:(nBpHST-2))
     paramsUB(0:(nBpHST-2)) = scaleFactorInv*zp(0:(nBpHST-2))
  ENDIF
     paramsLB((nBpHST+2*nStars):(nBpHST+3*nStars-1)) = scaleFactor2Inv*TMinStars 
     paramsUB((nBpHST+2*nStars):(nBpHST+3*nStars-1)) = scaleFactor2Inv*TMaxStars  

     paramsLB((nBpHST+3*nStars):(nBpHST+4*nStars-1)) = scaleFactor3Inv*GMinStars 
     paramsUB((nBpHST+3*nStars):(nBpHST+4*nStars-1)) = scaleFactor3Inv*GMaxStars  

     IF KEYWORD_SET(shift275) THEN BEGIN
        params(nBpHST+4*nStars) = scaleFactor275Inv*shift275
        paramsUB(nBpHST+4*nStars) = scaleFactor275Inv*shift275
        paramsLB(nBpHST+4*nStars) = -scaleFactor275Inv*shift275
     ENDIF ELSE BEGIN
        paramsUB(nBpHST+4*nStars) = scaleFactor275Inv*10 ; max ~100A shift in F275
        paramsLB(nBpHST+4*nStars) = -scaleFactor275Inv*10 ; max ~100A shift in F275
     ENDELSE

     IF KEYWORD_SET(FixedRv) THEN BEGIN
        params(nBPHST-1) = scaleFactorInv*FixedRv
        paramsLB(nBPHST-1) = scaleFactorInv*FixedRv
        paramsUB(nBPHST-1) = scaleFactorInv*FixedRv
     ENDIF ELSE BEGIN
        params(nBPHST-1)= scaleFactorInv*3.1
        paramsLB(nBPHST-1) = scaleFactorInv*2.1
        paramsUB(nBPHST-1) = scaleFactorInv*4.0
     ENDELSE

     paramsLB(nBpHST:(nBpHST+nStars-1)) = 0.0               ; extinction always >0

     status = -1
     
     paramsB = [paramsLB, paramsUB]
     CONSTRAINED_MIN, params, paramsB, [[0], [0]], 0, 'ExtinctFuncAll',status, REPORT='Allminreport.txt', EPSTOP=1.e-4
;     CONSTRAINED_MIN, params, paramsB, [[0], [0]], 0, 'ExtinctFuncAll',status, REPORT='Allminreport.txt', EPSTOP=10.
     ;; PRINTF, 3,'-------------------------'
     ;; PRINTF, 3, 'status: ', status
     ;; PRINTF, 3, params
     ;; PRINTF, 3, sampleHST - fitResult
     ;; PRINTF, 3, chisqRes
     ;; PRINTF, 3, '-------------------------'

     zp[0:(nBpHST-2)] = scaleFactor*params[0:(nBpHST-2)]
     zp[nBpHST-1] = -TOTAL(zp[0:(nBpHST-2)]) ; zpBand sums to zero

     zpStar = scaleFactor*params((nBpHST+nStars):(nBpHST+2*nStars-1))
     
     EBV = scalefactor*params[nBPHST:(nBpHST+nStars-1)]
     R = scaleFactor*params[nBpHST-1]

     Tstars = scalefactor2*params((nBpHST+2*nStars):(nBpHST+3*nStars-1))
     Gstars = scalefactor3*params((nBpHST+3*nStars):(nBpHST+4*nStars-1))

     FOR k = 1, nStars DO BEGIN
        idx = (k-1)*nBpHST
        resid = fitresult[idx:(idx+nBpHst-1)] - sampleHST[idx:(idx+nBpHst-1)]
;        PRINTF, 3, 'resid:', k, resid
     ENDFOR

     ;; PRINTF, 3, 'zp:', zp
     ;; PRINTF, 3, 'Av:', R*EBV
     ;; PRINTF, 3, 'temps:', Tstars
     ;; PRINTF, 3, 'zpStars:', zpStar



  RETURN, params
END
