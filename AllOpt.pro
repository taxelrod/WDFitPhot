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
  COMMON TopInfo,  nStars, nBp, nBpHST, bpData, bandList, bandDict, zp, Tstars, Gstars, zpStar, EBV, modelWl, modelFluxes, sampleHST

  COMMON ZpMinimizeInfo, fitResult, chisqRes, scaleFactor, fluxGD153, obsMagGD153f336w, obsMagGD153f475w, obsMagGD153f625w

  COMMON ScaleInfo, scaleFactor2, scaleFactor2Inv


  zpBand = DBLARR(nBpHST)
  zpBand(0:(nBpHST-2)) = scaleFactor*paramVec(0:(nBpHST-2))
  zpBand(nBpHST-1) = -TOTAL(zpBand(0:(nBpHST-2))) ; zpBand sums to zero
  R = paramVec(nBpHST-1)
  EBV = scaleFactor*paramVec(nBpHST:(nBpHST+nStars-1))
  zpStar = scaleFactor*paramVec((nBpHST+nStars):(nBpHST+2*nStars-1))
  Tstars = scaleFactor2*paramVec((nBpHST+2*nStars):(nBpHST+3*nStars-1))

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
        synMag(idx) = synphot2(modelWl, modelFluxes[n-1]*extinc, bp.wavelength, bp.throughput, 0)
        fitResult(idx) = zpBand[i-1] + zpStar[n-1] + synMag(idx)

     ENDFOR
  ENDFOR

; Stddev of GD153 mags as proxy for reddening
;
  bp336w = bpData[nBp+bandDict['F336W']]
  synMagGD153f336w = synphot2(modelWl, fluxGD153, bp336w.wavelength, bp336w.throughput, zpBand[bandDict['F336W']])

  bp475w = bpData[nBp+bandDict['F475W']]
  synMagGD153f475w = synphot2(modelWl, fluxGD153, bp475w.wavelength, bp475w.throughput, zpBand[bandDict['F475W']])

  bp625w = bpData[nBp+bandDict['F625W']]
  synMagGD153f625w = synphot2(modelWl, fluxGD153, bp625w.wavelength, bp625w.throughput, zpBand[bandDict['F625W']])

  GD153delta = DBLARR(3)
  GD153delta[0] = synMagGD153f336w - obsMagGD153f336w
  GD153delta[1] = synMagGD153f475w - obsMagGD153f475w
  GD153delta[2] = synMagGD153f625w - obsMagGD153f625w

  GD153red = TOTAL((GD153delta - MEAN(GD153delta))^2) / 0.01^2

  chisqRes = TOTAL((sampleHST - fitResult)^2) + GD153red
  PRINT, paramVec
  PRINT, 'chisq:', chisqRes, GD153red
  RETURN, chisqRes
END


FUNCTION AllOpt, FIXEDRV=fixedRv

;
;****************************************************************************
;* samples from multiD gaussian distribution and calls LSST_dump_func
;
;***************************************************************************
  COMMON TopInfo,  nStars, nBp, nBpHST, bpData, bandList, bandDict, zp, Tstars, Gstars, zpStar, EBV, modelWl, modelFluxes, sampleHST

  COMMON ZpMinimizeInfo, fitResult, chisqRes, scaleFactor, fluxGD153, obsMagGD153f336w, obsMagGD153f475w, obsMagGD153f625w

  COMMON ScaleInfo, scaleFactor2, scaleFactor2Inv

  
  nParams = nStars+nBpHST+nStars+nStars
  params = DBLARR(nParams)      ; Av for each star, zp offset for each band (sum to zero), R, zp for each star, temp for each star
  paramsUB = REPLICATE(1.0e30, nParams)
  paramsLB = REPLICATE(-1.0e30, nParams)

  scaleFactor = 1.0e3           ; to fool CONSTRAINED_MIN into doing the right stepsize
  scaleFactorInv = 1./scaleFactor

  scaleFactor2 = 1.0e6           
  scaleFactor2Inv = 1./scaleFactor2
;
; Fit for zeropoint, ebv and (maybe) R
;
     params(0:(nBpHST-2)) = scaleFactorInv*0 ; delta zp for each band
     params(nBpHST:(nBpHST+nStars-1)) = scaleFactorInv*EBV ; scaleFactor*EBV for each star
     params((nBpHST+nStars):(nBpHST+2*nStars-1)) = 0
     params((nBpHST+nStars):(nBpHST+2*nStars-1)) = scaleFactorInv*zpStar
     params((nBpHST+2*nStars):(nBpHST+3*nStars-1)) = scaleFactor2Inv*Tstars  ; temp for each star

     paramsLB((nBpHST+2*nStars):(nBpHST+3*nStars-1)) = scaleFactor2Inv*Tstars  ; Don't allow temp to vary
     paramsUB((nBpHST+2*nStars):(nBpHST+3*nStars-1)) = scaleFactor2Inv*Tstars  

     IF KEYWORD_SET(FixedRv) THEN BEGIN
        params(nBPHST-1) = FixedRv
        paramsLB(nBPHST-1) = FixedRv
        paramsUB(nBPHST-1) = FixedRv
     ENDIF ELSE BEGIN
        params(nBPHST-1)= 3.1
        paramsLB(nBPHST-1) = 2.1
        paramsUB(nBPHST-1) = 4.0
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
     R = params[nBpHST-1]

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
