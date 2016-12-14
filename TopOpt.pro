;
; This is the top level sequencer for the optimization
;


PRO TopOpt,starFile,starList,nSample,outFileName,FIXEDRV=fixedRv,PLOTSTARS=doPlots,DUMPFLUX=fluxFileName

  COMMON TopInfo,  nStars, nBp, nBpHST, bpData, bandList, bandDict, zp, Tstars, Gstars, ZpStars, EBV, modelWl, modelFluxes, sampleHST
  COMMON ZpMinimizeInfo, fitResult, chisqRes, scaleFactor, fluxGD153, obsMagGD153f275w, obsMagGD153f336w, obsMagGD153f475w, obsMagGD153f625w, obsMagGD153f775w, obsMagGD153f160w

  CLOSE,2
  OPENW,2,outFileName

  IF KEYWORD_SET(fluxFileName) THEN BEGIN
     CLOSE,3
     OPENW,3,fluxFileName
  ENDIF

 ; Call read_WDData once just to get the band list and the number
 ; of bands in use


  read_WDDATA,starFile,'',Tret,Gret,EBVret,HSTObs,HSTObsUnc,sigmaT,sigmaG,bandList

  nBpHST = N_ELEMENTS(bandList)
;
; dictionary to look up bands by name
;
  bandDict = DICTIONARY()
  bandListStr = ''
  FOR N = 0, nBpHST-1 DO BEGIN
     bandDict[bandList[N]] = N
     bandListStr += bandList[N] + ' '
  ENDFOR

  
; Get the central values of the T, G, EBV distributions, plus other
; needed per-star parameters

  nStars = N_ELEMENTS(starList)
  covarHST = DBLARR(nBpHST, nBpHST, nStars)
  covarTG = DBLARR(2, 2, nStars)
  HSTObsAll = DBLARR(nBpHST*nStars)
  SampleHST = DBLARR(nBpHST*nStars)
  T0 = DBLARR(nStars)
  G0 = DBLARR(nStars)
  EBV0 = DBLARR(nStars)
  T = DBLARR(nStars)
  G = DBLARR(nStars)
  
  zpBand = DBLARR(nBpHST)

; Set up the array of zeropoints in the order of the bands in the
; input photometry file.
; These are the infinite aperture zeropoints from MAST
; NOTE - order must be the same as in the file input to read_WDDATA!
;        Using bandList together with zeropointsHSTDict should ensure that


; Zeropoints on Cycle 22 Vega system
  zeropointsHSTDict = DICTIONARY('F275W', 22.70837, 'F336W', 23.53994, 'F475W', 25.53714,  'F625W', 25.37076,  'F775W', 24.45988, 'F160W', 24.70294)
;
; Zeropoints from latest MAST
;  zeropointsHSTDict = DICTIONARY('F275W', 22.548, 'F336W', 23.550, 'F475W', 25.309,  'F625W', 25.718,  'F775W', 25.485, 'F160W', 28.1875)
; Zeropoints as used for Cycle 20 paper
;  zeropointsHSTDict = DICTIONARY('F336W', 23.4836, 'F475W', 25.7783,  'F625W', 25.3783,  'F775W', 24.4747, 'F160W', 24.6949)

  zeropointsHST = DBLARR(nBpHST)
  FOR n = 0, nBpHST-1 DO BEGIN
     bandName = bandList[n]
     zeropointsHST[n] = zeropointsHSTDict[bandName]
  ENDFOR
;
; Now read in the per star info
;    
  FOR n = 1, nStars DO BEGIN
     read_WDDATA,starFile,starList[n-1],Tret,Gret,EBVret,HSTObs,HSTObsUnc,sigmaT,sigmaG,bandList
     T0[n-1] = Tret
     G0[n-1] = Gret
     EBV0[n-1] = EBVret
     idx = (n-1)*nBpHST
     HSTObsAll[idx:(idx+nBpHST-1)] = HSTObs - zeropointsHST; these are now instrumental mags
;
; Set up covariance matrices, diagonal for now
;
     covarTG[*, *, n-1] = DIAG_MATRIX([sigmaT^2, sigmaG^2])
     covarHST[*, *, n-1] = DIAG_MATRIX(HSTObsUnc^2)
  ENDFOR
;
; set up AB standard flux
;
  fluxAB = 3631.                ; Jy
  stdwl = 1000. + 10.*FINDGEN(2000)
  stdflux = fluxAB / (3.34e4 * stdwl^2)
;
; read in the SDSS bandpasses
;
  RESTORE, 'RunData/LSST_filter_template.sav'
  bandPassRoot = 'RunData/SDSS/'
  bandPasses = ['u','g','r','i','z']

  nBp = SIZE(bandPasses,/N_ELEMENTS)
  bpData = LIST()
  zeropoints = LIST()
  FOR n = 1, nBp DO BEGIN
     bandPassFile = bandPassRoot + bandPasses(n-1) + '.datx'
     PRINT, bandPassFile
     data=READ_ASCII(bandPassFile,TEMPLATE=tmpl)
     PRINT, size(data.wavelength,/N_ELEMENTS)
     bpData.add, data
                                ; calculate zeropoint as integral of standard flux over the bandpass
     zp = synphot2(stdwl, stdflux, data.wavelength, data.throughput, 0)
     zeropoints.add, zp
     print, 'check zero:', synphot2(stdwl, stdflux, data.wavelength, data.throughput, -zp)
  ENDFOR

  print, zeropoints
;
; set K0 from GD153 (see code CalcKFromGD153)
;
  K0 = -30.893

;
; read in the HST bandpasses. Note files must be like 'F336W.csv'
;
  RESTORE, 'RunData/HST_filter_template.sav'
  bandPassRoot = 'RunData/'
;

  FOR n = 1, nBpHST DO BEGIN
     bandPassFile = bandPassRoot + bandList(n-1) + '.csv'
     data=READ_ASCII(bandPassFile,TEMPLATE=tmplHST)
     bpData.add, data
     zeropoints.add, zeropointsHST(n-1)
  ENDFOR

  format1 = '($,A20,5(e12.4))'
  format2 = '($,e14.6)'

  
 
  
  bandMag = DBLARR(nBP+nBpHST)
  FOR n = 1, nSample DO BEGIN
     FOR k = 1, nStars DO BEGIN
        sample = mrandomn(undef, covarTG[*,*,k-1])
        ;; T[k-1] = sample(0) + T0(k-1)
        ;; G[k-1] = sample(1) + G0(k-1)
        T[k-1] = T0(k-1)
        G[k-1] = G0(k-1)
        idx = (k-1)*nBpHST
        sampleHST(idx:(idx+nBpHST-1)) = mrandomn(undef, covarHST[*,*,k-1]) + HSTObsAll(idx:(idx+nBpHST-1))
     ENDFOR
;
; initial guesses for ZpStar and EBV
;
     ZpStars = DBLARR(nStars)
     EBV = EBV0
     Tstars = T
     Gstars = G
     zp = DBLARR(nBpHST)
;
; Use GD153 as a reddening constraint; T and logg from Bohlin
;
     tempGD153 = 38686.
     loggGD153 = 7.66
     LSST_dump_func, tempGD153, loggGD153, 0, 0, 1.0, 0, modelWl,  fluxGD153

; cycle 20 (paper) system
     ;; obsMagGD153f336w = 11.348 - zeropointsHSTDict['F336W'] ; from GN table
     ;; obsMagGD153f475w = 13.182 - zeropointsHSTDict['F475W'] ; from GN table
     ;; obsMagGD153f625w = 13.455 - zeropointsHSTDict['F625W'] ; from GN table

; cycle 22 instrumental mags
     obsMagGD153f275w = -11.85745
     obsMagGD153f336w = -12.02242
     obsMagGD153f475w = -12.47516
     obsMagGD153f625w = -11.80634
     obsMagGD153f775w = -10.71429 
     obsMagGD153f160w = -10.39641
;
; OPTIMIZER CALL    
; 
     paramsAll = AllOpt(FIXEDRV=fixedRv)
;
;

;
; Dump zpBand
;
     zpBand(0:(nBpHST-2)) = scaleFactor*paramsAll(0:(nBpHST-2))
     zpBand(nBpHST-1) = -TOTAL(zpBand(0:(nBpHST-2))) ; zpBand sums to zero

     PRINTF, 2, bandListStr
     PRINTF, 2, FORMAT='("zpBand: ",10e12.4)', zpBand
     formatHead = '(A0," ",A0," ",A0," ",A0)'
     PRINTF,2,FORMAT=formatHead, '# id T Torig logg Av sigma u g r i z (observed)', bandListStr,  '(fit-observed)', bandListStr
;
; Dump GD153 mags
;
     bp275w = bpData[nBp+bandDict['F275W']] 
     synMagGD153f275w = synphot2(modelWl, fluxGD153, bp275w.wavelength, bp275w.throughput, 0) + zpBand[bandDict['F275W']]

     bp336w = bpData[nBp+bandDict['F336W']] 
     synMagGD153f336w = synphot2(modelWl, fluxGD153, bp336w.wavelength, bp336w.throughput, 0) + zpBand[bandDict['F336W']]

     bp475w = bpData[nBp+bandDict['F475W']]
     synMagGD153f475w = synphot2(modelWl, fluxGD153, bp475w.wavelength, bp475w.throughput, 0) + zpBand[bandDict['F475W']]

     bp625w = bpData[nBp+bandDict['F625W']]
     synMagGD153f625w = synphot2(modelWl, fluxGD153, bp625w.wavelength, bp625w.throughput, 0) + zpBand[bandDict['F625W']]

     bp775w = bpData[nBp+bandDict['F775W']] 
     synMagGD153f775w = synphot2(modelWl, fluxGD153, bp775w.wavelength, bp775w.throughput, 0) + zpBand[bandDict['F775W']]

     bp160w = bpData[nBp+bandDict['F160W']] 
     synMagGD153f160w = synphot2(modelWl, fluxGD153, bp160w.wavelength, bp160w.throughput, 0) + zpBand[bandDict['F160W']]

     PRINTF, 2 , 'GD153 syn:', synMagGD153f275w, synMagGD153f336w, synMagGD153f475w, synMagGD153f625w, synMagGD153f775w, synMagGD153f160w
     PRINTF, 2, 'GD153 obs:', obsMagGD153f275w, obsMagGD153f336w, obsMagGD153f475w, obsMagGD153f625w, obsMagGD153f775w, obsMagGD153f160w
     PRINTF, 2, 'GD153 syn-obs:', synMagGD153f275w - obsMagGD153f275w, synMagGD153f336w - obsMagGD153f336w, synMagGD153f475w - obsMagGD153f475w, synMagGD153f625w - obsMagGD153f625w, synMagGD153f775w - obsMagGD153f775w, synMagGD153f160w - obsMagGD153f160w


;
; Calculate extincted flux
;
     R = paramsAll(nBPHST-1)
     FOR k = 1, nStars DO BEGIN
        LSST_dump_func, Tstars[k-1], G0[k-1], 0, EBV[k-1], 1.0, 0, modelWl, modelFluxes[k-1]
        extincMag = ext_odonnell(modelWl,R)*EBV[k-1]*R
        extinc = 10^(-0.4*extincMag)
        fluxe = modelFluxes[k-1]*extinc
        idx = (k-1)*nBpHST
        sigma = STDDEV(fitResult(idx:(idx+nBpHST-1)) - sampleHST(idx:(idx+nBpHST-1)))
        PRINTF, 2, FORMAT=format1, starList[k-1], Tstars[k-1], T[k-1], G[k-1], R*EBV[k-1], sigma
        IF KEYWORD_SET(fluxFileName) THEN BEGIN
           fac = 10^(-0.4*(ZpStars[k-1] - K0))
           PRINTF, 3, '# ', starList[k-1]
           PRINTF, 3, '# wl flux'
           FOR l = 1, n_elements(modelWl) DO BEGIN
              PRINTF, 3, modelWl[l-1], fluxe[l-1]*fac
           ENDFOR
        ENDIF

;
; Calculate mags in all bands


        FOR j = 1, nBp DO BEGIN
           bp = bpData[j-1]
           bandMag[j-1] = synphot2(modelWl, fluxe, bp.wavelength, bp.throughput, -zeropoints(j-1)) + ZpStars[k-1] - K0
           PRINTF, 2, FORMAT=format2, bandMag[j-1]
        ENDFOR

        FOR j = 1, nBpHST DO BEGIN
           PRINTF, 2, FORMAT=format2, sampleHST[idx+j-1]
        ENDFOR

        FOR j = 1, nBpHST DO BEGIN
           PRINTF, 2, FORMAT=format2, fitResult[idx+j-1] - sampleHST[idx+j-1]
        ENDFOR

        ;;                         Calculate colors
        ;; FOR j = 1, nBp-1 DO BEGIN
        ;;    color = bandMag[j-1] - bandMag[j]
        ;;    PRINTF, 2, FORMAT=format2, color
        ;; ENDFOR
        PRINTF, 2, ' '
     ENDFOR

  ENDFOR


; plot result for each star
  IF KEYWORD_SET(doPlots) THEN BEGIN
     FOR k = 1, nStars DO BEGIN
        AV = 3.1*EBV[k-1]       ; sleaze alert
        idx = (k-1)*nBpHST
        wl =[336.,  475.,   625.,   775.,  1600.]
        plt=plot(1000./wl, fitResult(idx:(idx+nBpHST-1)),SYMBOL='+')
        plt=plot(1000./wl, sampleHST(idx:(idx+nBpHST-1)), SYMBOL='X',/OVERPLOT)
        plt.TITLE = starList(k-1) + ' SD =' + STRING(STDDEV(fitResult(idx:(idx+nBpHST-1)) - sampleHST(idx:(idx+nBpHST-1))), FORMAT='(F7.4)') + ' A_V=' + STRING(AV, FORMAT='(F7.4)')
        ax = plt.AXES
        ax[0].TITLE = '1/um'
        ax[1].TITLE = 'mag'

        plt2=plot(1000./wl, HSTObsAll(idx:(idx+nBpHST-1))-fitResult(idx:(idx+nBpHST-1)),SYMBOL='+')
        plt2.TITLE = starList(k-1) + ': residuals from fit'
        ax = plt2.AXES
        ax[0].TITLE = '1/um'
        ax[1].TITLE = 'mag'
     ENDFOR
  ENDIF

  CLOSE, 2
  CLOSE, 3
  RETURN
END
