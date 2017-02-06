;
; This is the top level sequencer for the optimization
;
; Note, if zplist is given, the reddening constraint star is ignored.
; Also, starList should have only a single star
;

PRO TopOpt,starFile,starList,nSample,outFileName,FIXEDRV=fixedRv,FIXEDZP=zplist,PLOTSTARS=doPlots,DUMPFLUX=fluxFileName

  COMMON TopInfo,  nStars, nBp, nBpHST, bpData, bandList, bandDict, zp, Tstars, TMinStars, TMaxStars, Gstars, GMinStars, GMaxStars, zpStar, EBV, modelWl, modelFluxes, sampleHST
  COMMON ZpMinimizeInfo, fitResult, chisqRes, scaleFactor, idRS, fluxRS, obsMagRS

  IF KEYWORD_SET(zplist) THEN BEGIN
     RSStar = ''
     nStars = N_ELEMENTS(starList)
     IF nStars NE 1 THEN BEGIN
        PRINT, 'When ZP is given, starList can have only one star'
        RETURN
     ENDIF
     
  ENDIF ELSE BEGIN
     RSStar = 'G191B2B'
  ENDELSE


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
  EBV0 = DBLARR(nStars)
  TStars = DBLARR(nStars)
  TNominal = DBLARR(nStars)
  TMinStars = DBLARR(nStars)
  TMaxStars = DBLARR(nStars)
  Torig = DBLARR(nStars)
  GStars = DBLARR(nStars)
  GNominal = DBLARR(nStars)
  GMinStars = DBLARR(nStars)
  GMaxStars = DBLARR(nStars)
  Gorig = DBLARR(nStars)
  zpStar = DBLARR(nStars)
  
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
     IF starList[n-1] eq RSStar THEN BEGIN
        idRS = n
        PRINTF, 2, 'IDRS=', idRS
     ENDIF

     TNominal[n-1] = Tret
     TMinStars[n-1] = Tret - 2*sigmaT
     TMaxStars[n-1] = Tret + 2*sigmaT
     Torig[n-1] = Tret
     GNominal[n-1] = Gret
     GMinStars[n-1] = Gret - 2*sigmaG
     GMaxStars[n-1] = Gret + 2*sigmaG
     Gorig[n-1] = Gret
     EBV0[n-1] = EBVret
     idx = (n-1)*nBpHST
;     HSTObsAll[idx:(idx+nBpHST-1)] = HSTObs - zeropointsHST; these
;                                              are now instrumental
;                                              mags- OLD convention for
;                                              input mags
     HSTObsAll[idx:(idx+nBpHST-1)] = HSTObs; these are instrumental mags as input
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
; read in the HST bandpasses. Note files names must be like 'F336W.csv'
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

  format1 = '($,A20,7(e12.4))'
  format2 = '($,e14.6)'
  
  bandMag = DBLARR(nBP+nBpHST)
  FOR n = 1, nSample DO BEGIN
     FOR k = 1, nStars DO BEGIN
        idx = (k-1)*nBpHST
        sampleHST(idx:(idx+nBpHST-1)) = mrandomn(undef, covarHST[*,*,k-1]) + HSTObsAll(idx:(idx+nBpHST-1))
        randomTG = mrandomn(undef, covarTG[*,*,k-1])
        TStars = TNominal + randomTG(0)
        GStars = GNominal + randomTG(1)
     ENDFOR
     
     Torig = Tstars
     Gorig = Gstars
     
;
; initial guesses for ZpStar and EBV
;
     ZpStars = DBLARR(nStars)
     EBV = EBV0

     
  IF KEYWORD_SET(zplist) THEN BEGIN
     zp = zplist
  ENDIF ELSE BEGIN
     zp = DBLARR(nBpHST)
  ENDELSE

;
; Use RSStar as a reddening constraint
;
  IF idRS NE !NULL THEN BEGIN
     idx = (idRS-1)*nBpHST
     obsMagRS = HSTObsAll[idx:(idx+nBpHST-1)]
     tempRS = Tstars[idRS-1]
     loggRS = Gstars[idRS-1]
     LSST_dump_func, tempRS, loggRS, 0, 0, 1.0, 0, modelWl,  fluxRS
  ENDIF

;
; OPTIMIZER CALL    
;
; ------------------------------------------------------------ 
     paramsAll = AllOpt(FIXEDRV=fixedRv,FIXEDZP=zplist)
; ------------------------------------------------------------ 
;
;

;
; Dump zpBand
;
     zpBand(0:(nBpHST-2)) = scaleFactor*paramsAll(0:(nBpHST-2))
     zpBand(nBpHST-1) = -TOTAL(zpBand(0:(nBpHST-2))) ; zpBand sums to zero

     PRINTF, 2, bandListStr
     PRINTF, 2, FORMAT='("zpBand: ",10e12.4)', zpBand
;
; DUMP F275w shift
;
     bp275wShift = scalefactor*paramsAll(nBpHST+4*nStars)
     PRINTF, 2, FORMAT='("bp275wShift: ",f0.4)', bp275wShift
;
; Dump RS mags
;
  IF idRS NE !NULL THEN BEGIN
     bp275w = bpData[nBp+bandDict['F275W']] 
     synMagRSf275w = synphot2(modelWl, fluxRS, bp275w.wavelength, SHIFT(bp275w.throughput, bp275wShift), 0) + zpBand[bandDict['F275W']]
     obsMagRSf275w = obsMagRS[bandDict['F275W']]

     bp336w = bpData[nBp+bandDict['F336W']] 
     synMagRSf336w = synphot2(modelWl, fluxRS, bp336w.wavelength, bp336w.throughput, 0) + zpBand[bandDict['F336W']]
     obsMagRSf336w = obsMagRS[bandDict['F336W']]

     bp475w = bpData[nBp+bandDict['F475W']]
     synMagRSf475w = synphot2(modelWl, fluxRS, bp475w.wavelength, bp475w.throughput, 0) + zpBand[bandDict['F475W']]
     obsMagRSf475w = obsMagRS[bandDict['F475W']]

     bp625w = bpData[nBp+bandDict['F625W']]
     synMagRSf625w = synphot2(modelWl, fluxRS, bp625w.wavelength, bp625w.throughput, 0) + zpBand[bandDict['F625W']]
     obsMagRSf625w = obsMagRS[bandDict['F625W']]

     bp775w = bpData[nBp+bandDict['F775W']] 
     synMagRSf775w = synphot2(modelWl, fluxRS, bp775w.wavelength, bp775w.throughput, 0) + zpBand[bandDict['F775W']]
     obsMagRSf775w = obsMagRS[bandDict['F775W']]

     bp160w = bpData[nBp+bandDict['F160W']] 
     synMagRSf160w = synphot2(modelWl, fluxRS, bp160w.wavelength, bp160w.throughput, 0) + zpBand[bandDict['F160W']]
     obsMagRSf160w = obsMagRS[bandDict['F160W']]

     PRINTF, 2 , 'RS syn:', synMagRSf275w, synMagRSf336w, synMagRSf475w, synMagRSf625w, synMagRSf775w, synMagRSf160w
     PRINTF, 2, 'RS obs:', obsMagRSf275w, obsMagRSf336w, obsMagRSf475w, obsMagRSf625w, obsMagRSf775w, obsMagRSf160w
     PRINTF, 2, 'RS syn-obs:', synMagRSf275w - obsMagRSf275w, synMagRSf336w - obsMagRSf336w, synMagRSf475w - obsMagRSf475w, synMagRSf625w - obsMagRSf625w, synMagRSf775w - obsMagRSf775w, synMagRSf160w - obsMagRSf160w
  ENDIF

     formatHead = '(A0," ",A0," ",A0," ",A0)'
     PRINTF,2,FORMAT=formatHead, '# id T Torig logg loggorig Av sigma u g r i z (observed)', bandListStr,  '(fit-observed)', bandListStr
;
; Calculate extincted flux
;
     R = scaleFactor*paramsAll(nBPHST-1)
     PRINTF,2,'Rv: ', R

     FOR k = 1, nStars DO BEGIN
        LSST_dump_func, Tstars[k-1], GStars[k-1], 0, EBV[k-1], 1.0, 0, modelWl, modelFluxes[k-1]
        extincMag = ext_odonnell(modelWl,R)*EBV[k-1]*R
        extinc = 10^(-0.4*extincMag)
        fluxe = modelFluxes[k-1]*extinc
        idx = (k-1)*nBpHST
        sigma = STDDEV(fitResult(idx:(idx+nBpHST-1)) - sampleHST(idx:(idx+nBpHST-1)))
        PRINTF, 2, FORMAT=format1, starList[k-1], Tstars[k-1], Torig[k-1], Gstars[k-1], Gorig[k-1], R*EBV[k-1], sigma
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
