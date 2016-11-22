# WDFitPhot

## Synopsis

Optimal fitting of HST photometry to WD models.  The star independent
band offsets and per-star Av are calculated to minimize the RMS
between model photometry and the observed HST photometry.

## Code Example

The top level routine is TopOpt in TopOpt.pro.  The calling sequence is

TopOpt,starFile,starList,nSample,outFileName,FIXEDRV=fixedRv,PLOTSTARS=doPlots,DUMPFLUX=fluxFileName

The calling parameter meanings are:

    starFile:	      Path to the data file for the observed stars.  An example is 'RunData/WD_cycle22_20_mean.dat'
    starList:	      IDL list of names of stars to process.  The names should all be present in starFile
    nSample:	      Number of samples to draw from a normal distribution of HST fluxes, temps and logg values. 
    		      The sigmas for the gaussians are specified in starFile.  Covariance is assumed diagonal for
		      now.
    outFileName:      Name of the output file
    fixedRv:	      If Rv is to be fixed (recommended for now), this is its value					 doPlots:          Set to some non-null value to produce plots for each star (broken for the moment)
    fluxFileName:     If non-null, dump model fluxes for each star to this file

An example call to TopOpt:
   TopOpt,'RunData/WD_cycle22_20_mean.dat',Starlist8,1,'../Nov16Runs/8Star_1.out',FIXEDRV=3.1