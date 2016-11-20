function synphot2,x,spc,pbx,pb,zp,plot=doplot,oplot=dooplot,pbsplinit=dopbsplinit $
                 ,allowneg=doallowneg, natspline=donatspline, pbergs=dopbergs, normit=donormit

;
;  x:   wavelength points for spectrum  
;  spc: spectrum (f_lambda)
;
;  pbx: wavelength points for passband (same units as x, e.g. Angstroms)
;  pb:  passband (energy or photon sensitivity/wavelength) 
;  zp: passband zero point
;
;  pbergs: set /pbergs if the passband response is energy sensitivity
;
;  allowneg: allow the passband to be negative (it's nonnegative by default)
;  natspline: allow first derivative at passband edge to be nonzero
;             (use only if there's a sharp break at the edge)
;
;  plot: set to see a plot
;  oplot: just overplot the passband (good for repeated calls)
;
;  pbsplinit: set this to some variable = spl_init(pbx,pb) if you are 
;             doing the same passband over and over to save time
;
;  returns mag
;
return, zp - 2.5 * alog10( synflux2(x,spc,pbx,pb,plot=plot,oplot=oplot $
                                   ,pbsplinit=pbsplinit,allowneg=allowneg $
                                   ,natspline=natspline, pbergs=pbergs, normit=donormit))

end
