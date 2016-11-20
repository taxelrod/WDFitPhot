function synflux2,x,spc,pbx,pb,zp,plot=doplot,oplot=dooplot,pbsplinit=dopbsplinit $
                 ,allowneg=doallowneg, natspline=donatspline, pbergs=dopbergs, normit=donormit

;
;  x:   wavelength points for spectrum  
;  spc: spectrum (f_lambda)
;
;  pbx: wavelength points for passband (same units as x, e.g. Angstroms)
;  pb:  passband (energy or photon sensitivity/wavelength)
;  zp:  zeropoint of passband [not used in synflux]
;
;  pbergs: set /pbergs if the passband response is energy sensitivity
;
;  allowneg: allow the passband to be negative (it's nonnegative by default)
;  natspline: allow first derivative at passband edge to be nonzero
;             (use only if there's a sharp break at the edge)
;  
;
;  plot: plot results
;  oplot: overplot just the passband (good for repeated calls)
;
;  pbsplinit: set this to some variable = spl_init(pbx,pb) if you are 
;             doing the same passband over and over to save time
;
;  we cubic spline interpolate the passband onto the spectrum wavelengths
;  maybe we will add some sort of extrapolation procedure later
;
;  x and pbx should be monotonically increasing
;
;  returns flux
;

pbphot = not keyword_set(pbergs)

nx = n_elements(x)
if (n_elements(spc) ne nx) then begin
    print,'synflux: x and spc have different sizes'
    stop
endif

np = n_elements(pbx)
if (n_elements(pb) ne np) then begin
    print,'synflux: pbx and pb have different sizes'
    stop
endif

if nx eq 1 or np eq 1 then begin
    print,'synflux: warning! 1-element array passed, returning 0'
    return,spc[0]-spc[0]
endif

diffx = x[1:nx-1]-x[0:nx-2]
diffp = pbx[1:np-1]-pbx[0:np-2]

if (min(diffx) le 0) or (min(diffp) le 0) then begin
    print,'synflux: passed non-increasing wavelength array'
    stop
endif

if x[0] gt pbx[0] then begin
    print,'synflux: spectrum doesn''t go blue enough for passband!'
    stop
endif

if x[nx-1] lt pbx[np-1] then begin
    print,'synflux: spectrum doesn''t go red enough for passband!'
    stop
endif

g = where(x ge pbx[0] and x le pbx[np-1])  ; overlap range

; set up spline if not already set up
if (keyword_set(pbsplinit) and (n_elements(pbsplinit) eq np)) then begin
    spl = pbsplinit
endif else if keyword_set(natspline) then begin
    spl = spl_init(pbx,pb)
endif else begin
    spl = spl_init(pbx,pb,yp0=0.0,ypn_1=0.0)
endelse
pbspl = spl_interp(pbx,pb,spl,x[g])  ; do it

if not keyword_set(allowneg) then pbspl = pbspl > 0

if keyword_set(plot) then plot,x,spc
if keyword_set(oplot) or keyword_set(plot) then $
  oplot,x[g],!y.crange[0]+pbspl/max(pbspl)*(!y.crange[1]-!y.crange[0])*0.9d,linestyle=2


; if the response is photon sensitivity we need to integrate
;    f_lam*response*lambda*dlambda
; otherwise we integrate f_lam*response*dlambda
if (pbphot) then pbspl *= x[g]

if keyword_set(donormit) then begin
   res = trapint(x[g],pbspl*spc[g])/trapint(x[g],pbspl)
endif else begin
   res = trapint(x[g],pbspl*spc[g])
endelse

return, res

end
