;+
; NAME:                ULY_TGM_EXTR
;
; PURPOSE:             Compute a stellar template spectrum using TGM model
;                      
; USAGE:               
;               (wl, flux) = uly_tgm_extr([<spectrum>], <pars>, MODEL_FILE=<file>)
;
; DESCRIPTION:
;      Compute a 1D stellar template spectrum by using TGM model at a given
;      Teff, Log g and [Fe/H]. Rebin it with the same sampling as <spectrum>
;      and with a compatible grid.
;
;      If the natural wavelength range of the model is wider than the one
;      of <spectrum> the returned spectrum has exactly the same sampling as 
;      <spectrum> (same range and same step). If instead the wavelength
;      range of the model is narrower, the output spectrum will be truncated
;      with respect to <spectrum> but its wavelength bins will match bins
;      of <spectrum>.
;
;      Note that although the model may be rebinned, its resolution is not 
;      matched to that of <spectrum>.
;
;      If <spectrum> is omitted, i. e. if only one argument is given, the
;      returned spectrum has the original sampling of the model.
;       
; ARGUMENT:
;      <spectrum>    filename or spectrum structure serving as a reference 
;                    for the velscale.
;
;      <pars>        array of parameters([Teff,Logg,[Fe/H]) used
;                    to evaluate the template spectrum.
;
; KEYWORD:              
;      MODEL_FILE = <file>
;         Filename of the model to use, by default:
;         uly_root + '/models/elodie32_flux_tgm.fits' 
; 
; RETURNS:    
;      A list (wl, flux)
;            
; EXAMPLE:
;      To compute a template spectrum at a given atmosphere parameters
;      according to uly_root+'/data/cflib_114642.fits's velscale.
;            model = uly_tgm_extr(uly_root+'/data/cflib_114642.fits', $
;                                 [6400., 4., 0.])
; 
; HISTORY:
;            Yue WU 2008/08/20 
;-
; CATEGORY:  ULY_TGM
;--------------------------------------------------------------------------
function uly_tgm_extr, spectrum, pars, MODEL_FILE=model_file, REBIN_COEF=rebin_coef

compile_opt idl2
on_error, 2

if n_params() eq 1 then begin ; a priori gave only <pars>
   pars = spectrum
   undefine, spectrum
endif

if file_test(model_file) ne 1 then begin
    message, /INFO, 'Error, model file does not exist'
    return, 0
endif

if n_elements(spectrum) eq 0 then begin
   spectrum = model_file
endif

if size(spectrum, /TYPE) ne 8 then begin
    spec = uly_spect_read(spectrum, /quiet)
endif else begin
    spec = uly_spect_alloc(SPECTRUM=spectrum)
endelse
if not uly_spect_get(spec, /VALID) then begin
  message, 'Input spectrum is invalid', /CONT
endif

if n_elements(pars) ne 3 then $
  message, '"Pars" should be a 3 elements array [Teff, Logg, [Fe/H]]'
if min(finite(pars[0])) eq 0 then message, 'Require a finite Teff input value'
if min(finite(pars[1])) eq 0 then message, 'Require a finite Logg input value'
if min(finite(pars[2])) eq 0 then message, 'Require a finite [Fe/H] input value'

pars[0]  = alog(pars[0])   ;Teff in log

cmp = uly_tgm(MODEL_FILE=model_file, REBIN_COEF=rebin_coef)
waverange = uly_spect_get(spec, /WAVERANGE)

cmp = uly_tgm_init(cmp, WAVERANGE=waverange, SAMPLING=spec.sampling, STEP=spec.step)

model = uly_tgm_eval(*cmp.eval_data, pars)

uly_spect_free, spec
heap_free, cmp

nWl = n_elements(model)
wl = cmp.start + cmp.step*findgen(nWl)
if (cmp.sampling eq 1) then begin
   wl = exp(wl)
endif
return, list(wl, model)
end
