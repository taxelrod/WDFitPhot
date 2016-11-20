function trapint,x,f

; Integrates f(x)dx with the trapezoidal rule
;
; x has to be monotonically increasing
;
; trapezoidal integration is the way to go on noisy data like spectra,
; it is also flux conserving.
;

n = n_elements(x)
if (n_elements(f) ne n) then begin
    message,'array sizes differ!'
endif

if n eq 1 then begin
    message,'warning 1-element arrays! returning 0',/info
    return,f[0]-f[0]        ; return same type as f
endif

diff = x[1:n-1]-x[0:n-2]
if min(diff) le 0 then begin
    message,'x must be montonically increasing'
endif

;
; the lines below are equivalent to:
;   res = f[0] - f[0]
;   for i=0,n-2 do res = res + (x[i+1]-x[i])*(f[i]+f[i+1])/2.0
;   return, res
;
; but using idl array manipulation to make sum and diff and multiply
; and total is much much faster. Integrating a 2 million element
; double array takes over 30s on mondatta with a loop, only about
; 1 or 2s this way.
;
; it's also faster than implementing the extended trapezoidal rule with
; a loop

sum = f[0:n-2]+f[1:n-1]
return, total(sum*diff)/2.0

end
