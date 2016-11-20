 FUNCTION FITZPATRICK,LAM,Rv
;******************************************************************************
; FITZPATRICK returns the normalized extinction A(LAM)/E(B-V) as function
; of wavelength (LAM) and the Rv parameter.  Constructed from the tables
; Fitzpatrick (2004).  Extinction Curves valid from 915 A to 1.6x10^6 A and
; for 2.1 < Rv < 5.5.
; 
; JBH 03/07/13
;******************************************************************************
; Requires    FITZPATRICK_3.1.sav   Extinction curve for Rv=3.1
;             FITZPATRICK_TOT.sav   All Extinction Curves for Rv
; Ref:  Fitzpatrick, E. L. 2004, in ASP Conf. Ser. 309, Astrophysics
; of Dust, ed. A. N. Witt, G. C. Clayton, & B. T. Draine (San
; Francisco: ASP),33. 
;******************************************************************************
; Note:  Extinction Factor is EF = 10^0.4*A where A = FITZPATRICK* E(B-V) 
; To apply Extinction Flux_red   = Flux/EF
; To deredden         Flux_dered = Flux*EF
;******************************************************************************
; Example compute the Absorption at 7000 A for Rv=3.1
; A = fitzpatrick(7000,3.1)
;******************************************************************************

; LAM    ........................... Wavelength (A), scalar or vector
; Rv     ........................... R parameter default R = 3.1




RESTORE,'/home/tsa/Dropbox/WD/GauthamCode/DATA/FITZPATRICK_TOT'            ; R,W,A (1099 points)

RR = 10*(Rv-2.1)
R2 = CEIL(RR)
R1 = FLOOR(RR)

MAXROW = (SIZE(A))(2)

IF R1 LT 0 THEN BEGIN
   R2 = 1
   R1 = 0
ENDIF

IF R2 EQ R1 THEN BEGIN
   R2 = R2+1
ENDIF

IF R2 GE MAXROW THEN BEGIN
   R2 = MAXROW - 1
   R1 = R2 - 1
ENDIF
   
A1 = A(*,R1)
A2 = A(*,R2)

FAV = A1 + (A2 - A1)*(RR - R1)/(R2 - R1)

F = INTERPOL(FAV,W,LAM)

FITZPATRICK = (Rv + F)

RETURN, FITZPATRICK

END





