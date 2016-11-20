;*****************************************************************************
;+
;*NAME:
;
;    TRIM     JULY 13,1981
;  
;*CLASS:
;
;    Spectral Data Reduction
;  
;*CATEGORY:
;  
;*PURPOSE:
;
;    Use the values of the wave vector to eliminate all points not within the
;    specified range (LOW to HIGH) of the first vector from the other supplied
;    vectors.
;  
;*CALLING SEQUENCE:
;
;    TRIM,LOW,HIGH,V1,V2,v3,v4,header=header
;  
;*PARAMETERS:
;
;    LOW	(REQ) (I) (0) (I L R D)
;		Lowest valid value of first supplied vector, V1.
;
;    HIGH	(REQ) (I) (0) (I L R D)
;		Highest valid value of first supplied vector, V1.
;
;    V1		(REQ) (I/O) (1) (I L R D)
;		Fisrt vector.  Values for LOW and HIGH based on.  V1 will be
;		trimmed accorfing to the LOW and HIGH values.
;
;    V2		(REQ) (I/O) (1) (I L R D)
;		Second vector to be trimmed.
;
;    V3...V6	(OPT) (I/O) (1) (I L R D)
;		Additional vectors to be trimmed.
;
;    HEADER	(KEY) (I/O) (1) (S)
;		The fits header.
;  
;*EXAMPLES:
;
;    To plot only the IUESIPS fluxes from 1000 to 2000 angstroms for a low
;    dispersion spectrum:
;
;         iuespec,imaget,h,wave,flux,eps
;         trim,1000,2000,wave,flux,eps
;         iueplot,h,wave,flux,eps
;   
;    To eliminate points with EPS values < -400:
;
;         trim,-400,max(eps),eps,wave,flux
;         iueplot,h,wave,flux,eps
;
;    For NEWSIPS,
;
;	  readmx,filename,main,wave,flux,flags,sigma,bkgrd,net
;	  trim,1000,2000,wave,flux,flags,sigma,bkgrd,net,header=main
;	  nsplot,main,wave,flux,flags,sigma
;  
;*SYSTEM VARIABLES USED:
;
;    none
;
;*INTERACTIVE INPUT:
;
;    none
;
;*SUBROUTINES CALLED:
;
;    PARCHECK
;    ADDPAR
;  
;*FILES USED:
;
;    none
;
;*SIDE EFFECTS:
;  
;*RESTRICTIONS:
;  
;*NOTES:
;
;    The original vectors V1, V2, V3, V4, V5, and V6 are changed to include
;    only the points between LOW and HIGH.
;
;    tested with IDL Version 2.1.2  (sunos sparc)     23 Jul 91
;    tested with IDL Version 2.1.0  (ultrix mipsel)   23 Jul 91
;    tested with IDL Version 2.1.2  (vms vax)         23 Jul 91
;  
;*PROCEDURE:
;
;    All points in the V1 vector less than LOW and greater than HIGH are
;    eliminated and the V1, V2, V3, V4, V5, and V6 vectors are reordered to
;    include only the remaining points.  If the HEADER keyword is set, a
;    HISTORY line is added about the use of trim.
;  
;*I_HELP nn:
;  
;*MODIFICATION HISTORY:
;
;    F.H. SCHIFFER 3RD  VERSION 0  13-JULY-1981
;    4-30-87 RWT add PARCHECK
;    7-19-91 GRA Converted code to lowercase; cleaned up; tested on
;                SUN, DEC, VAX; updated prolog
;	19 Nov 91  PJL  corrected typos in prolog
;	25 Jun 93  PJL  generalized;  made third vector optional;  added 3 more
;			optional vectors;  added fits header keyword
;	22 Sep 93  PJL  added else to npar case statement
;       28 Sep 93  LLT  fix minor typo in prolog
;-
;******************************************************************************
 pro trim,low,high,v1,v2,v3,v4,v5,v6,header=header
;
 npar = n_params(0)
 if (npar eq 0) then begin
    print,'TRIM,LOW,HIGH,V1,V2,v3,v4,v5,v6,header=header'
    retall 
 endif  ; npar
 parcheck,npar,[4,5,6,7,8],'TRIM'
;
 i = where(((v1 ge low) and (v1 le high)),ct)
;
;  no data remain
;
 if (ct le 0) then begin
    print,'No data available between ' + strtrim(low,2) + ' and ' +   $
       strtrim(high,2) + '.'
    print,'No points removed.'
    return
 endif  ; ct le 0
;
 v1 = v1(i)
 v2 = v2(i)
;
;  optional vectors
;
 case npar of
    5:  v3 = v3(i)
    6:  begin
           v3 = v3(i)
           v4 = v4(i)
        end  ; npar eq 6
    7:  begin
           v3 = v3(i)
           v4 = v4(i)
           v5 = v5(i)
        end  ; npar eq 7
    8:  begin
           v3 = v3(i)
           v4 = v4(i)
           v5 = v5(i)
           v6 = v6(i)
        end  ; npar eq 8
    else:
 endcase  ; npar
;
;  fits header
;
 if keyword_set(header) then addpar,header,'HISTORY',   $
    'TRIM used with LOW = ' + strtrim(low,2) + ' and HIGH = ' +    $
    strtrim(high,2) + '   ' + !stime,'','IUEDAC'
;
 return 
 end  ; trim
