;PROCEDURE MAGXY takes data from l.fname and ls.fname created by
;MAG and generates xy plots of energy time series and spectra.

;INPUT PARAMETERS
runname  = ' '		;MAG output suffix 
fname1   = ' '		;filename of timeseries data
fname2   = ' '      	;filename of spectra data
lmax     = 1 		;lmax value of grid, probably 32, 64, 96, etc
minc 	 = 1		;azimuthal folding symmetry, usually 1, 2, 4 or 6

PRINT,'  '
PRINT,'Enter l & ls files suffix name:' 
READ, runname
runname = STRTRIM(runname, 2)

fname1 = 'ls.'  + runname
fname2 = 'l.' + runname

PRINT,' '
PRINT,'Enter Harmonic truncation Lmax:' & read, lmax
PRINT,'Enter aximuthal symmetry minc:' 
READ, minc

;READ IN SPECTRAL DATA IN LS.'FNAME' FOR ns DIFFERENT NPRNT VALUES.
max_nprnt = 100

t_spectra = FLTARR(max_nprnt)			    	;time of spectra 
KE_L      = FLTARR(lmax + 1, max_nprnt)			;KE density as a function of l
ME_L      = FLTARR(lmax + 1, max_nprnt)			;ME density as a function of l
KE_m      = FLTARR(FIX(lmax/minc) + 1, max_nprnt)	;KE density as a function of m
ME_m      = FLTARR(FIX(lmax/minc) + 1, max_nprnt)	;ME density as a function of m

t_arr = 0.0
l1_arr = fltarr(lmax +1)
l2_arr = fltarr(lmax +1)
m1_arr = fltarr(FIX(lmax/minc) +1)
m2_arr = fltarr(FIX(lmax/minc) +1)


arrhagay=fltarr(5)

OPENR, 1, fname1
ns = -1
num = ' '
WHILE (NOT EOF(1)) DO BEGIN
	ns = ns + 1
	IF (ns+1 gt max_nprnt) THEN PRINT,'INCREASE max_nprnt VALUE.'

	READF, 1, t_arr, l1_arr
	READF, 1, t_arr, l2_arr
	READF, 1, t_arr, m1_arr
	READF, 1, t_arr, m2_arr

	t_spectra(ns)     = t_arr
  	KE_l(*, ns)       = l1_arr          	
  	ME_l(*, ns)       = l2_arr          	
  	KE_m(*, ns)       = m1_arr          	
  	ME_m(*, ns)       = m2_arr          	
ENDWHILE


;READ TIMESERIES DATA FROM L.'FNAME' FOR n TIMESTEPS
OPENR, 2, fname2
n = 0
bigarr = FLTARR(11, 5000)
arr    = FLTARR(11)
WHILE (NOT EOF(2)) DO BEGIN
	n = n + 1
	READF, 2, arr
	bigarr(*, n-1) = arr

	readf,2,arrhagay
ENDWHILE
tseries = FLTARR(11, n)
tseries = bigarr(*, 0:n-1)

time       = tseries(0,*)
ke	   = tseries(1,*)
pol_ke     = tseries(2,*)
me         = tseries(3,*)
pol_me     = tseries(4,*)
axi_tor_ke = tseries(5,*)
axi_pol_ke = tseries(6,*)
axi_pol_me = tseries(7,*)
axi_tor_me = tseries(8,*)
top_nu 	   = tseries(9,*)
bot_nu	   = tseries(10,*)

CLOSE, 1
CLOSE, 2


;CALCULATE MAGNETIC ENERGY SCALE FACTOR FOR PLOTTING
mefac=1.
if max(me) eq 0 then mefac=0
if max(me) gt 0 and max(ke) gt 0 then mefac=max(ke)/max(me)


;#######################################################################
; Here begins plotting material:
;#######################################################################
; Define window size and plotting frames
  XDIM=18.0
  YDIM=24.0
  IGIFPS=0

  LCOLOR = 0 
  LOADCT, LCOLOR 		;loads B&W color table


;** WINDOW SIZE IN PIXELS (ALSO SIZE OF PS FILE)
  XWINDOW=500
  YWINDOW=575
  SCWINDOW=1.25
  XWIND=XWINDOW*SCWINDOW
  YWIND=YWINDOW*SCWINDOW
  CSZ=1.00 & CSB=1.50 & CSS=0.720       ; CHARACTER SIZES
  SZF=1.0                               ; CHARACTER SIZE FACTOR

;*** Normalized coordinates
  XY1=[ 1.5/XDIM,14./YDIM, 8.5/XDIM,22./YDIM] & XT1=5./XDIM & YT1=22.2/YDIM
	XP1 = 6.5/XDIM  & YP1 = 20.5/YDIM
  XY2=[ 10.5/XDIM,14./YDIM,17.5/XDIM,22./YDIM] & XT2=14./XDIM & YT2=22.2/YDIM
	XP2 = 15.5/XDIM  & YP2 = 20.5/YDIM
  XY3=[ 1.5/XDIM, 3.5/YDIM, 8.5/XDIM,11.5/YDIM] & XT3=5./XDIM & YT3=11.7/YDIM
  XY4=[ 10.5/XDIM, 3.5/YDIM,17.5/XDIM,11.5/YDIM] & XT4=14./XDIM & YT4=11.7/YDIM

;*************************************************************************
;**************  THE PLOTTING MENU STARTS HERE **************************
;*************************************************************************
LABEL0: IF IGIFPS EQ 1 THEN PRINT,"USE --GOTO, LABELOUT-- STATEMENT HERE."
        PRINT,'  ' & PRINT,'  ' & PRINT,'  ' & PRINT,'  ' & PRINT,'  ' 
        PRINT, "ENTER OPTION: "
	PRINT, ' '
        PRINT, "EXIT PROCEDURE = -1"
	PRINT, "TIMESERIES PLOTS = 1"
	PRINT, "SPECTRAL PLOTS = 2"
        PRINT, "CHANGE WINDOW SIZE=11"
        READ,IOPTION

        IF (IOPTION GE 1 AND IOPTION LE 2) THEN BEGIN
           WINDOW,0,xsize=XWIND,ysize=YWIND,title='GRAPHIC WINDOW'
           !P.CHARTHICK=1.5
           !P.FONT=0
        ENDIF

        CASE IOPTION OF
       -1:    GOTO, LABEL99
        0:    GOTO, LABEL0
        1:    GOTO, LABEL1
        2:    GOTO, LABEL2
        11:   GOTO, LABEL11          
        ELSE: GOTO, LABEL0
        ENDCASE
 
;####################################################################

;11111111111111111111111111111111111111111111111111111111111111111111111
LABEL1: &


;ERASE
!P.MULTI = [0,2,2,0,0]

!P.POSITION = XY1
PLOT,  time, ke, LINESTYLE = 0, YTITLE='Energy Density', XTITLE='Time'
OPLOT, time, mefac*me, LINESTYLE = 1 
XYOUTS,XT1,YT1,'Mean KE(solid) & ME(dashed) Density',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,/NORMAL
XYOUTS,XT1,YT1-0.2,'ME scale='+strtrim(string(mefac))

!P.POSITION = XY2
PLOT,  time, axi_pol_ke, LINESTYLE = 0, YTITLE='Energy Density', XTITLE='Time'
OPLOT, time, mefac*axi_pol_me, LINESTYLE = 1
XYOUTS,XT2,YT2,'Axisymmetric Poloidal KE & ME Density',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,/NORMAL

!P.POSITION = XY3
PLOT,  time, axi_tor_ke, LINESTYLE = 0, YTITLE='Energy Density', XTITLE='Time'
OPLOT, time, mefac*axi_tor_me, LINESTYLE = 1
XYOUTS,XT3,YT3,'Axisymmetric Toroidal KE & ME Density',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,/NORMAL

!P.POSITION = XY4
PLOT,  time, bot_nu, LINESTYLE = 1, YTITLE='Nu', XTITLE='Time'
OPLOT, time, top_nu, LINESTYLE = 0
XYOUTS,XT4,YT4,'Top & Bottom(dashed) Nusselt number',CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,/NORMAL

XYOUTS,0.5,22.7/YDIM,RUNNAME,/NORMAL,CHARSIZE=1.5*CSZ*SZF,ALIGNMENT=0.5

;POSTSCRIPT COPY
PRINT, 'SAVE AS POSTSCRIPT =1; NO SAVE =0:'
READ, IFPS
IF (IFPS EQ 1) THEN BEGIN
OUTFILE=''
PRINT, 'ENTER .PS OUTFILE NAME:'
READ,OUTFILE
SET_PLOT, 'PS'
DEVICE,FILENAME=OUTFILE
DEVICE,XSIZE=XDIM,YSIZE=YDIM,XOFFSET=1.8,YOFFSET=2.0

;ERASE
!P.MULTI = [0,2,2,0,0]

!P.POSITION = XY1
PLOT,  time, ke, LINESTYLE = 0, YTITLE='Energy Density', XTITLE='Time'
OPLOT, time, mefac*me, LINESTYLE = 1 
XYOUTS,XT1,YT1,'Mean KE(solid) & ME(dashed) Density',ALIGNMENT=0.5,/NORMAL
XYOUTS,XT1,YT1-0.2,'ME scale='+strtrim(string(mefac))

!P.POSITION = XY2
PLOT,  time, axi_pol_ke, LINESTYLE = 0, YTITLE='Energy Density', XTITLE='Time'
OPLOT, time, mefac*axi_pol_me, LINESTYLE = 1
XYOUTS,XT2,YT2,'Axisymmetric Poloidal KE & ME Density',ALIGNMENT=0.5,/NORMAL

!P.POSITION = XY3
PLOT,  time, axi_tor_ke, LINESTYLE = 0, YTITLE='Energy Density', XTITLE='Time'
OPLOT, time, mefac*axi_tor_me, LINESTYLE = 1
XYOUTS,XT3,YT3,'Axisymmetric Toroidal KE & ME Density',ALIGNMENT=0.5,/NORMAL

!P.POSITION = XY4
PLOT,  time, bot_nu, LINESTYLE = 1, YTITLE='Nu', XTITLE='Time'
OPLOT, time, top_nu, LINESTYLE = 0
XYOUTS,XT4,YT4,'Top & Bottom(dashed) Nusselt number',ALIGNMENT=0.5,/NORMAL

XYOUTS,0.5,22.7/YDIM,RUNNAME,/NORMAL,CHARSIZE=CSZ*SZF,ALIGNMENT=0.5

DEVICE,/CLOSE
SET_PLOT,'X'
ENDIF

;!P.MULTI = 0 
GOTO, LABEL0
;111111111111111111111111111111111111111111111111111111111111111111111111

;222222222222222222222222222222222222222222222222222222222222222222222222
LABEL2: PRINT,' ' & PRINT,' '
index = ' '
PRINT,'Number of spectra =', ns + 1
PRINT,'Enter number (-1 = return):'
READ, index
IF (index eq -1) THEN GOTO, LABEL0

str_time = 'Time='+ STRTRIM(STRING(t_spectra(index-1)), 2)
;ERASE


;PLOT TO SCREEN

!P.MULTI = [0,2,2,0,0]
!P.POSITION = XY1
PLOT,  KE_l(*,index-1), LINESTYLE = 0, YTITLE='Energy Density', XTITLE='l Value'
OPLOT, mefac*ME_l(*,index-1), LINESTYLE = 1 
XYOUTS, XT1, YT1, 'KE(solid) ME(dots)', CHARSIZE=1/25*CSZ*SZF, ALIGNMENT=0.5, /NORMAL
XYOUTS, XP1, YP1, str_time, CHARSIZE=1.25*CSZ*SZF, ALIGNMENT=0.5, /NORMAL

!P.POSITION = XY2
PLOT,  KE_m(*,index-1), LINESTYLE = 0, YTITLE='Energy Density', XTITLE='m Value'
OPLOT, mefac*ME_m(*,index-1), LINESTYLE = 1
XYOUTS,XT2,YT2,'ME scale='+strtrim(string(mefac)),CHARSIZE=CSZ*SZF,ALIGNMENT=0.5,/NORMAL

XYOUTS,0.5,22.7/YDIM,RUNNAME,/NORMAL,CHARSIZE=1.5*CSZ*SZF,ALIGNMENT=0.5

;POSTSCRIPT COPY
PRINT, 'SAVE AS POSTSCRIPT =1; NO SAVE =0:'
READ, IFPS
IF (IFPS EQ 1) THEN BEGIN
OUTFILE=''
PRINT, 'ENTER .PS OUTFILE NAME:'
READ,OUTFILE
SET_PLOT, 'PS'
DEVICE,FILENAME=OUTFILE
DEVICE,XSIZE=XDIM,YSIZE=YDIM,XOFFSET=1.8,YOFFSET=2.0

!P.MULTI = [0,2,2,0,0]
!P.POSITION = XY1
PLOT,  KE_l(*,index-1), LINESTYLE = 0, YTITLE='Energy Density', XTITLE='l Value'
OPLOT, mefac*ME_l(*,index-1), LINESTYLE = 1 
XYOUTS, XT1, YT1, 'KE(solid) ME(dots)', ALIGNMENT=0.5, /NORMAL
XYOUTS, XP1, YP1, str_time,  ALIGNMENT=0.5, /NORMAL

!P.POSITION = XY2
PLOT,  KE_m(*,index-1), LINESTYLE = 0, YTITLE='Energy Density', XTITLE='m Value'
OPLOT, mefac*ME_m(*,index-1), LINESTYLE = 1
XYOUTS,XT2,YT2,'ME scale='+strtrim(string(mefac)),ALIGNMENT=0.5,/NORMAL

XYOUTS,0.5,22.7/YDIM,RUNNAME,/NORMAL,CHARSIZE=CSZ*SZF,ALIGNMENT=0.5

DEVICE,/CLOSE
SET_PLOT,'X'
ENDIF


;!P.MULTI = 0 
GOTO, LABEL2
;22222222222222222222222222222222222222222222222222222222222222222222222

;11 1111 11  11 11 1111 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11  11 
LABEL11:  PRINT, FORMAT='("WINDOW SCALE FACTOR?  CURRENT=",F7.3)',SCWINDOW
          READ, SCWINDOW
          XWIND=XWINDOW*SCWINDOW & YWIND=YWINDOW*SCWINDOW
          CSZ=SCWINDOW*0.8 & CSB=1.12*SCWINDOW & CSS=0.60*SCWINDOW
          GOTO, LABEL0
;11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 

;***********************************************************************


;99 99 99 99 99 99 99 99 99 99 99 99 99 99 99 99 99 99 99 99 99 99 99 99
LABEL99: !P.MULTI = 0
	 SET_PLOT, 'X'	
			

END
