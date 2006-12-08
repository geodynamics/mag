;Procedure magxy
;
;INPUT PARAMETERS
runname  = ' '		;MAG output suffix 
fname   = ' '      	;filename of spectra data
lmax     = 1 		;lmax value of grid, probably 32, 42, 64, 96, etc
mmax 	 = 1		;m-max value of grid
fname=' '               ;ls-file name

set_plot,'x'
device, retain=2, decomposed=0


PRINT,'  '
PRINT,'Enter ls file name:' 
READ, fname

PRINT,' '
PRINT,'Enter Degree truncation Lmax:' & read, lmax
PRINT,'Enter Order truncation Mmax:' & read, mmax

;READ IN SPECTRAL DATA IN LS.'FNAME' FOR ns DIFFERENT NPRNT VALUES.
max_nprnt = 250

time      = FLTARR(max_nprnt)			;time of spectra 
KE_L      = FLTARR(lmax + 1, max_nprnt)		;KE density as a function of l
ME_L      = FLTARR(lmax + 1, max_nprnt)		;ME density as a function of l
KE_m      = FLTARR(mmax + 1, max_nprnt)		;KE density as a function of m
ME_m      = FLTARR(mmax + 1, max_nprnt)		;ME density as a function of m

t_arr = 0.0
l1_arr = fltarr(lmax +1)
l2_arr = fltarr(lmax +1)
m1_arr = fltarr(mmax +1)
m2_arr = fltarr(mmax +1)

OPENR, 1, fname
ns = -1
num = ' '
WHILE (NOT EOF(1)) DO BEGIN
	ns = ns + 1
	IF (ns+1 gt max_nprnt) THEN PRINT,'INCREASE max_nprnt VALUE.'

	READF, 1, t_arr, l1_arr
	READF, 1, t_arr, l2_arr
	READF, 1, t_arr, m1_arr
	READF, 1, t_arr, m2_arr

	time(ns)          = t_arr
  	KE_l(*, ns)       = l1_arr          	
  	ME_l(*, ns)       = l2_arr          	
  	KE_m(*, ns)       = m1_arr          	
  	ME_m(*, ns)       = m2_arr          	
ENDWHILE

CLOSE, 1
;;;;;;END OF DATA READ;;;;;;;;;;;

;make trimed spectral arrays
lp=lmax
mp=mmax
KEL=KE_l(0:lp,0:ns)
MEL=ME_l(0:lp,0:ns)
KEM=KE_m(0:mp,0:ns)
MEM=ME_m(0:mp,0:ns)

KEavel    = FLTARR(lp+1, ns)	 ;ave KE density as a function of l
MEavel    = FLTARR(lp+1)     ;ave ME density as a function of l
KEavem    = FLTARR(mp+1)     ;ave KE density as a function of m
MEavem    = FLTARR(mp+1)     ;ave ME density as a function of m

KEstdl    = FLTARR(lp+1)	;sd KE density as a function of l
MEstdl    = FLTARR(lp+1)	;sd ME density as a function of l
KEstdm    = FLTARR(mp+1)	;sd KE density as a function of m
MEstdm    = FLTARR(mp+1)	;sd ME density as a function of m

;Calculate mean and standard deviation spectra

for ln=0,lp do begin
    result=moment(MEL(ln,*),sdev=sd) 
    MEavel(ln)=result(0) 
    MEstdl(ln)=sd
endfor	

for ln=0,lp do begin
    result=moment(KEL(ln,*),sdev=sd) 
    KEavel(ln)=result(0) 
    KEstdl(ln)=sd
endfor	

for mn=0,mp do begin
    result=moment(MEM(mn,*),sdev=sd) 
    MEavem(mn)=result(0) 
    MEstdm(mn)=sd
endfor	

for mn=0,mp do begin
    result=moment(KEM(mn,*),sdev=sd) 
    KEavem(mn)=result(0) 
    KEstdm(mn)=sd
endfor	
;end of mean and standard deviation spectra calculation


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
	PRINT, "TIME ave SPECTRA = 1"
	PRINT, "SPECTRA vs TIME = 2"
        PRINT, "CHANGE WINDOW SIZE=11"
	PRINT, "PRINT PS = 21"
        READ,IOPTION

        IF (IOPTION GE 1 AND IOPTION LE 2) THEN BEGIN
           WINDOW,0,xsize=XWIND,ysize=YWIND,title='GRAPHIC WINDOW'
           !P.CHARTHICK=1.5
           !P.FONT=0
        ENDIF

        CASE IOPTION OF
       -1:    GOTO, LABEL98
        0:    GOTO, LABEL0
        1:    GOTO, LABEL1
        2:    GOTO, LABEL2
        11:   GOTO, LABEL11 
	21:   GOTO, LABEL21         
        ELSE: GOTO, LABEL0
        ENDCASE
 
;####################################################################

;11111111111111111111111111111111111111111111111111111111111111111111111
LABEL1: IPAGE=1 
ERASE
lend=lp & mend=mp
PRINT, 'ENTER END POINTS lend, mend:'
read,lend,mend

;PLOT TO SCREEN

!P.MULTI = [0,2,2,0]

PLOT,  KEavel+KEstdl, LINESTYLE = 1, YTITLE='Power', XTITLE='Degree, l',$
	 TITLE="KE",XRANGE=[0,lend]
OPLOT, KEavel-KEstdl, LINESTYLE = 1 
OPLOT, KEavel, LINESTYLE = 0 

PLOT,  MEavel+MEstdl, LINESTYLE = 1, YTITLE='Power', XTITLE='Degree, l',$
	 TITLE="ME",XRANGE=[0,lend]
OPLOT, MEavel-MEstdl, LINESTYLE = 1 
OPLOT, MEavel, LINESTYLE = 0 


PLOT,  KEavem+KEstdm, LINESTYLE = 1, YTITLE='Power', XTITLE='Order, m', $
	TITLE="KE",XRANGE=[0,mend]
OPLOT, KEavem-KEstdm, LINESTYLE = 1 
OPLOT, KEavem, LINESTYLE = 0 

PLOT,  MEavem+MEstdm, LINESTYLE = 1, YTITLE='Power', XTITLE='Order, m', $
	TITLE="ME",XRANGE=[0,mend]
OPLOT, MEavem-MEstdm, LINESTYLE = 1 
OPLOT, MEavem, LINESTYLE = 0 


;!P.MULTI = 0 
GOTO, LABEL99
;1111111111111111111111111111111111111111111111111111111111111111111111111111


;222222222222222222222222222222222222222222222222222222222222222222222
LABEL2: IPAGE=2
ERASE

PRINT, 'CHANGE TIME SCALE'
tscale=1
PRINT, ' '
PRINT, 'ENTER TIME SCALE FACTOR (def=1):' 
read,tscal
;calculate start, end times
tstart = tscal*time(0) & tend=tscal*time(ns)
tstart=string(tstart)  & tend=string(tend)
sscale = "Time on, off = " + tstart + ' , ' + tend

lend=lp & mend=mp
PRINT, 'ENTER END POINTS lend, mend:'
read,lend,mend

xl=findgen(lp+1)
xm=findgen(mp+1)
y=findgen(ns+1)

;transpose and smooth spectra
kelt=transpose(smooth(KEL,3,/edge_truncate))
melt=transpose(smooth(MEL,3,/edge_truncate))
kemt=transpose(smooth(KEM,3,/edge_truncate))
memt=transpose(smooth(MEM,3,/edge_truncate))

;PLOT TO SCREEN

!P.MULTI = [0,1,4,0]

LOADCT,39

klevels=max(kelt)*findgen(20)/20
CONTOUR,kelt,y,xl,TITLE="KE(l)",levels=klevels,$
	xstyle=1,ystyle=1,yrange=[0,lend],/fill,ytitle='Degree, l'

mlevels=max(melt)*findgen(20)/20
CONTOUR,melt,y,xl,TITLE="ME(l)",levels=mlevels,$
	xstyle=1,ystyle=1,yrange=[0,lend],/fill,ytitle='Degree, l'

klevelm=max(kemt)*findgen(20)/20
CONTOUR,kemt,y,xm,TITLE="KE(m)",levels=klevelm,$
	xstyle=1,ystyle=1,yrange=[0,mend],/fill,ytitle='Order m'

mlevelm=max(memt)*findgen(20)/20
CONTOUR,memt,y,xm,TITLE="ME(m)",levels=mlevelm,$
	xstyle=1,ystyle=1,yrange=[0,mend],/fill,ytitle='Order m',$
	xtitle=sscale


;!P.MULTI = 0 
GOTO, LABEL99
;22222222222222222222222222222222222222222222222222222222222222

;11 1111 11  11 11 1111 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11  11 
LABEL11:  PRINT, FORMAT='("WINDOW SCALE FACTOR?  CURRENT=",F7.3)',SCWINDOW
          READ, SCWINDOW
          XWIND=XWINDOW*SCWINDOW & YWIND=YWINDOW*SCWINDOW
          CSZ=SCWINDOW*0.8 & CSB=1.12*SCWINDOW & CSS=0.60*SCWINDOW
          GOTO, LABEL0
;11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 11 

;***********************************************************************

;21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21
LABEL21: &
OUTFILE= ' '
PRINT,'ENTER FULL .PS FILE NAME:'
READ, OUTFILE 
SET_PLOT,'PS'
DEVICE,FILENAME=OUTFILE,/COLOR 
DEVICE,XSIZE=XDIM,YSIZE=YDIM,XOFFSET=1.8,YOFFSET=2.0 

       CASE IPAGE OF
        1:    GOTO, LABEL1
        2:    GOTO, LABEL2
       ELSE: GOTO, LABEL0
        ENDCASE

;21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21 21

;99 99 99 99 99 99 99 99 99 99 99 99 99 99 99 99 99 99 99 99 99 99 99 99
LABEL99: IF IOPTION EQ 21 THEN DEVICE,/CLOSE 
	 !P.MULTI = 0
	 SET_PLOT, 'X'	
	 GOTO,LABEL0	
;99 99 99 99 99 99 99 99 99 99 99 99 99 99 99 99 99 99 99 99			

LABEL98: PRINT,'EXIT MAGXY'

END
