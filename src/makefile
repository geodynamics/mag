###################################################################
#$URL: $
#$Rev: $
#$Date: $
###################################################################

all: magx

 magx  : nl.o prep.o amhd.o legtf.o rderiv.o stor.o prnt.o ludc.o dtchck.o kei.o mei.o pbar.o gquad.o random.o chebtf.o chebi.o rffti.o fact.o cftrig.o fourtf.o fax.o fftrig.o sgesl.o sgefa.o fft99a.o fft99b.o vpassm.o wpassm.o spherictf.o filter.o graphout.o moveout.o movaout.o movmout.o cmbcoeff.o spectrum.o mapdata.o copydat.o stopiteration.o 
	g77 nl.o prep.o amhd.o legtf.o rderiv.o stor.o prnt.o ludc.o dtchck.o kei.o mei.o pbar.o gquad.o random.o chebtf.o chebi.o rffti.o fact.o cftrig.o fourtf.o fax.o fftrig.o sgesl.o sgefa.o fft99a.o fft99b.o vpassm.o wpassm.o spherictf.o filter.o graphout.o moveout.o movaout.o movmout.o cmbcoeff.o spectrum.o mapdata.o copydat.o stopiteration.o -o magx

nl.o	: nl.f param.f com3.f com8.f
	g77   -O2 -funroll-all-loops  -c -C  nl.f   
prep.o	: prep.f param.f com1.f com2.f com3.f com4.f com5.f com6.f com7.f com8.f
	g77   -O2 -funroll-all-loops  -c -C  prep.f   
amhd.o	: amhd.f param.f com1.f com2.f com3.f com4.f com5.f com6.f com7.f com8.f
	g77   -O2 -funroll-all-loops  -c -C  amhd.f   
legtf.o	: legtf.f param.f com4.f com5.f com6.f
	g77   -O2 -funroll-all-loops  -c -C  legtf.f   
rderiv.o	: rderiv.f param.f com7.f 
	g77   -O2 -funroll-all-loops  -c -C  rderiv.f
dtchck.o	: dtchck.f
	g77   -O2 -funroll-all-loops  -c -C  dtchck.f
kei.o	: kei.f param.f com1.f com4.f com5.f
	g77   -O2 -funroll-all-loops  -c -C  kei.f
mei.o	: mei.f param.f com1.f com4.f com5.f
	g77   -O2 -funroll-all-loops  -c -C  mei.f
pbar.o	: pbar.f
	g77   -O2 -funroll-all-loops  -c -C  pbar.f   
ludc.o	: ludc.f param.f com1.f com2.f com3.f com4.f com7.f
	g77   -O2 -funroll-all-loops  -c -C  ludc.f   
gquad.o	: gquad.f
	g77   -O2 -funroll-all-loops  -c -C  gquad.f   
random.o	: random.f
	g77   -O2 -funroll-all-loops  -c -C  random.f   
chebtf.o	: chebtf.f
	g77   -O2 -funroll-all-loops  -c -C  chebtf.f  
chebi.o	: chebi.f
	g77   -O2 -funroll-all-loops  -c -C  chebi.f   
rffti.o	: rffti.f
	g77   -O2 -funroll-all-loops  -c -C  rffti.f   
fact.o	: fact.f
	g77   -O2 -funroll-all-loops  -c -C  fact.f   
cftrig.o	: cftrig.f
	g77   -O2 -funroll-all-loops  -c -C  cftrig.f  
fourtf.o	: fourtf.f
	g77   -O2 -funroll-all-loops  -c -C  fourtf.f  
fax.o	: fax.f
	g77   -O2 -funroll-all-loops  -c -C  fax.f   
fftrig.o	: fftrig.f
	g77   -O2 -funroll-all-loops  -c -C  fftrig.f  
sgesl.o	: sgesl.f
	g77   -O2 -funroll-all-loops  -c -C  sgesl.f   
sgefa.o	: sgefa.f
	g77   -O2 -funroll-all-loops  -c -C  sgefa.f   
fft99a.o	: fft99a.f
	g77   -O2 -funroll-all-loops  -c -C  fft99a.f  
fft99b.o	: fft99b.f
	g77   -O2 -funroll-all-loops  -c -C  fft99b.f  
vpassm.o	: vpassm.f
	g77   -O2 -funroll-all-loops  -c -C  vpassm.f  
wpassm.o	: wpassm.f
	g77   -O2 -funroll-all-loops  -c -C  wpassm.f  
spherictf.o	: spherictf.f param.f com4.f
	g77   -O2 -funroll-all-loops  -c -C  spherictf.f  
filter.o	: filter.f param.f com4.f
	g77   -O2 -funroll-all-loops  -c -C  filter.f  
graphout.o	: graphout.f param.f com1.f com3.f com4.f com5.f com6.f com8.f
	g77   -O2 -funroll-all-loops  -c -C  graphout.f  
moveout.o	: moveout.f param.f com1.f com3.f com4.f com6.f com8.f
	g77   -O2 -funroll-all-loops  -c -C  moveout.f  
movaout.o	: movaout.f param.f com1.f com3.f com4.f com6.f com8.f
	g77   -O2 -funroll-all-loops  -c -C  movaout.f  
movmout.o	: movmout.f param.f com1.f com3.f com4.f com6.f com8.f
	g77   -O2 -funroll-all-loops  -c -C  movmout.f  
cmbcoeff.o	: cmbcoeff.f param.f com1.f com3.f com5.f com8.f
	g77   -O2 -funroll-all-loops  -c -C  cmbcoeff.f  
spectrum.o	: spectrum.f param.f com1.f com3.f com4.f com5.f com8.f
	g77   -O2 -funroll-all-loops  -c -C  spectrum.f  
stor.o	: stor.f param.f com1.f com3.f com4.f com5.f 
	g77   -O2 -funroll-all-loops  -c -C  stor.f  
prnt.o	: prnt.f param.f com1.f com3.f com4.f com5.f com8.f
	g77   -O2 -funroll-all-loops  -c -C  prnt.f  
mapdata.o	: mapdata.f param.f  com3.f  com5.f
	g77   -O2 -funroll-all-loops  -c -C  mapdata.f  
copydat.o	: copydat.f
	g77   -O2 -funroll-all-loops  -c -C  copydat.f  
stopiteration.o	: stopiteration.f
	g77   -O2 -funroll-all-loops  -c -C stopiteration.f

clean: 
		rm -f *.o magx*
