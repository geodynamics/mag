c formulas in unscramble.f were taken from the prep.f and other mag files
c it is used to unscramble the cc-files, converting the cc-file
c output to arrays in harmonic degree l and order m. 
c It defines arrays al(lm) and am(lm) for lm=1,nlma, 
c which contain the l and m  values we need to convert the cc-file contents


c input lmax, minc, and nlma  from cc file header HERE:
      lmax=??
      minc=??
      nlma=??

c define intermediate indices  
      mmax=(lmax/minc)*minc
      nmaf=mmax+1
      nlaf=lmax+1

c define the unscramble array mclm(lm) 	
      lm=0
      do 35 mc=1,nmaf,minc
      do 31 lc=mc,nlaf
       lm=lm+1
       mclm(lm)=mc
   31 continue
   35 continue

c define la, ma arrays
     do 36 lm=1,nlma
      ma(lm)=mclm(lm)-1
c define al in three terms
      tl1=lm+ma(lm)-1
      tl2=ma(lm)*(ma(lm)-minc)/(2.*minc)
      tl3=-ma(lm)*(lmax+1)/minc
      la(lm)=tl1+tl2+tl3
c PRINT lm, la(lm), ma(lm) HERE
  36 continue

c if this works, then read the contents of the cc-file block in pairs, 
c and assign the new indices
  
    do 37 lm=1,nlma
     READ (cc-file) c1,c2
c assign new indices
     l=la(lm)
     m=ma(lm)
     anew(l,m)=c1
     bnew(l,m)=c2
  37 continue


     