      subroutine cmbcoeff(kc)
c
c  called in amhd
c  supplied by urc
c
c  output of complex spherical harmonic coefficients
c  for the poloidal magnetic potential b at the outer
c  boundary and of the poloidal velocity potential w,
c  its radial derivative dw, and the toroidal velocity
c  potential z at the radial level kc
c
      include 'param.f'
      include 'com1.f'
      include 'com3.f'
      include 'com5.f'
      include 'com8.f'
c
c  write header
c
      write(21,2100) nlma,lmax,minc,r(1),r(kc),time/tscale
 2100 format(/, 4x,"nlma=",i3,4x,"lmax=",i3,4x,"minc=",i3,4x,
     $   "r(1)=",f7.4,4x,"r(kc)=",f7.4,4x,"time/tscale=",f9.6)
c
c  write data
c
      write(21,2101) (b(lm,1),lm=1,nlma)
      write(21,2102) (w(lm,kc),lm=1,nlma)
      write(21,2102) (dw(lm,kc),lm=1,nlma)
      write(21,2102) (z(lm,kc),lm=1,nlma)
 2101 format(256(1X,f9.5))
 2102 format(256(1X,f9.3))
c
      return
      end
