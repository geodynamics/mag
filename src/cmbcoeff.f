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
      write(21) nlma,lmax,minc,r(1),r(kc),time/tscale
c
c  write data
c
      write(21) (b(lm,1),lm=1,nlma)
      write(21) (w(lm,kc),lm=1,nlma)
      write(21) (dw(lm,kc),lm=1,nlma)
      write(21) (z(lm,kc),lm=1,nlma)
c
      return
      end
