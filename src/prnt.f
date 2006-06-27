      subroutine prnt
c
c---------------------------------------------------------------
      include 'param.f'
      include 'com1.f'
      include 'com3.f'
      include 'com4.f'
      include 'com5.f'
      include 'com8.f'
c---------------------------------------------------------------
c
c urc: modified to find absolute maximum in terms of velocity/field
c
      absmax=0.
      kcmax=0
      lmmax=0
      do 6501 kc=1,nn
      do 6501 lm=2,nlma
      l=nint(ql(lm,4))
      a1=l*(l+1)*abs(w(lm,kc))*qk(kc,1)
      if(a1 .gt. absmax) then       
       absmax=a1
       kcmax=kc
       lmmax=lm
      endif
 6501 continue
      m=mclm(lmmax)-1
      l=nint(ql(lmmax,4))
      rad=r(kcmax)-r(nn)
      write(6,6500) m,l,rad,absmax/vscale
 6500 format(/,4x,"w",4x,"m =",i3,4x,"l =",i3,4x,"r =",f7.4,
     $   4x,"absmax =",f9.3)
c
      absmax=0.
      kcmax=0
      lmmax=0
      do 6601 kc=1,nn
      do 6601 lm=2,nlma
      l=nint(ql(lm,4))
      a1=l*abs(z(lm,kc))*qk(kc,3)
      if(a1 .gt. absmax) then        
       absmax=a1
       kcmax=kc
       lmmax=lm
      endif
 6601 continue
      m=mclm(lmmax)-1
      l=nint(ql(lmmax,4))
      rad=r(kcmax)-r(nn)
      write(6,6600) m,l,rad,absmax/vscale
 6600 format(4x,"z",4x,"m =",i3,4x,"l =",i3,4x,"r =",f7.4,
     $   4x,"absmax =",f9.3)
c
      absmax=0.
      kcmax=0
      lmmax=0
      do 6701 kc=1,nn
      do 6701 lm=2,nlma
      a1=abs(s(lm,kc))
      if(a1 .gt. absmax) then         
       absmax=a1
       kcmax=kc
       lmmax=lm
      endif
 6701 continue
      m=mclm(lmmax)-1
      l=nint(ql(lmmax,4))
      rad=r(kcmax)-r(nn)
      write(6,6700) m,l,rad,absmax
 6700 format(4x,"s",4x,"m =",i3,4x,"l =",i3,4x,"r =",f7.4,
     $   4x,"absmax =",f9.5)
c
      absmax=0.
      kcmax=0
      lmmax=0
      do 6801 kc=1,nn
      do 6801 lm=2,nlma
      a1=abs(p(lm,kc))
      if(a1 .gt. absmax) then       
       absmax=a1
       kcmax=kc
       lmmax=lm
      endif
 6801 continue
      m=mclm(lmmax)-1
      l=nint(ql(lmmax,4))
      rad=r(kcmax)-r(nn)
      write(6,6800) m,l,rad,absmax/pscale
 6800 format(4x,"p",4x,"m =",i3,4x,"l =",i3,4x,"r =",f7.4,
     $   4x,"absmax =",f9.4)
c
      absmax=0.
      kcmax=0
      lmmax=0
      do 691 kc=1,nn
      do 691 lm=2,nlma
      l=nint(ql(lm,4))
      a1=l*(l+1)*abs(b(lm,kc))*qk(kc,1)
      if(a1 .gt. absmax) then       
       absmax=a1
       kcmax=kc
       lmmax=lm
      endif
  691 continue
      m=mclm(lmmax)-1
      l=nint(ql(lmmax,4))
      rad=r(kcmax)-r(nn)
      write(6,690) m,l,rad,absmax
  690 format(4x,"b",4x,"m =",i3,4x,"l =",i3,4x,"r =",f7.4,
     $   4x,"absmax =",f9.5)
c
      absmax=0.
      kcmax=0
      lmmax=0
      do 701 kc=1,nn
      do 701 lm=2,nlma
      l=nint(ql(lm,4))
      a1=l*abs(aj(lm,kc))*qk(kc,3)
      if(a1 .gt. absmax) then        
       absmax=a1
       kcmax=kc
       lmmax=lm
      endif
  701 continue
      m=mclm(lmmax)-1
      l=nint(ql(lmmax,4))
      rad=r(kcmax)-r(nn)
      write(6,700) m,l,rad,absmax
  700 format(4x,"j",4x,"m =",i3,4x,"l =",i3,4x,"r =",f7.4,
     $   4x,"absmax =",f9.5)
c
      return
      end
