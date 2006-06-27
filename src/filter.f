      subroutine filter(ain,aout,kc,alfilt,nfilt,dipfilt)
c
c    -filter depending on harmonic degree l acting on radial
c    -level kc of field ain(lm,kc)
c    -super-Gaussian filter  F = exp (-[l/alfilt]^nfilt)
c    -or cos-taper   F = 0.5(1-sin(pi*[l-nfilt]/|alfilt|)
c    -if lfilt<0
c    -result written into aout(lm)
c
c---------------------------------------------------------------
      include 'param.f'
      include 'com4.f'
c
      dimension filt(lmax)
      complex 
     $ ain(nlma,nnp1),aout(nlma)
c
      if (alfilt.gt.0.0) then
      do 200 l=1,lmax
         filt(l)=exp(-((real(l)/alfilt)**nfilt))
  200 continue
      else
      do 201 l=1,lmax
         arg=(l-nfilt)/alfilt
         if(arg.gt.0.5) then 
          filt(l)=1.0
         else if(arg.lt.-0.5) then
          filt(l)=0.0
         else
          filt(l)=0.5*(1.+sin(pi*arg))
         endif
  201 continue
      endif
c
      do 203 lm=nlma,2,-1
         l=nint(ql(lm,4))
         aout(lm)=ain(lm,kc)*filt(l)
  203 continue
      aout(1)=0.
      aout(2)=dipfilt*aout(2)
c
      return
      end
