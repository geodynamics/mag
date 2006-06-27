      subroutine fftrig(trigs,n,mode)
c
c     called in prep
c
      dimension trigs(*)
      pi=2.0*asin(1.0)
      imode=iabs(mode)
      nn=n
      if (imode.gt.1.and.imode.lt.6) nn=n/2
      del=(pi+pi)/float(nn)
      l=nn+nn
      do 10 i=1,l,2
      angle=0.5*float(i-1)*del
      trigs(i)=cos(angle)
      trigs(i+1)=sin(angle)
   10 continue
      if (imode.eq.1) return
      if (imode.eq.8) return
      del=0.5*del
      nh=(nn+1)/2
      l=nh+nh
      la=nn+nn
      do 20 i=1,l,2
      angle=0.5*float(i-1)*del
      trigs(la+i)=cos(angle)
      trigs(la+i+1)=sin(angle)
   20 continue
      if (imode.le.3) return
      del=0.5*del
      la=la+nn
      if (mode.eq.5) go to 40
      do 30 i=2,nn
      angle=float(i-1)*del
      trigs(la+i)=2.0*sin(angle)
   30 continue
      return
   40 continue
      del=0.5*del
      do 50 i=2,n
      angle=float(i-1)*del
      trigs(la+i)=sin(angle)
   50 continue
      return
      end
