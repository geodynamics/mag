      subroutine cftrig(n,trigs)
c
c     called in chebi
c
      dimension trigs(*)
      pi=2.0*asin(1.0)
      del=(pi+pi)/float(n)
      l=n+n
      do 10 i=1,l,2
      angle=0.5*float(i-1)*del
      trigs(i)=cos(angle)
      trigs(i+1)=sin(angle)
   10 continue
      return
      end
