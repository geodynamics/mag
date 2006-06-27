      subroutine rffti(n,wsave)
c
c     called in chebi
c
      dimension wsave(*)
c
      ns2 = n/2
      nqm = (ns2-1)/2
      tpi = 8.*atan(1.)
      dt = tpi/float(n)
      dc = cos(dt)
      ds = sin(dt)
      wsave(1) = dc
      wsave(ns2-1) = ds
      if (nqm .lt. 2) return
      do 101 k=2,nqm
      kc = ns2-k
      wsave(k) = dc*wsave(k-1)-ds*wsave(kc+1)
      wsave(kc) = ds*wsave(k-1)+dc*wsave(kc+1)
  101 continue
      return
      end
