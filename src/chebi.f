      subroutine chebi(n,wsave,trigs,ifax,k2k)
c
c     initialization for sub chebtf
c
c     called in prep
c
      dimension wsave(*),trigs(*),ifax(*),k2k(*)
c
      if(n .le. 3) stop '41'
      nm1=n-1
      if(mod(nm1,4) .ne. 0) stop '42'
      np1=n+1
      ns2 = nm1/2
      ns2m = ns2-1
      iw1 = ns2+1
      nm2=n-2
      k2k(1)=1
      k2k(2)=2
c
      do 10 k=3,nm2,2
      k2k(k)=np1-k
      k2k(k+1)=k2k(k)+1
   10 continue
c
      pi = 4.*atan(1.)
      dt = pi/float(nm1)
      dcs = cos(dt)
      dss = sin(dt)
      wsave(1) = dss+dss
      wck = dcs+dcs
      if (ns2m .lt. 2) go to 102
      do 101 k=2,ns2m
      wsave(k) = dss*wck+dcs*wsave(k-1)
      wck = dcs*wck-dss*wsave(k-1)
  101 continue
  102 call rffti (nm1,wsave(iw1))
      call fact(ns2,ifax)
      k=ifax(1)
      if((k .lt. 1) .or. (ifax(k+1) .gt. 5)) stop '43'
      call cftrig(ns2,trigs)
c
      return
      end
