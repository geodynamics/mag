      subroutine gquad(l,root,w)
c
c  gquad (linked with pbar) finds the l roots (in theta)
c  and gaussian weights associated with
c  the legendre polynomial of degree l > 1
c
c  called in prep
c
      dimension root(l),w(l)
c
      pi=4.*atan(1.)
      del=pi/float(4*l)
      l1=l+1
      co=float(2*l+3)/float(l1**2)
      p2=1.
      t2=-del
      l2=l/2
      k=1
c
      do 10 i=1,l2
   20 t1=t2
      t2=t1+del
      theta=t2
      call pbar(theta,l,0,p)
      p1=p2
      p2=p
      if((k*p2) .gt. 0.) go to 20
      k=-k
   40 s=(t2-t1)/(p2-p1)
      t1=t2
      t2=t2-s*p2
      theta=t2
      call pbar(theta,l,0,p)
      p1=p2
      p2=p
      if(abs(p) .le. 1.e-10) go to 30
      if(p2 .eq. p1) then
c        write(6,*) 'sub gquad: zero = ',p,' at i = ',i
         go to 30
      endif
      go to 40
   30 root(i)=theta
      call pbar(theta,l1,0,p)
      w(i)=co*(sin(theta)/p)**2
   10 continue
c
      l22=2*l2
      if(l22 .eq. l) go to 70
      l2=l2+1
      theta=pi/2.
      root(l2)=theta
      call pbar(theta,l1,0,p)
      w(l2)=co/p**2
   70 continue
c
      l3=l2+1
      do 50 i=l3,l
      root(i)=pi-root(l-i+1)
      w(i)=w(l-i+1)
   50 continue
c
      return
      end
