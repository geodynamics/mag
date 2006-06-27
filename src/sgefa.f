      subroutine sgefa(a,ia,n,ip,info)
c
c  *******************************************************************
c     like the linpack routine
c
c     lu decomposes the real matrix a(n,n) via gaussian elimination
c
c     called in ludc and prep
c  *******************************************************************
c
      dimension a(ia,n),ip(n)
c
      if(n .le. 1) stop '45'
      info=0
      nm1=n-1
c
      do 50 k=1,nm1
      kp1=k+1
      l=k
      do 60 i=kp1,n
      if(abs(a(i,k)) .gt. abs(a(l,k))) l=i
   60 continue
      ip(k)=l
      if(a(l,k) .eq. 0.) go to 40
      if(l .eq. k) go to 10
      do 70 i=1,n
      t=a(k,i)
      a(k,i)=a(l,i)
      a(l,i)=t
   70 continue
   10 continue
      t=1./a(k,k)
      do 80 i=kp1,n
      a(i,k)=t*a(i,k)
   80 continue
      do 30 j=kp1,n
      do 90 i=kp1,n
      a(i,j)=a(i,j)-a(k,j)*a(i,k)
   90 continue
   30 continue
      go to 50
   40 continue
      info=k
   50 continue
c
      ip(n)=n
      if(a(n,n) .eq. 0.) info=n
c
      return
      end
