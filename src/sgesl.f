      subroutine sgesl(a,ia,n,ip,b,ijob)
c
c  ************************************************************
c     like the linpack routine
c
c     solves  a * x = b  via lu decomposition
c
c     sub sgefa must be called once first to initialize a and ip
c
c     n is the order of the linear matrix equation
c
c     ia .ge. n .gt. 1
c
c     on return, the solution vector x is stored in b
c
c     called in amhd and prep
c  ************************************************************
c
      dimension a(ia,*),b(*),ip(*)
c
      np1=n+1
      nm1=n-1
c urc if(ijob .ne. 0) stop '44'                                    
c
c     permute vector b
c
      do 1 k=1,nm1
      m=ip(k)
      if(m .eq. k) go to 1
      c=b(m)
      b(m)=b(k)
      b(k)=c
    1 continue
c
c     solve  l * y = b
c
      do 2 k=1,nm1
      kp1=k+1
      do 2 i=kp1,n
      b(i)=b(i)-b(k)*a(i,k)
    2 continue
c
c     solve  u * x = y
c
      do 3 kb=1,nm1
      k=np1-kb
      b(k)=b(k)/a(k,k)
      km1=k-1
      do 3 i=1,km1
      b(i)=b(i)-b(k)*a(i,k)
    3 continue
      b(1)=b(1)/a(1,1)
c
      return
      end
