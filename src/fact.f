      subroutine fact(n,ifax)
c
c     factorization routine that first extracts all factors of 4
c
c     called in chebi
c
      dimension ifax(*)
      if (n.gt.1) go to 10
      ifax(1) = 0
      if (n.lt.1) ifax(1) = -99
      return
   10 nn=n
      k=1
c     test for factors of 4
   20 if (mod(nn,4).ne.0) go to 30
      k=k+1
      ifax(k)=4
      nn=nn/4
      if (nn.eq.1) go to 80
      go to 20
c     test for extra factor of 2
   30 if (mod(nn,2).ne.0) go to 40
      k=k+1
      ifax(k)=2
      nn=nn/2
      if (nn.eq.1) go to 80
c     test for factors of 3
   40 if (mod(nn,3).ne.0) go to 50
      k=k+1
      ifax(k)=3
      nn=nn/3
      if (nn.eq.1) go to 80
      go to 40
c     now find remaining factors
   50 l=5
      max = sqrt(float(nn))
      inc=2
c     inc alternately takes on values 2 and 4
   60 if (mod(nn,l).ne.0) go to 70
      k=k+1
      ifax(k)=l
      nn=nn/l
      if (nn.eq.1) go to 80
      go to 60
   70 if (l.gt.max) go to 75
      l=l+inc
      inc=6-inc
      go to 60
   75 k = k+1
      ifax(k) = nn
   80 ifax(1)=k-1
c     ifax(1) now contains number of factors
      return
      end
