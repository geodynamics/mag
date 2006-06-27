      subroutine pbar(the,l,m,p)
c
c  pbar calculates the value of the normalized associated
c  legendre function of the first kind, of degree l,
c  of order m, for the real argument cos(the), and returns
c  it in the variable p
c  0 .le. m .le. l
c
c  called in gquad and prep
c
      s=sin(the)
      c=cos(the)
      p=1./sqrt(2.)
      if(m .eq. 0) go to 22
      do 20 i=1,m
      p=sqrt(float(2*i+1)/float(2*i))*s*p
   20 continue
   22 continue
      if(l .eq. m) return
      p1=1.
      m1=m+1
      do 30 j=m1,l
      p2=p1
      p1=p
      p=2.*sqrt((float(j**2)-0.25)/float(j**2-m**2))*c*p1-
     $sqrt(float((2*j+1)*(j-m-1)*(j+m-1))/
     $float((2*j-3)*(j-m)*(j+m)))*p2
   30 continue
      return
      end
