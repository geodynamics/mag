      subroutine rderiv(od,anorm,f,df)
c
c  calculates derivative of Chebychev polynomia
c  called in prep
c
c---------------------------------------------------------------
      include 'param.f'
      include 'com7.f'
c
      dimension f(*),df(*),g(nn)
c---------------------------------------------------------------
c
      do 10 nc=1,nn
      df(nc)=0.5*(f(1)+f(nn)*cheb(nc,nn))
   10 continue
      do 11 kc=2,nn1
      do 11 nc=1,nn
      df(nc)=df(nc)+f(kc)*cheb(nc,kc)
   11 continue
      do 12 nc=1,nn
      df(nc)=anorm*df(nc)
   12 continue
c
      g(nn)=0.
      g(nn1)=float(nn1)*df(nn)
      do 20 n=nn2,1,-1
      g(n)=g(n+2)+float(2*n)*df(n+1)
   20 continue
c
      do 30 kc=1,nn
      df(kc)=0.5*g(nn)*cheb(nn,kc)
   30 continue
      do 31 nc=nn1,2,-1
      do 31 kc=1,nn
      df(kc)=df(kc)+g(nc)*cheb(nc,kc)
   31 continue
      do 32 kc=1,nn
      df(kc)=2.*anorm*(df(kc)+0.5*g(1))
   32 continue
c
      return
      end
