      subroutine kei(envp,envt,adrke,amcke)
c
c  calculates total kinetic energy  = 1/2 Integral (v^2 dV)
c
c---------------------------------------------------------------
      include 'param.f'
      include 'com1.f'
      include 'com4.f'
      include 'com5.f'
c
      complex c
      cabssq(c)=real(c)**2+aimag(c)**2
c---------------------------------------------------------------
c
      do 20 kc=1,nn
         rvap(kc)=0.
         rvat(kc)=0.
         rvb(kc)=0.
         rvc(kc)=0.
   20 continue
c
      do 30 lm=nlma,2,-1
         do 31 kc=1,nn
            rvap(kc)=rvap(kc)+
     $         ql(lm,6)*(ql(lm,3)*qk(kc,1)*cabssq(w(lm,kc))
     $         +4.*cabssq(dw(lm,kc)))
            rvat(kc)=rvat(kc)+ql(lm,6)*cabssq(z(lm,kc))
   31    continue
   30 continue
c
      do 35 lm=nlaf,2,-1
         do 36 kc=1,nn
            rvb(kc)=rvb(kc)+ql(lm,6)*cabssq(z(lm,kc))
            rvc(kc)=rvc(kc)+
     $         ql(lm,6)*(ql(lm,3)*qk(kc,1)*cabssq(w(lm,kc))+
     $         4.*cabssq(dw(lm,kc)))
   36    continue
   35 continue
c
      call chebtf(1,ns2,1,nn1,ns2,rvap,rvap,rvap,wsave,
     $work,work(2,1),work(2,1),work(2,1),trigsc,ifaxc,k2k)
      call chebtf(1,ns2,1,nn1,ns2,rvat,rvat,rvat,wsave,
     $work,work(2,1),work(2,1),work(2,1),trigsc,ifaxc,k2k)
      call chebtf(1,ns2,1,nn1,ns2,rvb,rvb,rvb,wsave,
     $work,work(2,1),work(2,1),work(2,1),trigsc,ifaxc,k2k)
      call chebtf(1,ns2,1,nn1,ns2,rvc,rvc,rvc,wsave,
     $work,work(2,1),work(2,1),work(2,1),trigsc,ifaxc,k2k)
c
      envp=0.
      envt=0.
      adrke=0.
      amcke=0.
      rvap(1)=0.5*rvap(1)
      rvap(nn)=0.5*rvap(nn)
      rvat(1)=0.5*rvat(1)
      rvat(nn)=0.5*rvat(nn)
      rvb(1)=0.5*rvb(1)
      rvb(nn)=0.5*rvb(nn)
      rvc(1)=0.5*rvc(1)
      rvc(nn)=0.5*rvc(nn)
      do 50 ncb=1,nn,2
         nc=nnp1-ncb
         envp=envp+rvap(nc)*qn(nc,3)
         envt=envt+rvat(nc)*qn(nc,3)
         adrke=adrke+rvb(nc)*qn(nc,3)
         amcke=amcke+rvc(nc)*qn(nc,3)
   50 continue
      envp=0.5*anorm*envp
      envt=0.5*anorm*envt
      adrke=0.5*anorm*adrke
      amcke=0.5*anorm*amcke
c
      return
      end
