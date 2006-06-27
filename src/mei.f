      subroutine mei(enbp,enbt,apome,atome)
c
c  calculates total magnetic energy  = 1/2 Integral(B^2 dV)
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
     $         ql(lm,6)*(ql(lm,3)*qk(kc,1)*cabssq(b(lm,kc))
     $         +4.*cabssq(db(lm,kc)))
            rvat(kc)=rvat(kc)+ql(lm,6)*cabssq(aj(lm,kc))
   31    continue
   30 continue
c
      do 35 lm=nlaf,2,-1
         do 36 kc=1,nn
            rvb(kc)=rvb(kc)+ql(lm,6)*cabssq(aj(lm,kc))
            rvc(kc)=rvc(kc)+
     $         ql(lm,6)*(ql(lm,3)*qk(kc,1)*cabssq(b(lm,kc))+
     $         4.*cabssq(db(lm,kc)))
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
      enbp=0.
      enbt=0.
      atome=0.
      apome=0.
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
         enbp=enbp+rvap(nc)*qn(nc,3)
         enbt=enbt+rvat(nc)*qn(nc,3)
         atome=atome+rvb(nc)*qn(nc,3)
         apome=apome+rvc(nc)*qn(nc,3)
   50 continue
      enbp= 0.5*oekpm*anorm*enbp
      enbt= 0.5*oekpm*anorm*enbt
      atome=0.5*oekpm*anorm*atome
      apome=0.5*oekpm*anorm*apome
c
      return
      end
