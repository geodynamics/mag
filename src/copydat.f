      subroutine copydat(dat,dato,nlma,nlmao,lmax,lmaxo,mmax,mmaxo,
     $  minc,minco,nr,nrold,nr1,nr32,r,rold)
c
c  copy data into correct locations (old grid > new grid)
c
c  called in mapdata
c
      logical ifold
      complex dat(nlma,nr1),dato(nlmao,*)
      dimension r(nr),rold(nr32)
      dimension fcmin(129)
      integer kcmin(129)
c
      nn=nr1-1
      if(nrold.ne.nr) then
        kcmin(1)=1
        fcmin(1)=1.
        kcmin(nn)=nrold+nn-nr-1
        fcmin(nn)=0.
        do kc=2,nn-1
          kco=kcmin(kc-1)
   10     if(r(kc).gt.rold(kco+1)) then
           kcmin(kc)=kco
           fcmin(kc)=(rold(kco+1)-r(kc))/(rold(kco+1)-rold(kco))
          else
           kco=kco+1
           if(kco.gt.nrold) stop '55'
           go to 10
          endif
        enddo
      endif
c
      do m=0,mmax,minc
        mm=m*(m-minc)/(2*minc)
        if(m.gt.mmaxo .or. mod(m,minco).ne.0) then
          ifold=.false.
        else
          ifold=.true.
          mmo=m*(m-minco)/(2*minco)
        endif
c
        do l=m,lmax
          lm=m*(lmax+1)/minc-mm+l-m+1
          if( ifold .and. l.le.lmaxo ) then
            lmo=m*(lmaxo+1)/minco-mmo+l-m+1
            if(nr.eq.nrold) then
             do kc=1,nn
              dat(lm,kc)=dato(lmo,kc)
             enddo
            else
             do kc=1,nn
              kco=kcmin(kc)
              dat(lm,kc)=fcmin(kc)*dato(lmo,kco)
     $           +(1.-fcmin(kc))*dato(lmo,kco+1)
             enddo
             if(nr.gt.nn) dat(lm,nr)=dato(lmo,nrold)
            endif
          else
            do kc=1,nr
              dat(lm,kc)=(0.0,0.0)
            enddo
          endif
        enddo
      enddo
c
      return
      end
