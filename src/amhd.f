      subroutine amhd
c
c---------------------------------------------------------------
      include 'param.f'
      include 'com1.f'
      include 'com2.f'
      include 'com3.f'
      include 'com4.f'
      include 'com5.f'
      include 'com6.f'
      include 'com7.f'
      include 'com8.f'
c
      complex flms1(nlmpa),flms2(nlmpa),flms3(nlmpa),
     $flmw1(nlmpa),flmw2(nlmpa),flmw3(nlmpa),
     $flmb1(nlmpa),flmb2(nlmpa),flmb3(nlmpa),flmkb(nlma,nnp1),
     $ddj(nlma,nnp1),
     $flmks(nlma,nnp1),dflmk(nlma,nnp1),
     $ds(nlma,nnp1),dds(nlma,nnp1)
      complex carg2	
c
      dimension vr(nrp,ni),vt(nrp,ni),vp(nrp,ni),
     $dvrdr(nrp,ni),dvtdr(nrp,ni),dvpdr(nrp,ni),
     $cvr(nrp,ni),dvrdt(nrp,ni),dvrdp(nrp,ni),
     $dvtdp(nrp,ni),dvpdp(nrp,ni),
     $wnlr1(nrp,ni),wnlr2(nrp,ni),wnlr3(nrp,ni),
     $br(nrp,ni),bt(nrp,ni),bp(nrp,ni),
     $cbr(nrp,ni),cbt(nrp,ni),cbp(nrp,ni),
     $bnlr1(nrp,ni),bnlr2(nrp,ni),bnlr3(nrp,ni),
     $sr(nrp,ni),snlr1(nrp,ni),snlr2(nrp,ni),snlr3(nrp,ni),
     $frl(nn),fim(nn),frl2(nnx2),fim2(nnx2)
      dimension vrpoint(nnp1),vppoint(nnp1),vtpoint(nnp1)
c
      equivalence (vr,vrc),(vt,vtc),(vp,vpc),
     $(dvrdr,dvrdrc),
     $(dvtdr,dvtdrc),
     $(dvpdr,dvpdrc),
     $(cvr,cvrc),
     $(dvrdt,dvrdtc),
     $(dvrdp,dvrdpc),
     $(dvtdp,dvtdpc),
     $(dvpdp,dvpdpc),
     $(wnlr1,wnlc1),(wnlr2,wnlc2),(wnlr3,wnlc3),
     $(br,brc),(bt,btc),(bp,bpc),
     $(cbr,cbrc),(cbt,cbtc),(cbp,cbpc),
     $(bnlr1,bnlc1),(bnlr2,bnlc2),(bnlr3,bnlc3),
     $(dj,flmkb),(ddw,ddj),
     $(sr,sc),(snlr1,snlc1),(snlr2,snlc2),(snlr3,snlc3),
     $(ddz,ds,flmks),(dddw,dds,dflmk)

c
      cabssq(carg2)=real(carg2)**2+aimag(carg2)**2
c
c
      external stopiteration
c---------------------------------------------------------------
c
      isignal=30
      istop=0
c
      kcour=0
      kc0=0
      newdt=0
      kel=0
      toplum=0.
      botlum=0.
      vmax=0.
      courmax=.0
      couhmax=.0
      kcrmax=0
      kchmax=0
      reyn=0.
      bmax=0.
      alfrmax=0.
      alfhmax=0.
      karmax=0
      kahmax=0
      elsa=0.
      oz=0.
c
c **********************************************************************
c     integrate for nstep time steps.
c     nonlinear and coriolis terms via an explicit adams-bashforth method.
c     linear terms via an implicit method.
c **********************************************************************
c
c  ***  loop written explicitly:     do 1001 istep=1,nstep,2
c
      istep=-1
  100 istep=istep+2
c
      if(dtmin .eq. 0.) go to 1000
c
      call signal(isignal,stopiteration)
      if(istop.gt.0) istep=nstep
c
      do 1002 inel=1,2
c
      if(inel .eq. 1) then
         jnel=2
      else
         jnel=1
      endif
      kstep=kstep+1
      lcour=0               
      if((icour.gt.0).and.(mod(kstep+icour-1,icour).eq.0
     $ .or. kcour.gt.0)) then
       lcour=1
      dtr=100.*dtmax
      dth2=dtr*dtr
      endif
      time=time+dt
      if(istep .eq. nstep) kel=inel
      if((iprnt.eq.nprnt) .and. (kel.eq.2))
     $ call graphout(kc0,ngform)
      if(newdt .ne. 0) then
         call ludc
         newdt=0
      endif
      w2=-dt/(2.*dtold)
      w1=1.-w2
      dtold=dt
c
c ********************************************************************
c     advection of entropy and momentum,
c     coriolis, centrifugal, and lorentz forces,
c     induction of magnetic field,
c     joule, viscous and internal heating.
c ********************************************************************
c
c
      do 200 kc=1,nn
c
c    -legendre transform from (k,l,m) to (k,i,m)
c
      call legtf(kc)
c
c    -fourier transform from (k,i,m) to (k,i,j)
c
      call fourtf(sc,work,trigsf,ifaxf,1,nrp,nja,ni,1)
      call fourtf(vrc,work,trigsf,ifaxf,1,nrp,nja,ni,1)
      call fourtf(vtc,work,trigsf,ifaxf,1,nrp,nja,ni,1)
      call fourtf(vpc,work,trigsf,ifaxf,1,nrp,nja,ni,1)
      call fourtf(cvrc,work,trigsf,ifaxf,1,nrp,nja,ni,1)
      call fourtf(dvrdrc,work,trigsf,ifaxf,1,nrp,nja,ni,1)
      call fourtf(dvtdrc,work,trigsf,ifaxf,1,nrp,nja,ni,1)
      call fourtf(dvpdrc,work,trigsf,ifaxf,1,nrp,nja,ni,1)
      call fourtf(dvrdtc,work,trigsf,ifaxf,1,nrp,nja,ni,1)
      call fourtf(dvrdpc,work,trigsf,ifaxf,1,nrp,nja,ni,1)
      call fourtf(dvtdpc,work,trigsf,ifaxf,1,nrp,nja,ni,1)
      call fourtf(dvpdpc,work,trigsf,ifaxf,1,nrp,nja,ni,1)
      call fourtf(brc,work,trigsf,ifaxf,1,nrp,nja,ni,1)
      call fourtf(btc,work,trigsf,ifaxf,1,nrp,nja,ni,1)
      call fourtf(bpc,work,trigsf,ifaxf,1,nrp,nja,ni,1)
      call fourtf(cbrc,work,trigsf,ifaxf,1,nrp,nja,ni,1)
      call fourtf(cbtc,work,trigsf,ifaxf,1,nrp,nja,ni,1)
      call fourtf(cbpc,work,trigsf,ifaxf,1,nrp,nja,ni,1)
c
c    -courant condition check
c
c  urc modification: instead of the full Alfven velocity
c      a modified Alfven velocity is employed that takes
c      viscous and Joule damping into account. Different              
c      Courant factors are used for the fluid velocity and
c      the such modified Alfven velocity
c
      valri2=(0.5*(1.+opm)/delxr(kc))**2
      valhi2=(0.5*(1.+opm))**2/delxh2(kc)
c
      if(lcour .gt. 0) then  
         vrmax=0.
         vh2max=0.
c
         do 210 ic=1,ni
            do 211 jc=1,nja
               vflr=abs(vr(jc,ic)*qk(kc,1))
               valr2=oekpm*(br(jc,ic)*qk(kc,1))**2
               valr=valr2/sqrt(valr2+valri2)
               vrmax=max(vrmax,
     $            courfac*vflr
     $            +alffac*valr                      
     $            )
               vflh2=
     $            (vt(jc,ic)*vt(jc,ic)+vp(jc,ic)*vp(jc,ic))*
     $            qk(kc,1)*qi(ic,1)
               valh2=
     $            +(bt(jc,ic)*bt(jc,ic)+bp(jc,ic)*bp(jc,ic))*
     $            oekpm*qk(kc,1)*qi(ic,1)
               valh2m=valh2*valh2/(valh2+valhi2)
               vh2max=max(vh2max,
     $           courfac*courfac*vflh2
     $           +alffac*alffac*valh2m
     $           )
  211       continue
  210    continue
         if(vrmax .ne. 0.) dtr=min(dtr,delxr(kc)/vrmax)
         if(vh2max .ne. 0.) dth2=min(dth2,delxh2(kc)/vh2max)
      endif
c
c    -max velocity and magnetic field;
c
      if(kel .gt. 0) then
         vmax1=0.
         courmax1=.0
         couhmax1=.0
         bmax1=0.
         alfrmax1=0.
         alfhmax1=0.
         do 230 ic=1,ni
            do 231 jc=1,nja
               vflr=vr(jc,ic)*qk(kc,1)
               vflh2=
     $            (vt(jc,ic)*vt(jc,ic)+vp(jc,ic)*vp(jc,ic))*
     $            qk(kc,1)*qi(ic,1)
               vmax1=max(vmax1,vflr*vflr+vflh2)
               br2=(br(jc,ic)*qk(kc,1))**2
               bh2=
     $            +(bt(jc,ic)*bt(jc,ic)+bp(jc,ic)*bp(jc,ic))*
     $            qk(kc,1)*qi(ic,1)
               bmax1=max(bmax1,br2+bh2)
               courmax1=max(courmax1,abs(vflr))
               couhmax1=max(couhmax1,vflh2)
               valr2=br2*oekpm
               valh2=bh2*oekpm
               alfrmax1=max(alfrmax1,valr2*valr2/(valr2+valri2))
               alfhmax1=max(alfhmax1,valh2*valh2/(valh2+valhi2))
  231       continue
  230    continue
         vmax=max(vmax,vmax1)
         bmax=max(bmax,bmax1)
c
         courmax1=courmax1/delxr(kc)
         if(courmax1.gt.courmax) then
           courmax=courmax1
           kcrmax=kc
         endif
         couhmax1=sqrt(couhmax1/delxh2(kc))
         if(couhmax1.gt.couhmax) then
           couhmax=couhmax1
           kchmax=kc
         endif
         alfrmax1=sqrt(alfrmax1)/delxr(kc)
         if(alfrmax1.gt.alfrmax) then
           alfrmax=alfrmax1
           karmax=kc
         endif
         alfhmax1=sqrt(alfhmax1/delxh2(kc))
         if(alfhmax1.gt.alfhmax) then
           alfhmax=alfhmax1
           kahmax=kc
         endif
      endif
c
      if(nplog.ne.0) then 
        if(mod(kstep,nplog).eq.0) then
          vrpoint(kc)=vr(1,ni/2)*qk(kc,1)/vscale
          vppoint(kc)=vp(1,ni/2)*qk(kc,3)/(qi(ni/2,3)*vscale)
          vtpoint(kc)=vp(1,ni/4)*qk(kc,3)/(qi(ni/4,3)*vscale)
        endif
      endif
c
c urc :  graphics output      
c
      if((mod(kc-1,ngrad).eq.0) .and. (iprnt.eq.nprnt) .and.
     $ (kel.eq.2)) call graphout(kc,ngform)
c
c  urc:  movie output
c
      if(time/tscale.ge.tmovnext.and.imovct.le.iframes) then
        if(kc.eq.nn) then
         if(imovopt.ge.1000) then
          kcv=imovopt/1000
          call cmbcoeff(kcv)
         endif
         imovct=imovct+1
         tmovnext=tmovnext+tmovstep
        endif
        kvp=mod(imovopt,1000)/100
        if(mod(imovopt,10).ge.1) call moveout(kc)
        if(mod(imovopt,100).ge.10) call movaout(kc)
        if(mod(imovopt,1000).ge.100) call movmout(kc,kvp)
      endif
c
c  urc: decomposition of (br * u) at surface into poloidal
c       and toroidal parts
c
cc    if(kc.eq.1 .and. iprnt.eq. nprnt .and. kel.eq.2) then 
cc     do ic=1,ni
cc       do jc=1,nja
cc         bvp(jc,ic)=          vp(jc,ic)
cc         bvt(jc,ic)=          vt(jc,ic)
cccc       bvp(jc,ic)=br(jc,ic)*vp(jc,ic)
cccc       bvt(jc,ic)=br(jc,ic)*vt(jc,ic)
cc       enddo
cc     enddo
cc     call bvdecompose(bvp,bvt)
cc    endif
c
c    -quadratic products in physical space
c
      do  ic=1,ni
      do  jc=1,nja
      wnlr1(jc,ic)=0.            ! Inertia & Lorentz force, r-comp.
     $   -qk(kc,1)*(vr(jc,ic)*(dvrdr(jc,ic)-
     $   qk(kc,6)*vr(jc,ic))+qi(ic,1)*(vt(jc,ic)*
     $   (dvrdt(jc,ic)-r(kc)*vt(jc,ic))+vp(jc,ic)*
     $   (dvrdp(jc,ic)-r(kc)*vp(jc,ic))))
     $   +oekpm*qi(ic,1)*(cbt(jc,ic)*bp(jc,ic)-
     $   cbp(jc,ic)*bt(jc,ic))
c
      wnlr2(jc,ic)=0.            ! Inertia & Lorentz force, t-comp.
     $   +qk(kc,4)*qi(ic,1)*(vr(jc,ic)*
     $   (-dvtdr(jc,ic))+vt(jc,ic)*(qi(ic,2)*
     $   vt(jc,ic)+dvpdp(jc,ic)+dvrdr(jc,ic))+
     $   vp(jc,ic)*(qi(ic,2)*vp(jc,ic)-dvtdp(jc,ic))
     $   +oekpm*(cbp(jc,ic)*
     $   br(jc,ic)-cbr(jc,ic)*bp(jc,ic)))
c
      wnlr3(jc,ic)=0.            ! Inertia & Lorentz force, p-comp.
     $   +qk(kc,4)*qi(ic,1)*(vr(jc,ic)*
     $   (-dvpdr(jc,ic))-vt(jc,ic)*
     $   (dvtdp(jc,ic)+cvr(jc,ic))-
     $   vp(jc,ic)*dvpdp(jc,ic))
     $   +oekpm*qk(kc,4)*qi(ic,1)*(cbr(jc,ic)*
     $   bt(jc,ic)-cbt(jc,ic)*br(jc,ic))
      enddo
      enddo
c
      do  ic=1,ni
      do  jc=1,nja
      snlr1(jc,ic)=vr(jc,ic)*sr(jc,ic)
      snlr2(jc,ic)=qk(kc,1)*qi(ic,1)*(vt(jc,ic)*sr(jc,ic))
      snlr3(jc,ic)=qk(kc,1)*qi(ic,1)*(vp(jc,ic)*sr(jc,ic))
      enddo
      enddo
c
      do  ic=1,ni
      do  jc=1,nja
      bnlr1(jc,ic)=qi(ic,1)*(vt(jc,ic)*
     $   bp(jc,ic)-vp(jc,ic)*bt(jc,ic))
      bnlr2(jc,ic)=qk(kc,4)*qi(ic,1)*(vp(jc,ic)*br(jc,ic)-
     $   vr(jc,ic)*bp(jc,ic))
      bnlr3(jc,ic)=qk(kc,4)*qi(ic,1)*(vr(jc,ic)*bt(jc,ic)-
     $   vt(jc,ic)*br(jc,ic))
      enddo
      enddo
c
      do 208 jc=nja+1,nja+2
         do 209 ic=1,ni
            wnlr1(jc,ic)=0.
            wnlr2(jc,ic)=0.
            wnlr3(jc,ic)=0.
            snlr1(jc,ic)=0.
            snlr2(jc,ic)=0.
            snlr3(jc,ic)=0.
            bnlr1(jc,ic)=0.
            bnlr2(jc,ic)=0.
            bnlr3(jc,ic)=0.
  209    continue
  208 continue
c
c    -fourier transform from (k,i,j) to (k,i,m)
c
      call fourtf(wnlr1,work,trigsf,ifaxf,1,nrp,nja,ni,-1)
      call fourtf(wnlr2,work,trigsf,ifaxf,1,nrp,nja,ni,-1)
      call fourtf(wnlr3,work,trigsf,ifaxf,1,nrp,nja,ni,-1)
      call fourtf(snlr1,work,trigsf,ifaxf,1,nrp,nja,ni,-1)
      call fourtf(snlr2,work,trigsf,ifaxf,1,nrp,nja,ni,-1)
      call fourtf(snlr3,work,trigsf,ifaxf,1,nrp,nja,ni,-1)
      call fourtf(bnlr1,work,trigsf,ifaxf,1,nrp,nja,ni,-1)
      call fourtf(bnlr2,work,trigsf,ifaxf,1,nrp,nja,ni,-1)
      call fourtf(bnlr3,work,trigsf,ifaxf,1,nrp,nja,ni,-1)
c
c    -legendre transform from (k,i,m) to (k,l,m)
c
      do 212 lmp=1,nlmpa
         flmw1(lmp)=0.
         flmw2(lmp)=0.
         flmw3(lmp)=0.
         flms1(lmp)=0.
         flms2(lmp)=0.
         flms3(lmp)=0.
         flmb1(lmp)=0.
         flmb2(lmp)=0.
         flmb3(lmp)=0.
  212 continue
      do 213 ic=1,ni
         do 214 lmp=1,nlmpa
            mca=mcalmp(lmp)
            flmw1(lmp)=flmw1(lmp)+wnlc1(mca,ic)*aleg2(lmp,ic)
            flmw2(lmp)=flmw2(lmp)+wnlc2(mca,ic)*aleg2(lmp,ic)
            flmw3(lmp)=flmw3(lmp)+wnlc3(mca,ic)*aleg2(lmp,ic)
            flms1(lmp)=flms1(lmp)+snlc1(mca,ic)*aleg2(lmp,ic)
            flms2(lmp)=flms2(lmp)+snlc2(mca,ic)*aleg2(lmp,ic)
            flms3(lmp)=flms3(lmp)+snlc3(mca,ic)*aleg2(lmp,ic)
            flmb1(lmp)=flmb1(lmp)+bnlc1(mca,ic)*aleg2(lmp,ic)
            flmb2(lmp)=flmb2(lmp)+bnlc2(mca,ic)*aleg2(lmp,ic)
            flmb3(lmp)=flmb3(lmp)+bnlc3(mca,ic)*aleg2(lmp,ic)
  214    continue
  213 continue
c
      dwdt(1,kc,jnel)=real(flmw1(1))*qk(kc,1)
      dsdt(1,kc,jnel)=epsc0
      flmks(1,kc)=flms1(1)
      flmkb(1,kc)=0.
c
      do lm=2,nlma
         lma1=min(lm+1,nlma)
         lmp=lm+(mclm(lm)-1)/minc
         dwdt(lm,kc,jnel)=
     $      flmw1(lmp)*qk(kc,1)+
     $      2.*oek*qk(kc,3)*(2.*cmplx(0.,ql(lm,13))*dw(lm,kc)+
     $      ql(lm,14)*z(lma1,kc)-ql(lm,15)*z(lm-1,kc))
         dzdt(lm,kc,jnel)=
     $      (ql(lm,7)*flmw3(lmp-1)-ql(lm,8)*flmw3(lmp+1)-
     $      cmplx(0.0,ql(lm,10))*flmw2(lmp))+
     $      2.*oek*qk(kc,1)*
     $      (cmplx(0.0,ql(lm,13))*z(lm,kc)+
     $      ql(lm,16)*(2.*dw(lma1,kc)+
     $      ql(lm,17)*qk(kc,3)*w(lma1,kc))+
     $      ql(lm,18)*(2.*dw(lm-1,kc)-
     $      ql(lm,19)*qk(kc,3)*w(lm-1,kc)))
c
         dpdt(lm,kc,jnel)=
     $      (ql(lm,7)*flmw2(lmp-1)-ql(lm,8)*flmw2(lmp+1)+
     $      cmplx(0.0,ql(lm,10))*flmw3(lmp))+
     $      2.*oek*qk(kc,1)*(-cmplx(0.0,ql(lm,13))
     $      *(2.*dw(lm,kc)+ql(lm,20)*qk(kc,3)*w(lm,kc))+
     $      ql(lm,16)*z(lma1,kc)+
     $      ql(lm,18)*z(lm-1,kc))
         dsdt(lm,kc,jnel)=-(ql(lm,7)*flms2(lmp-1)-ql(lm,8)*flms2(lmp+1)+
     $    cmplx(0.0,ql(lm,10))*flms3(lmp))
         flmks(lm,kc)=flms1(lmp)
         dbdt(lm,kc,jnel)=
     $      (ql(lm,7)*flmb3(lmp-1)-ql(lm,8)*flmb3(lmp+1)-
     $      cmplx(0.0,ql(lm,10))*flmb2(lmp))
         djdt(lm,kc,jnel)=ql(lm,3)*qk(kc,4)*flmb1(lmp)
         flmkb(lm,kc)=(ql(lm,7)*flmb2(lmp-1)-
     $      ql(lm,8)*flmb2(lmp+1)+
     $      cmplx(0.0,ql(lm,10))*flmb3(lmp))*(r(kc)*r(kc))
      enddo    
c
  200 continue
c
c    -radial derivatives of nl terms
c
      call chebtf(lot,ns2,lot,nnp1,nps2,
     $flmks,flmks,flmks,wsave,work,
     $work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
      do 249 lm=1,nlma
         dflmk(lm,nn)=0.
         dflmk(lm,nn1)=float(nn1)*flmks(lm,nn)
  249 continue
      do 240 n=nn2,1,-1
         do 240 lm=1,nlma
           dflmk(lm,n)=dflmk(lm,n+2)+float(2*n)*flmks(lm,n+1)
  240 continue
      call chebtf(lot,ns2,lot,nnp1,nps2,
     $dflmk,dflmk,dflmk,wsave,work,
     $work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
      do 242 kc=1,nn
         do 243 lm=1,nlma
            dsdt(lm,kc,jnel)=dsdt(lm,kc,jnel)-
     $         2.*qk(kc,1)*dflmk(lm,kc)
  243    continue
  242 continue
      call chebtf(lot,ns2,lot,nnp1,nps2,
     $flmkb,flmkb,flmkb,wsave,work,
     $work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
      do 219 lm=1,nlma
         dflmk(lm,nn)=0.
         dflmk(lm,nn1)=float(nn1)*flmkb(lm,nn)
  219 continue
      do 220 n=nn2,1,-1
         do 220 lm=1,nlma
           dflmk(lm,n)=dflmk(lm,n+2)+float(2*n)*flmkb(lm,n+1)
  220 continue
      call chebtf(lot,ns2,lot,nnp1,nps2,
     $dflmk,dflmk,dflmk,wsave,work,
     $work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
      do 222 kc=1,nn
         do 223 lm=2,nlma
            djdt(lm,kc,jnel)=djdt(lm,kc,jnel)+
     $         2.*qk(kc,1)*dflmk(lm,kc)
  223    continue
  222 continue
c
c    -time-step check and change if needed
c
      if(lcour .gt. 0) then   
         dth=sqrt(dth2)
         call dtchck(kstep,newdt,dt,dtnew,
     $      dtmin,dtmax,dtr,dth,ifirst,kcour)
      else
         dtnew=dt
      endif
c
      w2new=-dtnew/(2.*dt)
      coex=(1.-alpha)/w2new
c
c **************************************
c     update magnetic field
c **************************************
c
      if(ifbfrz) then
c
      call chebtf(lot,ns2,lot,nnp1,nps2,
     $b,b,b,wsave,work,
     $work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
      call chebtf(lot,ns2,lot,nnp1,nps2,
     $aj,aj,aj,wsave,work,
     $work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
c
      else
c
      do 730 lm=2,nlaf
         l=lm-1
         do kc=1,nn
            frl(kc)=
     $         (w1*real(dbdt(lm,kc,jnel))+
     $         w2*real(dbdt(lm,kc,inel)))+
     $         oodt*ql(lm,3)*qk(kc,1)*real(b(lm,kc))
            frl2(kc)=
     $         (w1*real(djdt(lm,kc,jnel))+
     $         w2*real(djdt(lm,kc,inel)))+
     $         oodt*ql(lm,3)*qk(kc,1)*real(aj(lm,kc))
         enddo
         frl(1)=0.
         frl(nn)=0.
         if(kbotb .eq. 2) frl(nn1)=.0
         frl2(1)=0.
         frl2(nn)=0.
         if(l.eq.2.and.imagcon.gt.0.and.imagcon.ne.12) then 
          frl2(nn)=bpeakbot
          frl2( 1)=bpeaktop
         endif
         if(l.eq.1.and.imagcon.eq.12) then
           frl2(nn)=bpeakbot
           frl2( 1)=bpeaktop
         endif
         if(l.eq.1.and.imagcon.lt. 0) frl(nn)= bpeakbot
         call sgesl(bmat(1,1,l),nn,nn,ib(1,l),frl,0)
         call sgesl(ajmat(1,1,l),nn,nn,ij(1,l),frl2,0)
c                
         do nc=1,nnaf
            b(lm,nc)=frl(nc)*bscl
            aj(lm,nc)=frl2(nc)*bscl
         enddo
  730 continue
      do 731 lm=nlaf+1,nlma
         l=nint(ql(lm,4))
         do kc=1,nn
            frl(kc)=
     $         (w1*real(dbdt(lm,kc,jnel))+
     $         w2*real(dbdt(lm,kc,inel)))+
     $         oodt*ql(lm,3)*qk(kc,1)*real(b(lm,kc))
            fim(kc)=
     $         (w1*aimag(dbdt(lm,kc,jnel))+
     $         w2*aimag(dbdt(lm,kc,inel)))+
     $         oodt*ql(lm,3)*qk(kc,1)*aimag(b(lm,kc))
            frl2(kc)=
     $         (w1*real(djdt(lm,kc,jnel))+
     $         w2*real(djdt(lm,kc,inel)))+
     $         oodt*ql(lm,3)*qk(kc,1)*real(aj(lm,kc))
            fim2(kc)=
     $         (w1*aimag(djdt(lm,kc,jnel))+
     $         w2*aimag(djdt(lm,kc,inel)))+
     $         oodt*ql(lm,3)*qk(kc,1)*aimag(aj(lm,kc))
         enddo
         frl(1)=0.
         fim(1)=0.
         frl(nn)=0.
         fim(nn)=0.
         if(kbotb .eq. 2) then
            frl(nn1)=0.
            fim(nn1)=0.
         endif
         frl2(1)=0.
         fim2(1)=0.
         frl2(nn)=0.
         fim2(nn)=0.
         call sgesl(bmat(1,1,l),nn,nn,ib(1,l),frl,0)
         call sgesl(bmat(1,1,l),nn,nn,ib(1,l),fim,0)
         call sgesl(ajmat(1,1,l),nn,nn,ij(1,l),frl2,0)
         call sgesl(ajmat(1,1,l),nn,nn,ij(1,l),fim2,0)
         do nc=1,nnaf
            b(lm,nc)=cmplx(frl(nc),fim(nc))*bscl
            aj(lm,nc)=cmplx(frl2(nc),fim2(nc))*bscl
         enddo
  731 continue
      if(nnaf .lt. nn) then
         do nc=nnaf+1,nn
            do lm=2,nlma
               b(lm,nc)=0.
               aj(lm,nc)=0.
            enddo
         enddo
      endif
c
      endif
c
c    -radial derivs of b and j
c
      do lm=1,nlma
         db(lm,nn)=0.
         db(lm,nn1)=float(nn1)*b(lm,nn)
         dj(lm,nn)=0.
         dj(lm,nn1)=float(nn1)*aj(lm,nn)
         ddb(lm,nn)=0.
         ddb(lm,nn1)=0.
         ddj(lm,nn)=0.
         ddj(lm,nn1)=0.
      enddo
      do n=nn2,1,-1
         do lm=1,nlma
            db(lm,n)=db(lm,n+2)+float(2*n)*b(lm,n+1)
            dj(lm,n)=dj(lm,n+2)+float(2*n)*aj(lm,n+1)
            ddb(lm,n)=ddb(lm,n+2)+float(2*n)*db(lm,n+1)
            ddj(lm,n)=ddj(lm,n+2)+float(2*n)*dj(lm,n+1)
         enddo
      enddo
c
c    -chebyshev transform b, j and derivs from (n,l,m) to (k,l,m)
c
      call chebtf(lot,ns2,lot,nnp1,nps2,
     $b,b,b,wsave,work,
     $work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
      call chebtf(lot,ns2,lot,nnp1,nps2,
     $db,db,db,wsave,work,
     $work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
      call chebtf(lot,ns2,lot,nnp1,nps2,
     $ddb,ddb,ddb,wsave,work,
     $work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
      call chebtf(lot,ns2,lot,nnp1,nps2,
     $aj,aj,aj,wsave,work,
     $work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
      call chebtf(lot,ns2,lot,nnp1,nps2,
     $dj,dj,dj,wsave,work,
     $work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
      call chebtf(lot,ns2,lot,nnp1,nps2,
     $ddj,ddj,ddj,wsave,work,
     $work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
      do 755 kc=1,nnp1
         b(1,kc)=0.
         db(1,kc)=0.
         ddb(1,kc)=0.
         aj(1,kc)=0.
         dj(1,kc)=0.
         ddj(1,kc)=0.
  755 continue
c
c    -explicit parts of the linear terms in the b and j equations
c
      if(alpha .lt. 1.) then
         do 750 kc=1,nn
            do 751 lm=2,nlma
               dbdt(lm,kc,jnel)=dbdt(lm,kc,jnel)+coex*
     $            ql(lm,11)*opm*qk(kc,1)*(4.*ddb(lm,kc)-
     $            ql(lm,3)*qk(kc,1)*b(lm,kc))
               djdt(lm,kc,jnel)=djdt(lm,kc,jnel)+coex*
     $            ql(lm,11)*opm*qk(kc,1)*(4.*ddj(lm,kc)-
     $            ql(lm,3)*qk(kc,1)*aj(lm,kc))
  751       continue
  750    continue
      endif
c
c **************************************
c     update entropy
c **************************************
c
      if(ifsfrz) then
c
      call chebtf(lot,ns2,lot,nnp1,nps2,
     $s,s,s,wsave,work,
     $work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
c
      else
c
      do kc=1,nn
         frl(kc)=real(s(1,kc))*oodt+
     $      (w1*real(dsdt(1,kc,jnel))+
     $      w2*real(dsdt(1,kc,inel)))
      enddo
      frl(1)=real(tops(0,0))
      frl(nn)=real(bots(0,0))
      call sgesl(s0mat,nn,nn,is0,frl,0)
      do nc=1,nnaf
         s(1,nc)=frl(nc)*sscl
      enddo
      do 530 lm=2,nlaf
         l=lm-1
         do kc=1,nn
            frl(kc)=real(s(lm,kc))*oodt+
     $         (w1*real(dsdt(lm,kc,jnel))+
     $         w2*real(dsdt(lm,kc,inel)))
         enddo
         frl(1)=real(tops(l,0))
         frl(nn)=real(bots(l,0))
         call sgesl(smat(1,1,l),nn,nn,is(1,l),frl,0)
         do nc=1,nnaf
            s(lm,nc)=frl(nc)*sscl
         enddo
  530 continue
      do 531 lm=nlaf+1,nlma
         l=nint(ql(lm,4))
         m=mclm(lm)-1
         do kc=1,nn
           frl(kc)=real(s(lm,kc))*oodt+
     $        (w1*real(dsdt(lm,kc,jnel))+
     $        w2*real(dsdt(lm,kc,inel)))
           fim(kc)=aimag(s(lm,kc))*oodt+
     $        (w1*aimag(dsdt(lm,kc,jnel))+
     $        w2*aimag(dsdt(lm,kc,inel)))
         enddo
         frl(1)=real(tops(l,m))
         frl(nn)=real(bots(l,m))
         fim(1)=aimag(tops(l,m))
         fim(nn)=aimag(bots(l,m))
         call sgesl(smat(1,1,l),nn,nn,is(1,l),frl,0)
         call sgesl(smat(1,1,l),nn,nn,is(1,l),fim,0)
         do nc=1,nnaf
            s(lm,nc)=cmplx(frl(nc),fim(nc))*sscl
         enddo
  531 continue
      if(nnaf .lt. nn) then
         do nc=nnaf+1,nn
            do lm=1,nlma
               s(lm,nc)=0.
            enddo
         enddo
      endif
c
      endif
c
c *** radial derivs of s
c
      do lm=1,nlma
         ds(lm,nn)=0.
         ds(lm,nn1)=float(nn1)*s(lm,nn)
         dds(lm,nn)=0.
         dds(lm,nn1)=0.
      enddo
      do n=nn2,1,-1
         do lm=1,nlma
            ds(lm,n)=ds(lm,n+2)+float(2*n)*s(lm,n+1)
            dds(lm,n)=dds(lm,n+2)+float(2*n)*ds(lm,n+1)
         enddo
      enddo
c
c *** chebyshev transform s and derivs from (n,l,m) to (k,l,m)
c
      call chebtf(lot,ns2,lot,nnp1,nps2,
     $s,s,s,wsave,work,
     $work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
      call chebtf(lot,ns2,lot,nnp1,nps2,
     $ds,ds,ds,wsave,work,
     $work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
      call chebtf(lot,ns2,lot,nnp1,nps2,
     $dds,dds,dds,wsave,work,
     $work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
c
c    -explicit parts of the linear terms in the s equation
c
      if(alpha .lt. 1.) then
         do 545 kc=1,nn
            do 546 lm=1,nlma
               dsdt(lm,kc,jnel)=dsdt(lm,kc,jnel)+coex*
     $            opr*(4.*dds(lm,kc)+
     $            qk(kc,2)*ds(lm,kc)-ql(lm,3)*qk(kc,1)*s(lm,kc))
  546       continue
  545    continue
      endif
c
c *********************************************************
c    -diagnostics
c *********************************************************
c
c urc & plo modified
      if(mod(kstep,nlogstep).eq.0.or.kel.eq.2) then
         botlum=-4.*pi*y00*r(nn)**2*opr*
     $      (2.*real(ds(1,nn)))
         toplum=-4.*pi*y00*r(1)**2*opr*
     $      (2.*real(ds(1,1)))
      endif
c
c ***************************************
c     update velocity and pressure
c ***************************************
c
      if(ifvfrz) then
c
      call chebtf(lot,ns2,lot,nnp1,nps2,
     $w,w,w,wsave,work,
     $work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
      call chebtf(lot,ns2,lot,nnp1,nps2,
     $z,z,z,wsave,work,
     $work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
      call chebtf(lot,ns2,lot,nnp1,nps2,
     $p,p,p,wsave,work,
     $work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
c
      else
c
      do 630 lm=2,nlaf
         l=lm-1
         do kc=1,nn
            frl(kc)=oodt*ql(lm,3)*qk(kc,1)*real(z(lm,kc))+
     $         (w1*real(dzdt(lm,kc,jnel))+
     $         w2*real(dzdt(lm,kc,inel)))
            frl2(kc)=(oodt*ql(lm,3)*qk(kc,1)*real(w(lm,kc))+
     $         (rapr*alpha*grav(kc)*real(s(lm,kc))+
     $         (w1*real(dwdt(lm,kc,jnel))+
     $         w2*real(dwdt(lm,kc,inel)))))
            frl2(kc+nn)=-oodt*ql(lm,3)*qk(kc,1)*
     $         (2.*real(dw(lm,kc)))+
     $         (w1*real(dpdt(lm,kc,jnel))+
     $         w2*real(dpdt(lm,kc,inel)))
         enddo
         frl(1)=0.
         frl(nn)=0.
         frl2(1)=0.
         frl2(nn)=0.
         frl2(nnp1)=0.
         frl2(nnx2)=0.
         call sgesl(zmat(1,1,l),nn,nn,iz(1,l),frl,0)
         call sgesl(wpmat(1,1,l),nnx2,nnx2,iwp(1,l),frl2,0)
         do nc=1,nnaf
            z(lm,nc)=frl(nc)*zscl
            w(lm,nc)=frl2(nc)*wscl
            p(lm,nc)=frl2(nc+nn)*pscl
         enddo
  630 continue
      do 631 lm=nlaf+1,nlma
         l=nint(ql(lm,4))
         do kc=1,nn
            frl(kc)=oodt*ql(lm,3)*qk(kc,1)*real(z(lm,kc))+
     $         (w1*real(dzdt(lm,kc,jnel))+
     $         w2*real(dzdt(lm,kc,inel)))
            frl2(kc)=(oodt*ql(lm,3)*qk(kc,1)*real(w(lm,kc))+
     $         (rapr*alpha*grav(kc)*real(s(lm,kc))+
     $         (w1*real(dwdt(lm,kc,jnel))+
     $         w2*real(dwdt(lm,kc,inel)))))
            frl2(kc+nn)=-oodt*ql(lm,3)*qk(kc,1)*
     $         (2.*real(dw(lm,kc)))+
     $         (w1*real(dpdt(lm,kc,jnel))+
     $         w2*real(dpdt(lm,kc,inel)))
            fim(kc)=oodt*ql(lm,3)*qk(kc,1)*aimag(z(lm,kc))+
     $         (w1*aimag(dzdt(lm,kc,jnel))+
     $         w2*aimag(dzdt(lm,kc,inel)))
            fim2(kc)=(oodt*ql(lm,3)*qk(kc,1)*aimag(w(lm,kc))+
     $         (rapr*alpha*grav(kc)*aimag(s(lm,kc))+
     $         (w1*aimag(dwdt(lm,kc,jnel))+
     $         w2*aimag(dwdt(lm,kc,inel)))))
            fim2(kc+nn)=-oodt*ql(lm,3)*qk(kc,1)*
     $         (2.*aimag(dw(lm,kc)))+
     $         (w1*aimag(dpdt(lm,kc,jnel))+
     $         w2*aimag(dpdt(lm,kc,inel)))
         enddo
         frl(1)=0.
         frl(nn)=0.
         frl2(1)=0.
         frl2(nn)=0.
         frl2(nnp1)=0.
         frl2(nnx2)=0.
         fim(1)=0.
         fim(nn)=0.
         fim2(1)=0.
         fim2(nn)=0.
         fim2(nnp1)=0.
         fim2(nnx2)=0.
         call sgesl(zmat(1,1,l),nn,nn,iz(1,l),frl,0)
         call sgesl(zmat(1,1,l),nn,nn,iz(1,l),fim,0)
         call sgesl(wpmat(1,1,l),nnx2,nnx2,iwp(1,l),frl2,0)
         call sgesl(wpmat(1,1,l),nnx2,nnx2,iwp(1,l),fim2,0)
         do nc=1,nnaf
            z(lm,nc)=cmplx(frl(nc),fim(nc))*zscl
            w(lm,nc)=cmplx(frl2(nc),fim2(nc))*wscl
            p(lm,nc)=cmplx(frl2(nc+nn),fim2(nc+nn))*pscl
         enddo
  631 continue
      if(nnaf .lt. nn) then
         do nc=nnaf+1,nn
            do lm=2,nlma
               z(lm,nc)=0.
               w(lm,nc)=0.
               p(lm,nc)=0.
            enddo
         enddo
      endif
c
      endif
c
c    -radial derivs of w and z
c
      do lm=1,nlma
         dw(lm,nn)=0.
         dw(lm,nn1)=float(nn1)*w(lm,nn)
         ddw(lm,nn)=0.
         ddw(lm,nn1)=0.
         dddw(lm,nn)=0.
         dddw(lm,nn1)=0.
         dz(lm,nn)=0.
         dz(lm,nn1)=float(nn1)*z(lm,nn)
         ddz(lm,nn)=0.
         ddz(lm,nn1)=0.
      end do
      do n=nn2,1,-1
         do lm=1,nlma
            dw(lm,n)=dw(lm,n+2)+float(2*n)*w(lm,n+1)
            ddw(lm,n)=ddw(lm,n+2)+float(2*n)*dw(lm,n+1)
            dddw(lm,n)=dddw(lm,n+2)+float(2*n)*ddw(lm,n+1)
            dz(lm,n)=dz(lm,n+2)+float(2*n)*z(lm,n+1)
            ddz(lm,n)=ddz(lm,n+2)+float(2*n)*dz(lm,n+1)
         end do
      end do
c
c    -chebyshev transform w, z and derivs from (n,l,m) to (k,l,m)
c
      call chebtf(lot,ns2,lot,nnp1,nps2,
     $w,w,w,wsave,work,
     $work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
      call chebtf(lot,ns2,lot,nnp1,nps2,
     $dw,dw,dw,wsave,work,
     $work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
      call chebtf(lot,ns2,lot,nnp1,nps2,
     $ddw,ddw,ddw,wsave,work,
     $work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
      call chebtf(lot,ns2,lot,nnp1,nps2,
     $dddw,dddw,dddw,wsave,work,
     $work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
      call chebtf(lot,ns2,lot,nnp1,nps2,
     $z,z,z,wsave,work,
     $work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
      call chebtf(lot,ns2,lot,nnp1,nps2,
     $dz,dz,dz,wsave,work,
     $work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
      call chebtf(lot,ns2,lot,nnp1,nps2,
     $ddz,ddz,ddz,wsave,work,
     $work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
      do kc=1,nnp1
         w(1,kc)=0.
         dw(1,kc)=0.
         ddw(1,kc)=0.
         dddw(1,kc)=0.
         z(1,kc)=0.
         dz(1,kc)=0.
         ddz(1,kc)=0.
      end do
c
c    -update the l=0 part of pressure
c
      if(.not. ifvfrz) then
c
      do 652 kc=1,nn
         frl(kc)=real(dwdt(1,kc,jnel))+
     $      rapr*grav(kc)*real(s(1,kc))+
     $      oek*p00co*qk(kc,3)*real(z(2,kc))
  652 continue
      frl(nps2)=0.
c
      call sgesl(p0mat,nn,nn,ip0,frl,0)
      do 653 nc=1,nnaf
         p(1,nc)=frl(nc)
  653 continue
      if(nnaf .lt. nn) then
         do 654 nc=nnaf+1,nn
            p(1,nc)=0.
  654    continue
      endif
c
      endif
c
c    -radial derivs of p
c
      do lm=1,nlma
         dp(lm,nn)=0.
         dp(lm,nn1)=float(nn1)*p(lm,nn)
      end do
      do n=nn2,1,-1
         do lm=1,nlma
            dp(lm,n)=dp(lm,n+2)+float(2*n)*p(lm,n+1)
         end do
      end do
c
c    -chebyshev transform p and derivs from (n,l,m) to (k,l,m)
c
      call chebtf(lot,ns2,lot,nnp1,nps2,
     $p,p,p,wsave,work,
     $work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
      call chebtf(lot,ns2,lot,nnp1,nps2,
     $dp,dp,dp,wsave,work,
     $work(1,2),work(1,2),work(1,2),trigsc,ifaxc,k2k)
c
c    -explicit parts of the linear terms in the w and z equations
c
      if(alpha .lt. 1.) then
         do kc=1,nn
            do lm=2,nlma
               dwdt(lm,kc,jnel)=dwdt(lm,kc,jnel)+coex*
     $            (-2.*dp(lm,kc)+
     $            rapr*grav(kc)*s(lm,kc)+
     $            ql(lm,12)*qk(kc,1)*
     $            (4.*ddw(lm,kc)-
     $            (ql(lm,3)*qk(kc,1))*w(lm,kc)))
               dpdt(lm,kc,jnel)=dpdt(lm,kc,jnel)+coex*
     $            (ql(lm,3)*qk(kc,1)*p(lm,kc)+
     $            ql(lm,12)*qk(kc,1)*
     $            (-8.*dddw(lm,kc)+
     $            2.*(ql(lm,3)*qk(kc,1))*dw(lm,kc)-
     $            ql(lm,3)*qk(kc,5)*w(lm,kc)))
               dzdt(lm,kc,jnel)=dzdt(lm,kc,jnel)+coex*
     $            ql(lm,12)*qk(kc,1)*
     $            (4.*ddz(lm,kc)-
     $            (ql(lm,3)*qk(kc,1))*z(lm,kc))
            enddo
         enddo
      endif
c
c    -energies
c
c urc: distinguish between energy in toroidal and poloidal fields
      call kei(envp,envt,adrke,amcke)
      env=envp+envt
      call mei(enbp,enbt,apome,atome)
      enb=enbp+enbt
      ent=env+enb
c
c plo: extended l-file output
c
      topnuss=toplum/alum0
      botnuss=botlum/alum0
      vmean=sqrt(2.*env/ocorevol)/vscale
      bmean=sqrt(2.*enb/(ocorevol*oekpm))

c  plo: output dipole tilt, long, cmb axial and full dipole (rms) fields
      tiltdipole=0.0
      phidipole=0.0
      if(minc.eq.1) 
     $  tiltdipole=atan2(abs(b(lmax+2,1)),real(b(2,1)))
      if(minc.eq.1 .and. real(b(lmax+2,1)).ne.0) 	
     $   phidipole=atan2(-1.*aimag(b(lmax+2,1)),real(b(lmax+2,1))) 
      surface=4.*pi*radtop*radtop
      dipolax=sqrt(2.*ql(2,6)*(ql(2,3)*qk(1,1)*cabssq(b(2,1))
     $  +4.*cabssq(db(2,1))) /surface)
      dipole=dipolax
      if(minc.eq.1) then 
        lm=lmax+2
        dipole=sqrt(dipolax**2 +
     $   2.*ql(lm,6)*(ql(lm,3)*qk(1,1)*cabssq(b(lm,1))
     $  +4.*cabssq(db(lm,1))) /surface)
      endif

c  print to l-file
      if(mod(kstep,nlogstep).eq.0) write(15,'(f10.6,8f9.2,5f9.5,3f9.2)')
     $time/tscale,env/escale,envp/escale,enb/escale,enbp/escale,
     $adrke/escale,amcke/escale,apome/escale,atome/escale,
     $topnuss,botnuss,bmean,dipole,dipolax,
     $tiltdipole*180./pi,phidipole*180./pi,vmean
      if(nplog.gt.0) then
       if(mod(kstep,nplog).eq.0) write(17,'(f10.6,6(1x,f9.3))')
     $ time/tscale,
     $ vrpoint(nn/2),vppoint(nn/2),
     $ vrpoint(2*nn/3),vppoint(2*nn/3)
      endif
c
c    -change time-step
c
      dt=dtnew
      oodt=1./dt
c
 1002 continue
c
c ***** end of explicit loop:   1001 continue
      if(istep.lt.nstep) go to 100
c
 1000 continue
c
c
      vmax=sqrt(vmax)/vscale
      vmean=sqrt(2.*env/ocorevol)/vscale
c urc+2
      bmax=sqrt(bmax)
      bmean=sqrt(2.*enb/(ocorevol*oekpm))
      dth=sqrt(dth2)
      write(6,900) kstep,time/tscale
  900 format(/4x,"****",i6,1x,"steps",3x,
     $"time=",3x,f10.6," (visc.diff.time) ****")
      write(6,901) dt/tscale,dtr/tscale,dth/tscale
  901 format(/,2x,"dt =",f10.8,3x,"dtrmin =",f10.8,3x,
     $"dthmin =",f10.8)
      write(6,911) courmax*tscale,kcrmax,couhmax*tscale,
     $ kchmax,alfrmax*tscale,karmax,alfhmax*tscale,kahmax 
  911 format(2x,"cour= ",f7.0,i4,"  couh= ",f7.0,i4,
     $ "  alfr= ",f7.0,i4,"  alfh= ",f7.0,i4)
      write(6,902) ent/escale,env/escale,enb/escale
  902 format(2x,"ent = ",1pe10.3,1x,"env =",1pe10.3,1x,
     $"enb =",1pe10.3)
      write(6,907) vmax,vmean
  907 format(2x,"max/mean velocity =",f9.3,1x,f9.3)
      write(6,904) bmax,bmean
  904 format(2x,"max/mean field =",f9.4,1x,f9.4)
c     diflum=toplum-botlum
c     write(6,905) botlum,toplum,diflum
c 905 format(2x,"botlum =",1pe9.2,3x,"toplum =",1pe9.2,
c    $3x,"top-bot =",1pe9.2," ergs/s")
c urc + 7 
      topnusselt=toplum/alum0
      botnusselt=botlum/alum0
      write(6,909) topnusselt,botnusselt
  909 format(2x,"nusselt number top/bot =",2(1X,f8.5))
c
c    -print w, z, s, p, b, aj
c
      call prnt
c
c    -store current solution
c     note, this gets overwritten until istor increases
c
      if(ngform.gt.-1) call stor
      call spectrum(0)
      call spectrum(1)
c
c    -stop if dt too small
c
      if(dtmin .eq. 0.) stop '23'
c
      return
      end
