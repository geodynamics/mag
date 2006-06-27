      subroutine legtf(kc)
c
c    -legendre transform from (k,l,m) to (k,i,m)
c
c---------------------------------------------------------------
      include 'param.f'
      include 'com4.f'
      include 'com5.f'
      include 'com6.f'
c
      complex 
     $cs1(nlma),cs2(nlma),cs3(nlma),cs4(nlma),cs5(nlma)
     $,cs6(nlma),cs7(nlma),cs8(nlma),cs9(nlma),cs10(nlma)
     $,cs11(nlma),cs15(nlma),cs16(nlma)
     $,cs17(nlma),cs18(nlma),cs19(nlma),cs20(nlma)
     $,cs21(nlma),cs22(nlma),cs23(nlma),cs24(nlma)
      complex qim(nrp)
c
      do 200 mca=1,ncp
        m=(mca-1)*minc
        qim(mca)=ai*float(m)
  200 continue
c
      do 201 ic=1,ni
         do 202 mca=1,ncp
            sc(mca,ic)=0.
            vrc(mca,ic)=0.
            vtc(mca,ic)=0.
            vpc(mca,ic)=0.
            cvrc(mca,ic)=0.
            dvrdrc(mca,ic)=0.
            dvtdrc(mca,ic)=0.
            dvpdrc(mca,ic)=0.
            dvrdtc(mca,ic)=0.
            dvrdpc(mca,ic)=0.
            dvtdpc(mca,ic)=0.
            dvpdpc(mca,ic)=0.
            brc(mca,ic)=0.
            btc(mca,ic)=0.
            bpc(mca,ic)=0.
            cbrc(mca,ic)=0.
            cbtc(mca,ic)=0.
            cbpc(mca,ic)=0.
  202    continue
  201 continue
c
      do 203 lm=nlma,2,-1
         mca=mclma(lm)
         cs1(lm)=w(lm,kc)*ql(lm,3)
         cs2(lm)=2.*dw(lm,kc)
         cs3(lm)=z(lm,kc)*qim(mca)
         cs4(lm)=2.*dw(lm,kc)*qim(mca)
         cs5(lm)=-z(lm,kc)
         cs6(lm)=z(lm,kc)*ql(lm,3)
         cs7(lm)=2.*dw(lm,kc)*ql(lm,3)
         cs8(lm)=4.*ddw(lm,kc)
         cs9(lm)=2.*dz(lm,kc)*qim(mca)
         cs10(lm)=4.*ddw(lm,kc)*qim(mca)
         cs11(lm)=-2.*dz(lm,kc)
         cs15(lm)=b(lm,kc)*ql(lm,3)
         cs16(lm)=2.*db(lm,kc)
         cs17(lm)=aj(lm,kc)*qim(mca)
         cs18(lm)=2.*qim(mca)*db(lm,kc)
         cs19(lm)=-aj(lm,kc)
         cs20(lm)=aj(lm,kc)*ql(lm,3)
         cs21(lm)=2.*dj(lm,kc)
         cs22(lm)=qim(mca)*(ql(lm,3)*qk(kc,1)*b(lm,kc)-4.*ddb(lm,kc))
         cs23(lm)=2.*qim(mca)*dj(lm,kc)
         cs24(lm)=-ql(lm,3)*qk(kc,1)*b(lm,kc)+4.*ddb(lm,kc)
  203 continue
c
       do 204 ic=1,ni
         do 204 lm=nlma,2,-1
           mca=mclma(lm)
            sc(mca,ic)=sc(mca,ic)+s(lm,kc)*aleg1(lm,ic)
            vrc(mca,ic)=vrc(mca,ic)+cs1(lm)*aleg1(lm,ic)
            cvrc(mca,ic)=cvrc(mca,ic)+cs6(lm)*aleg1(lm,ic)
            dvrdrc(mca,ic)=dvrdrc(mca,ic)+cs7(lm)*aleg1(lm,ic)
            brc(mca,ic)=brc(mca,ic)+cs15(lm)*aleg1(lm,ic)
            cbrc(mca,ic)=cbrc(mca,ic)+cs20(lm)*aleg1(lm,ic)
  204  continue
       do 206 ic=1,ni
         do 206 lm=nlma,2,-1
           mca=mclma(lm)
            dvrdtc(mca,ic)=dvrdtc(mca,ic)+cs1(lm)*aleg3(lm,ic)
            vtc(mca,ic)=vtc(mca,ic)+(cs2(lm)*aleg3(lm,ic)+
     $         cs3(lm)*aleg1(lm,ic))
            vpc(mca,ic)=vpc(mca,ic)+(cs4(lm)*aleg1(lm,ic)+
     $         cs5(lm)*aleg3(lm,ic))
            dvtdrc(mca,ic)=dvtdrc(mca,ic)+(cs8(lm)*aleg3(lm,ic)+
     $         cs9(lm)*aleg1(lm,ic))
            dvpdrc(mca,ic)=dvpdrc(mca,ic)+(cs10(lm)*aleg1(lm,ic)+
     $         cs11(lm)*aleg3(lm,ic))
  206  continue
       do 208 ic=1,ni
         do 208 lm=nlma,2,-1
           mca=mclma(lm)
            btc(mca,ic)=btc(mca,ic)+(cs17(lm)*aleg1(lm,ic)+
     $         cs16(lm)*aleg3(lm,ic))
            bpc(mca,ic)=bpc(mca,ic)+(cs19(lm)*aleg3(lm,ic)+
     $         cs18(lm)*aleg1(lm,ic))
            cbtc(mca,ic)=cbtc(mca,ic)+(cs21(lm)*aleg3(lm,ic)+
     $         cs22(lm)*aleg1(lm,ic))
            cbpc(mca,ic)=cbpc(mca,ic)+(cs23(lm)*aleg1(lm,ic)+
     $         cs24(lm)*aleg3(lm,ic))
  208  continue
c
      do 215 ic=1,ni
         sc(1,ic)=sc(1,ic)+s(1,kc)*aleg1(1,ic)
  215 continue
c
      do 220 ic=1,ni
        do 220 mca=1,nmafa
            dvrdpc(mca,ic)=         qim(mca)*vrc(mca,ic)
            dvtdpc(mca,ic)=qi(ic,1)*qim(mca)*vtc(mca,ic)
            dvpdpc(mca,ic)=qi(ic,1)*qim(mca)*vpc(mca,ic)
  220 continue
c
      return
      end
