      subroutine spherictf(alm,aij)
c
c    -spherical harmonic transform from alm(l,m) to aij(i,j)
c
c---------------------------------------------------------------
      include 'param.f'
      include 'com4.f'
c
      complex 
     $cs1(nlma),alm(nlma),aij(ncp,ni)
c
      do 201 ic=1,ni
         do 202 mca=1,ncp
            aij(mca,ic)=0.
  202    continue
  201 continue
c
      do 203 lm=nlma,2,-1
         cs1(lm)=alm(lm)*ql(lm,3)
  203 continue
c
       do 204 ic=1,ni
         do 204 lm=nlma,2,-1
            mca=mclma(lm)
            aij(mca,ic)=aij(mca,ic)+cs1(lm)*aleg1(lm,ic)
  204  continue
c
      call fourtf(aij,work,trigsf,ifaxf,1,nrp,nja,ni,1)
c
      return
      end
