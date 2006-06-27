      subroutine mapdata(nnold,niold,njold,minco)
c
c---------------------------------------------------------------
c
c  map data from input file with different grid structure in the  
c  angular variables or different longitudinal symmetry (urc)
c
c  called in prep
c
c
      include 'param.f'
      parameter (nn32=3*nnp1/2)
      parameter (noldsize=3*nlma*nnp1)
      include 'com1.f'
      include 'com3.f'
      include 'com4.f'
      include 'com5.f'

      complex wo(noldsize),zo(noldsize),po(noldsize),so(noldsize)
      dimension rold(nn32)
      equivalence (wo,dw),(zo,dz),(po,db),(so,aleg1)
c
      if(mod(minco,minc).ne.0) 
     $  write(6,'('' Warning: Incompatible minc/minco= '',2i3)')
c     
      nnop1=nnold+1
      lmaxo= njold/3
      mmaxo= (lmaxo/minco) * minco
      nlmao= mmaxo*(lmaxo+1)/minco -mmaxo*(mmaxo-minco)/(2*minco)
     &                             +lmaxo-mmaxo+1
      ndata1=nlmao*nnop1
      ndata=nlmao*nnold
      if(nnold.gt.nn32) then
        write(6,'('' nnold='',i3,'' too large'')') nnold
        stop '54'
      endif
      if(ndata1.gt.3*nlma*nnp1) then
        write(6,'('' Old data set is too large '')')
        write(6,'('' New/Old Lmax= '',2I4,''  Mmax= '',2i4,  
     &   ''  Minc= '',2i3,'' Nr= '',2i3)')
     &  lmax,lmaxo,mmax,mmaxo,minc,minco,nn,nnold
        write(6,'('' Total old data > 3* new data !'')')
        stop '47'
      endif
      if(nn.ne.nnold) then
       pik=pi/float(nnold-1)
       do kco=1,nnold
         rold(kco)=radbot+0.5*(1.+cos(float(kco-1)*pik))
       enddo
      endif
c
      write(6,'(/'' Mapping data on different grid '')')
      write(6,'('' Old/New  Lmax= '',2I4,''  Mmax= '',2I4,
     $ ''  Minc= '',2I3,''  Nlma= '',2I5/)') lmaxo,lmax,mmaxo,
     $ mmax,minco,minc,nlmao,nlma
      if(nnold.ne.nn) write(6,'('' Old/New nn='',2i4)') nnold,nn
c
c  read and copy data for w,z,p,s
c
      read(8) (wo(i),i=1,ndata1),
     &        (zo(i),i=1,ndata1),
     &        (po(i),i=1,ndata1),
     &        (so(i),i=1,ndata1)
c
      call copydat(w,wo,nlma,nlmao,lmax,lmaxo,mmax,mmaxo,minc,minco,
     $ nnp1,nnop1,nnp1,nn32,r,rold)
      call copydat(z,zo,nlma,nlmao,lmax,lmaxo,mmax,mmaxo,minc,minco,
     $ nnp1,nnop1,nnp1,nn32,r,rold)
      call copydat(p,po,nlma,nlmao,lmax,lmaxo,mmax,mmaxo,minc,minco,
     $ nnp1,nnop1,nnp1,nn32,r,rold)
      call copydat(s,so,nlma,nlmao,lmax,lmaxo,mmax,mmaxo,minc,minco,
     $ nnp1,nnop1,nnp1,nn32,r,rold)
c
c  read and copy data for dsdt1,dwdt1,dzdt1,dpdt1
c
      read(8) (so(i),i=1,ndata),
     &        (wo(i),i=1,ndata),
     &        (zo(i),i=1,ndata),
     &        (po(i),i=1,ndata)
c
      call copydat(dsdt1,so,nlma,nlmao,lmax,lmaxo,mmax,mmaxo,minc,minco,
     $ nn,nnold,nnp1,nn32,r,rold)
      call copydat(dwdt1,wo,nlma,nlmao,lmax,lmaxo,mmax,mmaxo,minc,minco,
     $ nn,nnold,nnp1,nn32,r,rold)
      call copydat(dzdt1,zo,nlma,nlmao,lmax,lmaxo,mmax,mmaxo,minc,minco,
     $ nn,nnold,nnp1,nn32,r,rold)
      call copydat(dpdt1,po,nlma,nlmao,lmax,lmaxo,mmax,mmaxo,minc,minco,
     $ nn,nnold,nnp1,nn32,r,rold)
c
c  read and copy data for b, aj, dbdt1,djdt1
c
      if(init.ne.0) return
c
      read(8) (so(i),i=1,ndata1),
     &        (wo(i),i=1,ndata1),
     &        (zo(i),i=1,ndata),
     &        (po(i),i=1,ndata)
c
      call copydat(b,so,nlma,nlmao,lmax,lmaxo,mmax,mmaxo,minc,minco,
     $ nnp1,nnop1,nnp1,nn32,r,rold)
      call copydat(aj,wo,nlma,nlmao,lmax,lmaxo,mmax,mmaxo,minc,minco,
     $ nnp1,nnop1,nnp1,nn32,r,rold)
      call copydat(dbdt1,zo,nlma,nlmao,lmax,lmaxo,mmax,mmaxo,minc,minco,
     $ nn,nnold,nnp1,nn32,r,rold)
      call copydat(djdt1,po,nlma,nlmao,lmax,lmaxo,mmax,mmaxo,minc,minco,
     $ nn,nnold,nnp1,nn32,r,rold)
c
c  if starting from data file with longitudinal symmetry, add
c  weak non-axisymmetric dipole component
c
      if(minc.gt.1 .or. minco.eq.1 .or. tipdipole.eq.0.0) return
c
      do kc=1,nnp1
        b(lmax+2,kc)=0.707106781*tipdipole*b(2,kc)
      enddo
c
      return
      end 
